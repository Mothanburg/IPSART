int argmaxf(float *array, int len)
{
    int idx_max = 0;
    for (int i = 1; i < len; i++)
    {
        idx_max += isless(array[idx_max], array[i]) * (i - idx_max);
    }
    return idx_max;
}

__kernel void get_span(__global const float *c11, __global const float *c22, __global const float *c33,
                       __global float *span)
{
    int idx = get_global_id(0);
    span[idx] = c11[idx] + c22[idx] + c33[idx];
}

__kernel void filt_cij(__global const float *span, __global const float *cov_mat_ij, int height, int width,
                            __local float *shared_mem, int shared_height, int shared_width,
                            __constant float *flt_windows /* 7 * 7 * 8 */, int n_looks,
                            __global float *filted_cij)
{
    int wg_row = get_group_id(0);
    int wg_col = get_group_id(1);
    int wg_rows = get_num_groups(0);
    int wg_cols = get_num_groups(1);

    int local_row = get_local_id(0);
    int local_col = get_local_id(1);
    int wg_height = get_local_size(0);
    int wg_width = get_local_size(1);
    for (int c = local_col; c < shared_width; c += wg_width)
    {
        for (int r = local_row; r < shared_height; r += wg_height)
        {
            int row = min(max(wg_rows * wg_row + r - 3, 0), height - 1);
            int col = min(max(wg_cols * wg_col + c - 3, 0), width - 1);
            shared_mem[c * shared_height + r] = span[col * height + row];
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    float mean_mat[3][3];
    __attribute__((opencl_unroll_hint))
    for (int i = 0; i < 3; i++)
    {
        __attribute__((opencl_unroll_hint))
        for (int j = 0; j < 3; j++)
        {
            mean_mat[i][j] = 0.0;
            __attribute__((opencl_unroll_hint))
            for (int dc = 0; dc < 3; dc++)
            {
                __attribute__((opencl_unroll_hint))
                for (int dr = 0; dr < 3; dr++)
                {
                    int lr = local_row + i * 2 + dr;
                    int lc = local_col + j * 2 + dc;
                    mean_mat[i][j] += shared_mem[lc * shared_height + lr];
                }
            }
            // mean_mat[i][j] /= 9.0; This line can be unnecessary.
        }
    }

    float sums[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    sums[0] -= mean_mat[0][0];
    sums[0] += mean_mat[0][2];
    sums[0] -= mean_mat[1][0];
    sums[0] += mean_mat[1][2];
    sums[0] -= mean_mat[2][0];
    sums[0] += mean_mat[2][2];

    sums[1] += mean_mat[0][1];
    sums[1] += mean_mat[0][2];
    sums[1] -= mean_mat[1][1];
    sums[1] += mean_mat[1][2];
    sums[1] -= mean_mat[2][0];
    sums[1] -= mean_mat[2][1];

    sums[2] += mean_mat[0][0];
    sums[2] += mean_mat[0][1];
    sums[2] += mean_mat[0][2];
    sums[2] -= mean_mat[2][0];
    sums[2] -= mean_mat[2][1];
    sums[2] -= mean_mat[2][2];

    sums[3] += mean_mat[0][0];
    sums[3] += mean_mat[0][1];
    sums[3] += mean_mat[1][0];
    sums[3] -= mean_mat[1][2];
    sums[3] -= mean_mat[2][1];
    sums[3] -= mean_mat[2][2];

    sums[4] += mean_mat[0][0];
    sums[4] -= mean_mat[0][2];
    sums[4] += mean_mat[1][0];
    sums[4] -= mean_mat[1][2];
    sums[4] += mean_mat[2][0];
    sums[4] -= mean_mat[2][2];

    sums[5] -= mean_mat[0][1];
    sums[5] -= mean_mat[0][2];
    sums[5] += mean_mat[1][1];
    sums[5] -= mean_mat[1][2];
    sums[5] += mean_mat[2][0];
    sums[5] += mean_mat[2][1];

    sums[6] -= mean_mat[0][0];
    sums[6] -= mean_mat[0][1];
    sums[6] -= mean_mat[0][2];
    sums[6] += mean_mat[2][0];
    sums[6] += mean_mat[2][1];
    sums[6] += mean_mat[2][2];

    sums[7] -= mean_mat[0][0];
    sums[7] -= mean_mat[0][1];
    sums[7] -= mean_mat[1][0];
    sums[7] += mean_mat[1][2];
    sums[7] += mean_mat[2][1];
    sums[7] += mean_mat[2][2];

    int window_id = argmaxf(sums, 8);

    float mean_z = 0.0;
    __attribute__((opencl_unroll_hint))
    for (int j = 0; j < 7; j++)
    {
        __attribute__((opencl_unroll_hint))
        for (int i = 0; i < 7; i++)
        {
            mean_z += shared_mem[(local_col + j) * shared_height + local_row + i] *
                      flt_windows[49 * window_id + j * 7 + i];
        }
    }
    mean_z /= 28.0;

    float z_var = 0.0;
    __attribute__((opencl_unroll_hint))
    for (int j = 0; j < 7; j++)
    {
        __attribute__((opencl_unroll_hint))
        for (int i = 0; i < 7; i++)
        {
            z_var += pown(shared_mem[(local_col + j) * shared_height + local_row + i] *
                              flt_windows[49 * window_id + j * 7 + i] - mean_z,
                          2);
        }
    }
    z_var /= 28.0;

    float v_var = 1 / (float)n_looks;
    float x_var = (z_var - pown(mean_z, 2) * v_var) / (1.0 + v_var);
    float factor = x_var / z_var;
    uint no_nan = ~((uint)isnan(factor) * 0xFFFFFFFF) & *(uint*)&factor;
    factor = *(float *)&no_nan;

    barrier(CLK_LOCAL_MEM_FENCE);

    for (int c = local_col; c < shared_width; c += wg_width)
    {
        for (int r = local_row; r < shared_height; r += wg_height)
        {
            int row = min(max(wg_rows * wg_row + r - 3, 0), height - 1);
            int col = min(max(wg_cols * wg_col + c - 3, 0), width - 1);
            shared_mem[c * shared_height + r] = cov_mat_ij[col * height + row];
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    float mean_cij = 0.0;
    __attribute__((opencl_unroll_hint))
    for (int j = 0; j < 7; j++)
    {
        __attribute__((opencl_unroll_hint))
        for (int i = 0; i < 7; i++)
        {
            mean_cij += shared_mem[(local_col + j) * shared_height + local_row + i] *
                        flt_windows[49 * window_id + j * 7 + i];
        }
    }
    mean_cij /= 28.0;

    int g_row = get_global_id(0);
    int g_col = get_global_id(1);
    filted_cij[g_col * height + g_row] = mean_cij + factor * 
                                         (shared_mem[(local_col + 3) * shared_height + local_row + 3] - mean_cij);
}