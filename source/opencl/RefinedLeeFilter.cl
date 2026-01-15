// Convention: Leading dimension is dim 0, and data storage is row-major
#ifdef ENABLE_FP64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define dtype double
#else
#define dtype float
#endif

__kernel void
#ifdef MAT_SIZE_3X3
span_calc(__global dtype *span,
          __global const dtype *c11,
          __global const dtype *c22,
          __global const dtype *c33)
#else
span_calc(__global dtype *span,
          __global const dtype *c11,
          __global const dtype *c22)
#endif
{
    int idx = get_global_id(0) * get_global_size(1) + get_global_id(1);
#ifdef MAT_SIZE_3X3
    span[idx] = c11[idx] + c22[idx] + c33[idx];
#else
    span[idx] = c11[idx] + c22[idx];
#endif
}

// A 7x7 refined lee filter
__kernel void
page_filting(int rows, int cols,
             __global const dtype *span, // the SPAN of the PolSAR image
             int shared_rows, int shared_cols, __local dtype *shared_mem,
             __constant dtype *prewitt,   // the size of prewitt operator templates is 7x7x8
             int look_num,                // look number
             __global const dtype *input, // input Cij to be filted
             __global dtype *output)
{
    int wg_row = get_group_id(0);
    int wg_col = get_group_id(1);
    int enqd_wg_rows = get_enqueued_local_size(0);
    int enqd_wg_cols = get_enqueued_local_size(1);
    int wg_rows = get_local_size(0);
    int wg_cols = get_local_size(1);
    int local_row = get_local_id(0);
    int local_col = get_local_id(1);
    for (int r = local_row; r < shared_rows; r += wg_rows)
    {
        for (int c = local_col; c < shared_cols; c += wg_cols)
        {
            int row = min(max(wg_row * enqd_wg_rows + r - 3, 0), rows - 1);
            int col = min(max(wg_col * enqd_wg_cols + c - 3, 0), cols - 1);
            shared_mem[r * shared_cols + c] = span[row * cols + col];
        }
    }

    work_group_barrier(CLK_LOCAL_MEM_FENCE);

    dtype mean_mat[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    __attribute__((opencl_unroll_hint)) for (int i = 0; i < 3; i++)
    {
        __attribute__((opencl_unroll_hint)) for (int j = 0; j < 3; j++)
        {
            __attribute__((opencl_unroll_hint)) for (int dr = 0; dr < 3; dr++)
            {
                __attribute__((opencl_unroll_hint)) for (int dc = 0; dc < 3; dc++)
                {
                    int lr = local_row + i * 2 + dr;
                    int lc = local_col + j * 2 + dc;
                    mean_mat[i * 3 + j] += shared_mem[lr * shared_cols + lc];
                }
            }
        }
    }

    // Calculate Prewitt id

    // window id
    dtype sums[4];
    sums[0] = abs(mean_mat[2] + mean_mat[5] + mean_mat[8] - mean_mat[0] - mean_mat[3] - mean_mat[6]); 
    sums[1] = abs(mean_mat[1] + mean_mat[2] + mean_mat[5] - mean_mat[3] - mean_mat[6] - mean_mat[7]);
    sums[2] = abs(mean_mat[0] + mean_mat[1] + mean_mat[2] - mean_mat[6] - mean_mat[7] - mean_mat[8]);
    sums[3] = abs(mean_mat[0] + mean_mat[1] + mean_mat[3] - mean_mat[5] - mean_mat[7] - mean_mat[8]);
    //sums[4] = mean_mat[0] + mean_mat[3] + mean_mat[6] - mean_mat[2] - mean_mat[5] - mean_mat[8];
    //sums[5] = mean_mat[3] + mean_mat[6] + mean_mat[7] - mean_mat[1] - mean_mat[2] - mean_mat[5];
    //sums[6] = mean_mat[6] + mean_mat[7] + mean_mat[8] - mean_mat[0] - mean_mat[1] - mean_mat[2];
    //sums[7] = mean_mat[5] + mean_mat[7] + mean_mat[8] - mean_mat[0] - mean_mat[1] - mean_mat[3];

    int window_id = 0;
    __attribute__((opencl_unroll_hint)) for (int i = 1; i < 4; i++)
    {
        window_id += isless(sums[window_id], sums[i]) * (i - window_id);
    }

    dtype distance[8];
    distance[0] = abs(mean_mat[5] - mean_mat[4]); // m23 - m22
    distance[1] = abs(mean_mat[2] - mean_mat[4]); // m13 - m22
    distance[2] = abs(mean_mat[1] - mean_mat[4]); // m12 - m22
    distance[3] = abs(mean_mat[0] - mean_mat[4]); // m11 - m22
    distance[4] = abs(mean_mat[3] - mean_mat[4]); // m21 - m22
    distance[5] = abs(mean_mat[6] - mean_mat[4]); // m31 - m22
    distance[6] = abs(mean_mat[7] - mean_mat[4]); // m32 - m22
    distance[7] = abs(mean_mat[8] - mean_mat[4]); // m33 - m22

    prewitt_id = isgreater(distance[window_id], distance[window_id + 4]) * 4 + window_id;

    dtype z_mean = 0.0;
    __attribute__((opencl_unroll_hint)) for (int i = 0; i < 7; i++)
    {
        __attribute__((opencl_unroll_hint)) for (int j = 0; j < 7; j++)
        {
            z_mean += shared_mem[(local_row + i) * shared_cols + local_col + j] *
                      prewitt[49 * prewitt_id + i * 7 + j];
        }
    }
    z_mean /= 28.0;

    dtype z_var = 0.0;
    __attribute__((opencl_unroll_hint)) for (int i = 0; i < 7; i++)
    {
        __attribute__((opencl_unroll_hint)) for (int j = 0; j < 7; j++)
        {
            dtype z = shared_mem[(local_row + i) * shared_cols + local_col + j] *
                      prewitt[49 * prewitt_id + i * 7 + j];
            z_var += (z - z_mean) * (z - z_mean);
        }
    }
    z_var /= 28.0;

    dtype v_var = 1.0 / look_num;
    dtype x_var = (z_var - z_mean * z_mean * v_var) / (1.0 + v_var);
    dtype factor = (x_var + 1e-30) / (z_var + 1e-30); // avoid nan

    work_group_barrier(CLK_LOCAL_MEM_FENCE);

    for (int r = local_row; r < shared_rows; r += wg_rows)
    {
        for (int c = local_col; c < shared_cols; c += wg_cols)
        {
            int row = min(max(wg_row * enqd_wg_rows + r - 3, 0), rows - 1);
            int col = min(max(wg_col * enqd_wg_cols + c - 3, 0), cols - 1);
            shared_mem[r * shared_cols + c] = input[row * cols + col];
        }
    }

    work_group_barrier(CLK_LOCAL_MEM_FENCE);

    dtype mean_input = 0.0;
    __attribute__((opencl_unroll_hint)) for (int i = 0; i < 7; i++)
    {
        __attribute__((opencl_unroll_hint)) for (int j = 0; j < 7; j++)
        {
            mean_input += shared_mem[(local_row + i) * shared_cols + local_col + j] *
                          prewitt[49 * prewitt_id + i * 7 + j];
        }
    }
    mean_input /= 28.0;

    int g_row = get_global_id(0);
    int g_col = get_global_id(1);
    output[g_row * cols + g_col] = mean_input +
                                   factor * (shared_mem[(local_row + 3) * shared_cols + local_col + 3] - mean_input);
}