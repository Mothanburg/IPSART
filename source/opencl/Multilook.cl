// Convention: Leading dimension is dim 0, and data storage is row-major
#ifdef ENABLE_FP64
#pragma OPENCL EXTENSION cl_khr_fp64:enable
#define dtype double
#else
#define dtype float
#endif

__kernel void
multilook(int input_cols, __global const dtype *input,
          int output_cols, __global dtype *output,
          int row_look, int col_look)
{
    int row_out = get_global_id(0);
    int col_out = get_global_id(1);

    int row_in = row_out * row_look;
    int col_in = col_out * col_look;

    dtype sum = 0;
    for (int dr = 0; dr < row_look; dr++)
    {
        for (int dc = 0; dc < col_look; dc++)
        {
            sum += input[input_cols * (row_in + dr) + col_in + dc];
        }
    }

    output[output_cols * row_out + col_out] = sum / (row_look * col_look);
}