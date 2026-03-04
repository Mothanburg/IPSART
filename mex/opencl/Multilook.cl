// Convention: Leading dimension is dim 0, and data storage is col-major
#ifdef ENABLE_FP64
#pragma OPENCL EXTENSION cl_khr_fp64:enable
#define dtype double
#else
#define dtype float
#endif

__kernel void
multilook(int input_rows, __global const dtype *input,
          int output_rows, __global dtype *output,
          int row_look, int col_look, dtype inv_look)
{
    int row_out = get_global_id(0);
    int col_out = get_global_id(1);

    int row_in = row_out * row_look;
    int col_in = col_out * col_look;

    dtype sum = 0;
    for (int dc = 0; dc < col_look; dc++)
    {
        int col_offset = input_rows * (col_in + dc) + row_in;
        for (int dr = 0; dr < row_look; dr++)
        {
            sum += input[col_offset + dr];
        }
    }

    output[output_rows * col_out + row_out] = sum * inv_look;
}