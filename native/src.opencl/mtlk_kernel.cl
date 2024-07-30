// Convention: row-major
__kernel void
multilook(__global const float *image, int img_width,
          int row_look, int col_look,
          __global float *looked, int looked_width)
{
    // The global id is for looked image.
    int r_looked = get_global_id(0);
    int c_looked = get_global_id(1);

    int r0_img = r_looked * row_look;
    int c0_img = c_looked * col_look;

    float sum = 0.0;
    for (int dr = 0; dr < row_look; dr++)
    {
        for (int dc = 0; dc < col_look; dc++)
        {
            sum += image[img_width * (r0_img + dr) + c0_img + dc];
        }
    }

    looked[looked_width * r_looked + c_looked] = sum / (row_look * col_look);
}