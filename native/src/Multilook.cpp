#ifdef _MSC_VER
#define _SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING
#endif

#include "GaRS.h"

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

#include "mtlk_kernel.src.h"

using namespace std;

int Multilook(int height, int width, const float* image, int rowLook,
              int colLook, int outHeight, int outWidth, float* looked) {
  try {
    cl::Context context = cl::Context::getDefault();
    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    cl::CommandQueue queue(context);

    // Compile the program
    cl::Program program(context, src_mtlk_kernel);
    program.build(device, "-cl-std=CL2.0");

    // Create buffers
    cl::Buffer buf_image(context, image, image + (height * width), true);
    cl::Buffer buf_looked(context, CL_MEM_WRITE_ONLY,
                          sizeof(float) * outHeight * outWidth);

    // Prepare to execute the kernel
    cl::NDRange global_size(outHeight, outWidth);

    cl::Kernel krnl(program, "multilook");
    krnl.setArg(0, buf_image);
    krnl.setArg(1, width);
    krnl.setArg(2, rowLook);
    krnl.setArg(3, colLook);
    krnl.setArg(4, buf_looked);
    krnl.setArg(5, outWidth);
    queue.enqueueNDRangeKernel(krnl, cl::NullRange, global_size);

    queue.enqueueReadBuffer(buf_looked, false, 0,
                            sizeof(float) * outHeight * outWidth, looked);

    queue.finish();

  } catch (cl::Error& e) {
    return e.err();
  }
  return 0;
}