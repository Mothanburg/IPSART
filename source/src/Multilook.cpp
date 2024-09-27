#include "GaRS.h"

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

#include "Multilook.cl.h"

using namespace std;

template <typename TData>
static int multilook(int inRows, int inCols, const TData *input, int rowLook,
                     int colLook, int outRows, int outCols, TData *output) {
  try {
    cl::Context context = cl::Context::getDefault();
    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    cl::CommandQueue queue(context);

    // Compile the program
    cl::Program program(context, SRC_MULTILOOK);
    if constexpr (is_same_v<TData, float>) {
      program.build(device, "-cl-std=CL2.0");
    } else {
      /* TData is double */
      program.build(device, "-cl-std=CL2.0 -DENABLE_FP64");
    }

    // Create buffers
    cl::Buffer buf_input(context, input, input + (inRows * inCols), true);
    cl::Buffer buf_output(context, CL_MEM_WRITE_ONLY,
                          sizeof(TData) * outRows * outCols);

    // Prepare to execute the kernel
    cl::NDRange global_size(outRows, outCols);

    cl::Kernel krnl(program, "multilook");
    krnl.setArg(0, inCols);
    krnl.setArg(1, buf_input);
    krnl.setArg(2, outCols);
    krnl.setArg(3, buf_output);
    krnl.setArg(4, rowLook);
    krnl.setArg(5, colLook);
    queue.enqueueNDRangeKernel(krnl, cl::NullRange, global_size);

    queue.enqueueReadBuffer(buf_output, false, 0,
                            sizeof(TData) * outRows * outCols, output);

    queue.finish();
  } catch (cl::Error &e) {
    return e.err();
  }
  return 0;
}

int Multilookf(int inRows, int inCols, const float *input, int rowLook,
               int colLook, int outRows, int outCols, float *output) {
  return multilook(inRows, inCols, input, rowLook, colLook, outRows, outCols,
                   output);
}

int Multilookd(int inRows, int inCols, const double *input, int rowLook,
               int colLook, int outRows, int outCols, double *output) {
  return multilook(inRows, inCols, input, rowLook, colLook, outRows, outCols,
                   output);
}