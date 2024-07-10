#define _SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING

#include <array>
#include <numeric>

#include "GaRS.h"

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <cl/opencl.hpp>

#include "test_kernel.src.h"


int GaRSTestOpenCL() {
  constexpr auto TOTAL_LEN = 47;
  constexpr auto GROUP_LEN = 3;
  try {
    // Get the first device of the default platform
    cl::Context context = cl::Context::getDefault();
    auto device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    cl::CommandQueue queue(context);

    // Compile the program
    cl::Program program(context, src_test_kernel);
    program.build(device, "-cl-std=CL2.0");

    // Set data
    std::array<float, TOTAL_LEN> v1, v2;
    v1.fill(1.0f);
    v2.fill(-1.0f);

    cl::Buffer buf_v1(context, v1.begin(), v1.end(), true);
    cl::Buffer buf_v2(context, v2.begin(), v2.end(), true);

    float *result = new float[TOTAL_LEN];
    cl::Buffer buf_out(context, CL_MEM_WRITE_ONLY, sizeof(float) * TOTAL_LEN);

    // Create kernel
    cl::Kernel krnl(program, "vec_add");
    krnl.setArg(0, buf_v1);
    krnl.setArg(1, buf_v2);
    krnl.setArg(2, sizeof(float) * 3, nullptr);
    krnl.setArg(3, buf_out);

    // Execute the kernel
    queue.enqueueNDRangeKernel(krnl, cl::NullRange, cl::NDRange(TOTAL_LEN),
                               cl::NDRange(GROUP_LEN));
    queue.enqueueReadBuffer(buf_out, false, 0, sizeof(float) * TOTAL_LEN,
                            result);
    queue.finish();

    if (std::accumulate(result, result + TOTAL_LEN, 0.0f) != 0.0f) {
      return -100;
    }
  } catch (cl::Error &e) {
    return e.err();
  }
  return 0;
}