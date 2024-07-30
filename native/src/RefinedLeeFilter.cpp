#ifdef _MSC_VER
#define _SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING
#endif

#include <array>

#include "GaRS.h"
#include "utils.hpp"

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

#include "rlf_kernel.src.h"

using namespace std;

// clang-format off

// Prewitt operators (Row-major)
constexpr array<float, 7ull * 7 * 8> Prewitt{
  // Prewitt 1
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  // Prewitt 2
  1, 1, 1, 1, 1, 1, 1,
  0, 1, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 1, 1, 1,
  0, 0, 0, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 1,
  // Prewitt 3
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  // Prewitt 4
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 0,
  1, 1, 1, 1, 1, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0,
  // Prewitt 5
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  // Prewitt 6
  1, 0, 0, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0,
  1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 1, 0, 0,
  1, 1, 1, 1, 1, 1, 0,
  1, 1, 1, 1, 1, 1, 1,
  // Prewitt 7
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  // Prewitt 8
  0, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 1, 1,
  0, 0, 0, 0, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 1, 1,
  0, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1
};
// clang-format on

int RefinedLeeFilter3x3(int nLooks, int height, int width, const float *c11,
                        const float *c22, const float *c33, const float *c12r,
                        const float *c13r, const float *c23r, const float *c12i,
                        const float *c13i, const float *c23i, float *outC11,
                        float *outC22, float *outC33, float *outC12r,
                        float *outC13r, float *outC23r, float *outC12i,
                        float *outC13i, float *outC23i) {
  try {
    cl::Context context = cl::Context::getDefault();
    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    cl::CommandQueue queue(context);

    // Compile the program
    cl::Program program(context, src_rlf_kernel);
    program.build(device, "-cl-std=CL2.0");

    // Prepare to execute the kernel
    auto total_len = height * width;
    cl::NDRange global_size(height, width);

    auto max_size = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    auto group_height = int_sqrt(max_size);
    cl::NDRange group_size(group_height, group_height);

    // Create buffers
    cl::Buffer buf_c11(context, c11, c11 + total_len, true);
    cl::Buffer buf_c22(context, c22, c22 + total_len, true);
    cl::Buffer buf_c33(context, c33, c33 + total_len, true);
    cl::Buffer buf_c12r(context, c12r, c12r + total_len, true);
    cl::Buffer buf_c13r(context, c13r, c13r + total_len, true);
    cl::Buffer buf_c23r(context, c23r, c23r + total_len, true);
    cl::Buffer buf_c12i(context, c12i, c12i + total_len, true);
    cl::Buffer buf_c13i(context, c13i, c13i + total_len, true);
    cl::Buffer buf_c23i(context, c23i, c23i + total_len, true);
    cl::Buffer buf_span(context, CL_MEM_READ_WRITE, sizeof(float) * total_len);

    // We use the first device to calculate SPAN
    cl::Kernel krnl_getspan(program, "mat_add3");
    krnl_getspan.setArg(0, buf_span);
    krnl_getspan.setArg(1, buf_c11);
    krnl_getspan.setArg(2, buf_c22);
    krnl_getspan.setArg(3, buf_c33);
    queue.enqueueNDRangeKernel(krnl_getspan, cl::NullRange, global_size);

    // Now we can start filting process
    cl::Kernel krnl_rlf(program, "filt_cij");
    krnl_rlf.setArg(0, buf_span);
    krnl_rlf.setArg(1, height);
    krnl_rlf.setArg(2, width);

    auto shared_height = static_cast<int>(group_height) + 8;
    krnl_rlf.setArg(3, sizeof(float) * shared_height * shared_height, nullptr);
    krnl_rlf.setArg(4, shared_height);
    krnl_rlf.setArg(5, shared_height);

    cl::Buffer buf_prwt(context, Prewitt.begin(), Prewitt.end(), true);
    krnl_rlf.setArg(6, buf_prwt);
    krnl_rlf.setArg(7, nLooks);

    // For every channel

    cl::Buffer buf_oc11(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c11);
    krnl_rlf.setArg(9, buf_oc11);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc11, false, 0, sizeof(float) * total_len,
                            outC11);

    cl::Buffer buf_oc22(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c22);
    krnl_rlf.setArg(9, buf_oc22);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc22, false, 0, sizeof(float) * total_len,
                            outC22);

    cl::Buffer buf_oc33(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c33);
    krnl_rlf.setArg(9, buf_oc33);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc33, false, 0, sizeof(float) * total_len,
                            outC33);

    cl::Buffer buf_oc12r(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c12r);
    krnl_rlf.setArg(9, buf_oc12r);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc12r, false, 0, sizeof(float) * total_len,
                            outC12r);

    cl::Buffer buf_oc13r(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c13r);
    krnl_rlf.setArg(9, buf_oc13r);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc13r, false, 0, sizeof(float) * total_len,
                            outC13r);

    cl::Buffer buf_oc23r(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c23r);
    krnl_rlf.setArg(9, buf_oc23r);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc23r, false, 0, sizeof(float) * total_len,
                            outC23r);

    cl::Buffer buf_oc12i(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c12i);
    krnl_rlf.setArg(9, buf_oc12i);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc12i, false, 0, sizeof(float) * total_len,
                            outC12i);

    cl::Buffer buf_oc13i(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c13i);
    krnl_rlf.setArg(9, buf_oc13i);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc13i, false, 0, sizeof(float) * total_len,
                            outC13i);

    cl::Buffer buf_oc23i(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c23i);
    krnl_rlf.setArg(9, buf_oc23i);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc23i, false, 0, sizeof(float) * total_len,
                            outC23i);

    queue.finish();

  } catch (cl::Error &e) {
    return e.err();
  }
  return 0;
}

int RefinedLeeFilter2x2(int nLooks, int height, int width, const float *c11,
                        const float *c22, const float *c12r, const float *c12i,
                        float *outC11, float *outC22, float *outC12r,
                        float *outC12i) {
  try {
    cl::Context context = cl::Context::getDefault();
    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    cl::CommandQueue queue(context);

    // Compile the program
    cl::Program program(context, src_rlf_kernel);
    program.build(device, "-cl-std=CL2.0");

    // Prepare to execute the kernel
    auto total_len = height * width;
    cl::NDRange global_size(height, width);

    auto max_size = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    auto group_height = int_sqrt(max_size);
    cl::NDRange group_size(group_height, group_height);

    // Create buffers
    cl::Buffer buf_c11(context, c11, c11 + total_len, true);
    cl::Buffer buf_c22(context, c22, c22 + total_len, true);
    cl::Buffer buf_c12r(context, c12r, c12r + total_len, true);
    cl::Buffer buf_c12i(context, c12i, c12i + total_len, true);
    cl::Buffer buf_span(context, CL_MEM_READ_WRITE, sizeof(float) * total_len);

    // We use the first device to calculate SPAN
    cl::Kernel krnl_getspan(program, "mat_add2");
    krnl_getspan.setArg(0, buf_span);
    krnl_getspan.setArg(1, buf_c11);
    krnl_getspan.setArg(2, buf_c22);
    queue.enqueueNDRangeKernel(krnl_getspan, cl::NullRange, global_size);

    // Now we can start filting process
    cl::Kernel krnl_rlf(program, "filt_cij");
    krnl_rlf.setArg(0, buf_span);
    krnl_rlf.setArg(1, height);
    krnl_rlf.setArg(2, width);

    auto shared_height = static_cast<int>(group_height) + 8;
    krnl_rlf.setArg(3, sizeof(float) * shared_height * shared_height, nullptr);
    krnl_rlf.setArg(4, shared_height);
    krnl_rlf.setArg(5, shared_height);

    cl::Buffer buf_prwt(context, Prewitt.begin(), Prewitt.end(), true);
    krnl_rlf.setArg(6, buf_prwt);
    krnl_rlf.setArg(7, nLooks);

    // For every channel
    cl::Buffer buf_oc11(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c11);
    krnl_rlf.setArg(9, buf_oc11);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc11, false, 0, sizeof(float) * total_len,
                            outC11);

    cl::Buffer buf_oc22(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c22);
    krnl_rlf.setArg(9, buf_oc22);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc22, false, 0, sizeof(float) * total_len,
                            outC22);

    cl::Buffer buf_oc12r(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c12r);
    krnl_rlf.setArg(9, buf_oc12r);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc12r, false, 0, sizeof(float) * total_len,
                            outC12r);

    cl::Buffer buf_oc12i(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len);
    krnl_rlf.setArg(8, buf_c12i);
    krnl_rlf.setArg(9, buf_oc12i);
    queue.enqueueNDRangeKernel(krnl_rlf, cl::NullRange, global_size, group_size);
    queue.enqueueReadBuffer(buf_oc12i, false, 0, sizeof(float) * total_len,
                            outC12i);

    queue.finish();

  } catch (cl::Error &e) {
    return e.err();
  }
  return 0;
}