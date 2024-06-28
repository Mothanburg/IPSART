#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <span>
#include <vector>

#include "GaRS.h"
#include "utils.hpp"

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <cl/opencl.hpp>

#include "cl.rlf_kernel.h"

using namespace std;

// clang-format off

// Prewitt operators (Column Major)
constexpr array<float, 7ull * 7 * 8> prewitt_templates{
  // Prewitt 1
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  // Prewitt 2
  1, 0, 0, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0,
  1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 1, 0, 0,
  1, 1, 1, 1, 1, 1, 0,
  1, 1, 1, 1, 1, 1, 1,
  // Prewitt 3
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  // Prewitt 4
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 0,
  1, 1, 1, 1, 1, 0, 0,
  1, 1, 1, 1, 0, 0, 0,
  1, 1, 1, 0, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0,
  // Prewitt 5
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0,
  // Prewitt 6
  1, 1, 1, 1, 1, 1, 1,
  0, 1, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 0, 1, 1, 1,
  0, 0, 0, 0, 0, 1, 1,
  0, 0, 0, 0, 0, 0, 1,
  // Prewitt 7
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
  0, 0, 0, 1, 1, 1, 1,
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

static void rlf_impl(long height, long width, const PolMat3<const float> &c3,
                     const PolMat3<float> &out_c3) {
  cl::Platform plat = cl::Platform::getDefault();

  vector<cl::Device> devices;
  plat.getDevices(CL_DEVICE_TYPE_GPU, &devices);
  auto num_devices = devices.size();

  cl::Context context(devices);

  vector<cl::CommandQueue> queues(devices.size());

  for (const auto &device : devices) {
    queues.emplace_back(context, device,
                        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
  }

  // Compile the program
  cl::Program program(context, src_rlf_kernel);
  program.build(devices);

  // Prepare to execute the kernel
  auto total_len = height * width;
  vector<cl::NDRange> group_sizes(devices.size());
  for (const auto &device : devices) {
    auto max_size = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    auto size = int_sqrt(max_size);
    group_sizes.emplace_back(size, size, 1);
  }

  // We use the first device to calculate SPAN
  cl::Buffer buf_span(context, CL_MEM_READ_WRITE, sizeof(float) * total_len);
  cl::Buffer buf_c11(context, c3.m11.begin(), c3.m11.end(), true);
  cl::Buffer buf_c22(context, c3.m22.begin(), c3.m22.end(), true);
  cl::Buffer buf_c33(context, c3.m33.begin(), c3.m33.end(), true);

  cl::Kernel krnl_getspan(program, "get_span");
  krnl_getspan.setArg(0, buf_c11);
  krnl_getspan.setArg(1, buf_c22);
  krnl_getspan.setArg(2, buf_c33);
  krnl_getspan.setArg(3, buf_span);

  queues[0].enqueueNDRangeKernel(krnl_getspan, cl::NullRange, group_sizes[0]);
  queues[0].finish();

  // Now we can start filting process
  cl::Kernel krnl_rlf(program, "filt_cij");
  cl::Buffer buf_prwt(context, prewitt_templates.begin(),
                      prewitt_templates.end(), true);

  // Filt c11, c22, and c33 firstly -- reuse buffer
  for (auto i = 0; i < 3; i++) {
    auto queue_idx = num_devices % 3;
  }

  // cl::Device device = cl::Device::getDefault(&err_code);

  // cl::Context context(device, nullptr, nullptr, nullptr, &err_code);

  // cl::CommandQueue queue(context, device, cl::QueueProperties::None,
  // &err_code);

  // size_t total_len = static_cast<size_t>(height) * width;

  // float *span = new float[sizeof(float) * total_len];

  // create input buffer
  // cl::Buffer c11_buf(context, c11, c11 + total_len, true, false, &err_code);
  // cl::Buffer c22_buf(context, c22, c22 + total_len, true, false, &err_code);
  // cl::Buffer c33_buf(context, c33, c33 + total_len, true, false, &err_code);
  // cl::Buffer c12r_buf(context, c12_r, c12_r + total_len, true, false,
  //                     &err_code);
  // cl::Buffer c13r_buf(context, c13_r, c13_r + total_len, true, false,
  //                     &err_code);
  // cl::Buffer c23r_buf(context, c23_r, c23_r + total_len, true, false,
  //                     &err_code);
  // cl::Buffer c12i_buf(context, c12_i, c12_i + total_len, true, false,
  //                     &err_code);
  // cl::Buffer c13i_buf(context, c13_i, c13_i + total_len, true, false,
  //                     &err_code);
  // cl::Buffer c23i_buf(context, c23_i, c23_i + total_len, true, false,
  //                     &err_code);

  // // create output buffer
  // cl::Buffer out_c11_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                        outC11, &err_code);
  // cl::Buffer out_c22_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                        outC22, &err_code);
  // cl::Buffer out_c33_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                        outC33, &err_code);
  // cl::Buffer out_c12r_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                         outC12_r, &err_code);
  // cl::Buffer out_c13r_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                         outC13_r, &err_code);
  // cl::Buffer out_c23r_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                         outC23_r, &err_code);
  // cl::Buffer out_c12i_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                         outC12_i, &err_code);
  // cl::Buffer out_c13i_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                         outC13_i, &err_code);
  // cl::Buffer out_c23i_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) *
  // total_len,
  //                         outC23_i, &err_code);

  // cl::Program program(context, src_rlf_kernel, true, &err_code);

  // cl::Kernel kernel(program, "rlf_filt_part", &err_code);
}

void RefinedLeeFilter(long height, long width, const float *c11,
                      const float *c22, const float *c33, const float *c12_r,
                      const float *c13_r, const float *c23_r,
                      const float *c12_i, const float *c13_i,
                      const float *c23_i, float *outC11, float *outC22,
                      float *outC33, float *outC12_r, float *outC13_r,
                      float *outC23_r, float *outC12_i, float *outC13_i,
                      float *outC23_i) {
  long long total_len = static_cast<long long>(height) * width;
  PolMat3<const float> c3(total_len, c11, c22, c33, c12_r, c13_r, c23_r, c12_i,
                          c13_i, c23_i);
  PolMat3<float> out_c3(total_len, outC11, outC22, outC33, outC12_r, outC13_r,
                        outC23_r, outC12_i, outC13_i, outC23_i);
  rlf_impl(height, width, c3, out_c3);
}