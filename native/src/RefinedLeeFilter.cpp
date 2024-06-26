#include <algorithm>
#include <iostream>

#include "GaRS.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/Core>

#define CL_HPP_TARGET_OPENCL_VERSION 200
#include <cl/opencl.hpp>

#include "cl.rlf_kernel.h"

using Arrayf = Eigen::Array<float, 1, Eigen::Dynamic>;

void RefinedLeeFilter(long height, long width, const double *c11,
                      const double *c22, const double *c33, const double *c12_r,
                      const double *c13_r, const double *c23_r,
                      const double *c12_i, const double *c13_i,
                      const double *c23_i, double *outC11, double *outC22,
                      double *outC33, double *outC12_r, double *outC13_r,
                      double *outC23_r, double *outC12_i, double *outC13_i,
                      double *outC23_i) {
  cl_int err_code;

  cl::Device device = cl::Device::getDefault(&err_code);

  cl::Context context(device, nullptr, nullptr, nullptr, &err_code);

  cl::CommandQueue queue(context, device, cl::QueueProperties::None, &err_code);

  size_t total_len = static_cast<size_t>(height) * width;

  float *span = new float[sizeof(float) * total_len];
  Eigen::Map<Arrayf> map(span, total_len);

  // create input buffer
  cl::Buffer c11_buf(context, c11, c11 + total_len, true, false, &err_code);
  cl::Buffer c22_buf(context, c22, c22 + total_len, true, false, &err_code);
  cl::Buffer c33_buf(context, c33, c33 + total_len, true, false, &err_code);
  cl::Buffer c12r_buf(context, c12_r, c12_r + total_len, true, false,
                      &err_code);
  cl::Buffer c13r_buf(context, c13_r, c13_r + total_len, true, false,
                      &err_code);
  cl::Buffer c23r_buf(context, c23_r, c23_r + total_len, true, false,
                      &err_code);
  cl::Buffer c12i_buf(context, c12_i, c12_i + total_len, true, false,
                      &err_code);
  cl::Buffer c13i_buf(context, c13_i, c13_i + total_len, true, false,
                      &err_code);
  cl::Buffer c23i_buf(context, c23_i, c23_i + total_len, true, false,
                      &err_code);

  // create output buffer
  cl::Buffer out_c11_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                         outC11, &err_code);
  cl::Buffer out_c22_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                         outC22, &err_code);
  cl::Buffer out_c33_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                         outC33, &err_code);
  cl::Buffer out_c12r_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                          outC12_r, &err_code);
  cl::Buffer out_c13r_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                          outC13_r, &err_code);
  cl::Buffer out_c23r_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                          outC23_r, &err_code);
  cl::Buffer out_c12i_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                          outC12_i, &err_code);
  cl::Buffer out_c13i_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                          outC13_i, &err_code);
  cl::Buffer out_c23i_buf(context, CL_MEM_WRITE_ONLY, sizeof(float) * total_len,
                          outC23_i, &err_code);

  cl::Program program(context, src_rlf_kernel, true, &err_code);

  cl::Kernel kernel(program, "rlf_filt_part", &err_code);
}
