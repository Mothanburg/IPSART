#include <array>
#include <type_traits>

#include "GaRS.h"

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

#include "RefinedLeeFilter.cl.h"

using namespace std;

// clang-format off
// Prewitt operators (Row-major)
template <typename TData>
  requires is_same_v<TData, float> || is_same_v<TData, double>
constexpr array<TData, 7 * 7 * 8> Prewitt {
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

template <typename TData, int Dim>
static int refined_lee_filter(int lookNum, int rows, int cols,
                              const array<const TData *, Dim * Dim> &inputs,
                              const array<TData *, Dim * Dim> &outputs) {
  constexpr int MatSize = Dim * Dim;
  try {
    cl::Context context = cl::Context::getDefault();
    cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
    cl::CommandQueue queue(context);

    cl::Program program(context, SRC_REFINEDLEEFILTER);
    if constexpr (is_same_v<TData, float>) {
      if constexpr (Dim == 3) {
        program.build(device, "-cl-std=CL2.0 -DMAT_SIZE_3X3");
      } else {
        program.build(device, "-cl-std=CL2.0");
      }
    } else {
      /* TData is double */
      if constexpr (Dim == 3) {
        program.build(device, "-cl-std=CL2.0 -DENABLE_FP64 -DMAT_SIZE_3X3");
      } else {
        program.build(device, "-cl-std=CL2.0 -DENABLE_FP64");
      }
    }

    int total_len = rows * cols;
    cl::NDRange global_size(rows, cols);

    size_t max_size = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    auto group_height =
        static_cast<size_t>(floor(sqrt(static_cast<double>(max_size))));
    cl::NDRange group_size(group_height, group_height);

    array<cl::Buffer, Dim * Dim> buf_inputs, buf_outputs;
    for (auto idx = 0; idx < MatSize; idx++) {
      buf_inputs[idx] = std::move(
          cl::Buffer(context, inputs[idx], inputs[idx] + total_len, true));
      buf_outputs[idx] = std::move(
          cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(TData) * total_len));
    }

    cl::Buffer buf_span(context, CL_MEM_READ_WRITE, sizeof(TData) * total_len);
    cl::Kernel krnl_span_calc(program, "span_calc");
    krnl_span_calc.setArg(0, buf_span);
    krnl_span_calc.setArg(1, buf_inputs[0]);
    krnl_span_calc.setArg(2, buf_inputs[1]);
    if constexpr (Dim == 3) {
      krnl_span_calc.setArg(3, buf_inputs[2]);
    }
    queue.enqueueNDRangeKernel(krnl_span_calc, cl::NullRange, global_size,
                               group_size);

    cl::Kernel krnl_filter(program, "page_filting");
    krnl_filter.setArg(0, rows);
    krnl_filter.setArg(1, cols);
    krnl_filter.setArg(2, buf_span);

    auto shm_size = static_cast<int>(group_height) + 8;
    krnl_filter.setArg(3, shm_size);
    krnl_filter.setArg(4, shm_size);
    krnl_filter.setArg(5, sizeof(TData) * shm_size * shm_size, nullptr);

    cl::Buffer buf_prwt(context, Prewitt<TData>.begin(), Prewitt<TData>.end(),
                        true);
    krnl_filter.setArg(6, buf_prwt);
    krnl_filter.setArg(7, lookNum);

    for (auto idx = 0; idx < MatSize; idx++) {
      krnl_filter.setArg(8, buf_inputs[idx]);
      krnl_filter.setArg(9, buf_outputs[idx]);
      queue.enqueueNDRangeKernel(krnl_filter, cl::NullRange, global_size,
                                 group_size);
      queue.enqueueReadBuffer(buf_outputs[idx], false, 0,
                              sizeof(TData) * total_len, outputs[idx]);
    }

    queue.finish();
  } catch (cl::Error &e) {
    return e.err();
  }
  return 0;
}

int RefinedLeeFilting3f(int lookNum, int rows, int cols, const float *c11,
                        const float *c22, const float *c33, const float *c12r,
                        const float *c13r, const float *c23r, const float *c12i,
                        const float *c13i, const float *c23i, float *outC11,
                        float *outC22, float *outC33, float *outC12r,
                        float *outC13r, float *outC23r, float *outC12i,
                        float *outC13i, float *outC23i) {
  const array input{c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i};
  const array output{outC11,  outC22,  outC33,  outC12r, outC13r,
                     outC23r, outC12i, outC13i, outC23i};
  return refined_lee_filter<float, 3>(lookNum, rows, cols, input, output);
}

int RefinedLeeFilting3d(int lookNum, int rows, int cols, const double *c11,
                        const double *c22, const double *c33,
                        const double *c12r, const double *c13r,
                        const double *c23r, const double *c12i,
                        const double *c13i, const double *c23i, double *outC11,
                        double *outC22, double *outC33, double *outC12r,
                        double *outC13r, double *outC23r, double *outC12i,
                        double *outC13i, double *outC23i) {
  const array input{c11, c22, c33, c12r, c13r, c23r, c12i, c13i, c23i};
  const array output{outC11,  outC22,  outC33,  outC12r, outC13r,
                     outC23r, outC12i, outC13i, outC23i};
  return refined_lee_filter<double, 3>(lookNum, rows, cols, input, output);
}

int RefinedLeeFilting2f(int lookNum, int rows, int cols, const float *c11,
                        const float *c22, const float *c12r, const float *c12i,
                        float *outC11, float *outC22, float *outC12r,
                        float *outC12i) {
  const array input{c11, c22, c12r, c12i};
  const array output{outC11, outC22, outC12r, outC12i};
  return refined_lee_filter<float, 2>(lookNum, rows, cols, input, output);
}

int RefinedLeeFilting2d(int lookNum, int rows, int cols, const double *c11,
                        const double *c22, const double *c12r,
                        const double *c12i, double *outC11, double *outC22,
                        double *outC12r, double *outC12i) {
  const array input{c11, c22, c12r, c12i};
  const array output{outC11, outC22, outC12r, outC12i};
  return refined_lee_filter<double, 2>(lookNum, rows, cols, input, output);
}