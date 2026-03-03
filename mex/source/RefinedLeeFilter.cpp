#include "IPSART.h"
#include "util.hpp"

#include <array>
#include <cstdint>
#include <ranges>
#include <span>
#include <type_traits>
#include <vector>

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

#include "autogen/RefinedLeeFilter.cl.h"

using namespace std;
using namespace matlab;

// clang-format off
// Prewitt operators (Col-major)
template<typename Float>
constexpr std::array<Float, 7 * 7 * 8> PrewittMask {
// Prewitt 1
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0,
1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
// Prewitt 2
// 1, 1, 1, 1, 1, 1, 1,
// 0, 1, 1, 1, 1, 1, 1,
// 0, 0, 1, 1, 1, 1, 1,
// 0, 0, 0, 1, 1, 1, 1,
// 0, 0, 0, 0, 1, 1, 1,
// 0, 0, 0, 0, 0, 1, 1,
// 0, 0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0, 0,
1, 1, 0, 0, 0, 0, 0,
1, 1, 1, 0, 0, 0, 0,
1, 1, 1, 1, 0, 0, 0,
1, 1, 1, 1, 1, 0, 0,
1, 1, 1, 1, 1, 1, 0,
1, 1, 1, 1, 1, 1, 1,
// Prewitt 3
// 1, 1, 1, 1, 1, 1, 1,
// 1, 1, 1, 1, 1, 1, 1,
// 1, 1, 1, 1, 1, 1, 1,
// 1, 1, 1, 1, 1, 1, 1,
// 0, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, 0, 0, 0,
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
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0,
// Prewitt 6
// 1, 0, 0, 0, 0, 0, 0,
// 1, 1, 0, 0, 0, 0, 0,
// 1, 1, 1, 0, 0, 0, 0,
// 1, 1, 1, 1, 0, 0, 0,
// 1, 1, 1, 1, 1, 0, 0,
// 1, 1, 1, 1, 1, 1, 0,
// 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1,
0, 1, 1, 1, 1, 1, 1,
0, 0, 1, 1, 1, 1, 1,
0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 0, 1, 1, 1,
0, 0, 0, 0, 0, 1, 1,
0, 0, 0, 0, 0, 0, 1,
// Prewitt 7
// 0, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, 0, 0, 0,
// 0, 0, 0, 0, 0, 0, 0,
// 1, 1, 1, 1, 1, 1, 1,
// 1, 1, 1, 1, 1, 1, 1,
// 1, 1, 1, 1, 1, 1, 1,
// 1, 1, 1, 1, 1, 1, 1,
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

template <typename Float, int Dim>
  requires(Dim > 1 && Dim < 5 && is_floating_point_v<Float>)
static void
refined_lee_filter(const vector<data::TypedArray<Float>> &in_elements,
                   int look_num, vector<vector<Float>> &out_elements) {
  assert(in_elements.size() == out_elements.size());
  // initialize OpenCL
  static cl::Context context = cl::Context::getDefault();
  static cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
  static cl::CommandQueue queue(context, device);

  // Compile kernel program
  static cl::Program program = [&]<typename Ty, int D>() {
    cl::Program p(context, SRC_REFINEDLEEFILTER);
    if constexpr (std::is_same_v<Ty, float>) {
      if constexpr (D == 3) {
        p.build(device, "-cl-std=CL2.0 -DMAT_SIZE_3X3");
      } else {
        p.build(device, "-cl-std=CL2.0");
      }
    } else {
      /* Float is double */
      if constexpr (D == 3) {
        p.build(device, "-cl-std=CL2.0 -DENABLE_FP64 -DMAT_SIZE_3X3");
      } else {
        p.build(device, "-cl-std=CL2.0 -DENABLE_FP64");
      }
    }
    return p;
  }.template operator()<Float, Dim>();
  static cl::Kernel krnl_span_calc(program, "span_calc");
  static cl::Kernel krnl_filter(program, "page_filting");

  // set running parameters
  auto dim = in_elements[0].getDimensions();
  int rows = dim[0];
  int cols = dim[1];
  int total_len = rows * cols;
  cl::NDRange global_size(rows, cols);

  size_t max_size = device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
  size_t group_height = std::floor(std::sqrt(static_cast<double>(max_size)));
  cl::NDRange group_size(group_height, group_height);

  // set in/out buffers - only need min inputs for span_calc + 1 reusable output
  constexpr int MinInputs = (Dim == 2) ? 2 : 3;
  constexpr int BufCount = MinInputs + 1;
  std::array<cl::Buffer, BufCount> buf_inputs;
  cl::Buffer buf_output;

  // Initialize input buffers for span_calc
  for (auto idx = 0; idx < MinInputs; idx++) {
    auto ref_element = marray_to_span(in_elements[idx]);
    buf_inputs[idx] =
        cl::Buffer(context, ref_element.begin(), ref_element.end(), true);
  }
  // Reusable output buffer
  buf_output = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(Float) * total_len);

  // calculate span
  cl::Buffer buf_span(context, CL_MEM_READ_WRITE, sizeof(Float) * total_len);
  krnl_span_calc.setArg(0, buf_span);
  krnl_span_calc.setArg(1, buf_inputs[0]);
  krnl_span_calc.setArg(2, buf_inputs[1]);
  if constexpr (Dim == 3) {
    krnl_span_calc.setArg(3, buf_inputs[2]);
  }
  queue.enqueueNDRangeKernel(krnl_span_calc, cl::NullRange, global_size,
                             group_size);
  // No explicit finish() needed - in-order queue guarantees span_calc completes
  // before filter kernel starts

  // refined lee filting
  krnl_filter.setArg(0, rows);
  krnl_filter.setArg(1, cols);
  krnl_filter.setArg(2, buf_span);
  int shm_size = 8 + group_height;
  krnl_filter.setArg(3, shm_size);
  krnl_filter.setArg(4, shm_size);
  krnl_filter.setArg(5, sizeof(Float) * shm_size * shm_size * 2, nullptr);
  cl::Buffer buf_prwt(context, PrewittMask<Float>.begin(),
                      PrewittMask<Float>.end(), true);
  krnl_filter.setArg(6, buf_prwt);
  krnl_filter.setArg(7, look_num);

  // For each input element: upload to buffer, run kernel, read result
  const auto n_elements = in_elements.size();
  for (auto idx = 0; idx < n_elements; idx++) {
    // Reuse input buffer at index 0 for subsequent elements (after span_calc is done)
    if (idx >= MinInputs) {
      auto ref_element = marray_to_span(in_elements[idx]);
      // Copy data into existing buffer
      queue.enqueueWriteBuffer(buf_inputs[0], true, 0,
                               sizeof(Float) * total_len,
                               ref_element.data());
    } else {
      // First MinInputs elements already in buffers from initialization
    }
    // Reuse single output buffer
    krnl_filter.setArg(8, buf_inputs[idx < MinInputs ? idx : 0]);
    krnl_filter.setArg(9, buf_output);
    queue.enqueueNDRangeKernel(krnl_filter, cl::NullRange, global_size,
                               group_size);
    queue.enqueueReadBuffer(buf_output, true, 0,
                            sizeof(Float) * total_len,
                            out_elements[idx].data());
  }
  queue.finish();
}

vector<data::Array> RefinedLeeFilter(const vector<data::Array> &input,
                                     data::ArrayFactory &af) {
  auto in_num = input.size();
  auto in_type = input[0].getType();
  auto in_dims = input[0].getDimensions();
  assert(input.size() == 4 + 1 || input.size() == 9 + 1);
  assert(ranges::all_of(input | views::take(input.size() - 1), [&](auto &arr) {
    return arr.getDimensions() == in_dims && arr.getType() == in_type;
  }));
  auto numel = in_dims[0] * in_dims[1];

  auto execute = [&]<typename T>() {
    const int n_elements = static_cast<int>(input.size()) - 1;
    const int look_num = input.back()[0];

    vector<vector<T>> out_elements(n_elements, vector<T>(numel));
    auto in_view = input | views::take(n_elements);

    if (n_elements == 4) { // Dim = 2
      auto in_arr = in_view | ranges::to<vector<data::TypedArray<T>>>();
      refined_lee_filter<T, 2>(in_arr, look_num, out_elements);
    } else { // Dim = 3
      auto in_arr = in_view | ranges::to<vector<data::TypedArray<T>>>();
      refined_lee_filter<T, 3>(in_arr, look_num, out_elements);
    }

    vector<data::Array> result;
    result.reserve(n_elements);
    for (size_t idx = 0; idx < n_elements; idx++) {
      result.push_back(af.createArray(in_dims, out_elements[idx].begin(),
                                      out_elements[idx].end()));
      vector<T>().swap(out_elements[idx]); // free the memory by advance
    }
    return result;
  };

  if (in_type == data::ArrayType::SINGLE) {
    return execute.template operator()<float>();
  } else {
    return execute.template operator()<double>();
  }
}