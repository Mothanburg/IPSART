#include "IPSART.h"
#include "util.hpp"

#include <cmath>
#include <cstdint>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

#include "autogen/Multilook.cl.h"

using namespace std;
using namespace matlab;

template <typename Number> using TResult = tuple<vector<Number>, int, int>;

template <typename Number>
  requires std::is_same_v<Number, int32_t> || std::is_floating_point_v<Number>
static TResult<Number> multilook(const data::TypedArray<Number> &input,
                                 int row_look, int col_look) {
  // initialize OpenCL
  static cl::Context context = cl::Context::getDefault();
  static cl::Device device = context.getInfo<CL_CONTEXT_DEVICES>()[0];
  static cl::CommandQueue queue(context, device);

  // Compile the program
  static cl::Program program = [&]<typename Ty>() {
    cl::Program p(context, SRC_MULTILOOK);
    if constexpr (std::is_same_v<Ty, double>) {
      p.build(device, "-cl-std=CL2.0 -DENABLE_FP64");
    } else {
      p.build(device, "-cl-std=CL2.0");
    }
    return p;
  }.template operator()<Number>();

  // prepare kernel
  static cl::Kernel krnl(program, "multilook");

  // Calculate paramout dimensions
  auto in_dim = input.getDimensions();
  int in_rows = in_dim[0];
  int in_cols = in_dim[1];
  int out_rows = in_rows / row_look;
  int out_cols = in_cols / col_look;

  // create buffers
  auto ref_input = marray_to_span(input);
  cl::Buffer buf_input(context, ref_input.begin(), ref_input.end(), true);

  int out_len = out_rows * out_cols;
  vector<Number> result(out_len);
  cl::Buffer buf_output(context, CL_MEM_WRITE_ONLY, sizeof(Number) * out_len);

  // Prepare to execute the kernel
  cl::NDRange global_size(out_rows, out_cols);
  krnl.setArg(0, in_rows);
  krnl.setArg(1, buf_input);
  krnl.setArg(2, out_rows);
  krnl.setArg(3, buf_output);
  krnl.setArg(4, row_look);
  krnl.setArg(5, col_look);

  queue.enqueueNDRangeKernel(krnl, cl::NullRange, global_size);
  queue.enqueueReadBuffer(buf_output, false, 0, sizeof(Number) * out_len,
                          result.data());
  queue.finish();

  return {result, out_rows, out_cols};
}

vector<data::Array> Multilook(const vector<data::Array> &input,
                              data::ArrayFactory &af) {
  auto in_type = input[0].getType();
  assert(input.size() == 3);
  assert(in_type == data::ArrayType::INT32 ||
         in_type == data::ArrayType::SINGLE ||
         in_type == data::ArrayType::DOUBLE);
  assert(input[1].getNumberOfElements() == 1 &&
         input[1].getType() == data::ArrayType::INT32); // row look
  assert(input[2].getNumberOfElements() == 1 &&
         input[2].getType() == data::ArrayType::INT32); // col look
  const int32_t row_look = data::TypedArray<int32_t>(input[1])[0];
  const int32_t col_look = data::TypedArray<int32_t>(input[2])[0];

  auto execute = [&]<typename Ty> {
    data::TypedArray<Ty> in = input[0];
    auto [result, out_rows, out_cols] = multilook<Ty>(in, row_look, col_look);
    vector dim{static_cast<size_t>(out_rows), static_cast<size_t>(out_cols)};
    return vector<data::Array>{
        af.createArray(dim, result.begin(), result.end())};
  };

  if (in_type == data::ArrayType::INT32) {
    return execute.template operator()<int32_t>();
  } else if (in_type == data::ArrayType::SINGLE) {
    return execute.template operator()<float>();
  } else {
    return execute.template operator()<double>();
  }
}
