#include "IPSART.h"
#include "ocl.hpp"
#include "util.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <ranges>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

#include "autogen/RefinedLeeFilter.cl.h"

using namespace std;
using namespace matlab;

/// Prewitt 边缘检测算子 (8个方向, 7x7 窗口, 列优先存储)
/// 用于 Refined Lee 滤波器的边缘检测
template <typename Float>
inline constexpr std::array<Float, 7 * 7 * 8> PrewittMask{
    // clang-format off
// Prewitt 1 (水平向右)
0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0,
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1,
// Prewitt 2 (对角线 右下)
1, 0, 0, 0, 0, 0, 0, 
1, 1, 0, 0, 0, 0, 0,
1, 1, 1, 0, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 1, 0, 0, 
1, 1, 1, 1, 1, 1, 0, 
1, 1, 1, 1, 1, 1, 1,
// Prewitt 3 (垂直向下)
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0,
// Prewitt 4 (对角线 左下)
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 0, 
1, 1, 1, 1, 1, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 0, 0, 0, 0, 
1, 1, 0, 0, 0, 0, 0, 
1, 0, 0, 0, 0, 0, 0,
// Prewitt 5 (水平向左)
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 
0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0,
// Prewitt 6 (对角线 左上)
1, 1, 1, 1, 1, 1, 1, 
0, 1, 1, 1, 1, 1, 1,
0, 0, 1, 1, 1, 1, 1, 
0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 0, 1, 1, 1, 
0, 0, 0, 0, 0, 1, 1, 
0, 0, 0, 0, 0, 0, 1,
// Prewitt 7 (垂直向上)
0, 0, 0, 1, 1, 1, 1, 
0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 1, 1, 1, 1, 
0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 1, 1, 1, 1,
0, 0, 0, 1, 1, 1, 1,
// Prewitt 8 (对角线 右上)
0, 0, 0, 0, 0, 0, 1, 
0, 0, 0, 0, 0, 1, 1, 
0, 0, 0, 0, 1, 1, 1, 
0, 0, 0, 1, 1, 1, 1,
0, 0, 1, 1, 1, 1, 1, 
0, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1
//clang-format on
};


template <typename Float, int Dim>
  requires(Dim > 1 && Dim < 5 && is_floating_point_v<Float>)
static void
refined_lee_filter(const vector<data::TypedArray<Float>> &in_elements,
                   int look_num, vector<vector<Float>> &out_elements) {
  assert(in_elements.size() == out_elements.size());

  // 使用 OpenCLManager 获取资源
  auto &ocl = ipsart::ocl::OpenCLManager::instance();

  // 编译程序
  string build_opts = []<typename Ty, int D>() -> string {
    if constexpr (std::is_same_v<Ty, float>) {
      if constexpr (D == 3) {
        return "-cl-std=CL2.0 -DMAT_SIZE_3X3";
      } else {
        return "-cl-std=CL2.0";
      }
    } else {
      if constexpr (D == 3) {
        return "-cl-std=CL2.0 -DENABLE_FP64 -DMAT_SIZE_3X3";
      } else {
        return "-cl-std=CL2.0 -DENABLE_FP64";
      }
    }
  }.template operator()<Float, Dim>();
  cl::Program program =
      ocl.getProgram("RefinedLeeFilter", build_opts, SRC_REFINEDLEEFILTER);

  // 创建内核
  cl::Kernel krnl_span_calc(program, "span_calc");
  cl::Kernel krnl_filter(program, "page_filting");

  // 创建 Prewitt 算子的 buffer
  cl::Buffer buf_prwt(ocl.context(), PrewittMask<Float>.begin(),
                      PrewittMask<Float>.end(), true);

  // 设置工作组参数
  auto dim = in_elements[0].getDimensions();
  int rows = dim[0];
  int cols = dim[1];
  int total_len = rows * cols;
  cl::NDRange global_size(rows, cols);

  // 查询每个内核支持的最大工作组大小，取较小值以确保两者都兼容
  // 使用 CL_KERNEL_WORK_GROUP_SIZE（内核级）而非 CL_DEVICE_MAX_WORK_GROUP_SIZE（设备级）
  // 因为内核编译后的实际限制可能小于设备限制（如由于寄存器压力、共享内存占用等）
  size_t max_wg_span =
      krnl_span_calc.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(ocl.device());
  size_t max_wg_filter =
      krnl_filter.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(ocl.device());
  size_t max_wg_size = std::min(max_wg_span, max_wg_filter);

  // 取平方根得到 2D 工作组维度，确保不超过内核支持的最大值
  size_t group_height = static_cast<size_t>(std::floor(std::sqrt(static_cast<double>(max_wg_size))));
  if (group_height < 1)
    group_height = 1;
  cl::NDRange group_size(group_height, group_height);

  // 设置计算 span 时的输入 buffer
  // 在后面滤波时，我们复用计算 span 时的数据 buffer 填充矩阵元素
  array<cl::Buffer, Dim> buf_inputs;
  for (auto group = 0; group < Dim; group++) {
    auto ref_element = ipsart::marray_to_span(in_elements[group]);
    buf_inputs[group] =
        cl::Buffer(ocl.context(), ref_element.begin(), ref_element.end(), true);
  }

  // 计算 span
  cl::Buffer buf_span(ocl.context(), CL_MEM_READ_WRITE,
                      sizeof(Float) * total_len);
  krnl_span_calc.setArg(0, buf_span);
  krnl_span_calc.setArg(1, buf_inputs[0]);
  krnl_span_calc.setArg(2, buf_inputs[1]);
  if constexpr (Dim == 3) {
    krnl_span_calc.setArg(3, buf_inputs[2]);
  }
  cl::Event ev_span_calc;
  ocl.queue().enqueueNDRangeKernel(krnl_span_calc, cl::NullRange, global_size,
                                   group_size, nullptr, &ev_span_calc);

  // 准备滤波内核
  krnl_filter.setArg(0, rows);
  krnl_filter.setArg(1, cols);
  krnl_filter.setArg(2, buf_span);
  int shm_size = 8 + group_height;
  krnl_filter.setArg(3, shm_size);
  krnl_filter.setArg(4, shm_size);
  krnl_filter.setArg(5, sizeof(Float) * shm_size * shm_size * 2, nullptr);
  krnl_filter.setArg(6, buf_prwt);
  krnl_filter.setArg(7, look_num);

  // 设置存放滤波结果的 buffer，与之前创建的输入 buffer 对应
  array<cl::Buffer, Dim> buf_outputs;
  for (auto group = 0; group < Dim; group++) {
    buf_outputs[group] =
        cl::Buffer(ocl.context(), CL_MEM_WRITE_ONLY, sizeof(Float) * total_len);
  }

  // 开始滤波，分为 Dim 组，每组串行滤波 Dim 次
  for (auto group = 0; group < Dim; group++) {
    vector<cl::Event> filter_wait_list{ev_span_calc};

    krnl_filter.setArg(8, buf_inputs[group]);
    krnl_filter.setArg(9, buf_outputs[group]);

    // 推入内核、读取结果、读入待滤波元素
    for (auto progress = 0; progress < Dim; progress++) {
      cl::Event ev_filter;
      ocl.queue().enqueueNDRangeKernel(krnl_filter, cl::NullRange, global_size,
                                       group_size, &filter_wait_list,
                                       &ev_filter);

      vector<cl::Event> rw_wait_list = {ev_filter};
      cl::Event ev_read_buf;
      ocl.queue().enqueueReadBuffer(buf_outputs[group], false, 0,
                                    sizeof(Float) * total_len,
                                    out_elements[progress * Dim + group].data(),
                                    &rw_wait_list, &ev_read_buf);
      filter_wait_list = {ev_read_buf};
      if (progress < Dim - 1) {
        cl::Event ev_write_buf;
        auto ref_element =
            ipsart::marray_to_span(in_elements[(progress + 1) * Dim + group]);
        ocl.queue().enqueueWriteBuffer(
            buf_inputs[group], false, 0, sizeof(Float) * total_len,
            ref_element.data(), &rw_wait_list, &ev_write_buf);
        filter_wait_list.push_back(ev_write_buf);
      }
    }
  }

  ocl.queue().finish();
}

namespace ipsart{

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

} // namespace ipsart