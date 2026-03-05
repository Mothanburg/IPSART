#include "ocl.hpp"

#include <stdexcept>
#include <vector>

using namespace std;

namespace ipsart {
namespace ocl {

OpenCLManager &OpenCLManager::instance() {
  static OpenCLManager instance;
  return instance;
}

OpenCLManager::OpenCLManager()
    : device_(pickFastestDevice()), context_(device_),
      queue_(context_, device_) {
  // 默认创建顺序执行队列，若设备支持，则创建乱序队列
  auto queue_properties = device_.getInfo<CL_DEVICE_QUEUE_PROPERTIES>();
  if (queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE) {
    queue_ = cl::CommandQueue(context_, device_,
                              CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
  }
}

cl::Device OpenCLManager::pickFastestDevice() {
  vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);

  if (platforms.empty()) {
    throw runtime_error("No OpenCL platform found");
  }

  cl::Device best_device;

  // 寻找评分最高的设备
  double best_score = -1;
  for (auto &p : platforms) {
    vector<cl::Device> devices;
    p.getDevices(CL_DEVICE_TYPE_ALL, &devices);

    for (auto &d : devices) {
      auto type = d.getInfo<CL_DEVICE_TYPE>();
      double cu = d.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
      double mhz = d.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
      double memGB = d.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() / 1e9;

      // 评分机制：GPU > 计算单元数 * 频率 > 内存大小
      double score =
          (type & CL_DEVICE_TYPE_GPU ? 1e6 : 0) + cu * mhz + memGB * 1000;

      if (score > best_score) {
        best_score = score;
        best_device = d;
      }
    }
  }

  if (best_score < 0) {
    throw runtime_error("No OpenCL device found");
  }

  return best_device;
}

cl::Program OpenCLManager::getProgram(const string &function_name,
                                      const string &build_opts,
                                      const char *source) {
  string cache_key = function_name + build_opts;

  // 查询缓存
  auto it = program_cache_.find(cache_key);
  if (it != program_cache_.end()) {
    return it->second;
  }

  // 编译新程序
  cl::Program program(context_, source);
  try {
    program.build(device_, build_opts.c_str());
  } catch (const cl::Error &e) {
    string build_log = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_);
    throw runtime_error("OpenCL build failed: " + build_log);
  }

  // 加入缓存
  program_cache_.emplace(cache_key, std::move(program));
  return program_cache_.at(cache_key);
}

string OpenCLManager::getDeviceName() const {
  return device_.getInfo<CL_DEVICE_NAME>();
}

string OpenCLManager::getDeviceVersion() const {
  return device_.getInfo<CL_DEVICE_OPENCL_C_VERSION>();
}

} // namespace ocl
} // namespace ipsart
