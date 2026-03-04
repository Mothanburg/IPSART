#pragma once

#ifndef _IPSART_OCL_
#define _IPSART_OCL_

#include <memory>
#include <string>
#include <unordered_map>

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL/opencl.hpp>

namespace ipsart {
namespace ocl {

/**
 * @brief OpenCL 单例资源管理器
 * @details 统一管理设备、上下文、命令队列和程序编译
 */
class OpenCLManager {
public:
  /// 获取单例实例
  static OpenCLManager &instance();

  // 禁止拷贝和移动
  OpenCLManager(const OpenCLManager &) = delete;
  OpenCLManager &operator=(const OpenCLManager &) = delete;
  OpenCLManager(OpenCLManager &&) = delete;
  OpenCLManager &operator=(OpenCLManager &&) = delete;

  /// 获取 OpenCL 设备
  const cl::Device &device() const { return device_; }

  /// 获取 OpenCL 上下文
  cl::Context &context() { return context_; }

  /// 获取命令队列
  cl::CommandQueue &queue() { return queue_; }

  /**
   * @brief 获取或编译 OpenCL 程序
   * @param function_name 内核函数名称
   * @param build_opts 编译选项
   * @param source OpenCL kernel 源码
   * @return 编译后的 cl::Program
   */
  cl::Program getProgram(const std::string &function_name,
                         const std::string &build_opts, const char *source);

  /// 获取设备信息
  std::string getDeviceName() const;
  std::string getDeviceVersion() const;

private:
  OpenCLManager();

  /// 自动选择最快的 OpenCL 设备
  static cl::Device pickFastestDevice();

  cl::Device device_;
  cl::Context context_;
  cl::CommandQueue queue_;
  std::unordered_map<std::string, cl::Program> program_cache_;
};

/// 获取 OpenCL 设备名称
inline std::string getOpenCLDeviceName() {
  return OpenCLManager::instance().getDeviceName();
}

/// 获取 OpenCL 设备版本
inline std::string getOpenCLDeviceVersion() {
  return OpenCLManager::instance().getDeviceVersion();
}

} // namespace ocl
} // namespace ipsart

#endif
