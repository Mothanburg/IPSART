#include "IPSART.h"
#include "ocl.hpp"

#include <cassert>
#include <cstdint>
#include <exception>
#include <format>
#include <ranges>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include <mexAdapter.hpp>

#ifndef NDEBUG // Windows 调试陷阱
#include <intrin.h>
#define DEBUG_BREAK() __debugbreak()
#else
#define DEBUG_BREAK() ((void)0)
#endif

using namespace std;
using namespace matlab;

// 方法处理器类型定义
using MethodHandler = vector<data::Array> (*)(const vector<data::Array> &,
                                              data::ArrayFactory &);

// 方法注册表
static const unordered_map<u16string, MethodHandler> method_registry = {
    {u"CloudePottier", ipsart::CloudePottier},
    {u"Multilook", ipsart::Multilook},
    {u"RefinedLeeFilter", ipsart::RefinedLeeFilter},
    {u"Yamaguchi", ipsart::Yamaguchi},
    {u"G4U", ipsart::G4U}};

class MexFunction : public matlab::mex::Function {
private:
  vector<data::Array> dispatch(const data::MATLABString &method,
                               const vector<data::Array> &input) {
    // 查找并调用对应的处理函数
    auto it = method_registry.find(*method);
    if (it != method_registry.end()) {
      return it->second(input, this->af);
    }

    // 未找到处理函数
    throw runtime_error(
        "Unknown method when calling native lib, IPSART maybe broken.");
  }

public:
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) final {
    // 在 DEBUG 模式下，强制触发断点，不然有可能因为跳过符号加载而无法调试
    DEBUG_BREAK();

    // 获取方法名称
    assert(inputs[0].getType() == data::ArrayType::MATLAB_STRING);
    const data::MATLABString method_name =
        data::TypedArray<data::MATLABString>(inputs[0])[0];
    assert(method_name.has_value());

    try {
      // 设置输入参数并分发调用
      vector<data::Array> params =
          inputs | views::drop(1) | ranges::to<vector>();
      auto results = dispatch(method_name, params);
      ranges::move(results, outputs.begin());
    } catch (exception &e) {
      // 捕获异常，并将错误信息传递回 MATLAB
      auto matlab = this->getEngine();
      string u8_method_name(from_range,
                            *method_name | views::transform([](char16_t c) {
                              return static_cast<char>(c);
                            }));
      string msg =
          format("Error occurred in {}.\nException type: {}\nMessage: {}\n",
                 u8_method_name, typeid(e).name(), e.what());
      // 如果是 OpenCL 错误，附加错误码和设备信息
      if (cl::Error *cl_e = dynamic_cast<cl::Error *>(&e)) {
        msg += format("OpenCL error code: {}\n", cl_e->err());
        msg += format("OpenCL device: {}, version: {}\n",
                      ipsart::ocl::getOpenCLDeviceName(),
                      ipsart::ocl::getOpenCLDeviceVersion());
      }
      matlab->feval(u"error", 0,
                    vector<data::Array>({this->af.createScalar(msg)}));
    }
#if !defined(NDEBUG)
    catch (...) {
      // 在 DEBUG 模式下，开启 /EHsa
      // 编译后，捕获底层异常（如访问违规、断言失败等）
      auto matlab = this->getEngine();
      string u8_method_name(from_range,
                            *method_name | views::transform([](char16_t c) {
                              return static_cast<char>(c);
                            }));
      string msg = format("Underlying error occurred in {}.\n", u8_method_name);
      matlab->feval(u"error", 0,
                    vector<data::Array>({this->af.createScalar(msg)}));
    }
#endif
  }

private:
  data::ArrayFactory af;
};