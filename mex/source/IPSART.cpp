#include "IPSART.h"

#include <cassert>
#include <cstdint>
#include <exception>
#include <format>
#include <ranges>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include <mexAdapter.hpp>

#ifndef NDEBUG // Windows debug trap
#include <intrin.h>
#define DEBUG_BREAK() __debugbreak()
#else
#define DEBUG_BREAK() ((void)0)
#endif

using namespace std;
using namespace matlab;

class MexFunction : public matlab::mex::Function {
private:
  vector<data::Array> dispatch(const data::MATLABString &method,
                               const vector<data::Array> &input) {
    if (*method == u"CloudePottier") {
      return CloudePottier(input, this->af);
    } else if (*method == u"Multilook") {
      return Multilook(input, this->af);
    } else if (*method == u"RefinedLeeFilter") {
      return RefinedLeeFilter(input, this->af);
    } else if (*method == u"Yamaguchi") {
      return Yamaguchi(input, this->af);
    } else if (*method == u"G4U") {
      return G4U(input, this->af);
    } else {
      throw runtime_error("Incorrect method name when calling native lib, "
                          "IPSART maybe broken.");
    }
  }

public:
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) final {
    DEBUG_BREAK();
    // get method name
    assert(inputs[0].getType() == data::ArrayType::MATLAB_STRING);
    const data::MATLABString method_name =
        data::TypedArray<data::MATLABString>(inputs[0])[0];
    assert(method_name.has_value());    

    try {
      // enter input arguments
      vector<data::Array> params =
          inputs | views::drop(1) | ranges::to<vector>();
      auto results = dispatch(method_name, params);
      ranges::move(results, outputs.begin());
    } catch (exception &e) {
      auto matlab = this->getEngine();
      string u8_method_name(from_range,
                            *method_name | views::transform([](char16_t c) {
                              return static_cast<char>(c);
                            }));
      string msg =
          format("Error occurred in {}.\n  Exception type: {}\n  Message: {}\n",
                 u8_method_name, typeid(e).name(), e.what());
      matlab->feval(u"error", 0,
                    vector<data::Array>({this->af.createScalar(msg)}));
    }
#if !defined(NDEBUG)
    catch (...) {
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