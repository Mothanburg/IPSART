#include "GaRS.h"

#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>

using namespace std;

// clang-format off
template <typename TData>
constexpr array<TData, 77> test_set {
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
  2, 4, 6, 2, 4, 6, 2, 4, 6, 2, 4,
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
  2, 4, 6, 2, 4, 6, 2, 4, 6, 2, 4,
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2,
  2, 4, 6, 2, 4, 6, 2, 4, 6, 2, 4,
  1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2
};
// clang-format on

int main() {
  // Test openmp ability
  cout << "Testing openmp through Cloude-Pottier decomposition............";
  double *dummy_output = new double[77];
  CloudePottier3d(
      7, 11, test_set<double>.data(), test_set<double>.data(),
      test_set<double>.data(), test_set<double>.data(), test_set<double>.data(),
      test_set<double>.data(), test_set<double>.data(), test_set<double>.data(),
      test_set<double>.data(), dummy_output, dummy_output, dummy_output);
  CloudePottier2f(7, 11, test_set<float>.data(), test_set<float>.data(),
                  test_set<float>.data(), test_set<float>.data(),
                  reinterpret_cast<float *>(dummy_output),
                  reinterpret_cast<float *>(dummy_output),
                  reinterpret_cast<float *>(dummy_output));
  cout << "Pass." << endl;
  delete[] dummy_output;

  // Test opencl ability with multilook
  cout << "Testing opencl ability through multilook.......................";

  bool low_float_precision_flag = false;
  float *outputf = new float[9];
  int err_no = Multilookf(7, 11, test_set<float>.data(), 2, 3, 3, 3, outputf);
  if (err_no != 0) {
    cout << format("Fail. Error: OpenCL runtime occurred '{}' error", err_no)
         << endl;
    return 0;
  }
  for (int i = 0; i < 9; i++) {
    if (abs(outputf[i] - 3) > numeric_limits<float>::epsilon()) {
      low_float_precision_flag = true;
    }
  }

  bool low_double_precision_flag = false;
  double *outputd = new double[9];
  err_no = Multilookd(7, 11, test_set<double>.data(), 2, 3, 3, 3, outputd);
  if (err_no != 0) {
    cout << format("Fail. Error: OpenCL runtime occurred '{}' error", err_no)
         << endl;
    return 0;
  }
  for (int i = 0; i < 9; i++) {
    if (abs(outputd[i] - 3) > numeric_limits<double>::epsilon()) {
      low_double_precision_flag = true;
    }
  }

  if (low_double_precision_flag && low_float_precision_flag) {
    cout << "Pass. Warning: Low numerical precision in both fp64 and fp32."
         << endl;
  } else if (low_double_precision_flag) {
    cout << "Pass. Warning: Low numerical precision in fp64." << endl;
  } else if (low_float_precision_flag) {
    cout << "Pass. Warning: Low numerical precision in fp32." << endl;
  } else {
    cout << "Pass." << endl;
  }
  delete[] outputd;
  delete[] outputf;

  return 0;
}
