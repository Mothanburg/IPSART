#include "GaRS.h"

int main() {
  constexpr long H = 25;
  constexpr long W = 25;
  constexpr unsigned long long LEN = static_cast<unsigned long long>(H * W);

  CloudePottierT3(H, W, new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN], new double[LEN]);
  return 0;
}