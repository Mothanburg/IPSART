#include "GaRS.h"

constexpr int H = 537;
constexpr int W = 379;
constexpr int LEN = H * W;

static float input[LEN] = {};
static float output[LEN] = {};

int main() {
  int code = GaRSTestOpenCL();
  code |= RefinedLeeFilter3x3(1, H, W, input, input, input, input, input, input,
                              input, input, input, output, output, output,
                              output, output, output, output, output, output);
  CloudePottier(H, W, input, input, input, input, input, input, input, input,
                input, output, output, output);
  Yamaguchi(H, W, input, input, input, input, input, input, input, input, input,
            output, output, output, output);
  G4U(H, W, input, input, input, input, input, input, input, input, input,
      output, output, output, output);
  code |= Multilook(H, W, input, 4, 4, H % 4, W % 4, output);
  return code;
}