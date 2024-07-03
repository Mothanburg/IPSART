#include <cmath>
#include <complex>
#include <numbers>

#include "GaRS.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/EigenValues>

#include "omp.h"

using namespace Eigen;
using namespace std::literals;
constexpr float sqrt2 = std::numbers::sqrt2_v<float>;

using complexf = std::complex<float>;

void Yamaguchi(long height, long width, const float *c11, const float *c22,
               const float *c33, const float *c12_r, const float *c13_r,
               const float *c23_r, const float *c12_i, const float *c13_i,
               const float *c23_i, float *outPs, float *outPd, float *outPv,
               float *outPh) {
  #pragma omp parallel for
  for (int row = 0; row < height; row++) {
    for (int col = 0; col < width; col++) {
      int idx = row * width + col;
      Array33cf c;
      c(0, 0) = complexf(c11[idx], 0.0);
      c(0, 1) = complexf(c12_r[idx], c12_i[idx]);
      c(0, 2) = complexf(c13_r[idx], c13_i[idx]);
      c(1, 0) = std::conj(c(0, 1));
      c(1, 1) = complexf(c22[idx], 0.0);
      c(1, 2) = complexf(c23_r[idx], c23_i[idx]);
      c(2, 0) = std::conj(c(0, 2));
      c(2, 1) = std::conj(c(1, 2));
      c(2, 2) = complexf(c33[idx], 0.0);

      // The form of helix scattering
      Array33cf ch;
      if (c(0, 1).imag() + c(1, 2).imag() > 0) {
        ch << 1.0, 1.0if * sqrt2, -1.0, -1.0if * sqrt2, 2.0, 1.0 * sqrt2, -1.0,
            -1.0if * sqrt2, 1.0;
        ch /= 4.0;
      } else {
        ch << 1.0, -1.0if * sqrt2, -1.0, 1.0if * sqrt2, 2.0, -1.0 * sqrt2, -1.0,
            1.0if * sqrt2, 1.0;
        ch /= 4.0;
      }
      float fh = 2.0f * std::abs(c(0, 1).imag() + c(1, 2).imag());

      // The form of volume scattering
      const float co_ratio =
          10.0f * std::log10(c(2, 2).real() / c(0, 0).real());
      Array33cf cv;
      float fv;
      if (co_ratio < -2.0) {
        cv << 8.0, 0.0, 2.0, 0.0, 4.0, 0.0, 2.0, 0.0, 3.0;
        cv /= 15.0;
        fv = 15.0f * (c(1, 1).real() - fh / 2.0f) / 4.0f;
      } else if (co_ratio < 2) {
        cv << 3.0, 0.0, 1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 3.0;
        cv /= 8.0;
        fv = 4.0f * (c(1, 1).real() - fh / 2.0f);
      } else {
        cv << 3.0, 0.0, 2.0, 0.0, 4.0, 0.0, 2.0, 0.0, 8.0;
        cv /= 15.0;
        fv = 15.0f * (c(1, 1).real() - fh / 2.0f) / 4.0f;
      }

      // Remove volume and helix scattering only when cross-pol scattering is
      // small
      if (c(1, 1).real() < c(0, 0).real() && c(1, 1).real() < c(2, 2).real()) {
        c = c - fh * ch - fv * cv;
      } else {
        fh = 0;
        fv = 0;
      }

      float a2, b2, fd, fs;
      if (c(0, 2).real() > 0) {
        a2 = 1;
        fd = ((c(2, 2) * c(0, 0) - c(0, 2) * c(2, 0)) /
              (c(2, 2) + c(0, 0) + c(0, 2) + c(2, 0)))
                 .real();
        fs = c(2, 2).real() - fd;
        complexf tmp = ((c(0, 2) + fd) / fs);
        b2 = std::abs(tmp * std::conj(tmp));
      } else {
        b2 = 1;
        fs = ((c(0, 2) * c(2, 0) - c(2, 2) * c(0, 0)) /
              (c(0, 2) + c(2, 0) - c(0, 0) - c(2, 2)))
                 .real();
        fd = c(2, 2).real() - fs;
        complexf tmp = ((c(0, 2) - fs) / fd);
        a2 = std::abs(tmp * std::conj(tmp));
      }

      outPs[idx] = std::abs(fs) * (1.0f + b2);
      outPd[idx] = std::abs(fd) * (1.0f + a2);
      outPh[idx] = std::abs(fh);
      outPv[idx] = std::abs(fv);
    }
  }
}