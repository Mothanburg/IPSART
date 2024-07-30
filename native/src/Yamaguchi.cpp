#include <cmath>
#include <complex>
#include <numbers>

#include "GaRS.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/Core>

using namespace Eigen;
using namespace std;
using namespace std::literals;
constexpr float sqrt2 = std::numbers::sqrt2_v<float>;

void Yamaguchi(int height, int width, const float *c11, const float *c22,
               const float *c33, const float *c12r, const float *c13r,
               const float *c23r, const float *c12i, const float *c13i,
               const float *c23i, float *outPs, float *outPd, float *outPv,
               float *outPh) {
  auto len = height * width;
# pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    Array33cf c;
    c << scomplex(c11[idx], 0.0f), scomplex(c12r[idx], c12i[idx]),
        scomplex(c13r[idx], c13i[idx]), scomplex(c12r[idx], -c12i[idx]),
        scomplex(c22[idx], 0.0f), scomplex(c23r[idx], c23i[idx]),
        scomplex(c13r[idx], -c13i[idx]), scomplex(c23r[idx], -c23i[idx]),
        scomplex(c33[idx], 0.0f);

    // The form of helix scattering
    Array33cf ch;
    if (c(0, 1).imag() + c(1, 2).imag() > 0) {
      ch << 1.0f, 1.0if * sqrt2, -1.0f, -1.0if * sqrt2, 2.0f, 1.0f * sqrt2,
          -1.0f, -1.0if * sqrt2, 1.0f;
      ch /= 4.0f;
    } else {
      ch << 1.0f, -1.0if * sqrt2, -1.0f, 1.0if * sqrt2, 2.0f, -1.0f * sqrt2,
          -1.0f, 1.0if * sqrt2, 1.0f;
      ch /= 4.0f;
    }
    float fh = 2.0f * abs(c(0, 1).imag() + c(1, 2).imag());

    // The form of volume scattering
    const float coratio = 10.0f * log10(c(2, 2).real() / c(0, 0).real());
    Array33cf cv;
    float fv;
    if (coratio < -2.0f) {
      cv << 8.0f, 0.0f, 2.0f, 0.0f, 4.0f, 0.0f, 2.0f, 0.0f, 3.0f;
      cv /= 15.0f;
      fv = 15.0f * (c(1, 1).real() - fh / 2.0f) / 4.0f;
    } else if (coratio < 2) {
      cv << 3.0f, 0.0f, 1.0f, 0.0f, 2.0f, 0.0f, 1.0f, 0.0f, 3.0f;
      cv /= 8.0f;
      fv = 4.0f * (c(1, 1).real() - fh / 2.0f);
    } else {
      cv << 3.0f, 0.0f, 2.0f, 0.0f, 4.0f, 0.0f, 2.0f, 0.0f, 8.0f;
      cv /= 15.0f;
      fv = 15.0f * (c(1, 1).real() - fh / 2.0f) / 4.0f;
    }

    // Remove volume and helix scattering only when cross-pol scattering is
    // small
    if (c(1, 1).real() < c(0, 0).real() && c(1, 1).real() < c(2, 2).real()) {
      c = c - fh * ch - fv * cv;
    } else {
      fh = 0.0f;
      fv = 0.0f;
    }

    float a2, b2, fd, fs;
    if (c(0, 2).real() > 0) {
      a2 = 1.0f;
      fd = ((c(2, 2) * c(0, 0) - c(0, 2) * c(2, 0)) /
            (c(2, 2) + c(0, 0) + c(0, 2) + c(2, 0)))
               .real();
      fs = c(2, 2).real() - fd;
      scomplex tmp = ((c(0, 2) + fd) / fs);
      b2 = abs(tmp * conj(tmp));
    } else {
      b2 = 1;
      fs = ((c(0, 2) * c(2, 0) - c(2, 2) * c(0, 0)) /
            (c(0, 2) + c(2, 0) - c(0, 0) - c(2, 2)))
               .real();
      fd = c(2, 2).real() - fs;
      scomplex tmp = ((c(0, 2) - fs) / fd);
      a2 = abs(tmp * conj(tmp));
    }

    outPs[idx] = abs(fs) * (1.0f + b2);
    outPd[idx] = abs(fd) * (1.0f + a2);
    outPh[idx] = abs(fh);
    outPv[idx] = abs(fv);
  }
}