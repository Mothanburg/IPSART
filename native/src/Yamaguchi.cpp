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
using std::numbers::sqrt2;

using complexd = std::complex<double>;

GARS_API void Yamaguchi(long height, long width, const double *c11,
                        const double *c22, const double *c33,
                        const double *c12_r, const double *c13_r,
                        const double *c23_r, const double *c12_i,
                        const double *c13_i, const double *c23_i, double *outPs,
                        double *outPd, double *outPv, double *outPh) {
#pragma omp parallel for
  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      int idx = j * width + i;
      Array3cd c;
      c(0, 0) = complexd(c11[idx], 0.0);
      c(0, 1) = complexd(c12_r[idx], c12_i[idx]);
      c(0, 2) = complexd(c13_r[idx], c13_i[idx]);
      c(1, 0) = std::conj(c(0, 1));
      c(1, 1) = complexd(c22[idx], 0.0);
      c(1, 2) = complexd(c23_r[idx], c23_i[idx]);
      c(2, 0) = std::conj(c(0, 2));
      c(2, 1) = std::conj(c(1, 2));
      c(2, 2) = complexd(c33[idx], 0.0);

      // The form of helix scattering;
      Array3cd ch{c(0, 1).imag() + c(1, 2).imag() > 0
                      ? Array3cd(1.0, 1.0i * sqrt2, -1.0, -1.0i * sqrt2, 2.0,
                                 1.0 * sqrt2, -1.0, -1.0i * sqrt2, 1.0) / 4.0
                      : Array3cd(1.0, -1.0i * sqrt2, -1.0, 1.0i * sqrt2, 2.0,
                                 -1.0 * sqrt2, -1.0, 1.0i * sqrt2, 1.0) / 4.0};
      double fh = 2 * std::abs(c(0, 1).imag() + c(1, 2).imag());

      // The form of volume scattering
      const double co_ratio =
          10.0 * std::log10(c(2, 2).real() / c(0, 0).real());
      Array3cd cv;
      double fv;
      if (co_ratio < -2.0) {
        cv = Array3cd(8.0, 0.0, 2.0, 0.0, 4.0, 0.0, 2.0, 0.0, 3.0) / 15.0;
        fv = 15.0 * (c(1, 1).real() - fh / 2.0) / 4.0;
      } else if (co_ratio < 2) {
        cv = Array3cd(3.0, 0.0, 1.0, 0.0, 2.0, 0.0, 1.0, 0.0, 3.0) / 8.0;
        fv = 4.0 * (c(1, 1).real() - fh / 2.0);
      } else {
        cv = Array3cd(3.0, 0.0, 2.0, 0.0, 4.0, 0.0, 2.0, 0.0, 8.0) / 15.0;
        fv = 15.0 * (c(1, 1).real() - fh / 2.0) / 4.0;
      }

      // Remove volume and helix scattering only when cross-pol scattering is
      // small
      if (c(1, 1).real() < c(0, 0).real() && c(1, 1).real() < c(2, 2).real()) {
        c = c - fh * ch - fv * cv;
      } else {
        fh = 0;
        fv = 0;
      }

      double a2, b2, fd, fs;
      if (c(0, 2).real() > 0) {
        a2 = 1;
        fd = ((c(2, 2) * c(0, 0) - c(0, 2) * c(2, 0)) /
              (c(2, 2) + c(0, 0) + c(0, 2) + c(2, 0)))
                 .real();
        fs = c(2, 2).real() - fd;
        complexd tmp = ((c(0, 2) + fd) / fs);
        b2 = std::abs(tmp * std::conj(tmp));
      } else {
        b2 = 1;
        fs = ((c(0, 2) * c(2, 0) - c(2, 2) * c(0, 0)) /
              (c(0, 2) + c(2, 0) - c(0, 0) - c(2, 2)))
                 .real();
        fd = c(2, 2).real() - fs;
        complexd tmp = ((c(0, 2) - fs) / fd);
        a2 = std::abs(tmp * std::conj(tmp));
      }

      outPs[idx] = std::abs(fs) * (1.0 + b2);
      outPd[idx] = std::abs(fd) * (1.0 + a2);
      outPh[idx] = std::abs(fh);
      outPv[idx] = std::abs(fv);
    }
  }
}