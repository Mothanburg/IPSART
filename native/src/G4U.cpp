#include <cmath>
#include <complex>

#include "GaRS.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

int G4U(int height, int width, const float *t11, const float *t22,
        const float *t33, const float *t12r, const float *t13r,
        const float *t23r, const float *t12i, const float *t13i,
        const float *t23i, float *outPs, float *outPd, float *outPv,
        float *outPh) {
  auto len = height * width;
  volatile bool err_flag{false};

# pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    // For almost all cases, error would not happen
    // We may not need this 'if' clause
    // if (err_flag) [[unlikely]] {
    //   continue;
    // }
    Matrix3cf t0;
    t0 << scomplex(t11[idx]), scomplex(t12r[idx], t12i[idx]),
        scomplex(t13r[idx], t13i[idx]), scomplex(t12r[idx], -t12i[idx]),
        scomplex(t22[idx]), scomplex(t23r[idx], t23i[idx]),
        scomplex(t13r[idx], -t13i[idx]), scomplex(t23r[idx], -t23i[idx]),
        scomplex(t33[idx]);

    float two_theta =
        atanf(t0(1, 2).real() * 2.0f / (t0(1, 1).real() - t0(2, 2).real())) /
        2.0f;

    Matrix3cf t = t0;
    if (!isnan(two_theta)) {
      Matrix3f r;
      r << 1.0f, 0.0f, 0.0f, 0.0f, cosf(two_theta), sinf(two_theta), 0.0f,
          -sinf(two_theta), cosf(two_theta);
      t = r * t0 * r.transpose();
    }
    float tp = t.trace().real();

    float fs{0.0f}, fd{0.0f}, fv{0.0f}, fh = abs(t(1, 2).imag()) * 2.0f;

    float c1 = t(0, 0).real() - t(1, 1).real() + t(2, 2).real() * 7.0f / 8.0f +
               fh / 16.0f;
    if (c1 > 0.0) {
      scomplex c;
      float coratio =
          log10f((t(0, 0).real() + t(1, 1).real() - t(0, 1).real() * 2.0f) /
                 (t(0, 0).real() + t(1, 1).real() + t(0, 1).real() * 2.0f)) *
          10.0f;
      if (coratio < -2.0f) {
        fv = (t(2, 2).real() * 2.0f - fh) * 15.0f / 8.0f;
        if (fv < 0.0f) {
          fh = 0.0f;
          fv = t(2, 2).real() * 2.0f * 15.0f / 8.0f;
        }
        c = t(0, 1) + t(0, 2) - fv / 6.0f;
      } else if (coratio > -2.0f && coratio < 2.0f) {
        fv = (t(2, 2).real() * 2.0f - fh) * 2.0f;
        if (fv < 0.0f) {
          fh = 0.0f;
          fv = t(2, 2).real() * 2.0f * 2.0f;
        }
        c = t(0, 1) + t(0, 2);
      } else {
        fv = (t(2, 2).real() * 2.0f - fh) * 15.0f / 8.0f;
        if (fv < 0.0f) {
          fh = 0.0f;
          fv = t(2, 2).real() * 2.0f * 15.0f / 8.0f;
        }
        c = t(0, 1) + t(0, 2) + fv / 6.0f;
      }

      float s = t(0, 0).real() - fv / 2.0f;
      float d = tp - fv - fh - s;

      if ((fv + fh) >= tp) {
        outPs[idx] = 0.0f;
        outPd[idx] = 0.0f;
        outPv[idx] = tp - fh;
        outPh[idx] = fh;
        continue;
      } else {
        float c0 = t(0, 0).real() * 2.0f + fh - tp;
        if (c0 > 0.0f) {
          fs = s + abs(c) * abs(c) / s;
          fd = d - abs(c) * abs(c) / s;
        } else {
          fs = s - abs(c) * abs(c) / d;
          fd = d + abs(c) * abs(c) / d;
        }
      }

    } else {
      fv = (t(2, 2).real() * 2.0f - fh) * 15.0f / 16.0f;
      if (fv < 0.0f) {
        fh = 0.0f;
        fv = t(2, 2).real() * 2.0f * 15.0f / 16.0f;
      }
      float s = t(0, 0).real();
      float d = tp - fv - fh - s;
      scomplex c = t(0, 1) + t(0, 2);
      fs = s - abs(c) * abs(c) / d;
      fd = d + abs(c) * abs(c) / d;
    }

    if (fs >= 0.0f && fd >= 0.0f) {
      outPs[idx] = fs;
      outPd[idx] = fd;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else if (fs >= 0.0f && fd < 0.0f) {
      outPs[idx] = tp - fv - fh;
      outPd[idx] = 0;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else if (fs < 0.0f && fd >= 0.0f) {
      outPs[idx] = 0;
      outPd[idx] = tp - fh - fv;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else {
      // This case should not happen for regular data
      err_flag = true;
    }
  }
  return err_flag ? -1 : 0;
}