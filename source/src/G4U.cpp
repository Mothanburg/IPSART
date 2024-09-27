#include "GaRS.h"
#include "utils.hpp"

using namespace std;
using namespace Eigen;

template <typename TData, int Dim>
static int g4u(const PolMatView<TData, Dim> &pol_mat, TData *outPs,
               TData *outPd, TData *outPv, TData *outPh) {
  using TComplex = PolMatView<TData, Dim>::TComplex;
  using TMat = PolMatView<TData, Dim>::TMat;
  using TMatReal =
      Matrix<TData, TMat::RowsAtCompileTime, TMat::ColsAtCompileTime>;

  volatile bool err_flag{false};

  int len = pol_mat.Rows * pol_mat.Cols;
# pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    TMat t0 = pol_mat.at(idx);

    TData two_theta =
        atan((t0(1, 2).real() * 2) / (t0(1, 1).real() - t0(2, 2).real())) / 2;

    TMat t = t0;
    if (!isnan(two_theta)) {
      TMatReal r;
      r << 1, 0, 0, 0, cos(two_theta), sin(two_theta), 0, -sin(two_theta),
          cos(two_theta);
      t = r * t0 * r.transpose();
    }
    TData tp = t.trace().real();

    TData fs, fd, fv;
    TData fh = abs(t(1, 2).imag()) * 2;
    TData c1 =
        t(0, 0).real() - t(1, 1).real() + t(2, 2).real() * 7 / 8 + fh / 16;

    if (c1 > 0) {
      TComplex c;
      TData coratio =
          10 * log10((t(0, 0).real() + t(1, 1).real() - t(0, 1).real() * 2) /
                     (t(0, 0).real() + t(1, 1).real() + t(0, 1).real() * 2));

      if (coratio < -2) {
        fv = (t(2, 2).real() * 2 - fh) * 15 / 8;
        if (fv < 0) {
          fh = 0;
          fv = t(2, 2).real() * 2 * 15 / 8;
        }
        c = t(0, 1) + t(0, 2) - fv / 6;
      } else if (coratio > -2 && coratio < 2) {
        fv = (t(2, 2).real() * 2 - fh) * 2;
        if (fv < 0) {
          fh = 0;
          fv = t(2, 2).real() * 4;
        }
        c = t(0, 1) + t(0, 2);
      } else {
        fv = (t(2, 2).real() * 2 - fh) * 15 / 8;
        if (fv < 0) {
          fh = 0;
          fv = t(2, 2).real() * 2 * 15 / 8;
        }
        c = t(0, 1) + t(0, 2) + fv / 6;
      }

      TData s = t(0, 0).real() - fv / 2;
      TData d = tp - fv - fh - s;

      if ((fv + fh) >= tp) {
        outPs[idx] = 0;
        outPd[idx] = 0;
        outPv[idx] = tp - fh;
        outPh[idx] = fh;
        continue;
      } else {
        TData c0 = t(0, 0).real() * 2 + fh - tp;
        if (c0 > 0) {
          fs = s + abs(c) * abs(c) / s;
          fd = d - abs(c) * abs(c) / s;
        } else {
          fs = s - abs(c) * abs(c) / d;
          fd = d + abs(c) * abs(c) / d;
        }
      }
    } else {
      fv = (t(2, 2).real() * 2 - fh) * 15 / 16;
      if (fv < 0) {
        fh = 0;
        fv = t(2, 2).real() * 2 * 15 / 16;
      }
      TData s = t(0, 0).real();
      TData d = tp - fv - fh - s;
      TComplex c = t(0, 1) + t(0, 2);
      fs = s - abs(c) * abs(c) / d;
      fd = d + abs(c) * abs(c) / d;
    }

    if (fs >= 0 && fd >= 0) {
      outPs[idx] = fs;
      outPd[idx] = fd;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else if (fs >= 0 && fd < 0) {
      outPs[idx] = tp - fv - fh;
      outPd[idx] = 0;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else if (fs < 0 && fd >= 0) {
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

int G4Uf(int rows, int cols, const float *t11, const float *t22,
         const float *t33, const float *t12r, const float *t13r,
         const float *t23r, const float *t12i, const float *t13i,
         const float *t23i, float *outPs, float *outPd, float *outPv,
         float *outPh) {
  PolMatView<float, 3> mat(rows, cols, t11, t22, t33, t12r, t13r, t23r, t12i,
                           t13i, t23i);
  return g4u(mat, outPs, outPd, outPv, outPh);
}

int G4Ud(int rows, int cols, const double *t11, const double *t22,
         const double *t33, const double *t12r, const double *t13r,
         const double *t23r, const double *t12i, const double *t13i,
         const double *t23i, double *outPs, double *outPd, double *outPv,
         double *outPh) {
  PolMatView<double, 3> mat(rows, cols, t11, t22, t33, t12r, t13r, t23r, t12i,
                            t13i, t23i);
  return g4u(mat, outPs, outPd, outPv, outPh);
}