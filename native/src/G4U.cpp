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
#pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    TMat t0 = pol_mat.at(idx);

    TData two_theta = atan(t0(1, 2).real() * static_cast<TData>(2.0) /
                           (t0(1, 1).real() - t0(2, 2).real())) /
                      static_cast<TData>(2.0);

    TMat t = t0;
    if (!isnan(two_theta)) {
      TMatReal r;
      r << static_cast<TData>(1.0), static_cast<TData>(0.0),
          static_cast<TData>(0.0), static_cast<TData>(0.0), cos(two_theta),
          sin(two_theta), static_cast<TData>(0.0), -sin(two_theta),
          cos(two_theta);
      t = r * t0 * r.transpose();
    }
    TData tp = t.trace().real();

    TData fs, fd, fv, fh = abs(t(1, 2).imag()) * static_cast<TData>(2.0);

    TData c1 =
        t(0, 0).real() - t(1, 1).real() +
        t(2, 2).real() * static_cast<TData>(7.0) / static_cast<TData>(8.0) +
        fh / static_cast<TData>(16.0);

    if (c1 > static_cast<TData>(0.0)) {
      TComplex c;
      TData coratio = log10((t(0, 0).real() + t(1, 1).real() -
                             t(0, 1).real() * static_cast<TData>(2.0)) /
                            (t(0, 0).real() + t(1, 1).real() +
                             t(0, 1).real() * static_cast<TData>(2.0))) *
                      static_cast<TData>(10.0);

      if (coratio < static_cast<TData>(-2.0)) {
        fv = (t(2, 2).real() * static_cast<TData>(2.0) - fh) *
             static_cast<TData>(15.0) / static_cast<TData>(8.0);
        if (fv < static_cast<TData>(0.0)) {
          fh = static_cast<TData>(0.0);
          fv = t(2, 2).real() * static_cast<TData>(2.0) *
               static_cast<TData>(15.0) / static_cast<TData>(8.0);
        }
        c = t(0, 1) + t(0, 2) - fv / static_cast<TData>(6.0);
      } else if (coratio > static_cast<TData>(-2.0) &&
                 coratio < static_cast<TData>(2.0)) {
        fv = (t(2, 2).real() * static_cast<TData>(2.0) - fh) *
             static_cast<TData>(2.0);
        if (fv < static_cast<TData>(0.0)) {
          fh = static_cast<TData>(0.0);
          fv = t(2, 2).real() * static_cast<TData>(2.0) *
               static_cast<TData>(2.0);
        }
        c = t(0, 1) + t(0, 2);
      } else {
        fv = (t(2, 2).real() * static_cast<TData>(2.0) - fh) *
             static_cast<TData>(15.0) / static_cast<TData>(8.0);
        if (fv < static_cast<TData>(0.0)) {
          fh = static_cast<TData>(0.0);
          fv = t(2, 2).real() * static_cast<TData>(2.0) *
               static_cast<TData>(15.0) / static_cast<TData>(8.0);
        }
        c = t(0, 1) + t(0, 2) + fv / static_cast<TData>(6.0);
      }

      TData s = t(0, 0).real() - fv / static_cast<TData>(2.0);
      TData d = tp - fv - fh - s;

      if ((fv + fh) >= tp) {
        outPs[idx] = static_cast<TData>(0.0);
        outPd[idx] = static_cast<TData>(0.0);
        outPv[idx] = tp - fh;
        outPh[idx] = fh;
        continue;
      } else {
        TData c0 = t(0, 0).real() * static_cast<TData>(2.0) + fh - tp;
        if (c0 > static_cast<TData>(0.0)) {
          fs = s + abs(c) * abs(c) / s;
          fd = d - abs(c) * abs(c) / s;
        } else {
          fs = s - abs(c) * abs(c) / d;
          fd = d + abs(c) * abs(c) / d;
        }
      }
    } else {
      fv = (t(2, 2).real() * 2.0f - fh) * static_cast<TData>(15.0) /
           static_cast<TData>(16.0);
      if (fv < static_cast<TData>(0.0)) {
        fh = static_cast<TData>(0.0);
        fv = t(2, 2).real() * static_cast<TData>(2.0) *
             static_cast<TData>(15.0) / static_cast<TData>(16.0);
      }
      TData s = t(0, 0).real();
      TData d = tp - fv - fh - s;
      TComplex c = t(0, 1) + t(0, 2);
      fs = s - abs(c) * abs(c) / d;
      fd = d + abs(c) * abs(c) / d;
    }

    if (fs >= static_cast<TData>(0.0) && fd >= static_cast<TData>(0.0)) {
      outPs[idx] = fs;
      outPd[idx] = fd;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else if (fs >= static_cast<TData>(0.0) && fd < static_cast<TData>(0.0)) {
      outPs[idx] = tp - fv - fh;
      outPd[idx] = 0;
      outPv[idx] = fv;
      outPh[idx] = fh;
    } else if (fs < static_cast<TData>(0.0) && fd >= static_cast<TData>(0.0)) {
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