#include "GaRS.h"

#include "utils.hpp"

using namespace Eigen;
using namespace std;
using namespace std::literals;

template <typename TData, int Dim>
static void yamaguchi(const PolMatView<TData, Dim> &pol_mat, TData *outPs,
                      TData *outPd, TData *outPv, TData *outPh) {
  using TComplex = PolMatView<TData, Dim>::TComplex;
  using TMat = PolMatView<TData, Dim>::TMat;
  using TArray =
      Array<TComplex, TMat::RowsAtCompileTime, TMat::ColsAtCompileTime>;

  constexpr TData sqrt2 = std::numbers::sqrt2_v<TData>;

  int len = pol_mat.Rows * pol_mat.Cols;
# pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    TArray c = pol_mat.at(idx).array();

    // The form of helix scattering
    TArray ch;
    if (c(0, 1).imag() + c(1, 2).imag() > 0) {
      ch << static_cast<TData>(1.0), TComplex(0, sqrt2),
          static_cast<TData>(-1.0), TComplex(0, -sqrt2),
          static_cast<TData>(2.0), TComplex(0, sqrt2), static_cast<TData>(-1.0),
          TComplex(0, -sqrt2), static_cast<TData>(1.0);
      ch /= static_cast<TData>(4.0);
    } else {
      ch << static_cast<TData>(1.0), TComplex(0, -sqrt2),
          static_cast<TData>(-1.0), TComplex(0, sqrt2), static_cast<TData>(2.0),
          TComplex(0, -sqrt2), static_cast<TData>(-1.0), TComplex(0, sqrt2),
          static_cast<TData>(1.0);
      ch /= static_cast<TData>(4.0);
    }
    TData fh = static_cast<TData>(2.0) * abs(c(0, 1).imag() + c(1, 2).imag());

    // The form of volume scattering
    const TData coratio =
        static_cast<TData>(10.0) * log10(c(2, 2).real() / c(0, 0).real());
    TArray cv;
    TData fv;
    if (coratio < static_cast<TData>(-2.0)) {
      cv << static_cast<TData>(8.0), static_cast<TData>(0.0),
          static_cast<TData>(2.0), static_cast<TData>(0.0),
          static_cast<TData>(4.0), static_cast<TData>(0.0),
          static_cast<TData>(2.0), static_cast<TData>(0.0),
          static_cast<TData>(3.0);
      cv /= static_cast<TData>(15.0);
      fv = static_cast<TData>(15.0) *
           (c(1, 1).real() - fh / static_cast<TData>(2.0)) /
           static_cast<TData>(4.0);
    } else if (coratio < 2) {
      cv << 3.0f, 0.0f, 1.0f, 0.0f, 2.0f, 0.0f, 1.0f, 0.0f, 3.0f;
      cv /= static_cast<TData>(8.0);
      fv = static_cast<TData>(4.0) *
           (c(1, 1).real() - fh / static_cast<TData>(2.0));
    } else {
      cv << static_cast<TData>(3.0), static_cast<TData>(0.0),
          static_cast<TData>(2.0), static_cast<TData>(0.0),
          static_cast<TData>(4.0), static_cast<TData>(0.0),
          static_cast<TData>(2.0), static_cast<TData>(0.0),
          static_cast<TData>(8.0);
      cv /= static_cast<TData>(15.0);
      fv = static_cast<TData>(15.0) *
           (c(1, 1).real() - fh / static_cast<TData>(2.0)) /
           static_cast<TData>(4.0);
    }

    // Remove volume and helix scattering only when cross-pol scattering is
    // small
    if (c(1, 1).real() < c(0, 0).real() && c(1, 1).real() < c(2, 2).real()) {
      c = c - fh * ch - fv * cv;
    } else {
      fh = static_cast<TData>(0.0);
      fv = static_cast<TData>(0.0);
    }

    TData a2, b2, fd, fs;
    if (c(0, 2).real() > 0) {
      a2 = 1.0f;
      fd = ((c(2, 2) * c(0, 0) - c(0, 2) * c(2, 0)) /
            (c(2, 2) + c(0, 0) + c(0, 2) + c(2, 0)))
               .real();
      fs = c(2, 2).real() - fd;
      TComplex tmp = ((c(0, 2) + fd) / fs);
      b2 = abs(tmp * conj(tmp));
    } else {
      b2 = 1;
      fs = ((c(0, 2) * c(2, 0) - c(2, 2) * c(0, 0)) /
            (c(0, 2) + c(2, 0) - c(0, 0) - c(2, 2)))
               .real();
      fd = c(2, 2).real() - fs;
      TComplex tmp = ((c(0, 2) - fs) / fd);
      a2 = abs(tmp * conj(tmp));
    }

    outPs[idx] = abs(fs) * (static_cast<TData>(1.0) + b2);
    outPd[idx] = abs(fd) * (static_cast<TData>(1.0) + a2);
    outPh[idx] = abs(fh);
    outPv[idx] = abs(fv);
  }
}

void Yamaguchif(int rows, int cols, const float *c11, const float *c22,
                const float *c33, const float *c12r, const float *c13r,
                const float *c23r, const float *c12i, const float *c13i,
                const float *c23i, float *outPs, float *outPd, float *outPv,
                float *outPh) {
  PolMatView<float, 3> mat(rows, cols, c11, c22, c33, c12r, c13r, c23r, c12i,
                           c13i, c23i);
  yamaguchi(mat, outPs, outPd, outPv, outPh);
}

void Yamaguchid(int rows, int cols, const double *c11, const double *c22,
                const double *c33, const double *c12r, const double *c13r,
                const double *c23r, const double *c12i, const double *c13i,
                const double *c23i, double *outPs, double *outPd, double *outPv,
                double *outPh) {
  PolMatView<double, 3> mat(rows, cols, c11, c22, c33, c12r, c13r, c23r, c12i,
                            c13i, c23i);
  yamaguchi(mat, outPs, outPd, outPv, outPh);
}