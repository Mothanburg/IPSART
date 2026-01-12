#include "GaRS.h"

#include "utils.hpp"

using namespace std;
using namespace Eigen;

template <typename TData, int Dim>
static void cloude_pottier(const PolMatView<TData, Dim> &pol_mat, TData *outH,
                           TData *outAlpha, TData *outA) {
  using TMat = PolMatView<TData, Dim>::TMat;
  using TRowVec = Array<TData, Dim, 1>;

  constexpr TData PI = std::numbers::pi_v<TData>;
  constexpr auto inv_log3 = static_cast<TData>(1.0) / log(3.0);

  int len = pol_mat.Rows * pol_mat.Cols;
# pragma omp parallel for
  for (int idx = 0; idx < len; idx++) {
    const TMat t = pol_mat.at(idx);

    const SelfAdjointEigenSolver<TMat> eig(t);
    const TRowVec eig_vals = eig.eigenvalues().array().abs();
    const TMat eig_vecs = eig.eigenvectors();

    const TRowVec p = eig_vals / eig_vals.sum();
    const TRowVec alphas = eig_vecs.row(0).array().abs();

    if constexpr (Dim == 2) {
      outH[idx] =
          p.unaryExpr([](auto x) { return x != 0 ? -x * log2(x) : 0; }).sum();
      outA[idx] = abs(p(0) - p(1)) / p.sum();
    } else {
      outH[idx] =
          p.unaryExpr([](auto x) {
             return x != 0 ? -x * log(x) * inv_log3 : 0;
           }).sum();
      const TData p1 = p.maxCoeff();
      const TData p3 = p.minCoeff();
      const TData p_sum = p.sum();
      outA[idx] = (p_sum - p1 - 2 * p3) / (p_sum - p1);
    }

    outAlpha[idx] = 180 * (p * alphas.acos()).sum() / PI;
  }
}

void CloudePottier3f(int rows, int cols, const float *t11, const float *t22,
                     const float *t33, const float *t12r, const float *t13r,
                     const float *t23r, const float *t12i, const float *t13i,
                     const float *t23i, float *outH, float *outAlpha,
                     float *outA) {
  PolMatView<float, 3> mat(rows, cols, t11, t22, t33, t12r, t13r, t23r, t12i,
                           t13i, t23i);
  cloude_pottier(mat, outH, outAlpha, outA);
}

void CloudePottier3d(int rows, int cols, const double *t11, const double *t22,
                     const double *t33, const double *t12r, const double *t13r,
                     const double *t23r, const double *t12i, const double *t13i,
                     const double *t23i, double *outH, double *outAlpha,
                     double *outA) {
  PolMatView<double, 3> mat(rows, cols, t11, t22, t33, t12r, t13r, t23r, t12i,
                            t13i, t23i);
  cloude_pottier(mat, outH, outAlpha, outA);
}

void CloudePottier2f(int rows, int cols, const float *t11, const float *t22,
                     const float *t12r, const float *t12i, float *outH,
                     float *outAlpha, float *outA) {
  PolMatView<float, 2> mat(rows, cols, t11, t22, t12r, t12i);
  cloude_pottier(mat, outH, outAlpha, outA);
}

void CloudePottier2d(int rows, int cols, const double *t11, const double *t22,
                     const double *t12r, const double *t12i, double *outH,
                     double *outAlpha, double *outA) {
  PolMatView<double, 2> mat(rows, cols, t11, t22, t12r, t12i);
  cloude_pottier(mat, outH, outAlpha, outA);
}