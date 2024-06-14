#include "GaRS.h"

#include <cmath>
#include <complex>
#include <numbers>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/EigenValues>

#include "omp.h"

using complexd = std::complex<double>;

void RefinedLeeFilter(RawPolMat3 in, RawPolMat3 out, int nLooks) {}

void CloudePottierT3(long height, long width, double* m11, double* m22,
                     double* m33, double* m12_r, double* m13_r, double* m23_r,
                     double* m12_i, double* m13_i, double* m23_i, double* outH,
                     double* outAlpha, double* outA) {
  // #pragma omp parallel for
  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      int idx = j * width + i;
      Eigen::Matrix3cd t;
      t(0, 0) = complexd(m11[idx], 0.0);
      t(0, 1) = complexd(m12_r[idx], m12_i[idx]);
      t(0, 2) = complexd(m13_r[idx], m13_i[idx]);
      t(1, 0) = std::conj(t(0, 1));
      t(1, 1) = complexd(m22[idx], 0.0);
      t(1, 2) = complexd(m23_r[idx], m23_i[idx]);
      t(2, 0) = std::conj(t(0, 2));
      t(2, 1) = std::conj(t(1, 2));
      t(2, 2) = complexd(m33[idx], 0.0);

      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> eig(t);
      const Eigen::Array3cd eig_vals = eig.eigenvalues().array();
      const Eigen::Matrix3cd eig_vecs = eig.eigenvectors();

      const Eigen::Array3d p = eig_vals.real() / eig_vals.sum().real();
      outH[idx] = -p.unaryExpr([=](complexd x) {
                      return x.real() * std::log(x.real()) / std::log(3.0);
                    }).sum();

      const Eigen::Array3cd alphas = eig_vecs.row(0).array();
      outAlpha[idx] =
          (p * alphas.unaryExpr([](complexd x) {
            return std::acos(std::abs(x) / 180.0 * std::numbers::pi);
          })).sum();

      const double p1 = p.maxCoeff();
      const double p3 = p.minCoeff();
      const double p_sum = p.sum();
      outA[idx] = (p_sum - p1 - 2.0 * p3) / (p_sum - p1);
    }
  }
}