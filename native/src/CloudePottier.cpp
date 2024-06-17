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

using complexd = std::complex<double>;

void CloudePottier(long height, long width, const double *t11,
                     const double *t22, const double *t33, const double *t12_r,
                     const double *t13_r, const double *t23_r,
                     const double *t12_i, const double *t13_i,
                     const double *t23_i, double *outH, double *outAlpha,
                     double *outA) {
  #pragma omp parallel for
  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      int idx = j * width + i;
      Eigen::Matrix3cd t;
      t(0, 0) = complexd(t11[idx], 0.0);
      t(0, 1) = complexd(t12_r[idx], t12_i[idx]);
      t(0, 2) = complexd(t13_r[idx], t13_i[idx]);
      t(1, 0) = std::conj(t(0, 1));
      t(1, 1) = complexd(t22[idx], 0.0);
      t(1, 2) = complexd(t23_r[idx], t23_i[idx]);
      t(2, 0) = std::conj(t(0, 2));
      t(2, 1) = std::conj(t(1, 2));
      t(2, 2) = complexd(t33[idx], 0.0);

      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd> eig(t);
      const Eigen::Array3cd eig_vals = eig.eigenvalues().array();
      const Eigen::Matrix3cd eig_vecs = eig.eigenvectors();

      const Eigen::Array3d p = eig_vals.real() / eig_vals.sum().real();
      outH[idx] = -(p * p.log()).sum() / std::log(3.0);

      const Eigen::Array3d alphas = eig_vecs.row(0).array().abs();
      outAlpha[idx] = 180.0 * (p * alphas.acos()).sum() / std::numbers::pi;

      const double p1 = p.maxCoeff();
      const double p3 = p.minCoeff();
      const double p_sum = p.sum();
      outA[idx] = (p_sum - p1 - 2.0 * p3) / (p_sum - p1);
    }
  }
}
