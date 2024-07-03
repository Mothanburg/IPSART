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

using complexf = std::complex<float>;
constexpr float pi = std::numbers::pi_v<float>;

void CloudePottier(long height, long width, const float *t11, const float *t22,
                   const float *t33, const float *t12_r, const float *t13_r,
                   const float *t23_r, const float *t12_i, const float *t13_i,
                   const float *t23_i, float *outH, float *outAlpha,
                   float *outA) {
  #pragma omp parallel for
  for (int row = 0; row < height; row++) {
    for (int col = 0; col < width; col++) {
      int idx = row * width + col;
      Eigen::Matrix3cf t;
      t(0, 0) = complexf(t11[idx], 0.0);
      t(0, 1) = complexf(t12_r[idx], t12_i[idx]);
      t(0, 2) = complexf(t13_r[idx], t13_i[idx]);
      t(1, 0) = std::conj(t(0, 1));
      t(1, 1) = complexf(t22[idx], 0.0);
      t(1, 2) = complexf(t23_r[idx], t23_i[idx]);
      t(2, 0) = std::conj(t(0, 2));
      t(2, 1) = std::conj(t(1, 2));
      t(2, 2) = complexf(t33[idx], 0.0);

      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cf> eig(t);
      const Eigen::Array3cf eig_vals = eig.eigenvalues().array();
      const Eigen::Matrix3cf eig_vecs = eig.eigenvectors();

      const Eigen::Array3f p = eig_vals.real() / eig_vals.sum().real() + 1.0e-40f;
      outH[idx] = -(p * p.log()).sum() / std::log(3.0f);

      const Eigen::Array3f alphas = eig_vecs.row(0).array().abs();
      outAlpha[idx] = 180.0f * (p * alphas.acos()).sum() / pi;

      const float p1 = p.maxCoeff();
      const float p3 = p.minCoeff();
      const float p_sum = p.sum();
      outA[idx] = (p_sum - p1 - 2.0f * p3) / (p_sum - p1);
    }
  }
}
