#include <cmath>
#include <numbers>

#include "GaRS.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;

constexpr float pi = std::numbers::pi_v<float>;

void CloudePottier(int height, int width, const float *t11, const float *t22,
                   const float *t33, const float *t12r, const float *t13r,
                   const float *t23r, const float *t12i, const float *t13i,
                   const float *t23i, float *outH, float *outAlpha,
                   float *outA) {
  auto len = height * width;
# pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    Matrix3cf t;
    t << scomplex(t11[idx]), scomplex(t12r[idx], t12i[idx]),
        scomplex(t13r[idx], t13i[idx]), scomplex(t12r[idx], -t12i[idx]),
        scomplex(t22[idx]), scomplex(t23r[idx], t23i[idx]),
        scomplex(t13r[idx], -t13i[idx]), scomplex(t23r[idx], -t23i[idx]),
        scomplex(t33[idx]);

    const SelfAdjointEigenSolver<Matrix3cf> eig(t);
    const Array3f eig_vals = eig.eigenvalues().array().abs();
    const Matrix3cf eig_vecs = eig.eigenvectors();

    const Array3f p = eig_vals / eig_vals.sum();
    outH[idx] = p.unaryExpr([](auto x) {
                   return x != 0.0f ? -x * log(x) / log(3.0f) : 0.0f;
                 }).sum();

    const Array3f alphas = eig_vecs.row(0).array().abs();
    outAlpha[idx] = 180.0f * (p * alphas.acos()).sum() / pi;

    const float p1 = p.maxCoeff();
    const float p3 = p.minCoeff();
    const float p_sum = p.sum();
    outA[idx] = (p_sum - p1 - 2.0f * p3) / (p_sum - p1);
  }
}

void CloudePottierDP(int height, int width, const float *c11, const float *c22,
                     const float *c12r, const float *c12i, float *outH,
                     float *outAlpha, float *outA) {
  auto len = height * width;
# pragma omp parallel for
  for (auto idx = 0; idx < len; idx++) {
    Matrix2cf c;
    c << scomplex(c11[idx]), scomplex(c12r[idx], c12i[idx]),
        scomplex(c12r[idx], -c12i[idx]), scomplex(c22[idx]);

    const SelfAdjointEigenSolver<Matrix2cf> eig(c);
    const Array2f eig_vals = eig.eigenvalues().array().abs();
    const Matrix2cf eig_vecs = eig.eigenvectors();

    const Array2f p = eig_vals / eig_vals.sum();
    outH[idx] = p.unaryExpr([](auto x) {
                   return x != 0.0f ? -x * log2(x) : 0.0f;
                 }).sum();

    const Array2f alphas = eig_vecs.row(0).array().abs();
    outAlpha[idx] = 180.0f * (p * alphas.acos()).sum() / pi;

    outA[idx] = abs(p(0) - p(1)) / p.sum();
  }
}
