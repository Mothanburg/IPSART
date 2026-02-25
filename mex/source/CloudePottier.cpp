#include "IPSART.h"
#include "util.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <ranges>
#include <vector>

#include "omp.h"

#include <mex.hpp>

#include <Eigen/Dense>

using namespace std;
using namespace matlab;
using namespace Eigen;

template <typename Float, int Dim>
static void cloude_pottier(const PolMatView<Float, Dim> &T, vector<double> &H,
                           vector<double> &a, vector<double> &A) {
  using Matrix = Matrix<std::complex<double>, Dim, Dim>;
  using Row = Array<double, Dim, 1>;

  constexpr double PI = std::numbers::pi_v<double>;
  const double inv_log3 = 1.0 / log(3.0);

  int len = T.numel();
#pragma omp parallel for
  for (int idx = 0; idx < len; idx++) {
    const Matrix t = T.at(idx).template cast<std::complex<double>>();

    const SelfAdjointEigenSolver<Matrix> eig(t);
    const Row eig_vals = eig.eigenvalues().array().abs();
    const Matrix eig_vecs = eig.eigenvectors();

    const Row p = eig_vals / eig_vals.sum();
    const Row alphas = eig_vecs.row(0).array().abs();

    if constexpr (Dim == 2) {
      H[idx] = p.unaryExpr([](auto x) {
                  return x != 0.0 ? -x * log2(x) : 0.0;
                }).sum();
      A[idx] = abs(p(0) - p(1));
    } else {
      H[idx] = p.unaryExpr([&](auto x) {
                  return x != 0 ? -x * log(x) * inv_log3 : 0.0;
                }).sum();
      const double p1 = p.maxCoeff();
      const double p3 = p.minCoeff();
      const double p_sum = p.sum();
      A[idx] = (p_sum - p1 - 2.0 * p3) / (p_sum - p1);
    }
    a[idx] = 180.0 * (p * alphas.acos()).sum() / PI;
  }
}

vector<data::Array> CloudePottier(const vector<data::Array> &input,
                                  data::ArrayFactory &af) {
  auto in_num = input.size();
  auto in_type = input[0].getType();
  auto in_dims = input[0].getDimensions();
  assert(input.size() == 4 || input.size() == 9);
  assert(ranges::all_of(
      input | views::take(input.size() - 1), [&](const auto &arr) {
        return arr.getDimensions() == in_dims && arr.getType() == in_type;
      }));

  size_t numel = in_dims[0] * in_dims[1];
  vector<double> H(numel), a(numel), A(numel);

  auto execute = [&]<typename Ty>() {
    if (in_num == 4) {
      PolMatView<Ty, 2> T(input | ranges::to<vector<data::TypedArray<Ty>>>());
      cloude_pottier(T, H, a, A);
    } else {
      PolMatView<Ty, 3> T(input | ranges::to<vector<data::TypedArray<Ty>>>());
      cloude_pottier(T, H, a, A);
    }
  };

  if (in_type == data::ArrayType::SINGLE) {
    execute.template operator()<float>();
  } else {
    execute.template operator()<double>();
  }

  return {af.createArray(in_dims, H.begin(), H.end()),
          af.createArray(in_dims, a.begin(), a.end()),
          af.createArray(in_dims, A.begin(), A.end())};
}
