#include "IPSART.h"
#include "util.hpp"

#include <array>
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

template <typename Float>
static void yamaguchi(const PolMatView<Float, 3> &T, vector<double> &Ps,
                      vector<double> &Pd, vector<double> &Pv,
                      vector<double> &Ph) {
  using Complex = std::complex<double>;
  using Matrix = Eigen::Matrix<Complex, 3, 3>;
  using MatrixR = Eigen::Matrix<double, 3, 3>;

  constexpr double pi_4 = std::numbers::pi / 4.0;
  constexpr double pi_2 = std::numbers::pi / 2.0;

  int len = T.numel();
#pragma omp parallel for
  for (int idx = 0; idx < len; idx++) {
    Matrix t0 = T.at(idx).template cast<Complex>();

    double theta =
        atan2(2.0 * t0(1, 2).real(), t0(1, 1).real() - t0(2, 2).real()) / 4.0;
    if (theta < -pi_4) [[unlikely]] {
      theta += pi_2;
    } else if (theta > pi_4) [[unlikely]] {
      theta -= pi_2;
    } else if (isnan(theta)) [[unlikely]] {
      theta = 0.0;
    }
    double two_theta = theta * 2.0;

    MatrixR r = MatrixR::Zero();
    r(0, 0) = 1;
    r(1, 1) = cos(two_theta);
    r(1, 2) = sin(two_theta);
    r(2, 1) = -r(1, 2);
    r(2, 2) = r(1, 1);
    Matrix t = r * t0 * r.transpose();

    double t11 = t(0, 0).real();
    double t22 = t(1, 1).real();
    double t33 = t(2, 2).real();
    Complex t12 = t(0, 1);
    Complex t13 = t(0, 2);
    Complex t23 = t(1, 2);
    double tp = t11 + t22 + t33;

    double fs, fd, fv;
    double fh = abs(t23.imag()) * 2.0;

    double s, d;
    Complex c;
    if (t11 > t22) {
      double coratio = 10.0 * log10((t11 + t22 - t12.real() * 2.0) /
                                    (t11 + t22 + t12.real() * 2.0));
      if (coratio < -2.0) {
        fv = 15.0 * (t33 / 4.0 - fh / 8.0);
        if (fv < 0.0) {
          fh = 0.0;
          fv = t33 * 15.0 / 4.0;
        }
        s = t11 - fv / 2.0;
        d = t22 - 7.0 * fv / 30.0 - fh / 2.0;
        c = t12 - fv / 6.0;
      } else if (coratio < 2.0) {
        fv = 4.0 * t33 - 2.0 * fh;
        if (fv < 0.0) {
          fh = 0.0;
          fv = 4.0 * t33;
        }
        s = t11 - fv / 2.0;
        d = t22 - t33;
        c = t12;
      } else {
        fv = 15.0 * (t33 / 4.0 - fh / 8.0);
        if (fv < 0.0) {
          fh = 0.0;
          fv = 15.0 * t33 / 4.0;
        }
        s = t11 - fv / 2.0;
        d = t22 - 7.0 * fv / 30.0 - fh / 2.0;
        c = t12 + fv / 6.0;
      }

      if (fv + fh > tp) {
        fs = 0;
        fd = 0;
        fv = tp - fh;
        Ps[idx] = fs;
        Pd[idx] = fd;
        Pv[idx] = fv;
        Ph[idx] = fh;
        continue;
      } else {
        double c0 = t11 - t22 - t33 + fh;
        if (c0 > 0.0) {
          fs = s + norm(c) / s;
          fd = d - norm(c) / s;
        } else {
          fs = s - norm(c) / d;
          fd = d + norm(c) / d;
        }
      }
    } else {
      fv = 15.0 * (t33 - fh / 2.0) / 8.0;
      if (fv < 0.0) {
        fh = 0.0;
        fv = 15.0 * t33 / 8.0;
      }
      s = t11;
      d = t22 - 7.0 * fv / 15.0 - fh / 2.0;
      c = t12;
      fs = s - norm(c) / d;
      fd = d + norm(c) / d;
    }

    if (fs >= 0 && fd >= 0) {
      Ps[idx] = fs;
      Pd[idx] = fd;
      Pv[idx] = fv;
      Ph[idx] = fh;
    } else if (fs >= 0 && fd < 0) {
      Ps[idx] = tp - fv - fh;
      Pd[idx] = 0;
      Pv[idx] = fv;
      Ph[idx] = fh;
    } else if (fs < 0 && fd >= 0) {
      Ps[idx] = 0;
      Pd[idx] = tp - fh - fv;
      Pv[idx] = fv;
      Ph[idx] = fh;
    } else {
      // This case should not happen for regular data
      Ps[idx] = NAN;
      Pd[idx] = NAN;
      Pv[idx] = NAN;
      Ph[idx] = NAN;
    }
  }
}

vector<data::Array> Yamaguchi(const vector<data::Array> &input,
                              data::ArrayFactory &af) {
  auto in_type = input[0].getType();
  auto in_dims = input[0].getDimensions();
  assert(input.size() == 9);
  assert(ranges::all_of(
      input | views::take(input.size() - 1), [&](const auto &arr) {
        return arr.getDimensions() == in_dims && arr.getType() == in_type;
      }));

  size_t numel = in_dims[0] * in_dims[1];
  vector<double> Ps(numel), Pd(numel), Pv(numel), Ph(numel);

  if (in_type == data::ArrayType::SINGLE) {
    PolMatView<float, 3> T(input |
                           ranges::to<vector<data::TypedArray<float>>>());
    yamaguchi(T, Ps, Pd, Pv, Ph);
  } else {
    PolMatView<double, 3> T(input |
                            ranges::to<vector<data::TypedArray<double>>>());
    yamaguchi(T, Ps, Pd, Pv, Ph);
  }

  return {af.createArray(in_dims, Ps.begin(), Ps.end()),
          af.createArray(in_dims, Pd.begin(), Pd.end()),
          af.createArray(in_dims, Pv.begin(), Pv.end()),
          af.createArray(in_dims, Ph.begin(), Ph.end())};
}
