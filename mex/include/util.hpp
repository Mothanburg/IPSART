#pragma once

#ifndef _IPSART_UTIL_
#define _IPSART_UTIL_

#include <array>
#include <complex>
#include <cstddef>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <mex.hpp>

#include <Eigen/Dense>

namespace ipsart {
namespace detail {
template <typename Ty>
inline void check_contiguous(const matlab::data::TypedArray<Ty> &arr) {
#ifdef CHECK_CONTIGUOUS
  auto start = std::to_address(arr.begin());
  auto end = std::to_address(arr.end());
  auto size = arr.getNumberOfElements();
  if (end - start != static_cast<ptrdiff_t>(size)) {
    throw std::runtime_error("Input array is not contiguous in memory!");
  }
#endif
}
} // namespace detail

template <typename Ty>
inline std::span<Ty> marray_to_span(const matlab::data::TypedArray<Ty> &arr) {
  detail::check_contiguous(arr);
  auto src_ptr = const_cast<Ty *>(std::to_address(arr.begin()));
  auto size = arr.getNumberOfElements();
  return std::span(src_ptr, size);
}

template <typename Float, int Dim>
  requires(Dim > 1 && Dim < 5 && std::is_floating_point_v<Float>)
struct PolMatView {
  using mexarray = matlab::data::TypedArray<Float>;
  using Matrix = Eigen::Matrix<std::complex<Float>, Dim, Dim>;

  explicit PolMatView(const std::vector<mexarray> &elements) {
    assert(elements.size() == Dim * Dim);
    auto dims = elements[0].getDimensions();
    rows = dims[0];
    cols = dims[1];
    total_elements = static_cast<size_t>(rows) * cols;
    for (int i = 0; i < Dim * Dim; i++) {
      element_refs[i] = marray_to_span(elements[i]);
    }
  }

  size_t numel() const { return total_elements; }

  std::pair<int, int> dimensions() const { return {rows, cols}; }

  Matrix at(int idx) const {
    if constexpr (Dim == 2) {
      Matrix mat;
      mat(0, 0) = element_refs[0][idx];
      mat(1, 1) = element_refs[1][idx];
      mat(0, 1) = std::complex(element_refs[2][idx], element_refs[3][idx]);
      mat(1, 0) = std::conj(mat(0, 1));
      return mat;
    } else if constexpr (Dim == 3) {
      Matrix mat;
      mat(0, 0) = element_refs[0][idx];
      mat(1, 1) = element_refs[1][idx];
      mat(2, 2) = element_refs[2][idx];
      mat(0, 1) = std::complex(element_refs[3][idx], element_refs[6][idx]);
      mat(0, 2) = std::complex(element_refs[4][idx], element_refs[7][idx]);
      mat(1, 2) = std::complex(element_refs[5][idx], element_refs[8][idx]);
      mat(1, 0) = std::conj(mat(0, 1));
      mat(2, 0) = std::conj(mat(0, 2));
      mat(2, 1) = std::conj(mat(1, 2));
      return mat;
    } else /* Dim == 4 */ {
      Matrix mat;
      mat(0, 0) = element_refs[0][idx];
      mat(1, 1) = element_refs[1][idx];
      mat(2, 2) = element_refs[2][idx];
      mat(3, 3) = element_refs[3][idx];
      mat(0, 1) = std::complex(element_refs[4][idx], element_refs[10][idx]);
      mat(0, 2) = std::complex(element_refs[5][idx], element_refs[11][idx]);
      mat(0, 3) = std::complex(element_refs[6][idx], element_refs[12][idx]);
      mat(1, 2) = std::complex(element_refs[7][idx], element_refs[13][idx]);
      mat(1, 3) = std::complex(element_refs[8][idx], element_refs[14][idx]);
      mat(2, 3) = std::complex(element_refs[9][idx], element_refs[15][idx]);
      mat(1, 0) = std::conj(mat(0, 1));
      mat(2, 0) = std::conj(mat(0, 2));
      mat(3, 0) = std::conj(mat(3, 0));
      mat(2, 1) = std::conj(mat(1, 2));
      mat(3, 1) = std::conj(mat(1, 3));
      mat(3, 2) = std::conj(mat(2, 3));
      return mat;
    }
  }
  Matrix at(int row, int col) const {
    return this->at(rows * col + row); // default to col-major
  }

  std::array<std::span<Float>, Dim * Dim> element_refs;
  int rows, cols;
  size_t total_elements;
};

} // namespace ipsart




#endif