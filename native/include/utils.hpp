#pragma once

#ifndef GARS_UTIL
#define GARS_UTIL

#include "pch.h"

template <typename TData, int Dim>
struct PolMatView;

template <typename TData>
  requires std::is_same_v<TData, float> || std::is_same_v<TData, double>
struct PolMatView<TData, 2> {
 public:
  using TComplex = std::complex<TData>;
  using TMat = Eigen::Matrix2<TComplex>;

  PolMatView(int rows, int cols, const TData *t11, const TData *t22,
             const TData *t12r, const TData *t12i)
      : Rows(rows),
        Cols(cols),
        m11(t11, rows * cols),
        m22(t22, rows * cols),
        m12r(t12r, rows * cols),
        m12i(t12i, rows * cols) {}

  TMat at(int idx) const {
    TMat m;
    m << TComplex(m11[idx]), TComplex(m12r[idx], m12i[idx]),
        TComplex(m12r[idx], -m12i[idx]), TComplex(m22[idx]);
    return m;
  }

  TMat at(int row, int col) const { return this->at(Cols * row + col); }

  int Rows, Cols;
  std::span<const TData> m11, m22, m12r, m12i;
};

template <typename TData>
  requires std::is_same_v<TData, float> || std::is_same_v<TData, double>
struct PolMatView<TData, 3> {
  using TComplex = std::complex<TData>;
  using TMat = Eigen::Matrix3<TComplex>;

  PolMatView(int rows, int cols, const TData *t11, const TData *t22,
             const TData *t33, const TData *t12r, const TData *t13r,
             const TData *t23r, const TData *t12i, const TData *t13i,
             const TData *t23i)
      : Rows(rows),
        Cols(cols),
        m11(t11, rows * cols),
        m22(t22, rows * cols),
        m33(t33, rows * cols),
        m12r(t12r, rows * cols),
        m13r(t13r, rows * cols),
        m23r(t23r, rows * cols),
        m12i(t12i, rows * cols),
        m13i(t13i, rows * cols),
        m23i(t23i, rows * cols) {}

  TMat at(int idx) const {
    TMat m;
    m << TComplex(m11[idx]), TComplex(m12r[idx], m12i[idx]),
        TComplex(m13r[idx], m13i[idx]), TComplex(m12r[idx], -m12i[idx]),
        TComplex(m22[idx]), TComplex(m23r[idx], m23i[idx]),
        TComplex(m13r[idx], -m13i[idx]), TComplex(m23r[idx], -m23i[idx]),
        TComplex(m33[idx]);
    return m;
  }

  TMat at(int row, int col) const { return this->at(Cols * row + col); }

  int Rows, Cols;
  std::span<const TData> m11, m22, m33, m12r, m13r, m23r, m12i, m13i, m23i;
};

#endif