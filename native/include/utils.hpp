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
        m_11(t11),
        m_22(t22),
        m_12r(t12r),
        m_12i(t12i) {}

  TMat at(int idx) const {
    TMat m;
    m << TComplex(m_11[idx]), TComplex(m_12r[idx], m_12i[idx]),
        TComplex(m_12r[idx], -m_12i[idx]), TComplex(m_22[idx]);
    return m;
  }

  TMat at(int row, int col) const { return this->at(Cols * row + col); }

  int Rows, Cols;
  const TData *m_11, *m_22, *m_12r, *m_12i;
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
        m_11(t11),
        m_22(t22),
        m_33(t33),
        m_12r(t12r),
        m_13r(t13r),
        m_23r(t23r),
        m_12i(t12i),
        m_13i(t13i),
        m_23i(t23i) {}

  TMat at(int idx) const {
    TMat m;
    m << TComplex(m_11[idx]), TComplex(m_12r[idx], m_12i[idx]),
        TComplex(m_13r[idx], m_13i[idx]), TComplex(m_12r[idx], -m_12i[idx]),
        TComplex(m_22[idx]), TComplex(m_23r[idx], m_23i[idx]),
        TComplex(m_13r[idx], -m_13i[idx]), TComplex(m_23r[idx], -m_23i[idx]),
        TComplex(m_33[idx]);
    return m;
  }

  TMat at(int row, int col) const { return this->at(Cols * row + col); }

  int Rows, Cols;
  const TData *m_11, *m_22, *m_33, *m_12r, *m_13r, *m_23r, *m_12i, *m_13i, *m_23i;
};

#endif