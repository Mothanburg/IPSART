#pragma once

#ifndef GARS_PCH
#define GARS_PCH

#include <cmath>
#include <complex>
#include <concepts>
#include <numbers>
#include <span>
#include <type_traits>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <Eigen/Dense>

#endif