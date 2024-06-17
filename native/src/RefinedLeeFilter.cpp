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

using namespace Eigen;
using namespace std::literals;
using std::numbers::sqrt2;

using complexd = std::complex<double>;

void RefinedLeeFilterC3(long height, long width, const double *c11,
                      const double *c22, const double *c33, const double *c12_r,
                      const double *c13_r, const double *c23_r,
                      const double *c12_i, const double *c13_i,
                      const double *c23_i, double *outC11, double *outC22,
                      double *outC33, double *outC12_r, double *outC13_r,
                      double *outC23_r, double *outC12_i, double *outC13_i,
                      double *outC23_I) {
  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      int idx = j * width + i;


    }
  }
}
