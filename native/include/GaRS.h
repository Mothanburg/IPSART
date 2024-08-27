// clang-format off
#ifndef GARS
#define GARS

#ifdef _WIN32
# ifdef COMPILING_GARS
#   define GARS_API extern __declspec(dllexport)
# else
#   define GARS_API extern __declspec(dllimport)
# endif
#else
# define GARS_API extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------ BASIC FIlTERING ------------------------*/

// Multilook
GARS_API int Multilookf(int rowLook, int colLook,
                        int inRows, int inCols, const float *input,
                        int outRows, int outCols, float *output);

GARS_API int Multilookd(int rowLook, int colLook,
                        int inRows, int inCols, const double *input,
                        int outRows, int outCols, double *output);


/*------------------------ POLSAR FIlTERING ------------------------*/

// 7x7 Refined Lee Filter by OpenCL
GARS_API int RefinedLeeFilting3f(int lookNum,
                                 int rows, int cols,
                                 const float *c11, const float *c22, const float *c33,
                                 const float *c12r, const float *c13r, const float *c23r, 
                                 const float *c12i, const float *c13i, const float *c23i,
                                 float *outC11, float *outC22, float *outC33,
                                 float *outC12r, float *outC13r, float *outC23r,
                                 float *outC12i, float *outC13i, float *outC23i);

GARS_API int RefinedLeeFilting3d(int lookNum,
                                 int rows, int cols,
                                 const double *c11, const double *c22, const double *c33,
                                 const double *c12r, const double *c13r, const double *c23r,
                                 const double *c12i, const double *c13i, const double *c23i,
                                 double *outC11, double *outC22, double *outC33,
                                 double *outC12r, double *outC13r, double *outC23r,
                                 double *outC12i, double *outC13i, double *outC23i);

// Refined Lee Filter for 2x2 Covariance matrix
GARS_API int RefinedLeeFilting2f(int lookNum,
                                 int rows, int cols,
                                 const float *c11, const float *c22,
                                 const float *c12r, const float *c12i,
                                 float *outC11, float *outC22,
                                 float *outC12r, float *outC12i);

GARS_API int RefinedLeeFilting2d(int lookNum,
                                 int rows, int cols,
                                 const double *c11, const double *c22,
                                 const double *c12r, const double *c12i,
                                 double *outC11, double *outC22,
                                 double *outC12r, double *outC12i);


/*------------------------ POLSAR DECOMPOSITION ------------------------*/

// Cloude-Pottier H/a/A decomposition
GARS_API void CloudePottier3f(int rows, int cols,
                              const float *t11, const float *t22, const float *t33,
                              const float *t12r, const float *t13r, const float *t23r,
                              const float *t12i, const float *t13i, const float *t23i,
                              float *outH, float *outAlpha, float *outA);

GARS_API void CloudePottier3d(int rows, int cols,
                              const double *t11, const double *t22, const double *t33,
                              const double *t12r, const double *t13r, const double *t23r,
                              const double *t12i, const double *t13i, const double *t23i,
                              double *outH, double *outAlpha, double *outA);

// Cloude-Pottier H/a/A decomposition for dual-pol
GARS_API void CloudePottier2f(int rows, int cols,
                              const float *c11, const float *c22,
                              const float *c12r, const float *c12i,
                              float *outH, float *outAlpha, float *outA);

GARS_API void CloudePottier2d(int rows, int cols,
                              const double *c11, const double *c22,
                              const double *c12r, const double *c12i,
                              double *outH, double *outAlpha, double *outA);

// Yamaguchi four component decomposition
GARS_API void Yamaguchif(int rows, int cols,
                         const float *c11, const float *c22, const float *c33,
                         const float *c12r, const float *c13r, const float *c23r,
                         const float *c12i, const float *c13i, const float *c23i,
                         float *outPs, float *outPd, float *outPv, float *outPh);

GARS_API void Yamaguchid(int rows, int cols,
                         const double *c11, const double *c22, const double *c33,
                         const double *c12r, const double *c13r, const double *c23r,
                         const double *c12i, const double *c13i, const double *c23i,
                         double *outPs, double *outPd, double *outPv, double *outPh);

// General four component decomposition with unitary transformation T
GARS_API int G4Uf(int rows, int cols,
                  const float *t11, const float *t22, const float *t33,
                  const float *t12r, const float *t13r, const float *t23r,
                  const float *t12i, const float *t13i, const float *t23i,
                  float *outPs, float *outPd, float *outPv, float *outPh);

GARS_API int G4Ud(int rows, int cols,
                  const double *t11, const double *t22, const double *t33,
                  const double *t12r, const double *t13r, const double *t23r,
                  const double *t12i, const double *t13i, const double *t23i,
                  double *outPs, double *outPd, double *outPv, double *outPh);

#ifdef __cplusplus
}
#endif

#endif