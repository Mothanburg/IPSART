#ifndef GARS
#define GARS

// clang-format off
#ifdef _WIN32
# ifdef COMPILING_GARS
#   define GARS_API extern __declspec(dllexport)
# else
#   define GARS_API extern __declspec(dllimport)
# endif
#else
# define GARS_API extern
#endif
// clang-format on

#ifdef __cplusplus
extern "C" {
#endif

/*
 * GaRS library capability tests
 */

GARS_API int GaRSTestOpenCL();

/*
 * GaRS library algorithms
 */

// 7x7 Refined Lee Filter by OpenCL
GARS_API int RefinedLeeFilter3x3(
    int nLooks, int height, int width, const float *c11, const float *c22,
    const float *c33, const float *c12r, const float *c13r, const float *c23r,
    const float *c12i, const float *c13i, const float *c23i, float *outC11,
    float *outC22, float *outC33, float *outC12r, float *outC13r,
    float *outC23r, float *outC12i, float *outC13i, float *outC23i);

// Refined Lee Filter for 2x2 Covariance matrix
GARS_API int RefinedLeeFilter2x2(int nLooks, int height, int width,
                                 const float *c11, const float *c22,
                                 const float *c12r, const float *c12i,
                                 float *outC11, float *outC22, float *outC12r,
                                 float *outC12i);

// Multilook
GARS_API int Multilook(int height, int width, const float *image, int rowLook,
                       int colLook, int outHeight, int outWidth, float *looked);

// Cloude-Pottier H/a/A decomposition
GARS_API void CloudePottier(int height, int width, const float *t11,
                            const float *t22, const float *t33,
                            const float *t12r, const float *t13r,
                            const float *t23r, const float *t12i,
                            const float *t13i, const float *t23i, float *outH,
                            float *outAlpha, float *outA);

// Cloude-Pottier H/a/A decomposition for dual-pol
GARS_API void CloudePottierDP(int height, int width, const float *c11,
                              const float *c22, const float *c12r,
                              const float *c12i, float *outH, float *outAlpha,
                              float *outA);

// Yamaguchi four component decomposition
GARS_API void Yamaguchi(int height, int width, const float *c11,
                        const float *c22, const float *c33, const float *c12r,
                        const float *c13r, const float *c23r, const float *c12i,
                        const float *c13i, const float *c23i, float *outPs,
                        float *outPd, float *outPv, float *outPh);

// General four component decomposition with unitary transformation T
GARS_API int G4U(int height, int width, const float *t11, const float *t22,
                 const float *t33, const float *t12r, const float *t13r,
                 const float *t23r, const float *t12i, const float *t13i,
                 const float *t23i, float *outPs, float *outPd, float *outPv,
                 float *outPh);

#ifdef __cplusplus
}
#endif

#endif