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

// 7x7 RefinedLeeFilter by OpenCL
GARS_API int RefinedLeeFilter(long nLooks, long height, long width,
                               const float *c11, const float *c22,
                               const float *c33, const float *c12_r,
                               const float *c13_r, const float *c23_r,
                               const float *c12_i, const float *c13_i,
                               const float *c23_i, float *outC11, float *outC22,
                               float *outC33, float *outC12_r, float *outC13_r,
                               float *outC23_r, float *outC12_i,
                               float *outC13_i, float *outC23_i);

// Cloude-Pottier H/a/A decomposition
GARS_API void CloudePottier(long height, long width, const float *t11,
                            const float *t22, const float *t33,
                            const float *t12_r, const float *t13_r,
                            const float *t23_r, const float *t12_i,
                            const float *t13_i, const float *t23_i, float *outH,
                            float *outAlpha, float *outA);

// Yamaguchi four component decomposition
GARS_API void Yamaguchi(long height, long width, const float *c11,
                        const float *c22, const float *c33, const float *c12_r,
                        const float *c13_r, const float *c23_r,
                        const float *c12_i, const float *c13_i,
                        const float *c23_i, float *outPs, float *outPd,
                        float *outPv, float *outPh);

#ifdef __cplusplus
}
#endif

#endif