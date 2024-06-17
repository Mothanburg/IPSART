#ifndef TINYC
#define TINYC

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


GARS_API void RefinedLeeFilterC3(long height, long width, const double *c11,
                               const double *c22, const double *c33,
                               const double *c12_r, const double *c13_r,
                               const double *c23_r, const double *c12_i,
                               const double *c13_i, const double *c23_i,
                               double *outC11, double *outC22, double *outC33,
                               double *outC12_r, double *outC13_r,
                               double *outC23_r, double *outC12_i,
                               double *outC13_i, double *outC23_I);

GARS_API void CloudePottier(long height, long width, const double *t11,
                            const double *t22, const double *t33,
                            const double *t12_r, const double *t13_r,
                            const double *t23_r, const double *t12_i,
                            const double *t13_i, const double *t23_i,
                            double *outH, double *outAlpha, double *outA);

GARS_API void Yamaguchi(long height, long width, const double *c11,
                        const double *c22, const double *c33,
                        const double *c12_r, const double *c13_r,
                        const double *c23_r, const double *c12_i,
                        const double *c13_i, const double *c23_i, double *outPs,
                        double *outPd, double *outPv, double *outPh);

#ifdef __cplusplus
}
#endif

#endif