#ifndef TINYC
#define TINYC

#ifdef __cplusplus
extern "C" {
#endif

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

#include "cGaRS.h"

#ifdef __cplusplus
}
#endif

#endif