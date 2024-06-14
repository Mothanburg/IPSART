#ifndef GARS_API
#define GARS_API
#endif

struct RawPolMat3 {
  long height, width;
  double *m11, *m22, *m33;
  double *m12_r, *m13_r, *m23_r;
  double *m12_i, *m13_i, *m23_i;
};

GARS_API void RefinedLeeFilter(struct RawPolMat3 in, struct RawPolMat3 out,
                               int nLooks);

GARS_API void CloudePottierT3(long height, long width, double *m11, double *m22,
                              double *m33, double *m12_r, double *m13_r,
                              double *m23_r, double *m12_i, double *m13_i,
                              double *m23_i, double *outH, double *outAlpha,
                              double *outA);