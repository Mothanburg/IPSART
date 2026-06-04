## IPSART的MEX加速模块

部分算法在MATLAB中性能极差甚至不可用，因此基于mex接口开发了该模块。该模块依赖的第三方库包括Eigen（提供矩阵运算）和OpenCL（提供GPU加速），使用`../bootstrap/mex_compile.m`构建。