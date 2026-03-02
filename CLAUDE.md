# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

IPSART is a MATLAB-based library for InSAR (Interferometric Synthetic Aperture Radar) and PolSAR (Polarimetric SAR) image processing. The project combines pure MATLAB implementations with performance-critical C++ MEX functions.

## Initialization

Run `InitIPSART.m` in MATLAB to initialize the library:
```matlab
run('InitIPSART.m')
```
This adds all package directories to the MATLAB path and enables native MEX functions from the `bin` folder. It also starts a thread-based parallel pool.

## MATLAB Package Structure

The main packages (imported by InitIPSART) are:
- **Basic** - Core image processing utilities (Averagelook, HistStretch)
- **Data** - Data I/O and handling
- **SAR** - Generic SAR processing functions
- **PolSAR** - Polarimetric SAR decomposition and analysis
- **InSAR** - Interferometric SAR processing
- **Geometric** - Georeferencing and geometric transformations
- **Utils** - Utility functions (UnaryOp, ComposeOp, Hex2RGB, Hist2d)

## MEX Functions (C++ Acceleration)

Performance-critical functions are implemented in C++ and compiled as MATLAB MEX binaries located in `mex/source/`:

| MEX Function | Purpose |
|--------------|---------|
| CloudePottier | Cloude-Pottier decomposition |
| Yamaguchi | Yamaguchi decomposition |
| G4U | Four-component decomposition |
| Multilook | Multi-looking processing |
| RefinedLeeFilter | Refined Lee filter |

The MEX bridge is in `mex/source/IPSART.cpp` - it dispatches method calls based on the first input argument.

### Building MEX Functions

From MATLAB, run the bootstrap script:
```matlab
cd bootstrap
run mex_compile.m
```

This compiles all C++ sources in `mex/source/*.cpp` and outputs to the `bin` directory.

Requirements:
- MATLAB
- A C++ compiler on Windows (now only available for MSVC compiler)
- OpenCL SDK
- Eigen

The `mex/xmake.lua` file provides xmake configuration just for IntelliSense, it doesn't do anything on building.

### OpenCL Kernels

OpenCL kernel source files in `mex/opencl/*.cl` are embedded as header files during compilation via the bootstrap script.

## Development Notes

- MATLAB functions use the modern arguments block syntax for input validation
- C++ code uses C++latest, MATLAB Data API, and Eigen for matrix operations
- The project uses `parpool("Threads")` for parallel processing
- Build outputs go to `bin/` which is gitignored
- Edit MATLAB path configuration in `InitIPSART.m` if adding new packages
