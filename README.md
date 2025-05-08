TBCC and ZTCC Decoder Simulation

This repository contains simulation tools and decoder implementations for various convolutional codes, including:
- Zero-Terminated Convolutional Code (ZTCC)
- Tail-Biting Convolutional Codes (TBCC) with different lengths and puncturing patterns

## Code Configurations ##
Supported convolutional code configurations:
- (K=64, N=80) – ZTCC
- (K=64, N=128) – TBCC
- (K=64, N=256) – TBCC

Each configuration has its own header file containing constants:
- consts_K64N80.h
- consts_K64N128.h
- consts_K64N256.h

## Build Instructions ##
To compile for a specific configuration, define the appropriate macro during compilation:

# For (K=64, N=80) ZTCC
make CONFIG=K64N80

# For (K=64, N=128) TBCC
make CONFIG=K64N128

# For (K=64, N=256) TBCC
make CONFIG=K64N256

This ensures only one code configuration is active at a time, and the correct constants header is included via preprocessor logic in const.h.

To clean up build files:
make clean

# Running Simulations
After compiling, run the corresponding test using:
./run_test -t 8 (if CONFIG=K64N80 is used)
./run_test -t 4 (if CONFIG=K64N128 is used)
./run_test -t 0 (if CONFIG=K64N256 is used)

# To change EsN0 values
In run_test.cpp, modify the esno field in test_point contest[N_TEST] between lines 30-43.
Currently we set the EsN0 values to achieve FER=1e-3.

Notes
- Only one CONFIG=K64Nxx flag should be defined at a time.
- The decoder uses soft information (LLRs) and supports both ML and expurgated decoding.
- Requires C++17 or newer.
