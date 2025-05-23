#ifndef CONSTS_K64N256_H
#define CONSTS_K64N256_H

#include <vector>


/* --- Convolutional Code Parameters --- */
inline constexpr int k = 1;                    /* Number of input bits */
inline constexpr int n = 4;                    /* Number of output bits */
inline constexpr int K = 64;                   /* Number of input bits */
inline constexpr int N = 256;                  /* Number of output bits */
inline constexpr int V = 8;                    /* Number of memory elements */
inline constexpr int M = 12;                   /* Degree of CRC - 1 */
inline constexpr unsigned int CRC = 0x1005;    /* CRC polynomial */
inline const std::vector<int> numerators = {533, 575, 647, 711};  /* in octal */

inline constexpr int NUM_INFO_BITS = 64;       /* Number of information bits */
inline constexpr int NUM_CODED_SYMBOLS = 256;  /* Number of coded symbols */
inline constexpr int numPairs = 1 << (M + V);  /* Number of syndromes */
/* range of offsets in squared Euclidean distance */
inline constexpr int lowDist = 12;  
inline constexpr int highDist = 100;  // hamming distance 25

inline const std::vector<int> PUNCTURING_INDICES = {4, 10, 21, 24, 31, 37, 42, 48, 59, 62, 69, 75, 80, 86, 97, 100, 107, 113, 118, 124, 135, 138, 145, 151, 156, 162, 173, 176, 183, 189, 194, 200, 211, 214, 221, 227, 232, 238, 249, 252, 259, 265, 270, 276, 287, 290, 297, 303};

/* --- List Decoder Parameters --- */
inline constexpr int MAX_LISTSIZE = 16384;       /* Maximum list size */
inline constexpr char CODE_TYPE = 'T';     /* Select code type: {Z: zero-terminated CC, T: tail-biting CC} */
inline constexpr char DECODING_RULE = 'N';     /* Decoding rule: {P: projected, N: non-projected}*/
inline constexpr char ERROR_RUN_TYPE = 'U';    /* Accumulate to which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
inline constexpr int MAX_ERRORS = 20;           /* Maximum number of errors */
inline constexpr bool NOISELESS = false;        /* Noiseless simulation */
inline constexpr int LOGGING_ITERS = 1000;      /* Logging Interval*/
inline constexpr int BASE_SEED = 42;            /* Fixed base seed for simulation */

#endif