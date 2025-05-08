#ifndef CONSTS_K64N80_H
#define CONSTS_K64N80_H

#include <vector>


/* --- Convolutional Code Parameters --- */
inline constexpr int k = 1;                    /* Number of input bits */
inline constexpr int n = 1;                    /* Number of output bits */
inline constexpr int K = 64;                   /* Number of input bits */
inline constexpr int N = 80;                  /* Number of output bits */
inline constexpr int V = 8;                    /* Number of memory elements */
inline constexpr int M = 8;                   /* Degree of CRC - 1 */
inline constexpr unsigned int CRC = 0x18B;    /* CRC polynomial */
inline const std::vector<int> numerators = {205};  /* in octal */

inline constexpr int NUM_INFO_BITS = 64;       /* Number of information bits */
inline constexpr int NUM_CODED_SYMBOLS = 80;  /* Number of coded symbols */
inline constexpr int numPairs = 1 << (M + V);  /* Number of syndromes */
/* range of offsets in squared Euclidean distance */
inline constexpr int lowDist = 12;  
inline constexpr int highDist = 24;  // hamming distance 6

inline const std::vector<int> PUNCTURING_INDICES = {};

/* --- List Decoder Parameters --- */
inline constexpr int MAX_LISTSIZE = 256;       /* Maximum list size */
inline constexpr char CODE_TYPE = 'Z';     /* Select code type: {Z: zero-terminated CC, T: tail-biting CC} */
inline constexpr char DECODING_RULE = 'N';     /* Decoding rule: {P: projected, N: non-projected}*/
inline constexpr char ERROR_RUN_TYPE = 'U';    /* Accumulate to which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
inline constexpr int MAX_ERRORS = 20;           /* Maximum number of errors */
inline constexpr bool NOISELESS = false;        /* Noiseless simulation */
inline constexpr int LOGGING_ITERS = 1000;      /* Logging Interval*/
inline constexpr int BASE_SEED = 42;            /* Fixed base seed for simulation */

#endif