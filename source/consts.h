#ifndef CONST_H
#define CONST_H

#include <vector>

// initialize
/* --- Convolutional Code Parameters --- */
extern const int k;                    /* Number of input bits */
extern const int n;                    /* Number of output bits */
extern const int K;                   /* Number of input bits */
extern const int N;                  /* Number of output bits */
extern const int V;                    /* Number of memory elements */
extern const int M;                   /* Degree of CRC - 1 */
extern const unsigned int CRC;    /* CRC polynomial */
extern const std::vector<int> numerators;

extern const int NUM_INFO_BITS;       /* Number of information bits */
extern const int NUM_CODED_SYMBOLS;  /* Number of coded symbols */
extern const int numPairs;   /* Number of syndromes */

/* range of offsets in squared Euclidean distance */
extern const int lowDist;  
extern const int highDist;

/* to store offset tables */
using Row = std::vector<int>;
using RowList = std::vector<Row>;
extern RowList messageList;
extern RowList codewordList;

extern const std::vector<int> PUNCTURING_INDICES;

/* --- List Decoder Parameters --- */
extern const int MAX_LISTSIZE;       /* Maximum list size */
extern const char CODE_TYPE;     /* Select code type: {ZT: zero-terminated CC, TB: tail-biting CC} */
extern const char DECODING_RULE;     /* Decoding rule: {P: projected, N: non-projected}*/
extern const char ERROR_RUN_TYPE;    /* Accumulate to which type of error: {U: undetected, T: total}*/

/* --- Simulation Parameters --- */
extern const int MAX_ERRORS;           /* Maximum number of errors */
extern const bool NOISELESS;        /* Noiseless simulation */
extern const int LOGGING_ITERS;      /* Logging Interval*/
extern const int BASE_SEED;            /* Fixed base seed for simulation */


// import actual code-specific parameters
#if defined(CONFIG_K64N80)
#include "consts_K64N80.h"
#elif defined(CONFIG_K64N128)
#include "consts_K64N128.h"
#elif defined(CONFIG_K64N256)
#include "consts_K64N256.h"
#endif

#endif

