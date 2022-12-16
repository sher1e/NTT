#ifndef PARAMS_H
#define PARAMS_H

#ifndef SCHEME_N
#define SCHEME_N 512
#endif

#define SCHEME_NTT_LENGTH 128  // set NTT length

#if  (SCHEME_N == 512)
#define SCHEME_Q    3329
#define SCHEME_QINV -3327   // inverse_mod(p,2^16)
#define SCHEME_NTT_POLY 4   // set NTT polynomial number
#endif

//#define USE_SCHOOLBOOK_MULTIPLICATION 1
#if !(defined(USE_KARATSUBA_MULTIPLICATION) || defined(USE_SCHOOLBOOK_MULTIPLICATION))
#define USE_KARATSUBA_MULTIPLICATION 1
#endif

#endif
