#include "reduce.h"
#include "params.h"

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q,
*              where R=2^16 (see value of rlog)
*
* Arguments:   - int32_t a: input unsigned integer to be reduced; has to be in {-q*2^15,...,q*2^15}
*
* Returns:     signed integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int16_t montgomery_reduce(int32_t a) {
    int32_t t;
    int16_t u;

    u = a * SCHEME_QINV;
    t = (int32_t) u * SCHEME_Q;
    t = a - t;
    t >>= 16;
    return t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q,
*
* Arguments:   - int16_t a: input integer to be reduced; has to be in {-2^15,...,2^15}
*
* Returns:     integer in {0,...,q} congruent to a modulo q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
    int32_t t;
    const int32_t v = (1U << 26) / SCHEME_Q + 1;

    t = v * a;
    t >>= 26;
    t *= SCHEME_Q;
    return a - t;
}
