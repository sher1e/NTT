#include "inttypes.h"
#include "ntt.h"
#include "params.h"
#include "reduce.h"

/*************************************************
* Name:        ntt
*
* Description: Computes number-theoretic transform (NTT) of
*              a polynomial in place; inputs assumed to be in
*              normal order, output in bitreversed order
*
* Arguments:   - int16_t * a:           pointer to in/output polynomial
*              - const int16_t* gammas: pointer to input powers of root of unity gamma;
*                                       assumed to be in Montgomery domain
**************************************************/
void ntt(int16_t *a, const int16_t *gammas) {
    int k, len, start, j;
    int16_t G, temp;

    k = 1;
    for (len = SCHEME_NTT_LENGTH / 2; len >= 1; len >>= 1) {
        for (start = 0; start < SCHEME_NTT_LENGTH; start = j + len) {
            G = gammas[k++];
            for (j = start; j < start + len; ++j) {
                temp = montgomery_reduce(G * a[j + len]);
                a[j + len] = a[j] - temp; // Omit reduction (be lazy)
                a[j] = a[j] + temp; // Omit reduction (be lazy)
            }
        }
    }
}

/*************************************************
* Name:        invntt
*
* Description: Computes inverse number-theoretic transform (INTT) of
*              a polynomial in place; inputs assumed to be in
*              bitreversed order, output in normal order
*
* Arguments:   - int16_t * a:           pointer to in/output polynomial
*              - const int16_t* gammas: pointer to input powers of root of unity gamma;
*                                       assumed to be in Montgomery domain
**************************************************/
void invntt(int16_t *a, const int16_t *gammas) {
    int k, len, start, j;
    int16_t temp, G;
    const int16_t f = 1441; // mont^2/128

    k = 127;
    for (len = 1; len <= SCHEME_NTT_LENGTH / 2; len <<= 1) {
        for (start = 0; start < SCHEME_NTT_LENGTH; start = j + len) {
            G = -gammas[k--];
            for (j = start; j < start + len; ++j) {
                temp = a[j];
                a[j] = barrett_reduce(temp + a[j + len]);
                a[j + len] = montgomery_reduce(G * (temp - a[j + len]));
            }
        }
    }

    for (j = 0; j < SCHEME_NTT_LENGTH; j++)
        a[j] = montgomery_reduce(a[j] * f);
}

#if USE_SCHOOLBOOK_MULTIPLICATION

/*************************************************
* Name:        basemul
*
* Description: Multiply two polynomials in Zq[X]/(X^SCHEME_NTT_POLY-gamma) 
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void basemul(int16_t *r, const int16_t *a, const int16_t *b, int16_t gamma)
{
  int i,j;
  
  for(i=0;i<SCHEME_NTT_POLY;i++)
  {
    r[i*SCHEME_NTT_LENGTH] = 0;
    for(j=i+1;j<SCHEME_NTT_POLY;j++)
      r[i*SCHEME_NTT_LENGTH] += 
        montgomery_reduce(a[j*SCHEME_NTT_LENGTH] * b[(SCHEME_NTT_POLY+i-j)*SCHEME_NTT_LENGTH]);
    if (j!=i+1)
      r[i*SCHEME_NTT_LENGTH] = montgomery_reduce(r[i*SCHEME_NTT_LENGTH] * gamma);
    for(j=0;j<i+1;j++)
      r[i*SCHEME_NTT_LENGTH] += 
        montgomery_reduce(a[j*SCHEME_NTT_LENGTH] * b[(i-j)*SCHEME_NTT_LENGTH]);
  }
}

#elif USE_KARATSUBA_MULTIPLICATION

#define NTT_LENGTH(x) (x*SCHEME_NTT_LENGTH)

#define CALC_D(a, b, x, y, d) (montgomery_reduce((a[NTT_LENGTH(x)]+a[NTT_LENGTH(y)])* \
                                             (b[NTT_LENGTH(x)]+b[NTT_LENGTH(y)]))-d[x]-d[y])

/*************************************************
* Name:        basemul
*
* Description: Multiply two polynomials in Zq[X]/(X^SCHEME_NTT_POLY-gamma) 
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void basemul(int16_t *r, const int16_t *a, const int16_t *b, int16_t gamma) {
    int i;

    int16_t d[SCHEME_NTT_POLY];
    for (i = 0; i < SCHEME_NTT_POLY; i++)
        d[i] = montgomery_reduce(a[i * 128] * b[i * 128]);

#if   (SCHEME_N == 512)

    r[NTT_LENGTH(0)] = d[0] +
                       montgomery_reduce(
                               (CALC_D(a, b, 1, 3, d) + d[2]) * gamma);

    r[NTT_LENGTH(1)] = CALC_D(a, b, 0, 1, d) +
                       montgomery_reduce(
                               CALC_D(a, b, 2, 3, d) * gamma);

    r[NTT_LENGTH(2)] = barrett_reduce(CALC_D(a, b, 0, 2, d) + d[1] +
                                      montgomery_reduce(d[3] * gamma));

    r[NTT_LENGTH(3)] = barrett_reduce(CALC_D(a, b, 1, 2, d) +
                                      CALC_D(a, b, 0, 3, d));

#else
#error "SCHEME_N must be either 512, 768 or 1024"
#endif
}

#else
#error "USE_SCHOOLBOOK_MULTIPLICATION or USE_KARATSUBA_MULTIPLICATION must be defined."
#endif
