#include <string.h>
#include "poly.h"
#include "ntt.h"
#include "reduce.h"

extern int16_t zetas_128[];

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b) {
    int i;
    for (i = 0; i < SCHEME_N; i++)
        r->coeffs[i] = (a->coeffs[i] + b->coeffs[i]);
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b) {
    int i;
    for (i = 0; i < SCHEME_N; i++)
        r->coeffs[i] = (a->coeffs[i] - b->coeffs[i]);
}

/*************************************************
* Name:        poly_ntt
*
* Description: Forward NTT transform of a polynomial in place
*              Input is assumed to have coefficients in bitreversed order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r) {
    int i;
    for (i = 0; i < SCHEME_NTT_POLY; i++)
        ntt(r->coeffs + i * SCHEME_NTT_LENGTH, zetas_128);
}

/*************************************************
* Name:        poly_invntt
*
* Description: Inverse NTT transform of a polynomial in place
*              Input is assumed to have coefficients in normal order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_invntt(poly *r) {
    int i;
    for (i = 0; i < SCHEME_NTT_POLY; i++)
        invntt(r->coeffs + i * SCHEME_NTT_LENGTH, zetas_128);
}

/*************************************************
* Name:        poly_basemul
*
* Description: Multiply two polynomials in NTT domain.
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul(poly *r, const poly *a, const poly *b) {
    int i, k = SCHEME_NTT_LENGTH / 2;
    for (i = 0; i < SCHEME_NTT_LENGTH; i += 2) {
        basemul(r->coeffs + i,
                a->coeffs + i,
                b->coeffs + i,
                zetas_128[k]);

        basemul(r->coeffs + i + 1,
                a->coeffs + i + 1,
                b->coeffs + i + 1,
                -zetas_128[k++]);
    }
}
