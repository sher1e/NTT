#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct {
    int16_t coeffs[SCHEME_N];
} poly;

void poly_add(poly *r, const poly *a, const poly *b);

void poly_sub(poly *r, const poly *a, const poly *b);

void poly_ntt(poly *r);

void poly_invntt(poly *r);

void poly_basemul(poly *r, const poly *a, const poly *b);

#endif
