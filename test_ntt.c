#include <stdio.h>
#include <string.h>
#include "poly.h"

#define NTESTS 1000

void print_poly(char *s, poly *a, int len) {
    int i;
    printf("%s first %d elements: \n", s, len);
    for (i = 0; i < len; i++) {
        printf("%d, ", a->coeffs[i]);
        if (i % 16 == 15) printf("\n");
    }
    printf("\n\n");
}

void test_ntt_128() {
    poly ahat, shat, ahat_shat;

    for (int i = 0; i < SCHEME_N / 4; i++) {
        ahat.coeffs[i] = 1;
        shat.coeffs[i] = 1;
    }
    for (int i = SCHEME_N / 4; i < SCHEME_N; i++) {
        ahat.coeffs[i] = 0;
        shat.coeffs[i] = 0;
    }

    poly_ntt(&ahat);
    poly_ntt(&shat);
    poly_basemul(&ahat_shat, &shat, &ahat);
    poly_invntt(&ahat);
    poly_invntt(&ahat_shat);
    print_poly("ahat", &ahat, SCHEME_N);
    print_poly("ahat_shat", &ahat_shat, SCHEME_N);
}

int main() {
    test_ntt_128();
    return 0;
}

