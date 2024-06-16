#ifndef POLY_H
#define POLY_H

#include <arm_neon.h>
#include "params.h"
#include "poly.h"
#include "neon_batch_multiplication.h"
#include "neon_matrix_transpose.h"

#define POLY_N 1024
#define NTRU_N_PAD 720
#define NTRU_N_32 704

typedef struct {
    uint16_t coeffs[POLY_N];
} poly;


void poly_Rq_mul(poly *r, poly *a, poly *b);

#endif