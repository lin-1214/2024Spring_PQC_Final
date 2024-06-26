
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>

#include <arm_neon.h>

#include "params.h"
#include "poly.h"
#include "poly_NTT.h"
#include "neon_poly_mul.h"
#include "hal.h"

#define ITERATIONS 100000
#define NEON


uint64_t t0, t1;
uint64_t times[ITERATIONS];

static int cmp_uint64(const void *a, const void *b){
    return ((*((const uint64_t*)a)) - ((*((const uint64_t*)b))));
}


int main(void){

    poly a, b, c;
    // for (size_t i = 0; i < POLY_N; i++){
    //     printf("%d ", a[i]);
    // }

    for(size_t i = 0; i < ITERATIONS; i++){
        t0 = hal_get_time();
        
        #ifdef poly_Rq_mul_small
        poly_Rq_mul_small(&a, &b, &c);
        #endif
        
        #ifdef POLY_NTT
        poly_NTT(&a, &b, &c);
        #endif

        #ifdef NEON
        poly_neon(&a, &b, &c);
        #endif

        t1 = hal_get_time();
        times[i] = t1 - t0;
    }
    qsort(times, ITERATIONS, sizeof(uint64_t), cmp_uint64);
    printf("poly_mul: %ld\n", times[ITERATIONS >> 1]);

}






