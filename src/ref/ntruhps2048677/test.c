
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <assert.h>

#include <arm_neon.h>

#include "params.h"
#include "poly.h"
#include "poly_NTT.h"
#include "hal.h"

#define ITERATIONS 100
// #define ITERATIONS 1


static void schoolbook_core(uint16_t des[NTRU_N], const uint16_t src1[NTRU_N], const uint16_t src2[NTRU_N]){

    uint16_t tmp[NTRU_N * 2]; 

    for(size_t i = 0; i < NTRU_N * 2; i++){
        tmp[i] = 0;
    }   
    
    for(size_t i = 0; i < NTRU_N; i++){
        for(size_t j = 0; j < NTRU_N; j++){
            tmp[i + j] += src1[i] * src2[j]; 
        }
    }

    for(size_t i = 0; i < NTRU_N; i++){
        des[i] = tmp[i];
    }

    for(size_t i = NTRU_N; i < NTRU_N * 2; i++){
        des[i - NTRU_N] += tmp[i];
    }

}

static void schoolbook(poly *des, const poly *src1, const poly *src2){
    schoolbook_core(des->coeffs, (const uint16_t*)src1->coeffs, (const uint16_t*)src2->coeffs);
}

// schoolbook multilipcation to verify the result of poly_Rq_mul_small
int main(void){

    poly ref, src1, src2;
    poly res;

    for(size_t i = 0; i < ITERATIONS; i++){
        printf("ITERATION %d\n", i);
        for(size_t j = 0; j < NTRU_N; j++){
            src1.coeffs[j] = rand() % NTRU_Q;
            src2.coeffs[j] = rand() % 3;
            // src2.coeffs[j] = rand() % NTRU_Q;
            if(src2.coeffs[j] == 2){
                src2.coeffs[j] = NTRU_Q - 1;
            }
        }
        for(size_t j = NTRU_N; j < POLY_N; j++){
            src1.coeffs[j] = 0;
            src2.coeffs[j] = 0;
        }
        
        // printf("[*] src1: \n");
        // for (size_t j = 0; j < POLY_N; j++) {
        //     printf("%d ", src1.coeffs[j]);
        // }
        // printf("\n[*] src2: \n");
        // for (size_t j = 0; j < POLY_N; j++) {
        //     printf("%d ", src2.coeffs[j]);
        // }
        // printf("\n");

        schoolbook(&ref, &src1, &src2);
        poly_Rq_mul_small(&res, &src1, &src2);
        // poly_NTT(&res, &src1, &src2);

        // printf("[*] result of schoolbook: \n");
        for(size_t j = 0; j < NTRU_N; j++){
            assert(MODQ(ref.coeffs[j]) == MODQ(res.coeffs[j]));
            // printf("%d ", MODQ(ref.coeffs[j]));
        }
        // printf("\n");

        // printf("[*] result of poly_NTT: \n");
        // for(size_t j = 0; j < NTRU_N; j++){
        //     printf("%d ", MODQ(res.coeffs[j]));
        // }
        // printf("\n");

    }
    printf("poly_Rq_mul_small passed!\n");
    // printf("poly_NTT passed!\n");

}


