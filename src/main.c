#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "util.h"
#include "ntt.h"
#include "reduce.h"

#define LOOP 100

extern int ntt_3_layer_merged(int32_t a[256]);
extern uint64_t getcycles();
extern int ntt_1layer(int a[256]);

int main() {

    clock_t start = 0, end = 0, cycle = 0;
    int32_t a[256] = {0,};
    for(int i=0; i<256; i++){
        a[i] = i;
    }
    // srai, snez, lui -> 3cc

        ntt_asm_only_butterfly(a);
        ntt_asm_only_butterfly(a);
        ntt_asm_only_butterfly(a);
        ntt_asm_only_butterfly(a);
        ntt_asm_only_butterfly(a);

    for(int i = 0; i < LOOP; i++){
        start = getcycles();
        ntt_asm_only_butterfly(a);
        end = getcycles();
        cycle += (end - start);
    }
    printf("asm bufferfly: %d\n", cycle / LOOP);
    // for(int i = 0; i < 256; i++){
    //     printf("%08x ", a[i]);
    //     if((i+15)%16 == 0){
    //         printf("\n");
    //     }
    // }printf("\n\n\n");

    for(int i=0; i<256; i++){
        a[i] = i;
    }
    for(int i = 0; i < LOOP; i++){
        start = getcycles();
        PQCLEAN_DILITHIUM3_CLEAN_ntt(a);
        end = getcycles();
        cycle += (end - start);
    }
    printf("PQCLEAN: %d\n", cycle / LOOP);

    // for(int i = 0; i < 256; i++){
    //     printf("%08x ", a[i]);
    //     if((i+15)%16 == 0){
    //         printf("\n");
    //     }
    // }printf("\n\n\n");
    return 0;

}