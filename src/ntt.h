#ifndef NTT_H
#define NTT_H
#include "params.h"
#include <stdint.h>

void ntt(int32_t a[N]);
void intt(int32_t a[N]);

void invntt_tomont(int32_t a[N]);
void invntt_tomont_merged(int32_t a[N]);
void PQCLEAN_DILITHIUM2_CLEAN_ntt(int32_t a[N]);

#endif
