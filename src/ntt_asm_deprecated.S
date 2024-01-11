.section .text

# return value is stored in VALH
# Q and QINV must be stored in registers
.macro montgomery_reduce VALH, VALL, TMP1, TMP2, Q, QINV
    mul     \TMP1, \VALL, \QINV     # t <- a * QINV
    mulh    \TMP2, \TMP1, \Q
    mul     \TMP1, \TMP1, \Q
    sub     \VALH, \VALH, \TMP2
    sltu    \TMP1, \VALL, \TMP1
    add     \VALH, \VALH, \TMP1
.endm


# Cooley-Tucky Butterfly
.macro ct_butterfly_core C0, C1, THI, TLO, Q, QINV, TMP0, TMP1, ZETA
    mulh    \THI, \C1, \ZETA        # t3: high / t4: low
    mul     \TLO, \C1, \ZETA        # zeta * a[j + len]

    # montgomery reduction
    mul     \TMP0, \TLO, \QINV     # t <- a * QINV
    mulh    \TMP1, \TMP0, \Q
    mul     \TMP0, \TMP0, \Q
    sub     \THI,  \THI,  \TMP1
    sltu    \TMP0, \TLO,  \TMP0
    add     \THI,  \THI,  \TMP0

    sub     \C1, \C0, \THI
    add     \C0, \C0, \THI
.endm

.global ct_butterfly
.align 2
.type ct_butterfly, @function
ct_butterfly:
    # load coeffs
    slli    a1, a1, 2
    add     a1, a0, a1
    slli    a2, a2, 2
    add     a2, a0, a2
    lw      a6, 0(a1)
    lw      a7, 0(a2)

    ct_butterfly_core a6, a7, t3, t4, a3, a4, t5, t6, a5

    # store coeffs
    sw      a6, 0(a1)
    sw      a7, 0(a2)
    ret



.macro gs_butterfly_core A0, C0, C1, OFFSET0, OFFSET1, THI, TLO, Q, QINV, TMP0, TMP1, ZETA
    lw      \TMP0, \OFFSET0(\A0)    # t <- a[j]
    add     \C0, \TMP0, \C1
    sub     \TMP0, \TMP0, \C1

    mulh    \C1, \TMP0, \ZETA
    mul     \TLO, \TMP0, \ZETA
    montgomery_reduce \C1, \TLO, \TMP0, \TMP1, \Q, \QINV
.endm

# load coeffs into registers from array
.macro load_coeffs A0, C0, C1, C2, C3, C4, C5, C6, C7, OFF0, OFF1, OFF2, OFF3, OFF4, OFF5, OFF6, OFF7
    lw      \C0, \OFF0(\A0)
    lw      \C1, \OFF1(\A0)
    lw      \C2, \OFF2(\A0)
    lw      \C3, \OFF3(\A0)
    lw      \C4, \OFF4(\A0)
    lw      \C5, \OFF5(\A0)
    lw      \C6, \OFF6(\A0)
    lw      \C7, \OFF7(\A0)
.endm
 
# store coeffs into array from registers
.macro end_butterfly A0, C0, C1, C2, C3, C4, C5, C6, C7, OFF0, OFF1, OFF2, OFF3, OFF4, OFF5, OFF6, OFF7
    sw      \C0, \OFF0(\A0)
    sw      \C1, \OFF1(\A0)
    sw      \C2, \OFF2(\A0)
    sw      \C3, \OFF3(\A0)
    sw      \C4, \OFF4(\A0)
    sw      \C5, \OFF5(\A0)
    sw      \C6, \OFF6(\A0)
    sw      \C7, \OFF7(\A0)
.endm

.macro store_and_load_coeffs A0, C0, C1, C2, C3, C4, C5, C6, C7, OFF0, OFF1, OFF2, OFF3, OFF4, OFF5, OFF6, OFF7, OFFNEW0, OFFNEW1, OFFNEW2, OFFNEW3, OFFNEW4, OFFNEW5, OFFNEW6, OFFNEW7
    sw      \C0, \OFF0(\A0)
    sw      \C1, \OFF1(\A0)
    lw      \C0, \OFFNEW0(\A0)
    sw      \C2, \OFF2(\A0)
    lw      \C1, \OFFNEW1(\A0)
    sw      \C3, \OFF3(\A0)

    lw      \C2, \OFFNEW2(\A0)
    sw      \C4, \OFF4(\A0)
    lw      \C3, \OFFNEW3(\A0)

    sw      \C5, \OFF5(\A0)
    lw      \C4, \OFFNEW4(\A0)
    sw      \C6, \OFF6(\A0)
    lw      \C5, \OFFNEW5(\A0)
    
    sw      \C7, \OFF7(\A0)
    lw      \C6, \OFFNEW6(\A0)
    lw      \C7, \OFFNEW7(\A0)
.endm

# 3 layer merged CT Butterfly
# used to merge 1~3, 4~6 layer
.macro ct_butterfly_3layer A0, C0, C1, C2, C3, C4, C5, C6, C7, OFFSET0, OFFSET1, OFFSET2, OFFSET3, OFFSET4, OFFSET5, OFFSET6, OFFSET7 , THI, TLO, Q, QINV, TMP0, TMP1, Z0, Z1, Z2, Z3, Z4, Z5, Z6
    load_coeffs \A0, \C0, \C1, \C2, \C3, \C4, \C5, \C6, \C7, \OFFSET0, \OFFSET1, \OFFSET2, \OFFSET3, \OFFSET4, \OFFSET5, \OFFSET6, \OFFSET7
    # layer 1 / 4
    ct_butterfly_core \C0, \C4, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z0
    ct_butterfly_core \C1, \C5, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z0
    ct_butterfly_core \C2, \C6, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z0
    ct_butterfly_core \C3, \C7, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z0
    # layer 2 / 5
    ct_butterfly_core \C0, \C2, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z1
    ct_butterfly_core \C1, \C3, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z1
    ct_butterfly_core \C4, \C6, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z2
    ct_butterfly_core \C5, \C7, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z2
    # layer 3 / 6
    ct_butterfly_core \C0, \C1, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z3
    ct_butterfly_core \C2, \C3, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z4
    ct_butterfly_core \C4, \C5, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z5
    ct_butterfly_core \C6, \C7, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \Z6

    end_butterfly \A0, \C0, \C1, \C2, \C3, \C4, \C5, \C6, \C7, \OFFSET0, \OFFSET1, \OFFSET2, \OFFSET3, \OFFSET4, \OFFSET5, \OFFSET6, \OFFSET7
.endm

    # store_and_load_coeffs \A0, \C0, \C1, \C2, \C3, \C4, \C5, \C6, \C7, \OFFSET0, \OFFSET1, \OFFSET2, \OFFSET3, \OFFSET4, \OFFSET5, \OFFSET6, \OFFSET7, \OFFNEW0, \OFFNEW1, \OFFNEW2, \OFFNEW3, \OFFNEW4, \OFFNEW5, \OFFNEW6, \OFFNEW7 

# set zeta values
.macro set_zetas Z0, Z1, Z2, Z3, Z4, Z5, Z6
    li  a1, \Z0
    li  a2, \Z1
    li  a3, \Z2
    li  a4, \Z3
    li  a5, \Z4
    li  a6, \Z5
    li  a7, \Z6
.endm

.macro set_zetas_new ADDRZETA, OFFSET1, OFFSET2, OFFSET3, OFFSET4, OFFSET5, OFFSET6, OFFSET7
    lw  a1, \OFFSET1(\ADDRZETA)
    lw  a2, \OFFSET2(\ADDRZETA)
    lw  a3, \OFFSET3(\ADDRZETA)
    lw  a4, \OFFSET4(\ADDRZETA)
    lw  a5, \OFFSET5(\ADDRZETA)
    lw  a6, \OFFSET6(\ADDRZETA)
    lw  a7, \OFFSET7(\ADDRZETA)
.endm

# last 7~8 layer
.macro ntt_last_layer A0, C0, C1, C2, C3, OFFSET0, OFFSET1, OFFSET2, OFFSET3, THI, TLO, Q, QINV, TMP0, TMP1, ZREG0, ZREG1, ZREG2, Z0, Z1, Z2
    lw \C0, \OFFSET0(\A0)
    lw \C1, \OFFSET1(\A0)
    lw \C2, \OFFSET2(\A0)
    lw \C3, \OFFSET3(\A0)

    li \ZREG0, \Z0
    li \ZREG1, \Z1
    li \ZREG2, \Z2

    # layer 7
    ct_butterfly_core \C0, \C2, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \ZREG0
    ct_butterfly_core \C1, \C3, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \ZREG0
    # layer 8
    ct_butterfly_core \C0, \C1, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \ZREG1
    ct_butterfly_core \C2, \C3, \THI, \TLO, \Q, \QINV, \TMP0, \TMP1, \ZREG2

    sw \C0, \OFFSET0(\A0)
    sw \C1, \OFFSET1(\A0)
    sw \C2, \OFFSET2(\A0)
    sw \C3, \OFFSET3(\A0)
.endm

.macro ntt_3_layer_merged_new OFF0, OFF1, OFF2, OFF3, OFF4, OFF5, OFF6, OFF7
    // load coeffs
    lw      s0, \OFF0(a0)
    lw      s1, \OFF1(a0)
    lw      s2, \OFF2(a0)
    lw      s3, \OFF3(a0)
    lw      s4, \OFF4(a0)
    lw      s5, \OFF5(a0)
    lw      s6, \OFF6(a0)
    lw      s7, \OFF7(a0)

    // CT butterfly - 0, 4
    mulh    t2, s4, a2        # t3: high / t4: low
    mul     t3, s4, a2        # zeta * a[j + len]

    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0

    sub     s4, s0, t2
    add     s0, s0, t2

    // CT butterfly - 1, 5
    mulh    t2, s5, a2        # t3: high / t4: low
    mul     t3, s5, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s5, s1, t2
    add     s1, s1, t2

    // CT butterfly - 2, 6
    mulh    t2, s6, a2        # t3: high / t4: low
    mul     t3, s6, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s6, s2, t2
    add     s2, s2, t2

    // CT butterfly - 3, 7
    mulh    t2, s7, a2        # t3: high / t4: low
    mul     t3, s7, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s7, s3, t2
    add     s3, s3, t2

    // CT butterfly - 0, 2
    mulh    t2, s2, a2        # t3: high / t4: low
    mul     t3, s2, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s2, s0, t2
    add     s0, s0, t2

    // CT butterfly - 1, 3
    mulh    t2, s3, a2        # t3: high / t4: low
    mul     t3, s3, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s3, s1, t2
    add     s1, s1, t2

    // CT butterfly - 4, 6
    mulh    t2, s6, a2        # t3: high / t4: low
    mul     t3, s6, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s6, s4, t2
    add     s4, s4, t2

    // CT butterfly - 5, 7
    mulh    t2, s7, a2        # t3: high / t4: low
    mul     t3, s7, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s7, s5, t2
    add     s5, s5, t2

    // CT butterfly - 0, 1
    mulh    t2, s1, a2        # t3: high / t4: low
    mul     t3, s1, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s1, s0, t2
    add     s0, s0, t2

    // CT butterfly - 2, 3
    mulh    t2, s3, a2        # t3: high / t4: low
    mul     t3, s3, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s3, s2, t2
    add     s2, s2, t2

    // CT butterfly - 4,5
    mulh    t2, s5, a2        # t3: high / t4: low
    mul     t3, s5, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s5, s4, t2
    add     s4, s4, t2

    // CT butterfly - 6, 7
    mulh    t2, s7, a2        # t3: high / t4: low
    mul     t3, s7, a2        # zeta * a[j + len]
    mul     t0, t3, t6    
    mulh    t1, t0, t5
    mul     t0, t0, t5
    sub     t2, t2, t1
    sltu    t0, t3, t0
    add     t2, t2, t0
    sub     s7, s6, t2
    add     s6, s6, t2

    // store coeffs
    sw      s0, \OFF0(a0)
    sw      s1, \OFF1(a0)
    sw      s2, \OFF2(a0)
    sw      s3, \OFF3(a0)
    sw      s4, \OFF4(a0)
    sw      s5, \OFF5(a0)
    sw      s6, \OFF6(a0)
    sw      s7, \OFF7(a0)
.endm

# 3 layer merged NTT
# [Register usage]
# coeffs:             s2 ~ s9 (total 8)
# zetas:              a1 ~ a7 (total 7)
# mul res(high, low): t2,  t3 (for montgomery reduction)
# Q, QINV:            t5,  t6 (not changed)
# temporary:          t0,  t1 (for montgomery reduction)
# address of array:   a0      (passed by function argument)
# total:              22 registers
# Merged 1~3, 4~6, 7~8
# .align 4 
.global ntt_3_layer_merged
.align 2
.type ntt_3_layer_merged, @function
ntt_3_layer_merged:

    addi    sp, sp, -52
    sw      s0, 4(sp)
    sw      s1, 8(sp)
    sw      s2, 12(sp)
    sw      s3, 16(sp)
    sw      s4, 20(sp)
    sw      s5, 24(sp)
    sw      s6, 28(sp)
    sw      s7, 32(sp)
    sw      s8, 36(sp)
    sw      s9, 40(sp)
    sw      s10, 44(sp)
    sw      s11, 48(sp)

    li      t6, 0x03802001        # QINV
    li      t5, 0x007fe001        # Q
    la      t4, ZETA_TABLE
    
    # first 1~3 layer
    # set_zetas_new t4, 4, 8, 12, 16, 20, 24, 28
    set_zetas 25847, -2608894, -518909, 237124, -777960, -876248, 466468
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,   0, 128, 256, 384, 512, 640, 768,  896, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,   4, 132, 260, 388, 516, 644, 772,  900, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,   8, 136, 264, 392, 520, 648, 776,  904, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  12, 140, 268, 396, 524, 652, 780,  908, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  16, 144, 272, 400, 528, 656, 784,  912, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  20, 148, 276, 404, 532, 660, 788,  916, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  24, 152, 280, 408, 536, 664, 792,  920, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  28, 156, 284, 412, 540, 668, 796,  924, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  32, 160, 288, 416, 544, 672, 800,  928, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  36, 164, 292, 420, 548, 676, 804,  932, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  40, 168, 296, 424, 552, 680, 808,  936, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  44, 172, 300, 428, 556, 684, 812,  940, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  48, 176, 304, 432, 560, 688, 816,  944, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  52, 180, 308, 436, 564, 692, 820,  948, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  56, 184, 312, 440, 568, 696, 824,  952, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  60, 188, 316, 444, 572, 700, 828,  956, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  64, 192, 320, 448, 576, 704, 832,  960, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  68, 196, 324, 452, 580, 708, 836,  964, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  72, 200, 328, 456, 584, 712, 840,  968, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  76, 204, 332, 460, 588, 716, 844,  972, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  80, 208, 336, 464, 592, 720, 848,  976, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  84, 212, 340, 468, 596, 724, 852,  980, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  88, 216, 344, 472, 600, 728, 856,  984, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  92, 220, 348, 476, 604, 732, 860,  988, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  96, 224, 352, 480, 608, 736, 864,  992, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 100, 228, 356, 484, 612, 740, 868,  996, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 104, 232, 360, 488, 616, 744, 872, 1000, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 108, 236, 364, 492, 620, 748, 876, 1004, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 112, 240, 368, 496, 624, 752, 880, 1008, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 116, 244, 372, 500, 628, 756, 884, 1012, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 120, 248, 376, 504, 632, 760, 888, 1016, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 124, 252, 380, 508, 636, 764, 892, 1020, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    
    # 4~6 layer
    set_zetas_new t4, 32, 36, 40, 44, 48, 52, 56
    # set_zetas 1826347, 2725464, 1024112, 2706023, 95776, 3077325, 3530437
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  0, 16, 32, 48, 64, 80,  96, 112, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  4, 20, 36, 52, 68, 84, 100, 116, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9,  8, 24, 40, 56, 72, 88, 104, 120, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 12, 28, 44, 60, 76, 92, 108, 124, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7

    set_zetas_new t4, 60, 64, 68, 72, 76, 80, 84
    # set_zetas 2353451, -1079900, 3585928, -1661693, -3592148, -2537516, 3915439
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 128, 144, 160, 176, 192, 208, 224, 240, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 132, 148, 164, 180, 196, 212, 228, 244, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 136, 152, 168, 184, 200, 216, 232, 248, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 140, 156, 172, 188, 204, 220, 236, 252, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    
    set_zetas_new t4, 88, 92, 96, 100, 104, 108, 112
    # set_zetas -359251, -549488, -1119584, -3861115, -3043716, 3574422, -2867647
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 256, 272, 288, 304, 320, 336, 352, 368, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 260, 276, 292, 308, 324, 340, 356, 372, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 264, 280, 296, 312, 328, 344, 360, 376, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 268, 284, 300, 316, 332, 348, 364, 380, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7

    set_zetas_new t4, 116, 120, 124, 128, 132, 136, 140
    # set_zetas -2091905, 2619752, -2108549, 3539968, -300467, 2348700, -539299
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 384, 400, 416, 432, 448, 464, 480, 496, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 388, 404, 420, 436, 452, 468, 484, 500, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 392, 408, 424, 440, 456, 472, 488, 504, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 396, 412, 428, 444, 460, 476, 492, 508, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7

    set_zetas_new t4, 144, 148, 152, 156, 160, 164, 168
    # set_zetas 3119733, -2118186, -3859737, -1699267, -1643818, 3505694, -3821735
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 512, 528, 544, 560, 576, 592, 608, 624, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 516, 532, 548, 564, 580, 596, 612, 628, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 520, 536, 552, 568, 584, 600, 616, 632, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 524, 540, 556, 572, 588, 604, 620, 636, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
   
    set_zetas_new t4, 172, 176, 180, 184, 188, 192, 196
    # set_zetas -2884855, -1399561, -3277672, 3507263, -2140649, -1600420, 3699596
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 640, 656, 672, 688, 704, 720, 736, 752, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 644, 660, 676, 692, 708, 724, 740, 756, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 648, 664, 680, 696, 712, 728, 744, 760, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 652, 668, 684, 700, 716, 732, 748, 764, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7

    set_zetas_new t4, 200, 204, 208, 212, 216, 220, 224
    # et_zetas 3111497, 1757237, -19422, 811944, 531354, 954230, 3881043
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 768, 784, 800, 816, 832, 848, 864, 880, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 772, 788, 804, 820, 836, 852, 868, 884, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 776, 792, 808, 824, 840, 856, 872, 888, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 780, 796, 812, 828, 844, 860, 876, 892, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7

    set_zetas_new t4, 228, 232, 236, 240, 244, 248, 252
    # set_zetas 2680103, 4010497, 280005, 3900724, -2556880, 2071892, -2797779
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 896, 912, 928, 944, 960, 976,  992, 1008, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 900, 916, 932, 948, 964, 980,  996, 1012, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 904, 920, 936, 952, 968, 984, 1000, 1016, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7
    ct_butterfly_3layer a0, s2, s3, s4, s5, s6, s7, s8, s9, 908, 924, 940, 956, 972, 988, 1004, 1020, t2, t3, t5, t6, t0, t1, a1, a2, a3, a4, a5, a6, a7

    # 7~8 layer
    # ntt_last_layer a0, s2, s3, s4, s5,    0,    4,    8,   12, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3930395,  2091667, 3407706
    # ntt_last_layer a0, s2, s3, s4, s5,   16,   20,   24,   28, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1528703,  2316500, 3817976
    # ntt_last_layer a0, s2, s3, s4, s5,   32,   36,   40,   44, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3677745, -3342478, 2244091
    # ntt_last_layer a0, s2, s3, s4, s5,   48,   52,   56,   60, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3041255, -2446433, -3562462
    # ntt_last_layer a0, s2, s3, s4, s5,   64,   68,   72,   76, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1452451,   266997, 2434439
    # ntt_last_layer a0, s2, s3, s4, s5,   80,   84,   88,   92, t2, t3, t5, t6, t0, t1, a1, a2, a3,  3475950, -1235728, 3513181
    # ntt_last_layer a0, s2, s3, s4, s5,   96,  100,  104,  108, t2, t3, t5, t6, t0, t1, a1, a2, a3,  2176455, -3520352, -3759364
    # ntt_last_layer a0, s2, s3, s4, s5,  112,  116,  120,  124, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1585221, -1197226, -3193378
    # ntt_last_layer a0, s2, s3, s4, s5,  128,  132,  136,  140, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1257611,   900702, 1859098
    # ntt_last_layer a0, s2, s3, s4, s5,  144,  148,  152,  156, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1939314,   909542, 819034
    # ntt_last_layer a0, s2, s3, s4, s5,  160,  164,  168,  172, t2, t3, t5, t6, t0, t1, a1, a2, a3, -4083598,   495491, -1613174
    # ntt_last_layer a0, s2, s3, s4, s5,  176,  180,  184,  188, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1000202,   -43260, -522500
    # ntt_last_layer a0, s2, s3, s4, s5,  192,  196,  200,  204, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3190144,  -655327, -3122442
    # ntt_last_layer a0, s2, s3, s4, s5,  208,  212,  216,  220, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3157330,  2031748, 3207046
    # ntt_last_layer a0, s2, s3, s4, s5,  224,  228,  232,  236, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3632928, -3556995, -525098
    # ntt_last_layer a0, s2, s3, s4, s5,  240,  244,  248,  252, t2, t3, t5, t6, t0, t1, a1, a2, a3,   126922,  -768622, -3595838
    # ntt_last_layer a0, s2, s3, s4, s5,  256,  260,  264,  268, t2, t3, t5, t6, t0, t1, a1, a2, a3,  3412210,   342297, 286988
    # ntt_last_layer a0, s2, s3, s4, s5,  272,  276,  280,  284, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -983419, -2437823, 4108315
    # ntt_last_layer a0, s2, s3, s4, s5,  288,  292,  296,  300, t2, t3, t5, t6, t0, t1, a1, a2, a3,  2147896,  3437287, -3342277
    # ntt_last_layer a0, s2, s3, s4, s5,  304,  308,  312,  316, t2, t3, t5, t6, t0, t1, a1, a2, a3,  2715295,  1735879, 203044
    # ntt_last_layer a0, s2, s3, s4, s5,  320,  324,  328,  332, t2, t3, t5, t6, t0, t1, a1, a2, a3, -2967645,  2842341, 2691481
    # ntt_last_layer a0, s2, s3, s4, s5,  336,  340,  344,  348, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3693493, -2590150, 1265009
    # ntt_last_layer a0, s2, s3, s4, s5,  352,  356,  360,  364, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -411027,  4055324, 1247620
    # ntt_last_layer a0, s2, s3, s4, s5,  368,  372,  376,  380, t2, t3, t5, t6, t0, t1, a1, a2, a3, -2477047,  2486353, 1595974
    # ntt_last_layer a0, s2, s3, s4, s5,  384,  388,  392,  396, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -671102, -3767016, 1250494
    # ntt_last_layer a0, s2, s3, s4, s5,  400,  404,  408,  412, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1228525,  2635921, -3548272
    # ntt_last_layer a0, s2, s3, s4, s5,  416,  420,  424,  428, t2, t3, t5, t6, t0, t1, a1, a2, a3,   -22981, -2994039, 1869119
    # ntt_last_layer a0, s2, s3, s4, s5,  432,  436,  440,  444, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1308169,  1903435, -1050970
    # ntt_last_layer a0, s2, s3, s4, s5,  448,  452,  456,  460, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -381987, -1333058, 1237275
    # ntt_last_layer a0, s2, s3, s4, s5,  464,  468,  472,  476, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1349076, -3318210, -1430225
    # ntt_last_layer a0, s2, s3, s4, s5,  480,  484,  488,  492, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1852771,  -451100, 1312455
    # ntt_last_layer a0, s2, s3, s4, s5,  496,  500,  504,  508, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1430430,  3306115, -1962642
    # ntt_last_layer a0, s2, s3, s4, s5,  512,  516,  520,  524, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3343383, -1279661, 1917081
    # ntt_last_layer a0, s2, s3, s4, s5,  528,  532,  536,  540, t2, t3, t5, t6, t0, t1, a1, a2, a3,   264944, -2546312, -1374803
    # ntt_last_layer a0, s2, s3, s4, s5,  544,  548,  552,  556, t2, t3, t5, t6, t0, t1, a1, a2, a3,   508951,  1500165, 777191
    # ntt_last_layer a0, s2, s3, s4, s5,  560,  564,  568,  572, t2, t3, t5, t6, t0, t1, a1, a2, a3,  3097992,  2235880, 3406031
    # ntt_last_layer a0, s2, s3, s4, s5,  576,  580,  584,  588, t2, t3, t5, t6, t0, t1, a1, a2, a3,    44288,  -542412, -2831860
    # ntt_last_layer a0, s2, s3, s4, s5,  592,  596,  600,  604, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1100098, -1671176, -1846953
    # ntt_last_layer a0, s2, s3, s4, s5,  608,  612,  616,  620, t2, t3, t5, t6, t0, t1, a1, a2, a3,   904516, -2584293, -3724270
    # ntt_last_layer a0, s2, s3, s4, s5,  624,  628,  632,  636, t2, t3, t5, t6, t0, t1, a1, a2, a3,  3958618,   594136, -3776993
    # ntt_last_layer a0, s2, s3, s4, s5,  640,  644,  648,  652, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3724342, -2013608, 2432395
    # ntt_last_layer a0, s2, s3, s4, s5,  656,  660,  664,  668, t2, t3, t5, t6, t0, t1, a1, a2, a3,    -8578,  2454455, -164721
    # ntt_last_layer a0, s2, s3, s4, s5,  672,  676,  680,  684, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1653064,  1957272, 3369112
    # ntt_last_layer a0, s2, s3, s4, s5,  688,  692,  696,  700, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3249728,   185531, -1207385
    # ntt_last_layer a0, s2, s3, s4, s5,  704,  708,  712,  716, t2, t3, t5, t6, t0, t1, a1, a2, a3,  2389356, -3183426, 162844
    # ntt_last_layer a0, s2, s3, s4, s5,  720,  724,  728,  732, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -210977,  1616392, 3014001
    # ntt_last_layer a0, s2, s3, s4, s5,  736,  740,  744,  748, t2, t3, t5, t6, t0, t1, a1, a2, a3,   759969,   810149, 1652634
    # ntt_last_layer a0, s2, s3, s4, s5,  752,  756,  760,  764, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1316856, -3694233, -1799107
    # ntt_last_layer a0, s2, s3, s4, s5,  768,  772,  776,  780, t2, t3, t5, t6, t0, t1, a1, a2, a3,   189548, -3038916, 3523897
    # ntt_last_layer a0, s2, s3, s4, s5,  784,  788,  792,  796, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3553272,  3866901, 269760
    # ntt_last_layer a0, s2, s3, s4, s5,  800,  804,  808,  812, t2, t3, t5, t6, t0, t1, a1, a2, a3,  3159746,  2213111, -975884
    # ntt_last_layer a0, s2, s3, s4, s5,  816,  820,  824,  828, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1851402,  1717735, 472078
    # ntt_last_layer a0, s2, s3, s4, s5,  832,  836,  840,  844, t2, t3, t5, t6, t0, t1, a1, a2, a3, -2409325,  -426683, 1723600
    # ntt_last_layer a0, s2, s3, s4, s5,  848,  852,  856,  860, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -177440, -1803090, 1910376
    # ntt_last_layer a0, s2, s3, s4, s5,  864,  868,  872,  876, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1315589, -1667432, -1104333
    # ntt_last_layer a0, s2, s3, s4, s5,  880,  884,  888,  892, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1341330,  -260646, -3833893
    # ntt_last_layer a0, s2, s3, s4, s5,  896,  900,  904,  908, t2, t3, t5, t6, t0, t1, a1, a2, a3,  1285669, -2939036, -2235985
    # ntt_last_layer a0, s2, s3, s4, s5,  912,  916,  920,  924, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1584928,  -420899, -2286327
    # ntt_last_layer a0, s2, s3, s4, s5,  928,  932,  936,  940, t2, t3, t5, t6, t0, t1, a1, a2, a3,  -812732,   183443, -976891
    # ntt_last_layer a0, s2, s3, s4, s5,  944,  948,  952,  956, t2, t3, t5, t6, t0, t1, a1, a2, a3, -1439742,  1612842, -3545687
    # ntt_last_layer a0, s2, s3, s4, s5,  960,  964,  968,  972, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3019102,  -554416, 3919660
    # ntt_last_layer a0, s2, s3, s4, s5,  976,  980,  984,  988, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3881060,   -48306, -1362209
    # ntt_last_layer a0, s2, s3, s4, s5,  992,  996, 1000, 1004, t2, t3, t5, t6, t0, t1, a1, a2, a3, -3628969,  3937738, 1400424
    # ntt_last_layer a0, s2, s3, s4, s5, 1008, 1012, 1016, 1020, t2, t3, t5, t6, t0, t1, a1, a2, a3,  3839961,  -846154, 1976782
    

    lw      s0,  4(sp)
    lw      s1,  8(sp)
    lw      s2, 12(sp)
    lw      s3, 16(sp)
    lw      s4, 20(sp)
    lw      s5, 24(sp)
    lw      s6, 28(sp)
    lw      s7, 32(sp)
    lw      s8, 36(sp)
    lw      s9, 40(sp)
    lw      s10, 44(sp)
    lw      s11, 48(sp)
    addi    sp, sp, 52

    ret