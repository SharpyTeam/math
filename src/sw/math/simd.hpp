//
// Created by selya on 06.10.2019.
//

#ifndef SANDWICH_MATH_SIMD_HPP
#define SANDWICH_MATH_SIMD_HPP

#ifdef __clang__
#pragma clang attribute push (__attribute__((target("avx"))), apply_to=function)
#elif defined(__GNUC__)
#pragma GCC push_options
#pragma GCC target ("avx")
#pragma GCC optimize ("O3")
#pragma GCC tune ("skylake")
#endif

#include <immintrin.h>

namespace sw {
namespace math {
namespace simd {

/*inline void mat4d_mul_asm(double dst[16], const double a[16], const double b[16]) {
    asm volatile (
        "vmovupd (%%rax), %%ymm3\n"
        "vmovupd 32(%%rax), %%ymm4\n"
        "vmovupd 64(%%rax), %%ymm5\n"
        "vmovupd 96(%%rax), %%ymm6\n"

        // 0
        "vbroadcastsd (%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm3, %%ymm0\n"
        "vbroadcastsd 8(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm4, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 16(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm5, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 24(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm6, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"

        "vmovupd %%ymm0, (%%rcx)\n"

        // 1
        "vbroadcastsd 32(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm3, %%ymm0\n"
        "vbroadcastsd 40(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm4, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 48(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm5, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 56(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm6, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"

        "vmovupd %%ymm0, 32(%%rcx)\n"

        // 2
        "vbroadcastsd 64(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm3, %%ymm0\n"
        "vbroadcastsd 72(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm4, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 80(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm5, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 88(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm6, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"

        "vmovupd %%ymm0, 64(%%rcx)\n"

        // 3
        "vbroadcastsd 96(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm3, %%ymm0\n"
        "vbroadcastsd 104(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm4, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 112(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm5, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"
        "vbroadcastsd 120(%%rbx), %%ymm1\n"
        "vmulpd %%ymm1, %%ymm6, %%ymm1\n"
        "vaddpd %%ymm1, %%ymm0, %%ymm0\n"

        "vmovupd %%ymm0, 96(%%rcx)\n"

        :
        : "c" (dst), "a" (a), "b" (b)
        : "ymm0", "ymm1", "ymm2", "ymm3", "ymm4", "ymm5", "ymm6"
    );
}*/

inline void mat4d_mul(double dst[16], const double a[16], const double b[16]) {
    // 1. Store 4 columns of A matrix
    __m256d ac0 = _mm256_loadu_pd(a + 0);
    __m256d ac1 = _mm256_loadu_pd(a + 4);
    __m256d ac2 = _mm256_loadu_pd(a + 8);
    __m256d ac3 = _mm256_loadu_pd(a + 12);

    // 2. Calculate Out[k] (column k) (k in [0, 3]):
    // 3. Get A[i] (column 0), mul it by B[k][i] (i in [0, 3])
    // 4. Sum results
    // 5. Store resulting sum as Out[k]
    _mm256_storeu_pd(
            dst + 0,
            _mm256_add_pd(
                    _mm256_add_pd(
                            _mm256_mul_pd(ac0, _mm256_broadcast_sd(b + 0)),
                            _mm256_mul_pd(ac1, _mm256_broadcast_sd(b + 1))),
                    _mm256_add_pd(
                            _mm256_mul_pd(ac2, _mm256_broadcast_sd(b + 2)),
                            _mm256_mul_pd(ac3, _mm256_broadcast_sd(b + 3)))
            )
    );

    _mm256_storeu_pd(
            dst + 4,
            _mm256_add_pd(
                    _mm256_add_pd(
                            _mm256_mul_pd(ac0, _mm256_broadcast_sd(b + 4)),
                            _mm256_mul_pd(ac1, _mm256_broadcast_sd(b + 5))),
                    _mm256_add_pd(
                            _mm256_mul_pd(ac2, _mm256_broadcast_sd(b + 6)),
                            _mm256_mul_pd(ac3, _mm256_broadcast_sd(b + 7)))
            )
    );

    _mm256_storeu_pd(
            dst + 8,
            _mm256_add_pd(
                    _mm256_add_pd(
                            _mm256_mul_pd(ac0, _mm256_broadcast_sd(b + 8)),
                            _mm256_mul_pd(ac1, _mm256_broadcast_sd(b + 9))),
                    _mm256_add_pd(
                            _mm256_mul_pd(ac2, _mm256_broadcast_sd(b + 10)),
                            _mm256_mul_pd(ac3, _mm256_broadcast_sd(b + 11)))
            )
    );

    _mm256_storeu_pd(
            dst + 12,
            _mm256_add_pd(
                    _mm256_add_pd(
                            _mm256_mul_pd(ac0, _mm256_broadcast_sd(b + 12)),
                            _mm256_mul_pd(ac1, _mm256_broadcast_sd(b + 13))),
                    _mm256_add_pd(
                            _mm256_mul_pd(ac2, _mm256_broadcast_sd(b + 14)),
                            _mm256_mul_pd(ac3, _mm256_broadcast_sd(b + 15)))
            )
    );
}

} //namespace simd
} //namespace math
} //namespace sw

#ifdef __clang__
#pragma clang attribute pop
#elif __GNUC__
#pragma GCC pop_options
#endif

#endif //SANDWICH_MATH_SIMD_HPP
