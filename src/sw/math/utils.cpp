//
// Created by selya on 05.10.2019.
//

#include "utils.hpp"

#include <cstring>

#ifndef _WIN32
#include <cpuid.h>

#define _XCR_XFEATURE_ENABLED_MASK 0

uint64_t _xgetbv(unsigned int index) {
    uint32_t eax, edx;
    __asm__ __volatile__(
        "xgetbv;"
        : "=a" (eax), "=d"(edx)
        : "c" (index)
    );
    return ((uint64_t)edx << (uint64_t)32) | (uint64_t)eax;
}
#else
#include <immintrin.h>
#include <intrin.h>
#endif

namespace sw {
namespace math {
namespace utils {

const CPUFeatures CPUFeatures::instance;

CPUFeatures::CPUFeatures() noexcept {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    OS_x64 = true;
#endif

    uint32_t info[4];
    char name[13];
    CPUID(info, 0);
    std::memcpy(name + 0, &info[1], 4);
    std::memcpy(name + 4, &info[3], 4);
    std::memcpy(name + 8, &info[2], 4);
    name[12] = '\0';

    VENDOR_STRING = name;
    if (VENDOR_STRING == "GenuineIntel") {
        VENDOR_Intel = true;
    } else if (VENDOR_STRING == "AuthenticAMD") {
        VENDOR_AMD = true;
    }

    uint32_t number_of_ids = info[0];

    CPUID(info, 0x80000000u);
    uint32_t number_of_ext_ids = info[0];

    if (number_of_ids >= 0x00000001u) {
        CPUID(info, 0x00000001u);
    } else {
        info[0] = 0u;
        info[1] = 0u;
        info[2] = 0u;
        info[3] = 0u;
    }

    HW_MMX    = info[3] & (1u << 23u);
    HW_SSE    = info[3] & (1u << 25u);
    HW_SSE2   = info[3] & (1u << 26u);
    HW_SSE3   = info[2] & (1u << 0u);
    HW_SSSE3  = info[2] & (1u << 9u);
    HW_SSE41  = info[2] & (1u << 19u);
    HW_SSE42  = info[2] & (1u << 20u);
    HW_AES    = info[2] & (1u << 25u);
    HW_AVX    = info[2] & (1u << 28u);
    HW_FMA3   = info[2] & (1u << 12u);
    HW_RDRAND = info[2] & (1u << 30u);

    bool OS_USES_XSAVE_XRSTORE = info[2] & (1u << 27u);

    if (OS_USES_XSAVE_XRSTORE && HW_AVX) {
        uint64_t xcr_feature_mask = _xgetbv(_XCR_XFEATURE_ENABLED_MASK);
        OS_AVX = (xcr_feature_mask & 0x6u) == 0x6u;
        OS_AVX512 = (xcr_feature_mask & 0xe6u) == 0xe6u;
    }

    if (number_of_ids >= 0x00000007) {
        CPUID(info, 0x00000007);
    } else {
        info[0] = 0u;
        info[1] = 0u;
        info[2] = 0u;
        info[3] = 0u;
    }

    HW_AVX2         = info[1] & (1u <<  5u);
    HW_BMI1         = info[1] & (1u <<  3u);
    HW_BMI2         = info[1] & (1u <<  8u);
    HW_ADX          = info[1] & (1u << 19u);
    HW_MPX          = info[1] & (1u << 14u);
    HW_SHA          = info[1] & (1u << 29u);
    HW_PREFETCHWT1  = info[2] & (1u <<  0u);
    HW_AVX512_F     = info[1] & (1u << 16u);
    HW_AVX512_CD    = info[1] & (1u << 28u);
    HW_AVX512_PF    = info[1] & (1u << 26u);
    HW_AVX512_ER    = info[1] & (1u << 27u);
    HW_AVX512_VL    = info[1] & (1u << 31u);
    HW_AVX512_BW    = info[1] & (1u << 30u);
    HW_AVX512_DQ    = info[1] & (1u << 17u);
    HW_AVX512_IFMA  = info[1] & (1u << 21u);
    HW_AVX512_VBMI  = info[2] & (1u <<  1u);

    if (number_of_ext_ids >= 0x80000001) {
        CPUID(info, 0x80000001);
    } else {
        info[0] = 0u;
        info[1] = 0u;
        info[2] = 0u;
        info[3] = 0u;
    }

    HW_x64   = info[3] & (1u << 29u);
    HW_ABM   = info[2] & (1u <<  5u);
    HW_SSE4a = info[2] & (1u <<  6u);
    HW_FMA4  = info[2] & (1u << 16u);
    HW_XOP   = info[2] & (1u << 11u);
}

void CPUFeatures::CPUID(uint32_t *info, uint32_t index) {
#ifdef _WIN32
    __cpuidex((int *)info, index, 0);
#else
    __cpuid_count(index, 0, info[0], info[1], info[2], info[3]);
#endif
}

} //namespace utils
} //namespace math
} //namespace sw