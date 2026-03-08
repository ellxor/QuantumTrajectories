#pragma once
#define _GNU_SOURCE
#include <complex.h>
#include <stdatomic.h>

#define atomic(T) _Atomic(T)
#define swap(a,b) do { auto _tmp = a; a = b; b = _tmp; } while(0)
#define fill(begin, end, value) do { for (auto _p = begin; _p < end; ++_p) *_p = value; } while (0)

/*
 * Supplementary functions for C (based on C++ algorithm)
 */

static inline int min(int a, int b) { return (a < b) ? a : b; }
static inline int max(int a, int b) { return (a > b) ? a : b; }

/*
 * Supplementary functions for complex.h (based on C++ complex)
 */

static inline float cnormf(complex float c) {
	float r = crealf(c);
	float i = cimagf(c);
	return r*r + i*i;
}
