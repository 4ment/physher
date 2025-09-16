//
//  utilssse.c
//  physher
//
//  Created by Mathieu Fourment on 29/8/19.
//  Copyright © 2019 Mathieu Fourment. All rights reserved.
//

#include "utilssse.h"

#include <stddef.h>

#ifdef SSE3_ENABLED
#if defined(__aarch64__)
#include "neon2sse.h"
#else
#include <xmmintrin.h> // SSE
#include <pmmintrin.h> // SSE3
//#include <tmmintrin.h> // SSSE3
#endif
#endif

void add_vector_vector(double* out, const double* vec1, const double* vec2, size_t length, size_t chunk){
#ifdef SSE3_ENABLED
	if ((chunk & 1) == 0) {
		__m128d *pout = (__m128d*)out;
		__m128d *pvec1 = (__m128d*)vec1;
		__m128d *pvec2 = (__m128d*)vec2;
		for (size_t i = 0; i < length; i+=2) {
			*pout = _mm_add_pd(*pvec1, *pvec2);
			pout++;
			pvec1++;
			pvec2++;
		}
		return;
	}
#endif
	for (size_t i = 0; i < length; i++) {
		out[i] = vec1[i] + vec2[i];
	}
}

void mult_vector_vector(double* out, const double* vec1, const double* vec2, size_t length, size_t chunk){
#ifdef SSE3_ENABLED
	if ((chunk & 1) == 0) {
		__m128d *pout = (__m128d*)out;
		__m128d *pvec1 = (__m128d*)vec1;
		__m128d *pvec2 = (__m128d*)vec2;
		for (size_t i = 0; i < length; i+=2) {
			*pout = _mm_mul_pd(*pvec1, *pvec2);
			pout++;
			pvec1++;
			pvec2++;
		}
		return;
	}
#endif
	for (size_t i = 0; i < length; i++) {
		out[i] = vec1[i] * vec2[i];
	}
}

void mult_vector_vector_inplace(double* vec1, const double* vec2, size_t length, size_t chunk){
#ifdef SSE3_ENABLED
	if ((chunk & 1) == 0) {
		__m128d *pvec1 = (__m128d*)vec1;
		__m128d *pvec2 = (__m128d*)vec2;
		for (size_t i = 0; i < length; i+=2) {
			*pvec1 = _mm_mul_pd(*pvec1, *pvec2);
			pvec1++;
			pvec2++;
		}
		return;
	}
#endif
	for (size_t i = 0; i < length; i++) {
		vec1[i] *= vec2[i];
	}
}

void add_vector_vector_mult_vector(double* out, const double* vec1, const double* vec2, const double* vec3, size_t length, size_t chunk){
#ifdef SSE3_ENABLED
	if ((chunk & 1) == 0) {
		__m128d *pout = (__m128d*)out;
		__m128d *pvec1 = (__m128d*)vec1;
		__m128d *pvec2 = (__m128d*)vec2;
		__m128d *pvec3 = (__m128d*)vec3;
		for (size_t i = 0; i < length; i+=2) {
			*pout = _mm_mul_pd(_mm_add_pd(*pvec1, *pvec2), *pvec3);
			pout++;
			pvec1++;
			pvec2++;
			pvec3++;
		}
		return;
	}
#endif
	for (size_t i = 0; i < length; i++) {
		out[i] = (vec1[i] + vec2[i])*vec3[i];
	}
}
