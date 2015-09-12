//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include <cfloat>
#include <immintrin.h>
#include "scan_avx.hpp"

#define NSQ 8
#define BITSQ 8
#define NCENT 256

#define AVX_SZ 8

/** Base functions **/
__attribute__((always_inline))
static inline void simd_lookup_add_min(const uint8_t* const pqcode0,
		const uint8_t* const& pqcode1,
		const uint8_t* const& pqcode2,
		const uint8_t* const& pqcode3,
		const uint8_t* const& pqcode4,
		const uint8_t* const& pqcode5,
		const uint8_t* const& pqcode6,
		const uint8_t* const& pqcode7,
		const float*& dists,
		__m256& candidates,
		__m256& partial,
		float (&mcandidates)[AVX_SZ],
		__m256& min,
		float& minb,
		todo_binheap*& bh,
		long& binheap_op,
		unsigned long& i
		) {
	candidates = _mm256_set_ps(
			dists[pqcode7[0]],
			dists[pqcode6[0]],
			dists[pqcode5[0]],
			dists[pqcode4[0]],
			dists[pqcode3[0]],
			dists[pqcode2[0]],
			dists[pqcode1[0]],
			dists[pqcode0[0]]
			);
	// Such perf critical loop. Pls unroll
	for (unsigned j = 1; j < NSQ; ++j) {
		const float* const cdist = dists + j * NCENT;
		partial = _mm256_set_ps(
				cdist[pqcode7[j]],
				cdist[pqcode6[j]],
				cdist[pqcode5[j]],
				cdist[pqcode4[j]],
				cdist[pqcode3[j]],
				cdist[pqcode2[j]],
				cdist[pqcode1[j]],
				cdist[pqcode0[j]]
				);
		candidates = _mm256_add_ps(candidates, partial);
	}
	const __m256 cmp = _mm256_cmp_ps(candidates, min,_CMP_LT_OS);
	if (!_mm256_testz_ps (cmp, cmp)) {
		_mm256_store_ps(mcandidates, candidates);
		for(unsigned ii = 0; ii < AVX_SZ; ++ii) {
			if(mcandidates[ii] < minb) {
				bh.todo_add(i+ii, mcandidates[ii]);
				minb = bh.todo_peek();
				binheap_op++;
			}
		}
		min = _mm256_set1_ps(minb);
	}
}

void scan_bh_prefetch_avx(const char* partition, const float* dists,
		unsigned long n, pq_params pqp, todo_binheap* bh) {
	long binheap_op = 0;
	for (int t = 0; t < bh.todo_capacity(); ++t) {
		bh.todo_add(0, FLT_MAX - t);
	}
	float minb = bh.todo_peek();
	__m256 min = _mm256_set1_ps(minb);
	__m256 candidates = _mm256_setzero_ps();
	float mcandidates[AVX_SZ];
	__m256 partial;
	uint8_t* u8_buf = (uint8_t*) partition;
	for (unsigned long i = 0; i < n ; i += AVX_SZ) {
		const uint8_t* const pqcode0 = u8_buf + i * NSQ;
		const uint8_t* const pqcode1 = pqcode0 + NSQ;
		const uint8_t* const pqcode2 = pqcode0 + 2*NSQ;
		const uint8_t* const pqcode3 = pqcode0 + 3*NSQ;
		const uint8_t* const pqcode4 = pqcode0 + 4*NSQ;
		const uint8_t* const pqcode5 = pqcode0 + 5*NSQ;
		const uint8_t* const pqcode6 = pqcode0 + 6*NSQ;
		const uint8_t* const pqcode7 = pqcode0 + 7*NSQ;
		simd_lookup_add_min(pqcode0, pqcode1, pqcode2, pqcode3, pqcode4, pqcode5, pqcode6, pqcode7,
				dists, candidates, partial, mcandidates, min, minb, bh, binheap_op,
				i);
		_mm_prefetch(pqcode0 + 64*8, _MM_HINT_NTA);
	}
}

