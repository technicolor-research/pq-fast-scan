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
#include "scan_sse.hpp"
#include <cassert>

#define NSQ 8
#define NCENT 256
#define SSE_SZ 4

void scan_prefetch_sse(const char* partition, const float* dists,
		unsigned long n, pq_params pqp, binheap* bh) {
	long binheap_op = 0;
	for (int t = 0; t < bh->capacity(); ++t) {
		bh->push(0, FLT_MAX - t);
	}
	float minb = bh->max();
	__m128 min = _mm_set1_ps(minb);
	float mcandidates[SSE_SZ];
	__m128 partial;
	uint8_t* u8_buf = (uint8_t*) partition;
	for (unsigned long i = 0; i < n ; i += SSE_SZ) {
		const uint8_t* const pqcode0 = u8_buf + i * NSQ;
		const uint8_t* const pqcode1 = pqcode0 + NSQ;
		const uint8_t* const pqcode2 = pqcode0 + 2*NSQ;
		const uint8_t* const pqcode3 = pqcode0 + 3*NSQ;
		__m128 candidates = _mm_set_ps(
					dists[pqcode3[0]],
					dists[pqcode2[0]],
					dists[pqcode1[0]],
					dists[pqcode0[0]]
					);
			// Such perf critical loop. Pls unroll
			for (unsigned j = 1; j < NSQ; ++j) {
				const float* const cdist = dists + j * NCENT;
				partial = _mm_set_ps(
						cdist[pqcode3[j]],
						cdist[pqcode2[j]],
						cdist[pqcode1[j]],
						cdist[pqcode0[j]]
						);
				candidates = _mm_add_ps(candidates, partial);
			}
			const __m128 cmp = _mm_cmp_ps(candidates, min,_CMP_LT_OS);
			if (!_mm_testz_ps (cmp, cmp)) {
				_mm_store_ps(mcandidates, candidates);
				for(unsigned ii = 0; ii < SSE_SZ; ++ii) {
					if(mcandidates[ii] < minb) {
						bh->push(i+ii, mcandidates[ii]);
						minb = bh->max();
						binheap_op++;
					}
				}
				min = _mm_set1_ps(minb);
			}
		_mm_prefetch(pqcode0 + 64*4, _MM_HINT_NTA);
	}
}

#undef SSE_SZ

