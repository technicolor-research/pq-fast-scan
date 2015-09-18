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

#include "scan_gather.hpp"

#define NSQ 8
#define BITSQ 8
#define NCENT 256

#define AVX_SZ 8

void scan_bh_avx_gths(const char* partition, const float* dists, unsigned long n,
		pq_params pqp, binheap* bh) {
	for (int t = 0; t < bh->capacity(); ++t) {
		bh->push(0, FLT_MAX - t);
	}
	float minb = bh->max();
	__m256 min = _mm256_set1_ps(minb);
	float mcandidates[AVX_SZ];
	const __m256i mask_firstbyte = _mm256_set1_epi32(0xFF);
	const __m256i * u256_buf = (const __m256i *) partition;
	// Each iteration of the loop will process a block
	// 1 block = 32 pqcodes
	const unsigned long blk_n = n / 32;
	__m256i data[NSQ];
	// This loop processes a block
	// i: block number
	for (unsigned long i = 0; i < blk_n; i += 1) {
		// Load block data
		for (unsigned int j = 0; j < NSQ; j++) {
			const __m256i* __restrict__ pqcode = u256_buf + i * NSQ + j;
			data[j] = _mm256_load_si256(pqcode);  // Temporal load
		}
		// Prefetch next block. This has no significant impact.
		//_mm_prefetch(u256_buf + (i+1) * NSQ, _MM_HINT_NTA);
		// We handle 8 pqcodes/gather operation
		// 8 pqcodes/gather * 4 = 32 pqcodes / blocks
		for(unsigned s=0; s<4; s++) {
			// Gather and Add
			__m256 candidates = _mm256_i32gather_ps(dists,_mm256_and_si256(data[0],mask_firstbyte),4);
			for(unsigned int j=1;j< NSQ;j++) {
				const float *__restrict__ cdists = dists+j*NCENT;
				const __m256i value = _mm256_and_si256(data[j],mask_firstbyte);
				const __m256 partial = _mm256_i32gather_ps(cdists,value,4);
				candidates = _mm256_add_ps(candidates, partial);
			}
			// Check
			const __m256 cmp = _mm256_cmp_ps(candidates, min,_CMP_LT_OS);
			if (!_mm256_testz_ps (cmp, cmp)) {
				_mm256_store_ps(mcandidates, candidates);
				for(unsigned ii = 0; ii < AVX_SZ; ++ii) {
					if(mcandidates[ii] < minb) {
						const int index=(i*32)+(ii*4)+s;
						bh->push(index, mcandidates[ii]);
						minb = bh->max();
					}
				}
				min = _mm256_set1_ps(minb);
			}
			// Shift to next byte
			for(unsigned int j=0;j<NSQ;j++) {
				data[j] = _mm256_srli_epi32(data[j],8);
			}
		}
	}
}

void gather_check(const char* partition, const char* i32_partition,
		unsigned long n) {
	const __m256i * u256_buf = (const __m256i *) i32_partition;
	const unsigned long blk_n = n / (8 * 32);
	__m256i data[NSQ];
	char* data_mem;
	posix_memalign((void**) &data_mem, 32, 320);
	for (unsigned long i = 0; i < blk_n; ++i) {
		// Load
		for (unsigned int j = 0; j < NSQ; j++) {
			const __m256i* __restrict__ pqcode = u256_buf + i * NSQ + j;
			data[j] = _mm256_load_si256(pqcode);  // Temporal load
		}
		// Check
		for (unsigned int j = 0; j < NSQ; j++) {
			_mm256_store_si256((__m256i*) data_mem, data[j]);
			for(unsigned l=0; l<32; ++l) {
				if(partition[(i*32+l)*8+j] != data_mem[l]) {
					printf("i=%lu, j=%u, l=%u\n", i,j,l);
					exit(32);
				}
			}
		}
	}
	free(data_mem);
}

#undef AVX_SZ


