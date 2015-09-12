//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include <cstdint>
#include <cfloat>

#include <memory>
#include <limits> // I have no limits, but 32-bit ints do
#include <functional>

#include <immintrin.h>
#include <x86intrin.h>

#include "populate.hpp"
#include "benchmark.hpp"
#include "layout.hpp"

#define NSQ 8
#define NCENT 256

namespace Counters {
	long total_scan = 0;
	float quant_bound = 0;
};

FORCE_INLINE
static inline std::int8_t Q127(const float x, const float min, const float max) {
	const float ratio = ((x-min)/(max-min));
	const float ratio_saturated = ratio < 1.0f ? ratio : 1.0f;
	return (static_cast<std::int8_t>(ratio_saturated*127));
}

#define BITMASK(x) ((1 << (x)) - 1)

FORCE_INLINE
static inline float begin_scan(const std::uint8_t* partition, const unsigned* labels, const float* dists,
		unsigned long n, todo_binheap* bh) {
	for(int pqcode_i = 0; pqcode_i < bh.todo_capacity(); ++pqcode_i) {
		float candidate = 0;
		const std::uint8_t* const pqcode = partition + pqcode_i * NSQ;
		for(int comp_i = 0; comp_i < NSQ; ++comp_i) {
				candidate += dists[comp_i * NCENT + pqcode[comp_i]];
		}
		bh.todo_add(pqcode_i, candidate);
	}
	float bh_max = bh.todo_peek();
	for(unsigned pqcode_i = bh.todo_capacity(); pqcode_i < n; ++pqcode_i) {
		float candidate = 0;
		const std::uint8_t* const pqcode = partition + pqcode_i * NSQ;
		for(int comp_i = 0; comp_i < NSQ; ++comp_i) {
				candidate += dists[comp_i * NCENT + pqcode[comp_i]];
		}
		if(candidate < bh_max) {
			bh.todo_add(pqcode_i, candidate);
			bh_max = bh.todo_peek();
		}
	}
	return bh.todo_peek();
}

FORCE_INLINE
static inline float min4_small_table(float qmax, const float* dists, __m128i (&min4)[4]) {
	float groups_min[4][16];
	float qmin = FLT_MAX;
	for(int sq_i = 0; sq_i < 4; ++sq_i) {
		// Init group minimums
		for(int group_i = 0; group_i < 16; ++group_i) {
			groups_min[sq_i][group_i] = FLT_MAX;
		}
		// Find minimum for each group
		for(int cent_i = 0; cent_i < 256; ++cent_i) {
			const int group_i = cent_i % 16;
			const float candidate = dists[sq_i * NCENT + cent_i];
			if(candidate < groups_min[sq_i][group_i]) {
				groups_min[sq_i][group_i] = candidate;
				if(candidate < qmin) {
					qmin = candidate;
				}
			}
		}
	}
	for(int sq_i = 4; sq_i < 8; ++sq_i) {
		for(int cent_i = 0; cent_i < 256; ++cent_i) {
			const float candidate = dists[sq_i * NCENT + cent_i];
			if(candidate < qmin) {
				qmin = candidate;
			}
		}
	}
	// Quantize found minimums
	for(int sq_i = 0; sq_i < 4; ++sq_i) {
		min4[sq_i] = _mm_set_epi8(
				Q127(groups_min[sq_i][15], qmin, qmax),
				Q127(groups_min[sq_i][14], qmin, qmax),
				Q127(groups_min[sq_i][13], qmin, qmax),
				Q127(groups_min[sq_i][12], qmin, qmax),
				Q127(groups_min[sq_i][11], qmin, qmax),
				Q127(groups_min[sq_i][10], qmin, qmax),
				Q127(groups_min[sq_i][9], qmin, qmax),
				Q127(groups_min[sq_i][8], qmin, qmax),
				Q127(groups_min[sq_i][7], qmin, qmax),
				Q127(groups_min[sq_i][6], qmin, qmax),
				Q127(groups_min[sq_i][5], qmin, qmax),
				Q127(groups_min[sq_i][4], qmin, qmax),
				Q127(groups_min[sq_i][3], qmin, qmax),
				Q127(groups_min[sq_i][2], qmin, qmax),
				Q127(groups_min[sq_i][1], qmin, qmax),
				Q127(groups_min[sq_i][0], qmin, qmax)
				);

	}
	return qmin;
}

FORCE_INLINE
static inline void ft4_small_table(float qmax, const float* dists, __m128i (&ft4)[4][16], const float qmin) {
	for(int sq_i = 0; sq_i < 4; ++sq_i) {
		const float* const sq_dists = dists + (sq_i + 4) * NCENT;
		for(int h_cent_i = 0; h_cent_i < 16; ++h_cent_i) {
			const float* h_dists = sq_dists + h_cent_i * 16;
			ft4[sq_i][h_cent_i] = _mm_set_epi8(
					Q127(h_dists[15], qmin, qmax),
					Q127(h_dists[14], qmin, qmax),
					Q127(h_dists[13], qmin, qmax),
					Q127(h_dists[12], qmin, qmax),
					Q127(h_dists[11], qmin, qmax),
					Q127(h_dists[10], qmin, qmax),
					Q127(h_dists[9], qmin, qmax),
					Q127(h_dists[8], qmin, qmax),
					Q127(h_dists[7], qmin, qmax),
					Q127(h_dists[6], qmin, qmax),
					Q127(h_dists[5], qmin, qmax),
					Q127(h_dists[4], qmin, qmax),
					Q127(h_dists[3], qmin, qmax),
					Q127(h_dists[2], qmin, qmax),
					Q127(h_dists[1], qmin, qmax),
					Q127(h_dists[0], qmin, qmax)
					);
		}
	}
}

FORCE_INLINE
static inline float scan_pqcode_in_simd_block_1(int index_in_block,
		const std::uint8_t* simd_block, const std::uint8_t (&high)[4], const float* dists) {
	Counters::total_scan++;
	float candidate = 0;
	const int comp_block_size = 16;
	// Such perf critical loop. Pls unroll
	for (unsigned comp_i = 0; comp_i < 4; ++comp_i) {
		const std::uint8_t* const comp_block = simd_block + comp_i * comp_block_size;
		const float* const comp_dists = dists + comp_i * NCENT;
		const std::uint8_t comp = comp_block[index_in_block];
		candidate += comp_dists[comp];
	}
	// Components 4-5
	const std::uint8_t* const comp_block_45 = simd_block + 4 * comp_block_size;
	const std::uint8_t packed_45 = comp_block_45[index_in_block];
	const float* const comp_dists_4 = dists + 4 * NCENT;
	candidate += comp_dists_4[(packed_45 & 0xf) + high[0]];
	const float* const comp_dists_5 = dists + 5 * NCENT;
	candidate += comp_dists_5[((packed_45 & 0xf0) >> 4) + high[1]];

	// Component 6-7
	const std::uint8_t* comp_block_56 = simd_block + 5 * comp_block_size;
	const std::uint8_t packed_56 = comp_block_56[index_in_block];
	const float* const comp_dists_6 = dists + 6 * NCENT;
	candidate += comp_dists_6[(packed_56 & 0xf) + high[2]];
	const float* const comp_dists_7 = dists + 7 * NCENT;
	candidate += comp_dists_7[((packed_56 & 0xf0) >> 4) + high[3]];

	return candidate;
}

const std::uint64_t masktable[] = { 0x0, 0x0, 0x1, 0x100, 0x2, 0x200, 0x201,
		0x20100, 0x3, 0x300, 0x301, 0x30100, 0x302, 0x30200, 0x30201, 0x3020100,
		0x4, 0x400, 0x401, 0x40100, 0x402, 0x40200, 0x40201, 0x4020100, 0x403,
		0x40300, 0x40301, 0x4030100, 0x40302, 0x4030200, 0x4030201, 0x403020100,
		0x5, 0x500, 0x501, 0x50100, 0x502, 0x50200, 0x50201, 0x5020100, 0x503,
		0x50300, 0x50301, 0x5030100, 0x50302, 0x5030200, 0x5030201, 0x503020100,
		0x504, 0x50400, 0x50401, 0x5040100, 0x50402, 0x5040200, 0x5040201,
		0x504020100, 0x50403, 0x5040300, 0x5040301, 0x504030100, 0x5040302,
		0x504030200, 0x504030201, 0x50403020100, 0x6, 0x600, 0x601, 0x60100,
		0x602, 0x60200, 0x60201, 0x6020100, 0x603, 0x60300, 0x60301, 0x6030100,
		0x60302, 0x6030200, 0x6030201, 0x603020100, 0x604, 0x60400, 0x60401,
		0x6040100, 0x60402, 0x6040200, 0x6040201, 0x604020100, 0x60403,
		0x6040300, 0x6040301, 0x604030100, 0x6040302, 0x604030200, 0x604030201,
		0x60403020100, 0x605, 0x60500, 0x60501, 0x6050100, 0x60502, 0x6050200,
		0x6050201, 0x605020100, 0x60503, 0x6050300, 0x6050301, 0x605030100,
		0x6050302, 0x605030200, 0x605030201, 0x60503020100, 0x60504, 0x6050400,
		0x6050401, 0x605040100, 0x6050402, 0x605040200, 0x605040201,
		0x60504020100, 0x6050403, 0x605040300, 0x605040301, 0x60504030100,
		0x605040302, 0x60504030200, 0x60504030201, 0x6050403020100, 0x7, 0x700,
		0x701, 0x70100, 0x702, 0x70200, 0x70201, 0x7020100, 0x703, 0x70300,
		0x70301, 0x7030100, 0x70302, 0x7030200, 0x7030201, 0x703020100, 0x704,
		0x70400, 0x70401, 0x7040100, 0x70402, 0x7040200, 0x7040201, 0x704020100,
		0x70403, 0x7040300, 0x7040301, 0x704030100, 0x7040302, 0x704030200,
		0x704030201, 0x70403020100, 0x705, 0x70500, 0x70501, 0x7050100, 0x70502,
		0x7050200, 0x7050201, 0x705020100, 0x70503, 0x7050300, 0x7050301,
		0x705030100, 0x7050302, 0x705030200, 0x705030201, 0x70503020100,
		0x70504, 0x7050400, 0x7050401, 0x705040100, 0x7050402, 0x705040200,
		0x705040201, 0x70504020100, 0x7050403, 0x705040300, 0x705040301,
		0x70504030100, 0x705040302, 0x70504030200, 0x70504030201,
		0x7050403020100, 0x706, 0x70600, 0x70601, 0x7060100, 0x70602, 0x7060200,
		0x7060201, 0x706020100, 0x70603, 0x7060300, 0x7060301, 0x706030100,
		0x7060302, 0x706030200, 0x706030201, 0x70603020100, 0x70604, 0x7060400,
		0x7060401, 0x706040100, 0x7060402, 0x706040200, 0x706040201,
		0x70604020100, 0x7060403, 0x706040300, 0x706040301, 0x70604030100,
		0x706040302, 0x70604030200, 0x70604030201, 0x7060403020100, 0x70605,
		0x7060500, 0x7060501, 0x706050100, 0x7060502, 0x706050200, 0x706050201,
		0x70605020100, 0x7060503, 0x706050300, 0x706050301, 0x70605030100,
		0x706050302, 0x70605030200, 0x70605030201, 0x7060503020100, 0x7060504,
		0x706050400, 0x706050401, 0x70605040100, 0x706050402, 0x70605040200,
		0x70605040201, 0x7060504020100, 0x706050403, 0x70605040300,
		0x70605040301, 0x7060504030100, 0x70605040302, 0x7060504030200,
		0x7060504030201, 0x706050403020100 };

FORCE_INLINE
static inline void fast_scan_1(const std::uint8_t* partition, const unsigned* labels,
		const float* dists, const __m128i (&min4)[4], __m128i (&ft4)[4][16],
		const float qmin, const float qmax, todo_binheap* bh, unsigned scan_pqcode_count) {
	const unsigned simd_pqcode_count = 16;
	const int comp_block_size = 16;
	const unsigned simd_block_size = simd_pqcode_count * (4 * 1 + 4 * 0.5);
	const group_header* hdr;

	float bh_bound = qmax;
	__m128i bh_bound_quant = _mm_set1_epi8(Q127(bh_bound, qmin, bh_bound)); // CHK. Is 127

	for (;;) {
		// Parse group header
		hdr = reinterpret_cast<const group_header*>(partition);
		// Check if last group (All bits of size set to 1)
		if (hdr->size == std::numeric_limits<decltype(hdr->size)>::max()) {
			return;
		}
		partition += sizeof(*hdr);
		unsigned simd_block_count = (static_cast<unsigned>(hdr->size)
				+ simd_pqcode_count - 1) / simd_pqcode_count;
		// Load tables
		__m128i ft4_group[4];
		ft4_group[0] = ft4[0][hdr->values[0] >> 4];
		ft4_group[1] = ft4[1][hdr->values[1] >> 4];
		ft4_group[2] = ft4[2][hdr->values[2] >> 4];
		ft4_group[3] = ft4[3][hdr->values[3] >> 4];
		// Scan SIMD Blocks
		while (simd_block_count--) {
			const __m128i low_bits_mask = _mm_set_epi64x(0x0f0f0f0f0f0f0f0f,
					0x0f0f0f0f0f0f0f0f);

			// Component 0
			const __m128i comps_0 = _mm_loadu_si128(
					reinterpret_cast<const __m128i *>(partition));
			const __m128i masked_comps_0 = _mm_and_si128(comps_0, low_bits_mask);
			__m128i candidates = _mm_shuffle_epi8(min4[0], masked_comps_0);
			// Components 1..3
			for (int comp_i = 1; comp_i < 4; ++comp_i) {
				const __m128i comps = _mm_loadu_si128(
						reinterpret_cast<const __m128i *>(partition
								+ comp_i * comp_block_size));
				const __m128i masked_comps = _mm_and_si128(comps, low_bits_mask);
				const __m128i partial = _mm_shuffle_epi8(min4[comp_i], masked_comps);
				candidates = _mm_adds_epi8(candidates, partial);
			}
			// Components 4-5
			__m128i comps_45 = _mm_loadu_si128(
					reinterpret_cast<const __m128i *>(partition
							+ 4 * comp_block_size));
			const __m128i masked_comps_4 = _mm_and_si128(comps_45, low_bits_mask);
			const __m128i partial_4 = _mm_shuffle_epi8(ft4_group[0], masked_comps_4);
			candidates = _mm_adds_epi8(candidates, partial_4);

			comps_45 = _mm_srli_epi64(comps_45, 4);
			const __m128i masked_comps_5 = _mm_and_si128(comps_45, low_bits_mask);
			const __m128i partial_5 = _mm_shuffle_epi8(ft4_group[1], masked_comps_5);
			candidates = _mm_adds_epi8(candidates, partial_5);

			// Components 6-7
			__m128i comps_67 = _mm_loadu_si128(
					reinterpret_cast<const __m128i *>(partition
							+ 5 * comp_block_size));
			const __m128i masked_comps_6 = _mm_and_si128(comps_67, low_bits_mask);
			const __m128i partial_6 = _mm_shuffle_epi8(ft4_group[2], masked_comps_6);
			candidates = _mm_adds_epi8(candidates, partial_6);

			const __m128i comps_7 = _mm_srli_epi64(comps_67, 4);
			const __m128i masked_comp_7 = _mm_and_si128(comps_7, low_bits_mask);
			const __m128i partial_7 = _mm_shuffle_epi8(ft4_group[3], masked_comp_7);
			candidates = _mm_adds_epi8(candidates, partial_7);

			// Compare
			const __m128i compare = _mm_cmplt_epi8(candidates, bh_bound_quant);
			int cmp = _mm_movemask_epi8(compare);
			//std::uint64_t cmp_low = (_mm_cvtsi128_si64(compare));
			//std::uint64_t cmp_high = (_mm_extract_epi64(compare, 1));

			// Compute current block size
			int current_block_actual_size = 0;
			if(simd_block_count == 0) {
				current_block_actual_size = hdr->size % simd_pqcode_count;
				if(current_block_actual_size == 0) {
					current_block_actual_size = simd_pqcode_count;
				} else {
					/*__m128i mask;
					compute_simd_mask(current_block_actual_size, mask);
					compare = _mm_and_si128(compare, mask);*/
					/*
					std::uint64_t low_mask;
					std::uint64_t high_mask;
					compute_high_low_mask(current_block_actual_size, low_mask, high_mask);
					cmp_low = cmp_low & low_mask;
					cmp_high = cmp_high & high_mask;
					*/
					cmp = cmp & BITMASK(current_block_actual_size);
				}
			} else {
				current_block_actual_size = simd_pqcode_count;
			}

			if(cmp) {
				// Check low quadword
				const std::uint8_t cmp_low = cmp & 0xff;
				if (cmp_low) {
					/*const std::uint64_t low_possible_positions = 0x0706050403020100;
					const std::uint64_t match_positions = _pext_u64(
							low_possible_positions, cmp_low);*/
					const int match_count = _popcnt32(cmp_low);
					std::uint64_t match_pos = masktable[cmp_low];


					for (int i = 0; i < match_count; ++i) {
						const std::uint8_t pos = match_pos & 0xff;
						match_pos >>= 8;
						const float candidate = scan_pqcode_in_simd_block_1(pos,
								partition, hdr->values, dists);
						if (candidate < bh_bound) {
							bh.todo_add(labels[scan_pqcode_count + pos],
									candidate);
							bh_bound = bh.todo_peek();
							bh_bound_quant = _mm_set1_epi8(
									Q127(bh_bound, qmin, qmax));
						}
					}
				}

				// Check high quadword
				const std::uint8_t cmp_high = (cmp >> 8);
				if (cmp_high) {
					/*const std::uint64_t high_possible_positions = 0x0f0e0d0c0b0a0908;
					const std::uint64_t match_positions = _pext_u64(
							high_possible_positions, cmp_high);*/
					const int match_count = _popcnt32(cmp_high);
					std::uint64_t match_pos = masktable[cmp_high] + 0x0808080808080808;

					for (int i = 0; i < match_count; ++i) {
						const std::uint8_t pos = match_pos & 0xff;
						match_pos >>= 8;
						const float candidate = scan_pqcode_in_simd_block_1(pos,
								partition, hdr->values, dists);
						if (candidate < bh_bound) {
							bh.todo_add(labels[scan_pqcode_count + pos],
									candidate);
							bh_bound = bh.todo_peek();
							bh_bound_quant = _mm_set1_epi8(
									Q127(bh_bound, qmin, qmax));
						}
					}
				}
			}

			partition += simd_block_size;
			scan_pqcode_count += current_block_actual_size;
		}
	}
}

void scan_partition_1(const std::uint8_t* partition, const unsigned* labels,
		const float* dists, todo_binheap* bh) {
	// 0. Scan first keep pqcodes
	const std::uint32_t keep = *reinterpret_cast<const std::uint32_t*>(partition);
	partition += sizeof(keep);
	float qmax = begin_scan(partition, labels, dists, keep, bh);
	Counters::quant_bound = qmax;
	partition += keep * 8;
	// 1. Get small tables
	__m128i ft4[4][16];
	__m128i min4[4];
	float qmin = min4_small_table(qmax, dists, min4);
	ft4_small_table(qmax, dists, ft4, qmin);
	// 2. Fast scan (with lower bounds)
	fast_scan_1(partition, labels, dists, min4, ft4, qmin, qmax,
			bh, 0);
}
