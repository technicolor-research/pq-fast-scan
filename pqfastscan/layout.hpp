//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef LAYOUT_HPP_
#define LAYOUT_HPP_

#include <cstdint>
#include <vector>
#include <queue>
#include <memory>
#include <iostream>
#include "common.hpp"

/* Bits copy functions		*/

// Be careful when x >= 2**32
#define BITMASKU8(x) ((1U << (x)) - 1)

/**
 * Select low_sel bits from src starting from the lower bit and
 * copy the selected bits to dst at offset msb_off, starting from the most
 * significant bit.
 * The caller must ensure that low_sel + msb_off <= 8
 *
 * @param dst: Destination byte
 * @param src: Source byte
 * @param sel: Number of low bits to select from src
 * @param msb_off: Offset from msb in dst
 */
void byte_merge(std::uint8_t& dst, const std::uint8_t src,
		std::uint8_t sel, std::uint8_t msb_off);

/**
 * Copy an arbitrary number of *bits* from an arbitrary position (expressed
 * in bits) in the src buffer, to an abitrary position (expressed in *bits*) in
 * the dst buffer.
 * The caller must ensure that dst_offbits < 8 && src_offbits < 8. If you need
 * to have a bit offset > 8, set dst += bit offset / 8 and
 * dst_offbits = bit offset % 8 (respectively for src). low_sel may be >= 8.
 *
 * @param dst: Destination buffer
 * @param src: Source buffer
 * @param low_sel: Number of bits to copy from src to dst
 * @params src_offbits: [0,7] Offset in the first byte of the source buffer
 * @params dst_offbits: [0,7] Offset in the second byte of the destination
 * 		buffer
 */
void bitscpy(std::uint8_t* dst, const std::uint8_t* src, unsigned low_sel,
		std::uint8_t dst_offbits, std::uint8_t src_offbits = 0);

/* Components interleaving functions	 */

typedef std::vector<std::pair<unsigned, unsigned>> compblk_spec;
typedef std::vector<compblk_spec> interleave_spec;

inline unsigned simd_block_bitsize(interleave_spec spec, int simd_size) {
	unsigned simd_block_bitsize = 0;
	for(const auto& comp_spec: spec) {
		for(const auto& s: comp_spec) {
			simd_block_bitsize += s.second;
		}
	}
	return simd_block_bitsize * simd_size;
}

inline unsigned interleaved_partition_size(unsigned simd_block_bitsize,
		int pqcode_n, int inter_n) {
	// Number of simd_blocks (rounded up)
	const unsigned simd_block_count = (pqcode_n + inter_n - 1) / inter_n;
	const unsigned total_bits = simd_block_bitsize * simd_block_count;
	// Size in bytes
	unsigned total_size = (total_bits + 8 - 1) / 8;
	return total_size;
};

/**
 * Interleave a partition.
 *
 * @param src: Source buffer
 * @param dst: Destination buffer
 * @param pqcode_n: Number of pqcodes in the partition
 * @param inter_n: Number of pqcodes in an SIMD block
 * @param spec: Interleave specification
 */
void interleave_partition(std::uint8_t* dst, const std::uint8_t* src, int pqcode_n,
		int inter_n, interleave_spec spec, const pq_params& pqp);

/**
 * Interleave a component block.
 * PQ nx8 is asumed.
 *
 * @param dst: Destination buffer
 * @param bitoff: [0,7] Offset in bits in the first byte of the destination
 * 	buffer
 * @param src: Source buffer
 * @param n: Number of vectors in the component block
 * @param spec: Component block specification
 * @param pqp: PQ Parameters
 * @param src_n: Number of vectors in the source buffer. If src_n < n,
 * the component block will be padded using the last vector of the source
 * buffer.
 */
void interleave_comp_block(std::uint8_t*& dst, unsigned& bitoff,
		std::uint8_t* src, int n, const compblk_spec& spec, const pq_params& pqp,
		int src_n);

/* Grouping functions 		*/

template<typename T>
void apply_perms(unsigned* perms, T* arr, unsigned long n) {
	T* bck_arr;
	const unsigned long sz = n * sizeof(T);
	checked_aligned_alloc((void**) &bck_arr, sz);
	memcpy(bck_arr, arr, sz);
	for (unsigned i = 0; i < n; ++i) {
		if(perms[i] > n) {
			std::cout << i << " " << perms[i] << std::endl;
			exit(1);
		}
		arr[i] = bck_arr[perms[i]];
	}
	free(bck_arr);
}

template<typename T>
void apply_perms(T* dst, unsigned n, const unsigned* perms, int scale = 1, bool
		inverse_perms = false) {
	// Copy original array
	std::unique_ptr<T[]> src(new T[n*scale]);
	memcpy(src.get(), dst, n*sizeof(T)*scale);
	for(unsigned i = 0; i <n ; ++i) {
		// Select source and destination indexes
		unsigned dst_i;
		unsigned src_i;
		if(!inverse_perms) {
			dst_i = perms[i];
			src_i = i;
		} else {
			src_i = perms[i];
			dst_i = i;
		}
		// Copy data
		if(scale == 1) {
			dst[dst_i] = src[src_i];
		} else {
			memcpy(dst + dst_i*scale, src.get() + src_i*scale, scale);
		}
	}
}

void group_partition_simple(uint16_t* partition, unsigned long n,
		std::queue<int> grp_n_q, int grp_coord);

/* New grouping functions */
typedef std::queue<std::tuple<unsigned, unsigned>> grouping_spec;

void group_partition(std::uint8_t* partition, unsigned long n, unsigned* perms,
		std::vector<unsigned>& group_sizes, grouping_spec spec, const pq_params& pqp);

std::unique_ptr<unsigned[]> perms_alloc(unsigned n);

/* Full layout functions 	*/

struct group_header {
	std::uint32_t size;
	std::uint8_t values[4];
};

struct grouping_interleave_spec {
	grouping_spec group_spec;
	interleave_spec inter_spec;
	int simd_size;
};

extern const grouping_interleave_spec group_inter_s1;

template<typename T>
struct unique_buffer {
	size_t sz;
	std::unique_ptr<T[]> buf;
};

std::unique_ptr<std::uint8_t[]> group_interleave_partition(
		const std::uint8_t* const_partition, unsigned long n, unsigned* perms,
		const pq_params& pqp, const grouping_interleave_spec& spec);

unique_buffer<std::uint8_t> layout_partition(
		const std::uint8_t* const_partition, unsigned long n, unsigned* perms,
		const pq_params& pqp, const grouping_interleave_spec& spec,
		std::uint32_t keep);

#endif /* LAYOUT_HPP_ */
