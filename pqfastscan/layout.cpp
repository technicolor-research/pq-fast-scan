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
#include <cstring>
#include <vector>
#include <queue>
#include <memory>
#include <tuple>
#include <algorithm>
#include <iterator>
#include "layout.hpp"
#include "common.hpp"

/* Bits copy functions		*/

void byte_merge(std::uint8_t& dst, std::uint8_t src,
		std::uint8_t sel, std::uint8_t msb_off) {
	// Select (mask) bits from src and move them to the right position
	const std::uint8_t mask = BITMASKU8(sel);
	src = src & mask;
	std::uint8_t shift = 8 - (msb_off + sel);
	src = src << shift;
	// Clear bits in dst
	dst &= ~(mask << shift);
	// Merge
	dst |= src;
}

void bitscpy(std::uint8_t* dst, const std::uint8_t* src,
		unsigned select, std::uint8_t dst_bitoff, std::uint8_t src_bitoff) {
	while(select != 0) {
		// Number of bits to select from the first byte of src, and write
		// to the first byte of dst
		std::uint8_t select_for_byte = std::min(select, 8U - dst_bitoff);
		select_for_byte = std::min((unsigned) select_for_byte, 8U - src_bitoff);
		const std::uint8_t offseted_src = *src >> src_bitoff;
		byte_merge(*dst, offseted_src, select_for_byte, dst_bitoff);
		// Update dst offsets
		const std::uint8_t add_dst_bitoff = select_for_byte + dst_bitoff;
		dst += add_dst_bitoff / 8;
		dst_bitoff = add_dst_bitoff % 8;
		// Update src offsets
		src_bitoff += select_for_byte;
		src += src_bitoff / 8;
		src_bitoff = src_bitoff % 8;
		select -= select_for_byte;
	}
}

/* Components interleaving functions	 */

static inline bool is_44_comp_spec(const compblk_spec& spec) {
	return (spec.size() == 2 &&
			spec[0].second == 4 &&
			spec[1].second == 4);
}

static inline bool is_8_comp_spec(const compblk_spec& spec) {
	return (spec.size() == 1 && spec[0].second == 8);
}

void interleave_comp_block_8(std::uint8_t*& dst, unsigned& bitoff,
		const std::uint8_t* src, int n, const compblk_spec& spec, const pq_params& pqp,
		int src_n) {
	assert(pqp.bitsq % 8 == 0);
	const int bytes_per_comp = pqp.bitsq / 8;
	src_n = src_n == -1 ? n : src_n;
	const std::pair<unsigned, unsigned>& pair = spec[0];
	// Interleave vectors of the source buffer
	for (int pqcode_i = 0; pqcode_i < n; ++pqcode_i) {
		// PQ nx8 is assumed here
		const std::uint8_t* comp;
		if (pqcode_i < src_n) {
			comp = src + pqcode_i * pqp.nsq + pair.first * bytes_per_comp;
		} else {
			comp = src + (src_n - 1) * pqp.nsq + pair.first * bytes_per_comp;
		}
		*dst++ = *comp;
	}
}

void interleave_comp_block_44(std::uint8_t*& dst, unsigned& bitoff,
		const std::uint8_t* src, int n, const compblk_spec& spec, const pq_params& pqp,
		int src_n) {
	assert(pqp.bitsq % 8 == 0);
	const int bytes_per_comp = pqp.bitsq / 8;
	src_n = src_n == -1 ? n : src_n;
	const std::pair<unsigned, unsigned>& p1 = spec[0];
	const std::pair<unsigned, unsigned>& p2 = spec[1];
	// Interleave vectors of the source buffer
	for (int pqcode_i = 0; pqcode_i < n; ++pqcode_i) {
		// PQ nx8 is assumed here
		const std::uint8_t* pqcode;
		if (pqcode_i < src_n) {
			pqcode = src + pqcode_i * pqp.nsq;
		} else {
			pqcode = src + (src_n - 1) * pqp.nsq;
		}
		const std::uint8_t c1 = pqcode[p1.first * bytes_per_comp];
		const std::uint8_t c2 = pqcode[p2.first * bytes_per_comp];
		*dst++ = ((c1 & 0xf) << 4) + (c2 & 0xf);
	}
}

void interleave_comp_block_generic(std::uint8_t*& dst, unsigned& bitoff,
		const std::uint8_t* src, int n, const compblk_spec& spec, const pq_params& pqp,
		int src_n) {
	assert(pqp.bitsq % 8 == 0);
	const int bytes_per_comp = pqp.bitsq / 8;
	src_n = src_n == -1 ? n : src_n;
	// Interleave vectors of the source buffer
	for (int pqcode_i = 0; pqcode_i < n; ++pqcode_i) {
		for (auto& pair : spec) {
			// PQ nx8 is assumed here
			const std::uint8_t* comp;
			if (pqcode_i < src_n) {
				comp = src + pqcode_i * pqp.nsq + pair.first * bytes_per_comp;
			} else {
				comp = src + (src_n - 1) * pqp.nsq + pair.first * bytes_per_comp;
			}
			bitscpy(dst, comp, pair.second, bitoff, 0);
			bitoff += pair.second;
			dst += bitoff / 8;
			bitoff = bitoff % 8;
		}
	}
}

void interleave_comp_block(std::uint8_t*& dst, unsigned& bitoff,
		const std::uint8_t* src, int n, const compblk_spec& spec, const pq_params& pqp,
		int src_n) {
	if(is_8_comp_spec(spec)) {
		interleave_comp_block_8(dst, bitoff, src, n, spec, pqp, src_n);
	} else if(is_44_comp_spec(spec)) {
		interleave_comp_block_44(dst, bitoff, src, n, spec, pqp, src_n);
	} else {
		interleave_comp_block_generic(dst, bitoff, src, n, spec, pqp, src_n);
	}
}

void interleave_partition(std::uint8_t* dst, const std::uint8_t* src, int pqcode_n,
		int inter_n, interleave_spec spec, const pq_params& pqp) {
	assert(pqp.bitsq % 8 == 0);
	const int bytes_per_comp = pqp.bitsq / 8;
	// Number of SIMD blocks = pqcode_n / inter_n (rounded up)
	unsigned blk_n = (pqcode_n + inter_n - 1) / inter_n;
	for (unsigned blk_i = 0; blk_i < blk_n; ++blk_i) {
		for (auto& s : spec) {
			unsigned offbits = 0;
			interleave_comp_block(dst, offbits, src, inter_n, s, pqp, pqcode_n);
		}
		// PQ nx8 is assumed here
		src += inter_n * pqp.nsq * bytes_per_comp;
		pqcode_n -= inter_n;
	}
}

/* Grouping functions 		*/

void group_partition_simple_step(std::vector<std::vector<uint16_t> >& grouped_partition,
		uint16_t* partition, unsigned long n, int grp_n, int grp_coord) {
	grouped_partition.resize(grp_n);
	for (unsigned i = 0; i < n; ++i) {
		const uint16_t group_val = partition[i * 4 + grp_coord];
		const int group_card = 65536 / grp_n;
		const int group_id = group_val / group_card;
		grouped_partition[group_id].push_back(partition[i*4]);
		grouped_partition[group_id].push_back(partition[i*4 + 1]);
		grouped_partition[group_id].push_back(partition[i*4 + 2]);
		grouped_partition[group_id].push_back(partition[i*4 + 3]);
	}
}

void group_partition_simple(uint16_t* partition, unsigned long n,
		std::queue<int> grp_n_q, int grp_coord) {
	if (grp_n_q.empty())
		return;
	int grp_n = grp_n_q.front();
	grp_n_q.pop();
	std::vector<std::vector<uint16_t> > grouped_partition;
	group_partition_simple_step(grouped_partition, partition, n, grp_n, grp_coord);
	unsigned off = 0;
	for (auto& grp : grouped_partition) {
		//grp.shrink_to_fit();
		group_partition_simple(grp.data(), grp.size()/4, grp_n_q, grp_coord + 1);
		memcpy(partition + off, grp.data(), grp.size() * sizeof(uint16_t));
		off += grp.size();
	}
}

static inline std::uint64_t get_comp_value(const std::uint8_t* pqcode, int comp_i,
                        int bytes_per_comp) {
    const std::uint8_t* const comp_start = pqcode + comp_i * bytes_per_comp;
	std::uint64_t value = 0;
	memcpy(&value, comp_start, bytes_per_comp);
	return value;
}

typedef std::vector<std::vector<unsigned>> pqcode_id_group;

void group_pqcode_ids_generic(pqcode_id_group& grouped_pqcode_ids, const std::uint8_t* partition,
		unsigned long n, unsigned group_comp, unsigned group_count,
		const pq_params& pqp) {
	grouped_pqcode_ids.resize(group_count);
	assert(pqp.bitsq % 8 == 0);
	const int bytes_per_comp = pqp.bitsq / 8;
	const int values_per_group = (1 << pqp.bitsq) / group_count;
	for(unsigned long pqcode_i = 0; pqcode_i < n; ++pqcode_i) {
		const std::uint8_t* const pqcode = partition + pqcode_i * pqp.nsq * bytes_per_comp;
		std::uint64_t value = get_comp_value(pqcode, group_comp, bytes_per_comp);
		const int group_id = value / values_per_group;
		grouped_pqcode_ids[group_id].push_back(pqcode_i);
	}
}

template<typename Uint>
void group_pqcode_ids(pqcode_id_group& grouped_pqcode_ids, const std::uint8_t* partition,
		unsigned long n, unsigned group_comp, unsigned group_count,
		const pq_params& pqp) {
	grouped_pqcode_ids.resize(group_count);
	assert(pqp.bitsq % 8 == 0);
	const int bytes_per_comp = pqp.bitsq / 8;
	assert(bytes_per_comp == sizeof(Uint));
	const int values_per_group = (1 << pqp.bitsq) / group_count;
	for(unsigned long pqcode_i = 0; pqcode_i < n; ++pqcode_i) {
		const Uint* const pqcode = reinterpret_cast<const Uint*>(partition + pqcode_i * pqp.nsq * bytes_per_comp);
		const int group_id = pqcode[group_comp] / values_per_group;
		grouped_pqcode_ids[group_id].push_back(pqcode_i);
	}
}

decltype(&group_pqcode_ids_generic) group_pqcode_ids_func(const pq_params& pqp) {
	if(pqp.bitsq == 8) {
		return group_pqcode_ids<std::uint8_t>;
	} else if(pqp.bitsq == 16) {
		return group_pqcode_ids<std::uint16_t>;
	}
	return group_pqcode_ids_generic;
}

template<typename T>
std::unique_ptr<T[]> flatten_vectors(std::vector<std::vector<T>> vectors) {
	unsigned total_size = 0;
	std::for_each(vectors.begin(), vectors.end(),
			[&total_size](const std::vector<T>& v) {
				total_size += v.size();
			});
	std::unique_ptr<T[]> buffer(new T[total_size]);
	unsigned offset = 0;
	for(const auto& vec: vectors) {
		memcpy(buffer.get() + offset, vec.data(), vec.size() * sizeof(T));
		offset += vec.size();
	}
	return buffer;
}

void group_partition(std::uint8_t* partition, unsigned long n, unsigned* perms,
		std::vector<unsigned>& group_sizes, grouping_spec spec,
		const pq_params& pqp) {
	const int bytes_per_comp = pqp.bitsq / 8;
	const int bytes_per_pqcode = bytes_per_comp * pqp.nsq;
	auto group_pqcode_ids = group_pqcode_ids_func(pqp);
	// Perform one step of grouping
	unsigned group_comp;
	unsigned group_count;
	std::tie(group_comp, group_count) = spec.front();
	spec.pop();
	pqcode_id_group grouped_pqcode_ids;
	group_pqcode_ids(grouped_pqcode_ids, partition, n, group_comp, group_count, pqp);
	// Apply grouping
	std::unique_ptr<unsigned[]> flat_step_perms = flatten_vectors(
			grouped_pqcode_ids);
	apply_perms(partition, n, flat_step_perms.get(), bytes_per_comp * pqp.nsq,
			true);
	apply_perms(perms, n, flat_step_perms.get(), 1, true);
	// Check if "leaf", stop here if it is the case
	if (spec.empty()) {
		for (const auto& group : grouped_pqcode_ids) {
			group_sizes.push_back(group.size());
		}
		return;
	}
	// Otherwise, group "child" buffers
	unsigned offset = 0;
	for (const auto& group : grouped_pqcode_ids) {
		// It's vital that spec is **passed by copy** here
		group_partition(partition + offset * bytes_per_pqcode, group.size(), perms + offset,
				group_sizes, spec, pqp);
		offset += group.size();
	}
}

std::unique_ptr<std::uint8_t[]> group_partition_preserve(
		const std::uint8_t* const_partition, unsigned long n, unsigned* perms,
		std::vector<unsigned>& group_sizes, const grouping_spec& spec,
		const pq_params& pqp) {
	const int bytes_per_comp = pqp.bitsq / 8;
	const int bytes_per_pqcode = bytes_per_comp * pqp.nsq;

	// Copy partition
	std::unique_ptr<std::uint8_t[]> partition_ptr(
			new std::uint8_t[n * bytes_per_pqcode]);
	std::uint8_t* partition = partition_ptr.get();
	memcpy(partition, const_partition, n * bytes_per_pqcode);

	assert(!memcmp(partition, const_partition, n * bytes_per_pqcode));

	// Group partition
	group_partition(partition, n, perms, group_sizes, grouping_spec(spec), pqp);

	return partition_ptr;
}

std::unique_ptr<unsigned[]> perms_alloc(unsigned n) {
	std::unique_ptr<unsigned[]> vec_ids(new unsigned[n]);
	for (unsigned i = 0; i < n; ++i) {
		vec_ids[i] = i;
	}
	return vec_ids;
}

void interleave_all_groups(std::uint8_t* laidout_partition,
		const std::uint8_t* partition, const std::vector<unsigned>& group_sizes,
		const interleave_spec& inter_spec, const int simd_size,
		const pq_params& pqp) {
	const unsigned block_bitsize = simd_block_bitsize(inter_spec, simd_size);
	const int bytes_per_comp = pqp.bitsq / 8;
	const int bytes_per_pqcode = bytes_per_comp * pqp.nsq;
	// Interleave groups
	std::uint16_t values = 0;
	for (const auto& sz : group_sizes) {
		if (sz == 0) {
			values += 1;
			continue;
		}
		// Store header
		//std::cout.flush();
		group_header* hdr = reinterpret_cast<group_header*>(laidout_partition);
		assert(sz < std::numeric_limits<decltype(hdr->size)>::max());
		hdr->size = static_cast<decltype(hdr->size)>(sz);
		// Careful ! >> and << have higher priority than &
		hdr->values[3] = static_cast<std::uint8_t>(values & 0xf) << 4;
		hdr->values[2] = static_cast<std::uint8_t>(values & 0xf0);
		hdr->values[1] = static_cast<std::uint8_t>((values & 0xf00) >> 4);
		hdr->values[0] = static_cast<std::uint8_t>((values & 0xf000) >> 8);
		laidout_partition += sizeof(*hdr);
		// Store values
//		std::cout << "- Interleave partition " << values << " " << sz << std::endl;
//		std::cout << "laidout_partition=" << (void*) laidout_partition
//				<< " partition=" << (void*) partition << std::endl;
		interleave_partition(laidout_partition, partition, sz, simd_size, inter_spec,
				pqp);
		laidout_partition += interleaved_partition_size(block_bitsize, sz,
				simd_size);
		partition += sz * bytes_per_pqcode;
		values += 1;
	}
	// Store last header
	group_header* hdr = reinterpret_cast<group_header*>(laidout_partition);
	hdr->size = std::numeric_limits<decltype(hdr->size)>::max();
	hdr->values[0] = 0;
	hdr->values[1] = 0;
	hdr->values[2] = 0;
	hdr->values[3] = 0;
}

unsigned all_groups_interleaved_size(const std::vector<unsigned>& group_sizes,
		const interleave_spec& inter_spec, const int simd_size) {
	unsigned total_size = 0;
	const unsigned block_bitsize = simd_block_bitsize(inter_spec, simd_size);
	for (const auto& sz : group_sizes) {
		if (sz != 0) {
			total_size += interleaved_partition_size(block_bitsize, sz,
					simd_size);
			total_size += sizeof(group_header);
		}
	}
	// Room for last header (end of buffer marker)
	total_size += sizeof(group_header);
	return total_size;
}

const std::deque<std::tuple<unsigned, unsigned>> g1{ std::make_tuple(4, 16), std::make_tuple(5, 16),
		std::make_tuple(6, 16), std::make_tuple(7, 16) };
const grouping_spec gs1 = grouping_spec(g1);
const interleave_spec is1 { { { 0, 8 } }, { { 1, 8 } }, { { 2, 8 } },
		{ { 3, 8 } }, { { 5, 4 }, { 4, 4 } }, { { 7, 4 }, { 6, 4 } } };
const grouping_interleave_spec group_inter_s1 { gs1, is1, 16 };


std::unique_ptr<std::uint8_t[]> group_interleave_partition(
		const std::uint8_t* const_partition, unsigned long n, unsigned* perms,
		const pq_params& pqp, const grouping_interleave_spec& spec) {
	// Group partition
	std::vector<unsigned> group_sizes;

	std::unique_ptr<std::uint8_t[]> partition = group_partition_preserve(
			const_partition, n, perms, group_sizes, spec.group_spec,
			pqp);

	// Allocate final partition
	unsigned total_size = all_groups_interleaved_size(group_sizes, spec.inter_spec,
			spec.simd_size);
	std::unique_ptr<std::uint8_t[]> laidout_partition_ptr(
			new std::uint8_t[total_size]);

	// Interleave groups
	interleave_all_groups(laidout_partition_ptr.get(), partition.get(), group_sizes,
			spec.inter_spec, spec.simd_size, pqp);

	return laidout_partition_ptr;
}

unique_buffer<std::uint8_t> layout_partition(const std::uint8_t* const_partition, unsigned long n, unsigned* perms,
		const pq_params& pqp, const grouping_interleave_spec& spec, std::uint32_t keep) {
	const int bytes_per_comp = pqp.bitsq / 8;
	const int bytes_per_pqcode = bytes_per_comp * pqp.nsq;

	// Group partition
	std::vector<unsigned> group_sizes;
	std::unique_ptr<std::uint8_t[]> partition = group_partition_preserve(
			const_partition + keep * bytes_per_pqcode, n - keep, perms, group_sizes,
			spec.group_spec, pqp);

	// Allocate final partition
	unsigned total_size = all_groups_interleaved_size(group_sizes, spec.inter_spec,
			spec.simd_size);
	total_size += keep * bytes_per_pqcode;
	total_size += sizeof(keep);
	std::unique_ptr<std::uint8_t[]> laidout_partition_ptr(
			new std::uint8_t[total_size]);
	std::uint8_t* laidout_partition = laidout_partition_ptr.get();

	// Copy "kept" pqcodes
	*reinterpret_cast<std::uint32_t*>(laidout_partition) = keep;
	laidout_partition += sizeof(keep);
	std::copy(const_partition, const_partition + keep * bytes_per_pqcode, laidout_partition);
	laidout_partition += keep * bytes_per_pqcode;

	// Interleave remaining pqcodes
	interleave_all_groups(laidout_partition, partition.get(), group_sizes,
			spec.inter_spec, spec.simd_size, pqp);

	return unique_buffer<std::uint8_t> {total_size, std::move(laidout_partition_ptr)};
}
