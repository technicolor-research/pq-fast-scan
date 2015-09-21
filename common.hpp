//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <immintrin.h>

using namespace std;

enum color {txtrst = 0, bldblk, bldred};
extern const char* colors[];

void print_m128i(__m128i reg);
void print_uint64(std::uint64_t reg);
void print_buffer(const std::uint8_t* buf, unsigned size);

struct pq_params {
	int nsq;		// Number of sub-quantizers
	int bitsq;		// Bits per sub-quantizer
};

#define ALIGN_BYTES 4096
inline void checked_aligned_alloc(void** buf, unsigned long sz) {
	int err = posix_memalign((void**) buf, ALIGN_BYTES, sz);
	assert(err == 0);
}
;

inline void cprint(const char* str, color c) {
	cout << colors[c] << str << colors[txtrst];
}

// Returns number of centroids per sub-quantizer
inline int ncentsq(pq_params& pqp) {
	return pow(2, pqp.bitsq);
}

// Returns the size of a PQ pqcode
inline int pqcodesz(const pq_params& pqp) {
	return pqp.nsq * (pqp.bitsq/8);
}

// Returns the size of a PQ partition based on the number of pqcodes
inline unsigned long partitionsz(unsigned long n, const pq_params& pqp) {
	return n*pqcodesz(pqp);
}

// Returns the size of a sub-pqcode in a pqcode
inline unsigned long pqcodeslicesz(unsigned long n, pq_params& pqp, int s) {
	const int c = pqp.nsq / s;
	return c * (pqp.bitsq/8);
}

// Returns the size of a slice of a PQ sliced partition based on
// the number of pqcodes and the number of slices
inline unsigned long slicesz(unsigned long n, pq_params& pqp, int s) {
	return n * pqcodeslicesz(n, pqp, s);
}

inline void exitmsg(const char* msg, int code = 1) {
	cout << bldred << msg << txtrst << endl;
	exit(code);
}

#endif /* COMMON_HPP_ */
