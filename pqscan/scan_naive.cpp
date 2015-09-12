//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include "scan_naive.hpp"

#include <cfloat>
#include <cstdint>
#include <immintrin.h>

#define NSQ 8
#define NCENT 256

void scan_bh(const char* partition, const float* dists,
		unsigned long n, pq_params pqp, todo_binheap* bh) {
	long binheap_op = 0;
	for (int t = 0; t < bh.todo_capacity(); ++t) {
		bh.todo_add(0, FLT_MAX - t);
	}
	float min = bh.todo_peek();
	float candidate;
	uint8_t* u8_buf = (uint8_t*) partition;
	for (unsigned long i = 0; i < n; ++i) {
		uint8_t* pqcode = u8_buf + i * NSQ;
		candidate = 0;
		for (unsigned j = 0; j < NSQ; ++j) {
			candidate += dists[j * NCENT + pqcode[j]];
		}
		if (candidate < min) {
			bh.todo_add(i, candidate);
			min = bh.todo_peek();
			binheap_op++;
		}
	}
}

