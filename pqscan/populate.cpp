//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <array>
#include <queue>
#include <memory>
#include "populate.hpp"

using namespace std;

static mt19937 rnd_eng;
static uniform_int_distribution<int> int_dist;
static uniform_real_distribution<float> float_dist;

static void print_hist(unsigned long hist[256]) {
	cout << "Hist : ";
	for (int i = 0; i < 256; ++i) {
		cout << "[" << i << "] " << cout.flush();
	}
	cout << endl;
}

auto ulongcmp = [](const void* a, const void* b) {
	const unsigned long *ia = (const unsigned long *)a;
	const unsigned long *ib = (const unsigned long *)b;
	if(*ia < *ib) {
		return -1;
	} else if (*ia > *ib) {
		return 1;
	}
	return 0;
};

static void buffer_hist_check(char* buf, unsigned long sz) {
	unsigned long hist[256] = { 0 };
	unsigned char* ubuf = (unsigned char*) buf;
	while (sz--) {
		hist[*ubuf++]++;
	}
	unsigned long shist[256];
	memcpy(&shist, &hist, 256 * sizeof(unsigned long));
	qsort(&shist, 256, sizeof(unsigned long), ulongcmp);
	double ratio = (double) shist[255] / shist[0];
	if (ratio > 1) {
		print_hist(shist);
		exitmsg("Bad random data in partition");
	}
}

char* partition_new(unsigned long n, const pq_params& pqp) {
	char* partition;
	unsigned long sz = partitionsz(n, pqp);
	checked_aligned_alloc((void**) &partition, sz);
	return partition;
}

void partition_slice_populate(char* slice, unsigned long sz) {
	char* buf_end = slice + sz;
	const int s = sizeof(uint32_t);
	while (slice != buf_end) {
		*((uint32_t*) slice) = (uint32_t) rnd_eng();
		slice += s;
	}
}

void partition_populate(char* partition, unsigned long n, const pq_params& pqp) {
	unsigned long sz = partitionsz(n, pqp);
	partition_slice_populate(partition, sz);
}

void partition_interleave_pqcodes(char* i_partition, const char* partition,
		unsigned long n, pq_params& pqp, unsigned k) {
	// i : block number
	// j : component number
	// l : pqcode in block
	unsigned long blk_n = n / k;
	for (unsigned long i = 0; i < blk_n; ++i) {
		for (int j = 0; j < pqp.nsq; ++j) {
			for(unsigned l = 0; l < k; ++l) {
				const unsigned long block_off = i * pqp.nsq * k;
				const unsigned long comp_block_off = j * k;
				const unsigned long off = block_off + comp_block_off + l;
				const unsigned long pqcode_id = i * k + l;
				// Pas très hardware prefetecher compliant
				i_partition[off] = partition[pqcode_id * pqp.nsq + j];
			}
		}
	}
}

float* dists_new(pq_params& pqp) {
	float* table;
	int sz = pqp.nsq * ncentsq(pqp) * sizeof(float);
	checked_aligned_alloc((void**) &table, sz);
	return table;
}

void dists_populate(float* table, pq_params& pqp) {
	for (unsigned long i = 0; i < pqp.nsq * (unsigned long) ncentsq(pqp); i++) {
		table[i] = float_dist(rnd_eng);
	}
}

void partition_hist_check(char* buffer, unsigned long n, pq_params& pqp) {
	buffer_hist_check(buffer, partitionsz(n, pqp));
}


