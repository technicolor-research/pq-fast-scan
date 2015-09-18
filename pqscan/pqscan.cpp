//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include <functional>
#include "config.h"
#include "common.hpp"
#include "populate.hpp"
#include "benchmark.hpp"
#include "scan_avx.hpp"
#include "scan_gather.hpp"
#include "scan_naive.hpp"

#define NSQ 8
#define BITSQ 8
#define NCENT 256

using namespace std;
using namespace std::placeholders;

const unsigned long n = 25*1000*1000;
const int repeat = 30;

static const char* events[] = { "cycles", "instructions", "L1-dcache-loads",
		"L1-dcache-load-misses", "CYCLE_ACTIVITY:CYCLES_LDM_PENDING",
		"CYCLE_ACTIVITY:CYCLES_NO_EXECUTE",
		"CYCLE_ACTIVITY:STALLS_LDM_PENDING",
		"RESOURCE_STALLS:ANY",
		"UOPS_RETIRED:ANY"};
const int event_count = 9;

int main(int argc, char* argv[]) {
	auto bench_f = std::bind(perf_func_binheap_display, _1, _2, _3, _4, _5, events, event_count);

	// PQ Parameters
	pq_params pqp;
	pqp.nsq = NSQ;
	pqp.bitsq = BITSQ;

	// Linear partition
	char* partition = partition_new(n, pqp);
	float* dists = dists_new(pqp);
	partition_populate(partition, n, pqp);
	dists_populate(dists, pqp);

	// Interleaved partition for VGATHER
	// n+1 is for gather_shift2 as we may read some garbage at the end
	char* i32_partition = partition_new(n+1, pqp);
	// memset is also for gather_shift2 as we may read garbage at the end
	memset(i32_partition, 0, (n+1)*NSQ);
	partition_interleave_pqcodes(i32_partition, partition, n, pqp, 32);

	//
	const int k = 10;
	binheap* oracle = NULL;

	// Naive Scan
	oracle = bench_f(bind(scan_bh, partition, dists, n, pqp, _1), k, oracle,
			"Scan", repeat);

	// SSE Scan
	bench_f(bind(scan_prefetch_sse, partition, dists, n, pqp, _1),
			k, oracle, "Scan +sse +prefetch", repeat);

	// AVX Scan
	#ifdef AVX
	bench_f(bind(scan_bh_prefetch_avx, partition, dists, n, pqp, _1),
			k, oracle, "Scan +avx +prefetch", repeat);
	#endif

	// VGATHER Scan
	#ifdef AVX2

	bench_f(bind(scan_bh_avx_gths, i32_partition, dists, n, pqp, _1),
			k, oracle, "Scan +avx +vgather +prefetch", repeat);
	#endif

	free(partition);
	free(i32_partition);
	free(dists);
	delete oracle;

}
