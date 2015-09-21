//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef POPULATE_HPP_
#define POPULATE_HPP_

#include <cstdlib>
#include "common.hpp"
#include <queue>



// Allocate a new PQ partition comprising n pqcodes
char* partition_new(unsigned long n, const pq_params& pqp);

// Populate PQ partition with random data
void partition_populate(char* partition, unsigned long n, const pq_params& pqp);

// Allocate a new distance-to-centroids table
float* dists_new(pq_params& pqp);

// Populate a distance-to-centroids table
void dists_populate(float* table, pq_params& pqp);

void partition_hist_check(char* buffer, unsigned long n, pq_params& pqp);

// Interleaved partitions
void partition_interleave_pqcodes(char* i_partition, const char* partition,
		unsigned long n, pq_params& pqp, unsigned k);

#endif /* POPULATE_HPP_ */
