//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef SCAN_SSE_HPP_
#define SCAN_SSE_HPP_

#include "common.hpp"
#include "todo_binheap.hpp"

void scan_prefetch_sse(const char* partition, const float* dists,
		unsigned long n, pq_params pqp, todo_binheap* bh);

#endif /* SCAN_SSE_HPP_ */
