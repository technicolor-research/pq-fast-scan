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

#ifndef FASTSCAN_HPP_
#define FASTSCAN_HPP_

namespace Counters {
	extern long total_scan;
	extern float quant_bound;
};

void scan_partition_1(const std::uint8_t* partition, const unsigned* labels,
		const float* dists, todo_binheap* bh);

#endif /* FASTSCAN_HPP_ */
