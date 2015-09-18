//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include "binheap.hpp"

void binheap::push(int id, float value) {
	if (size_ == capacity_) {
		if (data_[0].second < value) {
			return;
		}
		std::pop_heap(data_.get(), data_.get() + size_, binheap::comparator);
		size_--;
	}
	data_[size_].first = id;
	data_[size_].second = value;
	size_++;
	std::push_heap(data_.get(), data_.get() + size_, binheap::comparator);
}

void binheap::sort(int ids[], float values[]) {
	// Sort tuples
	std::unique_ptr<TupleType[]> sorted_data(new TupleType[size_]);
	std::copy(data_.get(), data_.get() + size_, sorted_data.get());
	std::sort(sorted_data.get(), sorted_data.get() + size_, binheap::comparator);
	// Unpack tuples
	for(int i = 0; i < size_; ++i) {
		ids[i] = sorted_data[i].first;
		values[i] = sorted_data[i].second;
	}
}
