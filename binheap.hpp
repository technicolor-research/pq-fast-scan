//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef BINHEAP_HPP_
#define BINHEAP_HPP_

#include <memory>
#include <utility>
#include <algorithm>

class binheap {

private:
	typedef std::pair<int, float> TupleType;
	std::unique_ptr<TupleType[]> data_;
	int capacity_;
	int size_;

	static bool comparator(const TupleType a, const TupleType& b) {
		return a.second < b.second;
	}

public:
	binheap(int capacity) :
		capacity_(capacity), size_(0) {
		data_.reset(new std::pair<int, float>[capacity_+1]);
	}

	int capacity() const {
		return capacity_;
	}

	int size() const {
		return size_;
	}

	float max() const {
		return data_[0].second;
	}

	void push(int id, float value);

	void sort(int ids[], float values[]);
};

#endif /* BINHEAP_HPP_ */
