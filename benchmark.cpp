//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include <cfloat>
#include <cstdlib>
#include <string>
#include <functional>
#include "benchmark.hpp"
#include "todo_binheap.hpp"
#include "config.h"

using namespace std;
using namespace std::placeholders;

void display_binheap_elems(int* labels, float* dists, int k, int idiff) {
	cout << "[" << k << "] " << endl;
	for (int i = 0; i < k; ++i) {
		if(i != idiff) {
			cout << labels[i] << ":" << dists[i] << " | ";
		} else {
			cout << colors[bldred] << labels[i] << ":" << dists[i] << colors[txtrst]<< " | ";
		}
	}
	cout << endl;
}

void exit_display_binheap_elems(int* oracle_labels, float* oracle_dists,
		int oracle_k, int* labels, float* dists, int k, int i) {
	cout << "Oracle binheap and current binheap are different." << endl;
	cout << "Oracle binheap : ";
	display_binheap_elems(oracle_labels, oracle_dists, oracle_k, i);
	cout << "Current binheap : ";
	display_binheap_elems(labels, dists, k, i);
	exitmsg("Aborting.");
}

void binheap_sort_clean(todo_binheap* bh, int*& labels_out, float*& dists_out) {
	bh.todo_sort(labels_out, dists_out);
}

todo_binheap* binheap_oracle_check(todo_binheap* oracle_bh, todo_binheap* bh) {
	if (oracle_bh == NULL) {
		return bh;
	}
	int* labels = new int[bh.todo_size()];
	float* dists = new float[bh.todo_size()];
	binheap_sort_clean(bh, labels, dists);
	int* oracle_labels = new int[oracle_bh.todo_size()];
	float* oracle_dists = new float[oracle_bh.todo_size()];
	binheap_sort_clean(oracle_bh, oracle_labels, oracle_dists);
	if (oracle_bh.todo_size() != bh.todo_size()) {
		exit_display_binheap_elems(oracle_labels, oracle_dists, oracle_bh.todo_size(), labels, dists, bh.todo_size(), -1);
	}
	for (int i = 0; i < oracle_bh.todo_size(); ++i) {
		if (abs(oracle_dists[i] - dists[i]) > 0.1) {
			std::cout << abs(oracle_dists[i] - dists[i]) << std::endl;
			exit_display_binheap_elems(oracle_labels, oracle_dists, oracle_bh.todo_size(), labels, dists, bh.todo_size(), i);
		}
	}
	delete[] oracle_labels;
	delete[] oracle_dists;
	delete[] labels;
	delete[] dists;
	delete bh;
	return oracle_bh;
}

todo_binheap* time_func_binheap_display(binheap_scan_func func,  int k, todo_binheap* oracle_data,
		const char* desc, int repeat) {
	std::function<todo_binheap*(todo_binheap*, todo_binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<todo_binheap*()> setup_f = [k] () {
		return new todo_binheap(k);
	};
	return time_func_display(func, check_f, setup_f, oracle_data, desc, repeat);
}

todo_binheap* time_func_binheap(binheap_scan_func func, int k,
		todo_binheap* oracle_data, int repeat, unsigned long& time) {
	std::function<todo_binheap*(todo_binheap*, todo_binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<todo_binheap*()> setup_f = [k] () {
		return new todo_binheap(k);
	};
	return time_func(func, check_f, setup_f, oracle_data, time, repeat);
}

todo_binheap* perf_func_binheap(binheap_scan_func func, int k,
		todo_binheap* oracle_data,
		std::uint64_t event_values[],
		int repeat, const char* events[], int event_count) {
	std::function<todo_binheap*(todo_binheap*, todo_binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<todo_binheap*()> setup_f = [k] () {
		return new todo_binheap(k);
	};
	return perf_func(func, check_f, setup_f, oracle_data,
			event_values, repeat, events, event_count);
}

void simple_time_display(function<void()> func, const char* desc, int repeat) {
	std::function<void(void*)> f = [&func](void*) {
		func();
	};
	std::function<void*(void*, void*)> check_f = NULL;
	std::function<void*()> setup_f = NULL;
	void* null_ptr = NULL;
	time_func_display(f, check_f, setup_f, null_ptr, desc, repeat);
}
