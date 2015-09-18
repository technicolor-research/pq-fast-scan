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
#include "binheap.hpp"
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

binheap* binheap_oracle_check(binheap* oracle_bh, binheap* bh) {
	if (oracle_bh == NULL) {
		return bh;
	}
	int* labels = new int[bh->size()];
	float* dists = new float[bh->size()];
	bh->sort(labels, dists);
	int* oracle_labels = new int[oracle_bh->size()];
	float* oracle_dists = new float[oracle_bh->size()];
	oracle_bh->sort(oracle_labels, oracle_dists);
	if (oracle_bh->size() != bh->size()) {
		exit_display_binheap_elems(oracle_labels, oracle_dists, oracle_bh->size(), labels, dists, bh->size(), -1);
	}
	for (int i = 0; i < oracle_bh->size(); ++i) {
		if (abs(oracle_dists[i] - dists[i]) > 0.1) {
			std::cout << abs(oracle_dists[i] - dists[i]) << std::endl;
			exit_display_binheap_elems(oracle_labels, oracle_dists, oracle_bh->size(), labels, dists, bh->size(), i);
		}
	}
	delete[] oracle_labels;
	delete[] oracle_dists;
	delete[] labels;
	delete[] dists;
	delete bh;
	return oracle_bh;
}

binheap* time_func_binheap_display(binheap_scan_func func,  int k, binheap* oracle_data,
		const char* desc, int repeat) {
	std::function<binheap*(binheap*, binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<binheap*()> setup_f = [k] () {
		return new binheap(k);
	};
	return time_func_display(func, check_f, setup_f, oracle_data, desc, repeat);
}

binheap* time_func_binheap(binheap_scan_func func, int k,
		binheap* oracle_data, int repeat, unsigned long& time) {
	std::function<binheap*(binheap*, binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<binheap*()> setup_f = [k] () {
		return new binheap(k);
	};
	return time_func(func, check_f, setup_f, oracle_data, time, repeat);
}

binheap* perf_func_binheap(binheap_scan_func func, int k,
		binheap* oracle_data,
		std::uint64_t event_values[],
		int repeat, const char* events[], int event_count) {
	std::function<binheap*(binheap*, binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<binheap*()> setup_f = [k] () {
		return new binheap(k);
	};
	return perf_func(func, check_f, setup_f, oracle_data,
			event_values, repeat, events, event_count);
}

binheap* perf_func_binheap_display(binheap_scan_func func, int k, binheap* oracle_data,
		const char* desc, int repeat, const char* events[],
		int event_count) {
	std::function<binheap*(binheap*, binheap*)> check_f;
	check_f = binheap_oracle_check;
	std::function<binheap*()> setup_f = [k] () {
		return new binheap(k);
	};
	return perf_func_display(func, check_f, setup_f, oracle_data, desc, repeat, events, event_count);
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
