//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <sys/time.h>
#include <functional>
#include "common.hpp"
#include "perfevents.hpp"
#include "todo_binheap.hpp"
#include "config.h"

using namespace std;

inline std::uint64_t get_us_clock() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (std::uint64_t) tv.tv_sec * 1000000 + tv.tv_usec;
}

template<typename T>
T* time_func(std::function<void(T*)> bench_f, std::function<T*(T*, T*)> check_f,
		std::function<T*()> setup_f, T* oracle_data, unsigned long& time, int repeat) {
	T* out_data = NULL;
	time = 0;
	for (int r = 0; r < repeat; ++r) {
		if(setup_f != NULL) {
			out_data = setup_f();
		}
		unsigned long before_us = get_us_clock();
		bench_f(out_data);
		time += get_us_clock() - before_us;
		if (check_f != NULL) {
			out_data = check_f(oracle_data, out_data);
		}
	}
	time = time / repeat;
	return out_data;
}

template<typename T>
T* time_func_display(std::function<void(T*)> bench_f, std::function<T*(T*, T*)> check_f,
		std::function<T*()> setup_f, T* oracle_data, const char* desc, int repeat) {
	cout << desc << ": ";
	cout.flush();
	unsigned long time;
	T* out_data = time_func(bench_f, check_f, setup_f, oracle_data, time, repeat);
	cout << time << "us" << endl;
	return out_data;
}

typedef std::function<void(todo_binheap*)> binheap_scan_func;

todo_binheap* time_func_binheap_display(binheap_scan_func func,  int k, todo_binheap* oracle_data,
		const char* desc, int repeat);

todo_binheap* time_func_binheap(binheap_scan_func func, int k,
		todo_binheap* oracle_data, int repeat, unsigned long& time);

void simple_time_display(function<void()> func, const char* desc, int repeat);

struct perf_events {
	std::uint64_t cycles;
	std::uint64_t intructions;
	std::uint64_t l1d_load;
	std::uint64_t l1d_miss;
	std::uint64_t time;
};

//static const char* events[] = { "cycles", "instructions", "MEM_LOAD_UOPS_RETIRED:L1_HIT",
//		"MEM_LOAD_UOPS_RETIRED:L1_MISS"};
//const int event_count = 4;

template<typename T>
void divide_array(T* array, int size, int divisor) {
	for(int i = 0; i < size; ++i) {
		array[i] /= divisor;
	}
}

// libpfm and Linux perf_events - New Hardware Performance Counters interface
template<typename T>
T* perf_func(std::function<void(T*)> bench_f, std::function<T*(T*, T*)> check_f,
		std::function<T*()> setup_f, T* oracle_data,
		std::uint64_t event_values[], int repeat, const char* events[],
		int event_count) {
	T* out_data = NULL;
	event_values[0] = 0;
	int perf_events_fds[event_count];
	perf_events_open_by_name(perf_events_fds, event_count, events);
	for (int r = 0; r < repeat; ++r) {
		if (setup_f != NULL) {
			out_data = setup_f();
		}

		const unsigned long before_us = get_us_clock();

		perf_events_start();
		bench_f(out_data);
		perf_events_stop();

		event_values[0] += get_us_clock() - before_us;
		if (check_f != NULL) {
			out_data = check_f(oracle_data, out_data);
		}
	}
	//
	std::uint64_t raw[event_count][3];
	perf_events_read(perf_events_fds, event_count, event_values + 1, raw);
	divide_array(event_values + 1, event_count, repeat);
	event_values[0] = event_values[0] / repeat;
	perf_events_close(perf_events_fds, event_count);
	return out_data;
}

template<typename T>
T* perf_func_display(std::function<void(T*)> bench_f, std::function<T*(T*, T*)> check_f,
		std::function<T*()> setup_f, T* oracle_data, const char* desc, int repeat, const char* events[],
		int event_count) {
	cout << desc << ": ";
	cout.flush();
	std::uint64_t event_values[event_count + 1];
	T* out_data = perf_func(bench_f, check_f, setup_f, oracle_data, event_values, repeat, events, event_count);
	std::cout << "time=" <<  event_values[0] << ",";
	for(int event_id = 0; event_id < event_count - 1; ++event_id) {
		std::cout << events[event_id] << "=" << event_values[event_id + 1] << ",";
	}
	std::cout << events[event_count - 1] <<  "="  << event_values[event_count] << std::endl;
	return out_data;
}

todo_binheap* perf_func_binheap(binheap_scan_func func, int k,
		todo_binheap* oracle_data,
		std::uint64_t event_values[],
		int repeat, const char* events[], int event_count);

todo_binheap* perf_func_binheap_display(binheap_scan_func func, int k, todo_binheap* oracle_data,
		const char* desc, int repeat, const char* events[],
		int event_count);

#endif /* BENCHMARK_HPP_ */
