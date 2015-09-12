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
#include <cstdarg>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <functional>

#include <unistd.h>

#include "pqscan/scan_naive.hpp"
#include "benchmark.hpp"
#include "layout.hpp"
#include "common.hpp"
#include "fastscan.hpp"

#define NSQ 8
#define NCENT 256

using namespace std::placeholders;

// Benchmarking and display functions
template<typename T>
void display_csv(T array[], int size) {
	for(int i = 0; i < size - 1; ++i) {
		std::cout << array[i] << ",";
	}
	std::cout << array[size-1];
}

static const int repeat = 1;

typedef todo_binheap* (*benchmark_csv_func)(binheap_scan_func func, int k,
		todo_binheap* oracle_data);

static const char* events[] = {"cycles", "instructions", "L1-dcache-loads"};
const int event_count = 3;

static const char* csv_header_perf = "pq_us,pq_cycles,pq_instructions,pq_l1_loads,fast_pq_us,fast_pq_cycles,fast_pq_instructions,fast_pq_l1_loads";
todo_binheap* perf_func_csv(binheap_scan_func func, int k,
		todo_binheap* oracle_data) {
	int total_event_count = event_count + 1;
	std::uint64_t* event_values = new std::uint64_t[total_event_count];
	todo_binheap* ret = perf_func_binheap(func, k, oracle_data, event_values, repeat,
			events, event_count);
	display_csv(event_values, total_event_count);
	delete[] event_values;
	std::cout.flush();
	return ret;
}

static const char* csv_header_time = "pq_us,fast_pq_us";
todo_binheap* time_func_csv(binheap_scan_func func, int bh_size,
		todo_binheap* bh_oracle) {
	unsigned long us_time;
	todo_binheap* ret = time_func_binheap(func, bh_size, bh_oracle, repeat, us_time);
	std::cout << us_time;
	std::cout.flush();
	return ret;
}

// Parse command line
struct cmdargs {
	float keep_percent;
	unsigned bh_size;
	const char* partition_file;
	const char* query_set_file;
	const char* query_id_file;
	//
	benchmark_csv_func bench_func;
	const char* bench_header;
};

void usage(const char* progname) {
	std::cout << progname << ": [-k keep] [-b binheap_size] partition query_set query_ids" << std::endl;
	std::exit(1);
}

void parse_args(cmdargs& args, int argc, char* argv[]) {
	int opt;
	args.keep_percent = 0.01;
	args.bh_size = 100;
	args.bench_header = csv_header_time;
	args.bench_func = &time_func_csv;
	while ((opt = getopt(argc, argv, "k:b:p")) != -1) {
		switch (opt) {
		case 'k':
			args.keep_percent = std::atof(optarg);
			break;
		case 'b':
			args.bh_size = std::atoi(optarg);
			break;
		case 'p':
			args.bench_header = csv_header_perf;
			args.bench_func = perf_func_csv;
			break;
		default: /* '?' */
			usage(argv[0]);
		}
	}
	if(argc - optind < 3) {
		usage(argv[0]);
	}
	args.partition_file = argv[optind];
	args.query_set_file = argv[optind + 1];
	args.query_id_file = argv[optind + 2];
}

struct partition {
	unique_buffer<std::uint8_t> laidout_partition;
	std::unique_ptr<unsigned[]> permuted_labels;
	std::uint8_t* partition;
	int id;
	unsigned long n;
	unsigned keep;
};

struct query_set {
	std::vector<float> vectors;
	std::vector<unsigned> ids;
	std::vector<float> distance_tables;
	int dim;
};

void load_query_vectors(cmdargs& args,
		query_set& query) {

	// Load query vectors set
	long query_vecs_count;
	float* query_byte_vectors = todo_read_byte_vectors(args.query_set_file, &query_vecs_count, &query.dim);
	std::cerr << query.dim << " " << query_vecs_count << std::endl;

	// Load distance tables
	query.distance_tables = todo_load_distance_tables();

	// Copy selected vectors
	std::ifstream query_ids(args.query_id_file);
	int id;
	while (query_ids >> id) {
		if (id < 0 || id >= query_vecs_count) {
			std::cerr << "Query ids must be comprised between 0 and "
					<< query_vecs_count << std::endl;
			std::exit(1);
		}
		query.vectors.insert(query.vectors.end(), query_byte_vectors + id * query.dim,
				query_byte_vectors + (id + 1) * query.dim);
		query.vectors.push_back(id);
	}
	free(query_byte_vectors);
}

void process_query_vectors(cmdargs& args, partition& part, query_set& query) {
	int distance_table_size = NCENT * NSQ;
	int query_n = query.vectors.size() / query.dim;

	// Display header
	std::cout << "vec_id,partition_id,partition_n,bh_size,keep,";
	std::cout << args.bench_header << ",";
	std::cout << "fast_pq_scanned,quant_bound" << std::endl;

	// Process query vectors
	for(int i = 0; i < query_n; ++i) {
		float* dist_table = query.distance_tables.data() + i * distance_table_size;
		// "Global" counters
		Counters::total_scan = 0;
		Counters::quant_bound = 0;

		//
		std::cout << query.ids[i] << "," << part.id << "," << part.n << "," <<
				args.bh_size << "," << part.keep << ",";

		// Normal PQ Scan
		pq_params pqp {8, 8};
		todo_binheap* bh_oracle = nullptr;
		std::bind(pq_binheap_scan,
				reinterpret_cast<const char*>(part.partition),
				dist_table, part.n, pqp, _1);
		bh_oracle = args.bench_func(
				std::bind(scan_bh,
							reinterpret_cast<const char*>(part.partition),
							dist_table, part.n, pqp, _1), args.bh_size, bh_oracle);
		std::cout << ",";
		// Fast PQ Scan
		args.bench_func(
				std::bind(scan_partition_1, part.laidout_partition.buf.get(),
						part.permuted_labels.get(), dist_table, _1),
				args.bh_size, bh_oracle);

		std::cout <<  "," << Counters::total_scan << "," << Counters::quant_bound << std::endl;
	}
}

int main(int argc, char* argv[]) {
	// Parse arguments
	cmdargs args;
	parse_args(args, argc, argv);
	// Load partition
	partition& part;
	todo_load_partition(args, part);
	std::cerr << "Loaded partition"<< std::endl ;
	// Load query vectors
	query_set& query;
	load_query_vectors(args, query);
	std::cerr << "Loaded " << query.vectors.size() << " query vectors" << std::endl;
	process_query_vectors(args, part, query);
	free(ivf);
}
