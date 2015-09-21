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
#include <cstring>
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

typedef binheap* (*benchmark_csv_func)(binheap_scan_func func, int k,
		binheap* oracle_data);

static const char* events[] = {"cycles", "instructions", "L1-dcache-loads"};
const int event_count = 3;

static const char* csv_header_perf = "pq_us,pq_cycles,pq_instructions,pq_l1_loads,fast_pq_us,fast_pq_cycles,fast_pq_instructions,fast_pq_l1_loads";
binheap* perf_func_csv(binheap_scan_func func, int k,
		binheap* oracle_data) {
	int total_event_count = event_count + 1;
	std::uint64_t* event_values = new std::uint64_t[total_event_count];
	binheap* ret = perf_func_binheap(func, k, oracle_data, event_values, repeat,
			events, event_count);
	display_csv(event_values, total_event_count);
	delete[] event_values;
	std::cout.flush();
	return ret;
}

static const char* csv_header_time = "pq_us,fast_pq_us";
binheap* time_func_csv(binheap_scan_func func, int bh_size,
		binheap* bh_oracle) {
	unsigned long us_time;
	binheap* ret = time_func_binheap(func, bh_size, bh_oracle, repeat, us_time);
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
	const char* distance_tables_file;
	const char* query_id_file;
	//
	benchmark_csv_func bench_func;
	const char* bench_header;
};

void usage(const char* progname) {
	std::cerr << "Usage: " << progname << " [-p] [-k keep_ratio] [-b binheap_size]\n";
	std::cerr << "       " << std::string(std::strlen(progname), ' ') << " partition_file query_file tables_file ids_file\n";
	std::cerr << std::endl;
	std::cerr << "Benchmark PQ Fast Scan against PQ Scan and output results as CSV\n";
	std::cerr << std::endl;
	std::cerr << "Positional arguments:\n";
	std::cerr << "  partition_file\t Partition to scan (raw pqcodes)\n";
	std::cerr << "  query_file\t\t Full set of query vectors (.bvecs file format)\n";
	std::cerr << "  tables_file\t\t Full set of distance tables (.fvecs file\n";
	std::cerr << "  \t\t\t format)\n";
	std::cerr << "  ids_file\t\t IDs of vectors to load from query_file (plain\n";
	std::cerr << "  \t\t\t text, one ID per line)\n";
	std::cerr << std::endl;
	std::cerr << "Examples:\n";
	std::cerr << "  " << progname << " 100M1-partition-2.dat bigann_query.bvecs \\\n";
	std::cerr << "  \t bigann_distance_tables.fvecs 100M1-list-2.txt\n\n";
	std::cerr << "  " << progname << " -p -k 0.01 -b 200 100M1-partition-2.dat \\\n";
	std::cerr << "  \t bigann_query.bvecs bigann_distance_tables.fvecs 100M1-list-2.txt\n";
	std::cerr << std::endl;
	std::cerr << "Options:\n";
	std::cerr << "  -p\t\t\t Output performance counters (cycles,\n";
	std::cerr << "  \t\t\t intructions, L1-dcache-loads) in addition to\n";
	std::cerr << "  \t\t\t run times\n";
	std::cerr << "  -k keep_ratio\t\t Scan keep_ratio vectors to determine qmax\n";
	std::cerr << "  \t\t\t quantization bound (default: 0.005, i.e., 0.5%)\n";
	std::cerr << "  -b binheap_size\t Allocate binary heap of size binheap_size to\n";
	std::cerr << "  \t\t\t store nearest neighbors (default: 100)\n";
	std::exit(1);
}

void parse_args(cmdargs& args, int argc, char* argv[]) {
	int opt;
	args.keep_percent = 0.005;
	args.bh_size = 100;
	args.bench_header = csv_header_time;
	args.bench_func = &time_func_csv;
	while ((opt = getopt(argc, argv, "k:b:ph")) != -1) {
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
	if(argc - optind < 4) {
		usage(argv[0]);
	}
	args.partition_file = argv[optind];
	args.query_set_file = argv[optind + 1];
	args.distance_tables_file = argv[optind + 2];
	args.query_id_file = argv[optind + 3];
}

struct partition_data {
	unique_buffer<std::uint8_t> laidout_partition;
	std::unique_ptr<unsigned[]> permuted_labels;
	std::unique_ptr<std::uint8_t[]> partition;
	int id;
	unsigned long n;
	unsigned keep;
};

void check_fstream(const std::ifstream& file, const char* filename) {
	if (!file) {
		std::cerr << "Could not open " << filename << std::endl;
	}
}

int parse_number_before_extension(const char* filename) {
    filename = basename(filename);
    const char* rstart = std::strrchr(filename, '.');
    if(rstart == nullptr) {
        rstart = filename + std::strlen(filename) - 1;
    }
    // Iterate backwards until we find a number
    while(*rstart < '0' || *rstart > '9') {
        if(rstart == filename) {
            return -1;
        }
        rstart--;
    }
    // Iterate backwards until beginning of number
    while(rstart != filename
            && *(rstart -1 ) >= '0' && *(rstart -1) <= '9') {
        rstart--;
    }
    // Parse number
    return static_cast<int>(std::strtoull(rstart, nullptr, 10));
}

void print_positive_number(int num) {
	if(num < 0) {
		std::cout << '?';
	} else {
		std::cout << num;
	}
}
void check_same_id(cmdargs& args, int& partition_id) {
	partition_id = parse_number_before_extension(args.partition_file);
	int list_id = parse_number_before_extension(args.query_id_file);
	if(partition_id != list_id || partition_id == -1) {
		std::cerr << "Inconsistent partition_file and ids_file." << std::endl;
		std::cerr << "  partition_file: " << args.partition_file << " (last number=";
		print_positive_number(partition_id);
		std::cerr << ")" << std::endl;
		std::cerr << "  ids_file: " << args.query_id_file << " (last number=";
		print_positive_number(list_id);
		std::cerr << ")\n" << std::endl;
		std::cerr << "partition_file and ids_file MUST have the same last number." << std::endl;
		std::exit(1);
	}
}

void load_partition(cmdargs& args, partition_data& part) {
	// Check partition filename and query vectors ids filenames
	check_same_id(args, part.id);

	// Read original partition for PQ Scan
	std::ifstream infile(args.partition_file, std::ifstream::binary);
	check_fstream(infile, args.partition_file);
	std::int32_t pqcodes_count;
	infile.read(reinterpret_cast<char*>(&pqcodes_count), sizeof(pqcodes_count));
	part.n = pqcodes_count;
	part.partition.reset(new std::uint8_t[part.n * NSQ]);
	infile.read(
			reinterpret_cast<char*>(part.partition.get()), part.n * NSQ);

	// Layout partition for PQ Fast Scan
	pq_params pqp {8, 8};
	part.permuted_labels = perms_alloc(part.n);
	unsigned* labels = part.permuted_labels.get();
	part.keep = (unsigned) part.n * args.keep_percent;
	part.laidout_partition = layout_partition(part.partition.get(), part.n,
			labels, pqp, group_inter_s1, part.keep);
}

struct query_set {
	std::vector<float> vectors;
	std::vector<unsigned> ids;
	std::vector<float> distance_tables;
	int dimension;
};

template<typename T>
std::unique_ptr<T[]> read_vectors(const char* filename, long& number_vector,
		int& dimension) {
	long read_vectors = 0;
	std::ifstream infile(filename, std::ifstream::binary);
	check_fstream(infile, filename);

	// Read first dimension
	std::int32_t read_dimension;
	infile.read(reinterpret_cast<char*>(&read_dimension),
			sizeof(read_dimension));
	dimension = static_cast<int>(read_dimension);

	// Compute number of vectors in file
	infile.seekg(0, std::ifstream::end);
	number_vector = infile.tellg()
			/ (dimension * sizeof(T) + sizeof(read_dimension));
	infile.seekg(0, std::ifstream::beg);

	// Allocate buffer
	std::unique_ptr<T[]> vectors(new T[number_vector * dimension]);

	// Read all vectors
	while (read_vectors != number_vector) {
		// Read dimension
		infile.read(reinterpret_cast<char*>(&read_dimension),
				sizeof(read_dimension));
		// Check all dimensions are the same
		if (static_cast<int>(read_dimension) != dimension) {
			std::cerr << "Error while reading vectors from " << filename << "."
					<< std::endl;
			std::cerr << "Vector " << number_vector << " has " << read_dimension
					<< " dimensions while other vectors have " << dimension
					<< " dimensions" << std::endl;
			std::cerr << "All vectors must have the same number of "
					<< "dimensions" << std::endl;
			std::exit(1);
		}
		// Read vector data
		infile.read(
				reinterpret_cast<char*>(vectors.get() + read_vectors * dimension),
				sizeof(T) * dimension);
		read_vectors++;
	}
	return vectors;
}

void load_query_vectors(cmdargs& args,
		query_set& query) {

	// Load query vectors set
	long query_vectors_count;
	std::unique_ptr<char[]> query_vectors = read_vectors<char>(
			args.query_set_file, query_vectors_count, query.dimension);
	std::cerr << query.dimension << " " << query_vectors_count << std::endl;

	// Load distance tables
	long distance_tables_count;
	int tables_dimension;
	std::unique_ptr<float[]> distance_tables = read_vectors<float>(
			args.distance_tables_file, distance_tables_count, tables_dimension);

	// Checks
	if(distance_tables_count != query_vectors_count) {
		std::cerr << "Found " << distance_tables_count << " distance tables in "
				<< args.distance_tables_file << std::endl;
		std::cerr << "Found " << query_vectors_count << " query vectors in "
				<< args.query_set_file << std::endl;
		std::cerr << "There must be as many distance tables as query vectors." << std::endl;
		std::exit(1);
	}
	if(tables_dimension != 2048) {
		std::cerr << "Invalid distance tables file" << std::endl;
		std::exit(1);
	}

	// Copy selected vectors and distance tables
	std::ifstream query_ids(args.query_id_file);
	int id;
	while (query_ids >> id) {
		if (id < 0 || id >= query_vectors_count) {
			std::cerr << "Query ids must be comprised between 0 and "
					<< query_vectors_count << std::endl;
			std::exit(1);
		}
		query.vectors.insert(query.vectors.end(),
				query_vectors.get() + id * query.dimension,
				query_vectors.get() + (id + 1) * query.dimension);
		query.distance_tables.insert(query.distance_tables.end(),
				distance_tables.get() + id * tables_dimension,
				distance_tables.get() + (id + 1) * tables_dimension);
		query.ids.push_back(id);
	}
}

void process_query_vectors(cmdargs& args, partition_data& part,
		query_set& query) {
	int distance_table_size = NCENT * NSQ;
	int query_n = query.vectors.size() / query.dimension;

	// Display header
	std::cout << "vec_id,partition_id,partition_n,bh_size,keep,";
	std::cout << args.bench_header << ",";
	std::cout << "fast_pq_scanned,quant_bound" << std::endl;

	// Process query vectors
	for (int i = 0; i < query_n; ++i) {
		float* dist_table = query.distance_tables.data()
				+ i * distance_table_size;
		// "Global" counters
		Counters::total_scan = 0;
		Counters::quant_bound = 0;

		//
		std::cout << query.ids[i] << "," << part.id << "," << part.n << ","
				<< args.bh_size << "," << part.keep << ",";

		// Normal PQ Scan
		pq_params pqp { 8, 8 };
		binheap* bh_oracle = nullptr;
		bh_oracle = args.bench_func(
				std::bind(scan_bh,
						reinterpret_cast<const char*>(part.partition.get()),
						dist_table, part.n, pqp, _1), args.bh_size, bh_oracle);
		std::cout << ",";
		// Fast PQ Scan
		args.bench_func(
				std::bind(scan_partition_1, part.laidout_partition.buf.get(),
						part.permuted_labels.get(), dist_table, _1),
				args.bh_size, bh_oracle);

		std::cout << "," << Counters::total_scan << "," << Counters::quant_bound
				<< std::endl;
	}
}

int main(int argc, char* argv[]) {
	// Parse arguments
	cmdargs args;
	parse_args(args, argc, argv);
	// Load partition
	partition_data part;
	load_partition(args, part);
	std::cerr << "Loaded partition" << std::endl ;
	// Load query vectors
	query_set query;
	load_query_vectors(args, query);
	std::cerr << "Loaded " << query.vectors.size() / query.dimension << " query vectors" << std::endl;
	process_query_vectors(args, part, query);
}
