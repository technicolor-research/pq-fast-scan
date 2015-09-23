# PQ Fast Scan #

## Description ##

The PQ Fast Scan project is a C++11 implementation of fast vector scanning 
techniques for nearest neighbor search in large databases of high-dimensional
vectors.

PQ Fast Scan builds on Product Quantization (PQ), a widely used solution for 
nearest neighbor search. To find the nearest neighbor of a query vector, PQ 
computes the distance between the query vector and databases vectors using 
*lookup tables* (stored in the L1 cache). Thus, PQ Fast Scan performs many cache 
accesses, which limit its performance. L1 cache accesses do not parallelize well 
(maximum 2 concurrent accesses) and have a high latency (about 4 cycles). To 
avoid these issues, PQ Fast Scan uses lookup tables stored in SIMD registers, 
which can be queried in parallel (16 concurrent accesses) and with a low latency 
(1 cycle).

<p align="center">
    <img src="https://storage.sbg-1.runabove.io/v1/AUTH_9e87f0768a594dfe984fbd556ac7792b/pq-fast-scan/pq-fast-scan-overview.png">
</p>

This novel approach allows PQ Fast Scan to perform 4-6x faster than PQ Scan, 
while returning the exact same results.

<p align="center">
    <img src="https://storage.sbg-1.runabove.io/v1/AUTH_9e87f0768a594dfe984fbd556ac7792b/pq-fast-scan/pq-fast-scan-performance.png">
</p>

**Contact:**  
Nicolas Le Souarnec: nicolas.le-scouarnec technicolor.com  
Fabien André: fabien.andre technicolor.com  
Please replace the space by an at sign.

## License ##

The PQ Fast Scan project is made available under the Clear BSD license 
terms (See LICENSE file).

Copyright (c) 2015 – Thomson Licensing, SAS

## Building ##

### Requirements ###

Hardware:

* Processor supporting SSE2, SSE3, SSSE3 and POPCNT

Software:

* Reasonably recent Linux distribution
* [g++](https://gcc.gnu.org/) (**4.9.x recommended**)
* [CMake](http://www.cmake.org/) (2.8 or higher)
* [libpfm4](http://perfmon2.sourceforge.net/)

On Debian-based distributions, you can install these dependencies with the 
following command:

    $ sudo apt-get install build-essential gcc g++ libpfm4-dev cmake

PQ Fast Scan can be built with other g++ versions than g++ 4.9, but we 
benchmarked it with this version only. Obtaining best results with other g++ 
versions may require some adjustments to the source code. Therefore, if your 
distribution provides g++ 4.9, we recommend installing it:

    $ sudo apt-get install g++-4.9 gcc-4.9

Otherwise, you may either build PQ Fast Scan with another g++ version, or 
install g++ 4.9 from sources.
Installing g++ 4.9 from sources takes about 20 minutes on a reasonably recent 
system. You can do so with the following commands:

    $ sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev
    $ wget ftp://ftp.irisa.fr/pub/mirrors/gcc.gnu.org/gcc/\
    releases/gcc-4.9.3/gcc-4.9.3.tar.bz2
    $ tar xvf gcc-4.9.3.tar.bz2
    $ mdkir build-gcc-4.9; cd build-gcc-4.9
    $ ../gcc-4.9.3/configure --enable-languages=c,c++ --program-suffix=-4.9
    $ make -j8
    $ sudo make install
    $ which g++-4.9
    
Some quick tests with g++ 5.2 seemed to indicate that PQ Fast Scan performs 
well when compiled with this version. You may also be able to compile PQ Fast 
Scan with clang++ but we have not tested it.

### Building PQ Fast Scan ###

Once you have gathered all the requirements, building PQ Fast Scan is 
straightforward:

    $ git clone https://github.com/technicolor-research/pq-fast-scan.git
    $ mkdir build-pq-fast-scan
    $ cd build-pq-fast-scan

If you use g++ 4.9:

    $ CC=gcc-4.9 CXX=g++-4.9 cmake ../pq-fast-scan
    $ make

To use the default g++ version of your distribution:

    $ cmake ../pq-fast-scan
    $ make

## Usage ##

This project provides two executables :
    
* **pqscan**, to benchmark different PQ Scan implementations on *synthetic* data
* **pqfastscan**, to benchmark PQ Fast Scan against PQ Scan on *real* data

**pqscan** does not take any command-line argument. It outputs run times and 
performance counters for different PQ Scan implementations. You may need to run 
it as root to access the performance counters kernel API.

    $ sudo ./pqscan
    Scan: time=74758,cycles=269918483,instructions=800272344,L1-dcache-loads=
    400231701,L1-dcache-load-misses=3139139 [...]
    Scan +sse +prefetch: time=91214,cycles=329937021,instructions=787927217,
    L1-dcache-loads=462852107,L1-dcache-load-misses=3139943 [...]
    Scan +avx +prefetch: time=68490,cycles=231844613,instructions=618770051,
    L1-dcache-loads=409382281,L1-dcache-load-misses=3141246 [...]
    Scan +avx +vgather +prefetch: time=97841,cycles=331413774,
    instructions=188461772,L1-dcache-loads=228194484, [...]

**pqfastscan** requires four input files: a *partition* to scan, a *list*
of IDs of query vectors, a set of *query vectors* and a set of *distance*
*tables*.

We provide the following input files:

* *partitions* and *lists*: `100M1-8-partitions.tar.xz`  
    [Link 1](https://storage.sbg-1.runabove.io/v1/AUTH_9e87f0768a594dfe984fbd556ac7792b/pq-fast-scan/100M1-8-partitions.tar.xz)
* *query vectors*: `bigann_query.bvecs.gz`  
    [Corpus Texmex](http://corpus-texmex.irisa.fr/), download *ANN_SIFT1B Query Set*
* *distance tables*: `bigann_distance_tables.fvecs.xz`  
    [Link 1](https://storage.sbg-1.runabove.io/v1/AUTH_9e87f0768a594dfe984fbd556ac7792b/pq-fast-scan/bigann_distance_tables.fvecs.xz)

For more information about these datasets, see the Datasets section.

Decompress the datasets with the following commands:

    $ tar xvf 100M1-8-partitions.tar.xz
    $ gzip -d bigann_query.bvecs.gz
    $ xz -d bigann_distance_tables.fvecs.xz

You can then run pqfastscan, e.g. to scan partition 0:

    $ ./pqfastscan 100M1-8-partition/100M1-partition-0.dat bigann_query.bvecs\
    bigann_distance_tables.fvecs 100M1-8-partition/100M1-list-0.txt
    
    vec_id,partition_id,partition_n,bh_size,keep,pq_us,fast_pq_us, [...]
    1,0,25159451,100,125797,75392,12272,24404,93728.6
    4,0,25159451,100,125797,74754,17480,435694,96652.2
    [...]

You may replace 0 by another partition number [0-7]. Make sure you use the list 
file corresponding to the partition, e.g., `100M1-list-6.txt` with 
`100M1-partition-6.dat`.

pqfastscan outputs one line (comma separated values) for each query vector of
the list with the scan time for PQ Scan (`pq_us`) and PQ Fast Scan (`fast_pq_us`).

In addition to run times, you may also want to output performance counters.
To do so, use the `-p` option. Like pqscan, you may also need to run pqfastscan
as root when using performance counters:

    $ sudo ./pqfastscan -p 100M1-8-partition/100M1-partition-0.dat\
    bigann_query.bvecs bigann_distance_tables.fvecs\
    100M1-8-partition/100M1-list-0.txt
    
    vec_id,partition_id,partition_n,bh_size,keep,pq_us,pq_cycles, [...]
    1,0,25159451,100,125797,76247,272870815,805615127,402679306,12027, [...]
    4,0,25159451,100,125797,76000,273304611,805588015,402670181,18040, [...]
    [...]

For more information about pqfastscan options, see `./pqfastscan -h`.

## Datasets ##

To test pqfastscan, we provide:

* *partitions* and *lists*: `100M1-8-partitions.tar.xz`  
    [Link 1](https://storage.sbg-1.runabove.io/v1/AUTH_9e87f0768a594dfe984fbd556ac7792b/pq-fast-scan/100M1-8-partitions.tar.xz)
* *distance tables*: `bigann_distance_tables.fvecs.xz`  
    [Link 1](https://storage.sbg-1.runabove.io/v1/AUTH_9e87f0768a594dfe984fbd556ac7792b/pq-fast-scan/bigann_distance_tables.fvecs.xz)
    
Besides, we use the query set from the well-known ANN_SIFT1B dataset;

* *query vectors*: `bigann_query.bvecs.gz`  
    [Corpus Texmex](http://corpus-texmex.irisa.fr/), 
    download *ANN_SIFT1B Query Set*

We generated the *partitions*, *lists* and *distance tables* files we use from
the ANN_SIFT1B dataset.

You may re-generate these files yourself. To generate these files from the ANN_SIFT1B dataset,
you need an implementation of product quantization. The authors of the 
[seminal paper on product quantization](http://www.computer.org/csdl/trans/tp/2011/01/ttp2011010117-abs.html)
distribute 
[two implementations](http://people.rennes.inria.fr/Herve.Jegou/projects/ann.html)
of product quantization: a C implementation, and a matlab implementation.
We used the C implementation (libpq), which we obtained under a commercial
licence. You can obtain similar files using the matlab implementation,
freely available.

Once you have an implementation of product quantization, you can generate
*partitions*, *lists* and *distance tables* following these steps:

1. Extract the first 100 million vectors from the ANN_SIFT1B Base Set
(`bigann_base.bvecs`). We name this extracted dataset ANN_SIFT100M1.
2. Build an IVFADC (database) using the ANN_SIFT1B Learning Set 
(`bigann_learn.bvecs`). We used a *coarse quantizer* with *8 centroids* and a
*product quantizer* with *8 sub-quantizers* of *256 centroids* each.
3. Modify the assignement of sub-quantizer centroid indexes. To do so, use a
    [k-means variant which forces groups of same sizes](http://elki.dbs.ifi.lmu.de/wiki/Tutorial/SameSizeKMeans).
4. Add all vectors of ANN_SIFT100M1 to the IVFADC.
5. We name *partition* each Voronoi cell of the coarse quantizer. Save each 
partition [0-7] to a different file. See below for file format.
6. Using the coarse quantizer, map each query vector to its corresponding
partition (no multiple assignement, i.e. ma=1). Save the query vectors IDs 
corresponding to each partition in a different *list*. See below for file format.
7. Compute distance tables for each vector. Save *distance tables* sequentially
in a single file. See below for file format.

**partitions file format:** [Binary file]  
The number of pqcodes in the partition (32-bit int), followed by all pqcodes of 
the partition (8x8 bits per pqcode).  
*Filename:* The last number before the extension is parsed as the partition
number.

**distance tables file format:**    [Binary file]  
Each query vector has 8 distance tables of 256 floats associated with it.
The 8 distance tables of each query vector are stored as a vector of 2048
floats (8 x 256 floats). These vectors of 2048 floats are stored as an
[fvecs file](http://corpus-texmex.irisa.fr/).  
*Filename:* No restrictions.

**lists file format:**              [Text file]  
One vector ID per line.  
*Filename:* The last number before the extension is parsed as the partition
number. If the vector ID 42 is in the file xxxx-list-6.txt, it means that the
vector number 42 from the query set was mapped to the partition 6 by the coarse
quantizer.
