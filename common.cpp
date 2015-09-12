//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstdio>
#include "common.hpp"

const char* colors[] = {
	"\e[0m",
	"\e[1;30m",
	"\e[1;31m",
};

void print_buffer(const std::uint8_t* buf, unsigned size) {
	for(unsigned i = 0; i < size; ++i) {
		std::cout << std::setw(3) << (unsigned) buf[i] << " ";
		if(i % 16 == 15) {
			std::cout << std::endl;
		}
	}
	std::cout << std::dec << std::endl;
}

void print_m128i(__m128i reg) {
    std::uint8_t freg[16];
    _mm_store_si128(reinterpret_cast<__m128i*>(freg), reg);
    int i;
    for(i = 0; i < 16; ++i) {
        printf("%d ", freg[i]);
    }
    printf("\n");
}

void print_uint64(std::uint64_t reg) {
    std::uint8_t* freg=reinterpret_cast<std::uint8_t*>(&reg);
    int i;
    for(i = 0; i < 8; ++i) {
        printf("%d ", freg[i]);
    }
    printf("\n");
}
