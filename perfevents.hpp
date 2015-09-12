//
// Copyright (c) 2015 – Thomson Licensing, SAS
//
// The source code form of this open source project is subject to the terms of the
// Clear BSD license.
//
// You can redistribute it and/or modify it under the terms of the Clear BSD
// License (See LICENSE file).
//


#ifndef PERFEVENTS_HPP_
#define PERFEVENTS_HPP_

#include <cstdint>

void perf_events_open_by_name(int* fds, int n, const char** events);
void perf_events_read(int* fds, int n, std::uint64_t scaled[], std::uint64_t raw[][3]);

void perf_events_start();
void perf_events_stop();
void perf_events_reset(int* fds, int n);
void perf_events_close(int* fds, int n);

#endif /* PERFEVENTS_HPP_ */
