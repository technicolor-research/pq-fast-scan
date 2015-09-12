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
#include <cstdio>
#include <cstdint>
#include <cstdarg>
#include <cstring>

#include <err.h>
#include <signal.h>

#include <perfmon/pfmlib.h>
#include <perfmon/pfmlib_perf_event.h>

#include <iostream>

// Error management
static void pfm_success_or_die(int ret, const char* format, ...) {
	if (ret != PFM_SUCCESS) {
		va_list argptr;
		va_start(argptr, format);
		vfprintf(stderr, format, argptr);
		fprintf(stderr, "libpfm error: %s\n", pfm_strerror(ret));
		exit(1);
	}
}

static void not_error_or_die(int ret, int error_value, const char* format, ...) {
	if (ret == error_value) {
		va_list argptr;
		va_start(argptr, format);
		vfprintf(stderr, format, argptr);
		va_end(argptr);
		perror("");
		exit(1);
	}
}

static void equal_or_die(int ret, int good_value, const char* format, ...) {
	if(ret != good_value) {
		va_list argptr;
		va_start(argptr, format);
		vfprintf(stderr, format, argptr);
		va_end(argptr);
		perror("");
		exit(1);
	}
}

void perf_events_open_by_name(int* fds, int n, const char** events) {
	int ret = pfm_initialize();
	pfm_success_or_die(ret, "Could not initialize libpfm\n");

	pfm_perf_encode_arg_t arg;
	perf_event_attr attr;

	for (int i = 0; i < n; ++i) {
		memset(&attr , 0, sizeof(attr));
		memset(&arg, 0, sizeof(arg));
		arg.attr = &attr;
		int ret = pfm_get_os_event_encoding(events[i], PFM_PLM0 | PFM_PLM3,
				PFM_OS_PERF_EVENT_EXT, &arg);
		pfm_success_or_die(ret, "Could not get event %s\n", events[i]);
		// Disable event so that counting does not start just after
		// perf_event_open(...)
		attr.disabled = 1;
		attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED|PERF_FORMAT_TOTAL_TIME_RUNNING;
		fds[i] = perf_event_open(&attr, 0, -1, -1, 0);
		not_error_or_die(fds[i], -1, "Could not open event %d\n", i);
	}
}

std::uint64_t perf_scale(std::uint64_t values[3]) {
	uint64_t res = 0;

	if (!values[2] && !values[1] && values[0])
		warnx("WARNING: time_running = 0 = time_enabled, raw count not zero\n");

	if (values[2] > values[1])
		warnx("WARNING: time_running > time_enabled\n");

	if (values[2])
		res = (uint64_t) ((double) values[0] * values[1] / values[2]);
	return res;
}

void read_counter(int fd, std::uint64_t raw[3]) {
	const int size = 24;
	ssize_t ret = read(fd, raw, size);
	equal_or_die(ret, size, "Failed to read counter %d\n", fd);
}

void perf_events_read(int* fds, int n, std::uint64_t scaled[], std::uint64_t raw[][3]) {
	for(int i=0; i < n; ++i) {
		read_counter(fds[i], raw[i]);
		scaled[i] = perf_scale(raw[i]);
	}
}

void perf_events_start() {
	const int ret = prctl(PR_TASK_PERF_EVENTS_ENABLE);
	equal_or_die(ret, 0, "Could not enable perf events\n");
}

void perf_events_stop() {
	const int ret = prctl(PR_TASK_PERF_EVENTS_DISABLE);
	equal_or_die(ret, 0, "Could not disable perf events\n");
}

void perf_events_reset(int* fds, int n) {
	for(int i = 0; i < n; ++i) {
		ioctl(fds[i], PERF_EVENT_IOC_RESET, 0);
	}
}

void perf_events_close(int* fds, int n) {
	for(int i = 0; i < n; ++i) {
		close(fds[i]);
	}
}




