/*
 * Copyright (c) 2010-2012, Diego Fabregat-Traver and Paolo Bientinesi.
 * All rights reserved.
 *
 * This file is part of OmicABEL.
 * 
 * OmicABEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * OmicABEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with OmicABEL. If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * Coded by:
 *   Diego Fabregat-Traver (fabregat@aices.rwth-aachen.de)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>

#include "timing.h"

long get_diff_ms(struct timeval *start_time, struct timeval *end_time) {
    long seconds = end_time->tv_sec - start_time->tv_sec;
    long useconds = end_time->tv_usec - start_time->tv_usec;
    return ((seconds) * 1000 + useconds/1000.0) +0.5;
}

int read_clock(struct timeval *t)
{
    return gettimeofday(t, NULL);
}

int elapsed_time(struct timeval *start, struct timeval *end)
{
    return (int)(end->tv_sec  - start->tv_sec) * 1e6 + 
		   (int)(end->tv_usec - start->tv_usec);
}
