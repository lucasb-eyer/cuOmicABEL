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

#ifndef UTILS_H
#define UTILS_H

#define DOUBLE_NAN 0x7FF8000000000000

#include "GWAS.h"
/*#include "utils.h"*/

// Threading
void set_multi_threaded_BLAS( int nths );
void set_single_threaded_BLAS( void );

// Compute an approximation of the machine epsilon
// in double precision (within a factor of two)
double get_epsilon( void );

// Startup
void get_main_memory_size( size_t *totalMem, size_t *availMem );
void estimate_block_sizes( FGLS_config_t *cf, char var, int estimate_inc );

void print_timestamp( void );

// Sanity
void average( double *data, int n, int ncols, int threshold, const char *obj_type, char *obj_name, int namelength, int verbose );
void checkNoNans( size_t n, double *buff, const char* err_msg);

#endif // UTILS_H
