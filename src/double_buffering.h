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

#ifndef DOUBLE_BUFFERING_H
#define DOUBLE_BUFFERING_H

#define NUM_BUFFERS_PER_THREAD 2

#define MAX_AIO_SIZE (1L<<30)

#define DB_DOUBLE 0
#define DB_FLOAT  1

#include <aio.h>

#include "GWAS.h"

// Argument for some double buffering IO functions
// to specify whether we refer the IO buffer or
// the computation one
typedef enum {
	IO_BUFF,
	COMP_BUFF
} which_buffer;


// List of aiocb structs + additional info to keep track
// of how many aiocb's are needed, and how many aio operations
// have been issued
typedef struct {
	int n;
	const struct aiocb **aiocb_l;
	int issued;
	// auxiliar for the suspend (wait)
	const struct aiocb **aux_l;
} aiocb_list;

// Constructor and destructor
void aiocb_list_init( aiocb_list *b, int num_chunks, FILE *fp, double *buff );
void aiocb_list_destroy( aiocb_list *b );


// 
// Double buffering handler. One per each operand
//
typedef struct {
	// Single buffer size, and MAX_AIO_SIZE chunks
	size_t buff_size;
	int num_chunks;
	// Actual buffers
	double *buffs[NUM_BUFFERS_PER_THREAD];
	// AIO control blocks
	aiocb_list aiocb_l[NUM_BUFFERS_PER_THREAD];
	aiocb_list *comp_l;
	aiocb_list *io_l;
	// Pointer to problem configuration and data
	FGLS_config_t *cf;
} double_buffering;

// Constructor, destructor, view
void double_buffering_init( double_buffering *b, size_t block_size, FILE *fp, FGLS_config_t *c );
void double_buffering_destroy( double_buffering *b );
void double_buffering_view( double_buffering *b );
// Double-buffering swap, buffers query
void double_buffering_swap( double_buffering *b );
double * double_buffering_get_io_buffer(   double_buffering *b );
double * double_buffering_get_comp_buffer( double_buffering *b );

void double_buffering_read(  double_buffering *b, which_buffer wh_buff, size_t nbytes, off_t offset );
void double_buffering_write( double_buffering *b, which_buffer wh_buff, size_t nbytes, off_t offset );

// Loads
void double_buffering_read_Y(  double_buffering *b, which_buffer wh_buff, size_t from, size_t to );
void double_buffering_read_XR( double_buffering *b, which_buffer wh_buff, size_t from, size_t to );
// Stores
void double_buffering_write_B( double_buffering *b, which_buffer wh_buff, size_t from_i, size_t to_i, size_t from_j, size_t to_j );
/*void double_buffering_write_V( double_buffering *b, which_buffer wh_buff, size_t from_i, size_t to_i, size_t from_j, size_t to_j );*/
// Wait
void double_buffering_wait( double_buffering *b, which_buffer wh_buff );

#endif // DOUBLE_BUFFERING_H
