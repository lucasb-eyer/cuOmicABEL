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
#include <string.h>

#include <assert.h>

#include "GWAS.h"
#include "wrappers.h"
#include "double_buffering.h"

static size_t size_t_ceil( size_t num, size_t den );

// List of aiocb structs - Constructor
void aiocb_list_init( aiocb_list *b, int num_chunks, FILE *fp, double *buff )
{
	int i;

	b->n = num_chunks;
	b->issued = 0;

	b->aiocb_l = ( const struct aiocb ** ) fgls_malloc ( b->n * sizeof(struct aiocb *) );
	for ( i = 0; i < b->n; i++ )
	{
		struct aiocb *tmp_cb;
		tmp_cb = (struct aiocb *) fgls_malloc ( sizeof(struct aiocb) );
		// zero-out the structure and fill it up with "constant" data
		bzero( (char *)tmp_cb, sizeof(struct aiocb) );
		tmp_cb->aio_fildes = fileno( fp );
		tmp_cb->aio_buf = &buff[i*(size_t)MAX_AIO_SIZE/sizeof(double)];
		// assign to the const list
		b->aiocb_l[i] = tmp_cb;
	}
	// Will wait one by one. After all, I need them all to be done 
	// ( and usually only one call actually needed, block size is not that big )
	b->aux_l = ( const struct aiocb ** ) fgls_malloc ( sizeof(struct aiocb *) );
}

// List of aiocb structs - Destructor
void aiocb_list_destroy( aiocb_list *b )
{
	int i;
	struct aiocb *tmp_cb;

	// aiocb structs
	for ( i = 0; i < b->n; i++ )
	{
		tmp_cb = (struct aiocb *)b->aiocb_l[i];
		free( tmp_cb );
	}
	free( b->aiocb_l ); // WARNING ?
	// Reset
	b->n = 0;
	b->issued = 0;
	b->aiocb_l = NULL;
	free( b->aux_l ); // WARNING ?
	b->aux_l = NULL;
}

//
// Double buffering handler
//
// Constructor
void double_buffering_init( double_buffering *b, size_t buff_size, FILE *fp, FGLS_config_t *c )
{
	int i;
	
	// Single buffer size
	b->buff_size = buff_size;
	// Number of chunks into which an aio operation is split
	b->num_chunks = size_t_ceil( buff_size, MAX_AIO_SIZE );

	// Allocate double buffers
	for ( i = 0; i < NUM_BUFFERS_PER_THREAD; i++ )
		b->buffs[i] = ( double * ) fgls_malloc ( buff_size );

	// Create the necessary amount of aiocb lists,
	// one for io, one for comp,
	// and initialize them
	for ( i = 0; i < NUM_BUFFERS_PER_THREAD; i++ )
		aiocb_list_init( &b->aiocb_l[i], b->num_chunks, fp, b->buffs[i] );
	b->comp_l = &b->aiocb_l[0];
	b->io_l   = &b->aiocb_l[1];
	// Config
	b->cf = c;
}

// Destructor
void double_buffering_destroy( double_buffering *b )
{
	int i;

	// aiocb structs
	for ( i = 0; i < NUM_BUFFERS_PER_THREAD; i++ )
		aiocb_list_destroy( &b->aiocb_l[i] );
	// double buffering
	for ( i = 0; i < NUM_BUFFERS_PER_THREAD; i++ )
		free( b->buffs[i] );
}

// View - To be improved
void double_buffering_view( double_buffering *b )
{
	printf("Buffer size: %zu (%p, %p)\n", b->buff_size, b->buffs[0], b->buffs[1]);
	printf("Number of chunks: %d\n", b->num_chunks);

	printf("Waiting aiocb for:\n");
	printf("\tIssued: %d\n", b->comp_l->issued);
	printf("\tNBytes: %zu\n", b->comp_l->aiocb_l[0]->aio_nbytes);
	printf("\tOffset: %jd\n", b->comp_l->aiocb_l[0]->aio_offset);
	printf("\tFile: %d\n", b->comp_l->aiocb_l[0]->aio_fildes);
	fflush(stdout);
}

// Double-buffering swap, 
void double_buffering_swap( double_buffering *b )
{
	aiocb_list *tmp = b->comp_l;
	b->comp_l = b->io_l;
	b->io_l   = tmp;
}

// Buffers query
double * double_buffering_get_io_buffer(   double_buffering *b )
{
	return (double *)b->io_l->aiocb_l[0]->aio_buf;
}
double * double_buffering_get_comp_buffer( double_buffering *b )
{
	return (double *)b->comp_l->aiocb_l[0]->aio_buf;
}

// Loads
void double_buffering_read_Y(  double_buffering *b, which_buffer wh_buff, size_t from, size_t to )
{
	FGLS_config_t *cf = b->cf;

	// Which buffer
	/*aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;*/
	// Still buff_size available?
	size_t nbytes = (to - from + 1) * cf->n * sizeof(double);
	// Number of chunks to read
	/*int chunks_to_read = size_t_ceil( nbytes, (size_t)MAX_AIO_SIZE );*/
	// offset
	off_t offset = (off_t)from * cf->n * sizeof(double);

#if DEBUG
	printf("Reading Y: %zu bytes from %jd\n", nbytes, offset);
	fflush(stdout);
#endif
	double_buffering_read( b, wh_buff, nbytes, offset );
}

void double_buffering_read_XR( double_buffering *b, which_buffer wh_buff, size_t from, size_t to )
{
	FGLS_config_t *cf = b->cf;

	// Which buffer
	/*aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;*/
	// Still buff_size available?
	size_t nbytes = (to - from + 1) * cf->n * cf->wXR * sizeof(double);
	// Number of chunks to read
	/*int chunks_to_read = size_t_ceil( nbytes, (size_t)MAX_AIO_SIZE );*/
	// offset
	off_t offset = (off_t)from * cf->n * cf->wXR * sizeof(double);

#if DEBUG
	printf("Reading XR: %zu bytes from %jd\n", nbytes, offset);
	fflush(stdout);
#endif
	double_buffering_read( b, wh_buff, nbytes, offset );
}

// Stores
void double_buffering_write_B(  double_buffering *b, which_buffer wh_buff, size_t from_i, size_t to_i, size_t from_j, size_t to_j )
{
	FGLS_config_t *cf = b->cf;
	size_t size_one_b = cf->p + (cf->p*(cf->p+1)) / 2; // beta + half var-cov matrix

	// Which buffer
	/*aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;*/
	// Still buff_size available?
	size_t nbytes = (to_i - from_i + 1) * size_one_b *
	                (to_j - from_j + 1) * sizeof(double);
	// Number of chunks to write
	/*int chunks_to_write = size_t_ceil( nbytes, (size_t)MAX_AIO_SIZE );*/
	// offset
	off_t offset = ((off_t)from_j * cf->m * size_one_b + 
			        (off_t)from_i * size_one_b * (to_j - from_j + 1)) * sizeof(double);

#if DEBUG
	printf("Writing B: %zu bytes from %jd\n", nbytes, offset);
	fflush(stdout);
#endif
	double_buffering_write( b, wh_buff, nbytes, offset );
}

void double_buffering_write_V(  double_buffering *b, which_buffer wh_buff, size_t from_i, size_t to_i, size_t from_j, size_t to_j )
{
	FGLS_config_t *cf = b->cf;

	// Which buffer
	/*aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;*/
	// Still buff_size available?
	size_t nbytes = ((to_i - from_i + 1) * cf->p * cf->p *
	                 (to_j - from_j + 1)) * sizeof(double);
	// Number of chunks to write
	/*int chunks_to_write = size_t_ceil( nbytes, (size_t)MAX_AIO_SIZE );*/
	// offset
	off_t offset = ((off_t)from_j * cf->m * cf->p * cf->p + 
			        (off_t)from_i * cf->p * cf->p * (to_j - from_j + 1)) * sizeof(double);

#if DEBUG
	printf("Writing V: %zu bytes from %jd\n", nbytes, offset);
	fflush(stdout);
#endif
	double_buffering_write( b, wh_buff, nbytes, offset );
}

void double_buffering_read( double_buffering *b, which_buffer wh_buff, size_t nbytes, off_t offset )
{
	// Buffer
	aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;
	// Number of chunks to read
	int chunks_to_read = size_t_ceil( nbytes, (size_t)MAX_AIO_SIZE );

#if DEBUG
	printf("Reading %d chunks\n", chunks_to_read);
	fflush(stdout);
#endif

	int i;
	struct aiocb *tmp_cb;
	// First n-1 chunks are full ( MAX_AIO_SIZE )
	for ( i = 0; i < chunks_to_read - 1; i++ )
	{
		tmp_cb = (struct aiocb *)buff->aiocb_l[i];
		// fildes and buffer already set up in the initialization
		tmp_cb->aio_nbytes = (size_t)MAX_AIO_SIZE;
		tmp_cb->aio_offset = offset + (i * (off_t)MAX_AIO_SIZE);
		if ( aio_read( tmp_cb ) != 0 )
		{
			perror("[ERROR] Couldn't read");
			exit( EXIT_FAILURE );
		}
	}

	if (chunks_to_read)
	{
		i = chunks_to_read - 1;
		tmp_cb = (struct aiocb *)buff->aiocb_l[i];
		tmp_cb->aio_nbytes = nbytes - (i * (size_t)MAX_AIO_SIZE);
		tmp_cb->aio_offset = offset + (i * (off_t)MAX_AIO_SIZE);
		if ( aio_read( tmp_cb ) != 0 )
		{
			perror("[ERROR] Couldn't read from XR");
			exit( EXIT_FAILURE );
		}
	}

	buff->issued = chunks_to_read;
}

void double_buffering_write( double_buffering *b, which_buffer wh_buff, size_t nbytes, off_t offset )
{
	// Buffer
	aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;
	// Number of chunks to write
	int chunks_to_write = size_t_ceil( nbytes, (size_t)MAX_AIO_SIZE );

	int i;
	struct aiocb *tmp_cb;
	// First n-1 chunks are full ( MAX_AIO_SIZE )
	for ( i = 0; i < chunks_to_write - 1; i++ )
	{
		tmp_cb = (struct aiocb *)buff->aiocb_l[i];
		// fildes and buffer already set up in the initialization
		tmp_cb->aio_nbytes = (size_t)MAX_AIO_SIZE;
		tmp_cb->aio_offset = offset + (i * (size_t)MAX_AIO_SIZE);
		if ( aio_write( tmp_cb ) != 0 )
		{
			perror("[ERROR] Couldn't write");
			exit( EXIT_FAILURE );
		}
	}

	if (chunks_to_write)
	{
		i = chunks_to_write - 1;
		tmp_cb = (struct aiocb *)buff->aiocb_l[i];
		tmp_cb->aio_nbytes = nbytes - (i * (size_t)MAX_AIO_SIZE);
		tmp_cb->aio_offset = offset + (i * (off_t)MAX_AIO_SIZE);
		if ( aio_write( tmp_cb ) != 0 )
		{
			perror("[ERROR] Couldn't write");
			exit( EXIT_FAILURE );
		}
	}

	buff->issued = chunks_to_write;
}

// Wait for data to be read or written
void double_buffering_wait( double_buffering *b, which_buffer wh_buff )
{
	int i;
	ssize_t nbytes;
	// Which buffer
	aiocb_list *buff = wh_buff == IO_BUFF ? b->io_l : b->comp_l;

#if DEBUG
	printf("Waiting for %p\n%d operations issued\nLast chunk: %zu bytes\n", 
			buff, buff->issued, buff->aiocb_l[buff->issued-1]->aio_nbytes);
	fflush(stdout);
#endif
	for ( i = 0; i < buff->issued; i++ )
	{
		buff->aux_l[0] = (struct aiocb *)buff->aiocb_l[i];
#if DEBUG
	    printf("Waiting for chunk %d (of %d): %zu bytes\n", i, buff->issued, buff->aux_l[0]->aio_nbytes);
		printf("Start: "); print_timestamp();
		fflush(stdout);
#endif
		if ( aio_suspend( buff->aux_l, 1, NULL ) != 0 )
		{
			perror("aio_suspend error");
			exit( EXIT_FAILURE );
		}

		nbytes = aio_return( (struct aiocb *)buff->aux_l[0] );
		if ( nbytes != buff->aux_l[0]->aio_nbytes )
		{
			perror("aio_suspend error - data not read completely");
			printf("%zu read - %zu expected\n", nbytes, buff->aux_l[0]->aio_nbytes);
			exit( EXIT_FAILURE );
		}
#if DEBUG
		printf("End: "); print_timestamp();
		fflush(stdout);
#endif
	}
}


static size_t size_t_ceil( size_t num, size_t den )
{
	assert( den != 0 );
	return (num + ( den - 1 )) / den;
}
