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
#include <strings.h>

#include "GWAS.h"
#include "wrappers.h"

void sync_read(double *buffer, FILE* fp, size_t size, size_t start) 
{
    size_t out;

    fseek( fp, start*sizeof(double), SEEK_SET );
    out = fread( buffer, sizeof(double), size, fp );
    if ( out != size )
    {
        printf("[ERROR] "__FILE__ ": Read %d elements - %d expected\n", (int)out, (int)size);
    }

    return; 
}

void sync_write(double* buffer, FILE* fp, size_t size, size_t start) 
{
    size_t out;

    fseek(fp, start * sizeof(double), SEEK_SET);
    out = fwrite(buffer, sizeof(double), size, fp);
    if ( out != size ) 
    {
        printf("[ERROR] "__FILE__ ": Wrote %d elements - %d expected\n", (int)out, (int)size);
    }

    return; 
}

FILE * fgls_fopen( const char *path, const char *mode )
{
	FILE * f;
	char err[STR_BUFFER_SIZE];

	f = fopen( path, mode );
	if ( f == NULL )
	{
		snprintf(err, STR_BUFFER_SIZE, "Couldn't open %s", path);
		perror( err );
		exit( EXIT_FAILURE );
	}
	return f;
}

void load_file( const char *path, const char *mode, double *buff, size_t size )
{
	FILE *fp;

    fp = fgls_fopen( path, mode );
    sync_read( buff, fp, size, 0 );
    fclose( fp );
}

void * fgls_malloc_impl( const char* file, long line, size_t size )
{
    void *buf;

    if ( (buf = malloc( size )) == NULL ) {
        fprintf(stderr, "Couldn't allocate %ld bytes of memory in %s:%ld\n", size, file, line);
        exit(EXIT_FAILURE);
    }

    return buf;
}

/*void error_msg(char *file, int line, char *msg, int abort)*/
void error_msg(char *msg, int abort)
{
    /*fprintf(stderr, "[Error] %s(line %d): %s\n", file, line, msg);*/
    fprintf(stderr, "[Error]: %s\n", msg);
    if (abort)
        exit( EXIT_FAILURE );
}
