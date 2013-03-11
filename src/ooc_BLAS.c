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

/*#include <pthread.h>*/
#include <aio.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "GWAS.h"
#include "wrappers.h"
#include "utils.h"
#include "double_buffering.h"
#include "ooc_BLAS.h"

#if VAMPIR
  #include "vt_user.h"
#endif

#include <stdint.h>

/*
 * Out-of-core gemms:
 *   - Z' XR
 *   - Z' Y
 * Z is m x m
 * The other matrix is m x n
 */
void ooc_gemm( int m, int n, int ooc_b, double *Z, char *in, char *out, 
		int threshold, const char *obj_type, char *obj_name, int namelength )
{
	/* Files */
	FILE *fp_in  = fgls_fopen( in, "rb" );
	FILE *fp_out = fgls_fopen( out, "wb" );

    /* OOC Problem dimensions */
	/*size_t max_elems_per_buffer = 1L << 26; // 64MElems, 512 MBs*/
	/*max_elems_per_buffer = max_elems_per_buffer - max_elems_per_buffer % n;*/
	/*size_t num_cols_per_buff = max_elems_per_buffer / n;*/

    /* Asynchronous IO data structures */
	double *in_comp, *out_comp;
	double_buffering db_in, db_out; // B, C
	double_buffering_init( &db_in, ooc_b * m * sizeof(double),
			                fp_in, NULL ); // _fp, cf not needed in this case
	double_buffering_init( &db_out, ooc_b * m * sizeof(double),
			                fp_out, NULL ); // _fp, cf not needed in this case

    /* BLAS constants */
    double ONE  = 1.0;
    double ZERO = 0.0;

    /* Read first piece of "in" */
    double_buffering_read( &db_in, IO_BUFF,
                           MIN( (size_t)ooc_b * m, (size_t)m * n ) * sizeof(double), 0);
	double_buffering_swap( &db_in );

    int cur_n;
    int i;
    for ( i = 0; i < n; i += ooc_b ) 
    {
        /* Read next piece of "in" */
        size_t nbytes = i + ooc_b > n ? 1 : MIN( ooc_b * m, ( n - (size_t)( i + ooc_b ) ) * m ) * sizeof(double);
		off_t  offset = i + ooc_b > n ? 0 : (off_t)(i + ooc_b) * m * sizeof(double);
        double_buffering_read( &db_in, IO_BUFF, nbytes, offset );

        /* Wait for current piece of "in" */
#if VAMPIR
    VT_USER_START("OOC_GEMM_WAIT");
#endif
		double_buffering_wait( &db_in, COMP_BUFF );
#if VAMPIR
    VT_USER_END("OOC_GEMM_WAIT");
#endif

        /* Compute */
		in_comp  = double_buffering_get_comp_buffer( &db_in );
		out_comp = double_buffering_get_comp_buffer( &db_out );
        cur_n = MIN( ooc_b, (n - i) );
		/*printf("Compute\n");*/

		// Sanity check
		average( in_comp, m, cur_n, threshold, obj_type, &obj_name[i*namelength], namelength, 1 );
#if VAMPIR
    VT_USER_START("OOC_GEMM");
#endif
	/*printf("\nPRE: ");  print_timestamp(); fflush( stdout );*/
        dgemm_("T", "N", &m, &cur_n, &m, &ONE, Z, &m, in_comp, &m, &ZERO, out_comp, &m);
		/*printf("\nPOST: "); print_timestamp(); fflush( stdout );*/
#if VAMPIR
    VT_USER_END("OOC_GEMM");
#endif

        /* Wait until previous piece of "out" is written */
        if ( i > 0)
			double_buffering_wait( &db_out, IO_BUFF );

        /* Write current piece of "out" */
		double_buffering_write( &db_out, COMP_BUFF,
				                MIN( ooc_b * m, (size_t)(n - i) * m ) * sizeof(double),
                                (off_t)i * m * sizeof(double) );

        /* Swap buffers */
		double_buffering_swap( &db_in );
		double_buffering_swap( &db_out );
    }

    /* Wait for the remaining io calls issued */
	double_buffering_wait( &db_in,  COMP_BUFF );
	double_buffering_wait( &db_out, IO_BUFF );

	/* Clean-up */
	double_buffering_destroy( &db_in );
	double_buffering_destroy( &db_out );

	fclose( fp_in );
	fclose( fp_out );
}
