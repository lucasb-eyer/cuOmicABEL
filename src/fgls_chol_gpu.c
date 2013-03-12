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
 *   Lucas Beyer (lucas-b.eyer.be)
 */

#ifdef WITH_GPU
// Do nothing if there is no GPU support.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

#include <aio.h>

#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include <omp.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "GWAS.h"
#include "wrappers.h"
#include "timing.h"
#include "double_buffering.h"
#include "utils.h"
#include "fgls_chol_gpu.h"

#if VAMPIR
    #include "vt_user.h"
#endif

/*
 * Builds Phi as an SPD matrix, after the eigenvalues were fixed 
 * during the REML estimation
 */
void build_SPD_Phi( int n, double *eigVecs, double *eigVals, double *Phi );

/*
 * Cholesky-based GPU solution of the 
 *  sequence of Feasible Generalized Least-Squares problem
 *  in the context of GWAS:
 */
int fgls_chol_gpu( FGLS_config_t cf )
{
	int n = cf.n,
		   m = cf.m,
		   p = cf.p,
		   t = cf.t,
		   x_b = cf.x_b,
		   /*y_b = cf.y_b,*/
		   wXL = cf.wXL,
		   wXR = cf.wXR;
    /* In-core operands */
    double *Phi;
    double *M;
    double *ests;
    double *h2;
	double *res_sigma;
    double alpha;
    double beta;

    /* Out-of-core operands */
    double *Bij; // Auxiliary variables

    /* Reusable data thanks to constant XL */
    double *XL;
    double *XL_orig; // XL and a copy (XL is overwritten at every iteration of j)
    double *B_t;  // Top part of b ( in inv(S) b )
    double *V_tl; // Top-Left part of V

    /* BLAS / LAPACK constants */
    double ZERO = 0.0;
    double ONE = 1.0;
    int iONE = 1;
    /* LAPACK error value */
    int info;

    /* iterators and auxiliar vars */
    int ib, i, j, k, l; // size_t
    int nn = cf.n * cf.n; // size_t
	size_t size_one_b_record = p + (p*(p+1))/2;

	// Threading
	int id;
	double *tmpBs, *tmpVs; // Buffer with one B and one V per thread
	double *oneB, *oneV;   // Each thread pointer to its B and V

    fprintf(stdout, "\nHi from the GPU code!\n");

    if ( cf.y_b != 1 )
	{
        fprintf(stderr, "\n[Warning] y_b not used (set to 1)\n");
		cf.y_b = 1;
	}

    if ( cf.t > 1 ) {
        fprintf(stderr, "\nThe chol_gpu variant doesn't support multiple phenotypes. Use the chol or eigen variants.\n");
        exit(1); // Die a fast death.
    }

    /* Memory allocation */
    // In-core
	build_SPD_Phi( cf.n, cf.Z, cf.W, cf.Phi );
	Phi   = cf.Phi;
    M     = ( double * ) fgls_malloc ( (size_t)cf.n * cf.n * sizeof(double) );
    ests  = cf.ests;

	h2 = ests;
	res_sigma = &ests[2*cf.t];

    XL_orig = cf.XL;
    XL      = ( double * ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
    B_t  = ( double * ) fgls_malloc ( cf.wXL * sizeof(double) );
    V_tl = ( double * ) fgls_malloc ( cf.wXL * cf.wXL * sizeof(double) );

	// Temporary storage prior to copying in db_B
    tmpBs = ( double * ) fgls_malloc ( cf.p * cf.num_threads * sizeof(double) );
    tmpVs = ( double * ) fgls_malloc ( cf.p * cf.p * cf.num_threads * sizeof(double) );

    /* Files and pointers for out-of-core */
    double *XR_comp, *Y_comp, *B_comp;

    /* Asynchronous IO data structures */
	double_buffering db_XR, db_Y, db_B;
	double_buffering_init( &db_XR, (size_t)cf.n * cf.wXR * cf.x_b * sizeof(double),
			                cf.XR, &cf ); // _fp
	double_buffering_init( &db_Y, (size_t)cf.n * cf.y_b * sizeof(double),
			                cf.Y,  &cf );
	double_buffering_init( &db_B, (size_t)size_one_b_record * cf.x_b * cf.y_b * sizeof(double),
			                cf.B,  &cf );

#if VAMPIR
    VT_USER_START("READ_X");
#endif
    /* Read first block of XR's */
	double_buffering_read_XR( &db_XR, IO_BUFF, 0, (size_t)MIN( cf.x_b, cf.m ) - 1 );
	double_buffering_swap( &db_XR );
#if VAMPIR
    VT_USER_END("READ_X");
#endif
#if VAMPIR
    VT_USER_START("READ_Y");
#endif
    /* Read first Y */
	double_buffering_read_Y( &db_Y, IO_BUFF, 0, 0 );
	double_buffering_swap( &db_Y );
#if VAMPIR
    VT_USER_END("READ_Y");
#endif

    int iter = 0;
    for ( j = 0; j < t; j++ )
    {
        /* Set the number of threads for the multi-threaded BLAS */
		set_multi_threaded_BLAS( cf.num_threads );

#if VAMPIR
        VT_USER_START("READ_Y");
#endif
        /* Read next Y */
		size_t next_j = (j+1) >= t ? 0 : j+1;
		double_buffering_read_Y( &db_Y, IO_BUFF, next_j, next_j );
#if VAMPIR
        VT_USER_END("READ_Y");
#endif

#if VAMPIR
        VT_USER_START("COMP_J");
#endif
        /* M := sigma * ( h^2 Phi - (1 - h^2) I ) */
        memcpy( M, Phi, (size_t)n * n * sizeof(double) );
		alpha = res_sigma[j] * h2[j];
        beta  = res_sigma[j] * (1 - h2[j]);
        dscal_(&nn, &alpha, M, &iONE);
        for ( i = 0; i < n; i++ )
            M[i*n + i] = M[i*n + i] + beta;

        /* L * L' = M */
        dpotrf_(LOWER, &n, M, &n, &info);
        if (info != 0)
        {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "dpotrf(M) failed (info: %d)", info);
            error_msg(err, 1);
        }

        /* XL := inv(L) * XL */
        memcpy( XL, XL_orig, wXL * n * sizeof(double) );
        dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &wXL, &ONE, M, &n, XL, &n);

#if VAMPIR
        VT_USER_START("WAIT_Y");
#endif
        // Wait until current Y is available for computation
		double_buffering_wait( &db_Y, COMP_BUFF );
#if VAMPIR
        VT_USER_END("WAIT_Y");
#endif

        /* y := inv(L) * y */
		Y_comp = double_buffering_get_comp_buffer( &db_Y );
		// Sanity check
		average( Y_comp, n, 1, cf.threshold, "TRAIT", 
				&cf.Y_fvi->fvi_data[n*NAMELENGTH], NAMELENGTH, 0 );
        dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, M, &n, Y_comp, &iONE);

        /* B_t := XL' * y */
        dgemv_(TRANS, &n, &wXL, &ONE, XL, &n, Y_comp, &iONE, &ZERO, B_t, &iONE);

        /* V_tl := XL' * XL */
        dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, XL, &n, &ZERO, V_tl, &wXL);
#if VAMPIR
        VT_USER_END("COMP_J");
#endif
		/* Solve for x_b X's at once */
        for (ib = 0; ib < m; ib += x_b) 
        {
#if VAMPIR
            VT_USER_START("READ_X");
#endif
            /* Read next block of XR's */
			size_t next_x_from = ((size_t)ib + x_b) >= m ?  0 : (size_t)ib + x_b;
			size_t next_x_to   = ((size_t)ib + x_b) >= m ? MIN( (size_t)x_b, (size_t)m ) - 1 : 
				                                           next_x_from + MIN( (size_t)x_b, (size_t)m - next_x_from ) - 1;
			double_buffering_read_XR( &db_XR, IO_BUFF, next_x_from, next_x_to );
#if VAMPIR
            VT_USER_END("READ_X");
#endif

#if VAMPIR
            VT_USER_START("WAIT_X");
#endif
            /* Wait until current block of XR's is available for computation */
			double_buffering_wait( &db_XR, COMP_BUFF );
#if VAMPIR
            VT_USER_END("WAIT_X");
#endif

            /* Set the number of threads for the multi-threaded BLAS */
			set_multi_threaded_BLAS( cf.num_threads );

#if VAMPIR
            VT_USER_START("COMP_IB");
#endif
            /* XR := inv(L) XR */
			XR_comp = double_buffering_get_comp_buffer( &db_XR );
			// Auxiliar variables
            int x_inc = MIN(x_b, m - ib);
            int rhss  = wXR * x_inc;
			// Sanity check
			average( XR_comp, n, x_inc, cf.threshold, "SNP", 
					&cf.XR_fvi->fvi_data[(n+ib)*NAMELENGTH], NAMELENGTH, 1 );
			// Computation
            dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &rhss, &ONE, M, &n, XR_comp, &n);

#if VAMPIR
            VT_USER_END("COMP_IB");
#endif

#if CHOL_MIX_PARALLELISM
            /* Set the number of threads for the multi-threaded BLAS to 1.
             * The innermost loop is parallelized using OPENMP */
			set_single_threaded_BLAS();
#endif
#if VAMPIR
            VT_USER_START("COMP_I");
#endif
            B_comp = double_buffering_get_comp_buffer( &db_B );
#if CHOL_MIX_PARALLELISM
            #pragma omp parallel for private(Bij, oneB, oneV, i, k, info, id) schedule(static) num_threads(cf.num_threads)
#endif
            for (i = 0; i < x_inc; i++)
            {
				id = omp_get_thread_num();
				oneB = &tmpBs[ id * p ];
				oneV = &tmpVs[ id * p * p ];
                Bij = &B_comp[i * size_one_b_record];

                // Building B
                // Copy B_T
                memcpy(oneB, B_t, wXL * sizeof(double));
                // B_B := XR' * y
                dgemv_("T", 
                        &n, &wXR, 
                        &ONE, &XR_comp[i * wXR * n], &n, Y_comp, &iONE, 
                        &ZERO, &oneB[wXL], &iONE);

                // Building V
                // Copy V_TL
                for( k = 0; k < wXL; k++ )
                    dcopy_(&wXL, &V_tl[k*wXL], &iONE, &oneV[k*p], &iONE); // V_TL
                // V_BL := XR' * XL
                dgemm_("T", "N",
                        &wXR, &wXL, &n,
                        &ONE, &XR_comp[i * wXR * n], &n, XL, &n,
                        &ZERO, &oneV[wXL], &p); // V_BL
                // V_BR := XR' * XR
                dsyrk_("L", "T", 
                        &wXR, &n, 
                        &ONE, &XR_comp[i * wXR * n], &n, 
                        &ZERO, &oneV[wXL * p + wXL], &p); // V_BR

                // B := inv(V) * B
                dpotrf_(LOWER, &p, oneV, &p, &info);
                if (info != 0)
                {
					for ( k = 0; k < size_one_b_record; k++ )
						Bij[k] = 0.0/0.0; //nan("char-sequence");
					continue;
                }
                dtrsv_(LOWER, NO_TRANS, NON_UNIT, &p, oneV, &p, oneB, &iONE);
                dtrsv_(LOWER,    TRANS, NON_UNIT, &p, oneV, &p, oneB, &iONE);

                /* V := res_sigma * inv( X' inv(M) X) */
                dpotri_(LOWER, &p, oneV, &p, &info);
                if (info != 0)
                {
                    char err[STR_BUFFER_SIZE];
                    snprintf(err, STR_BUFFER_SIZE, "dpotri failed (info: %d)", info);
                    error_msg(err, 1);
                }

				// Copy output
				for ( k = 0; k < p; k++ )
					Bij[k] = (float) oneB[k];
                for ( k = 0; k < p; k++ )
                    Bij[p+k] = (float)sqrt(oneV[k*p+k]);
				int idx = 0;
				for ( k = 0; k < p-1; k++ ) // Cols of V
					for ( l = k+1; l < p; l++ ) // Rows of V
					{
						Bij[p+p+idx] = (float)oneV[k*p+l];
						idx++;
					}
#if 0
			  printf("Chi square: %.6f\n", ( (oneB[p-1] / Bij[p+p-1]) * (oneB[p-1] / Bij[p+p-1]) ) );
#endif
            }
#if VAMPIR
            VT_USER_END("COMP_I");
#endif

#if VAMPIR
            VT_USER_START("WAIT_BV");
#endif
            /* Wait until the previous blocks of B's and V's are written */
            if ( iter > 0)
                double_buffering_wait( &db_B, IO_BUFF );
#if VAMPIR
            VT_USER_END("WAIT_BV");
#endif
            /* Write current blocks of B's and V's */
#if VAMPIR
            VT_USER_START("WRITE_BV");
#endif
			double_buffering_write_B( &db_B, COMP_BUFF, ib, ib+x_inc - 1, j, j );
#if VAMPIR
            VT_USER_END("WRITE_BV");
#endif

            /* Swap buffers */
			double_buffering_swap( &db_XR );
			double_buffering_swap( &db_B  );
            iter++;
        }
        /* Swap buffers */
		double_buffering_swap( &db_Y );
    }

#if VAMPIR
    VT_USER_START("WAIT_ALL");
#endif
    /* Wait for the remaining IO operations issued */
	double_buffering_wait( &db_XR, COMP_BUFF );
	double_buffering_wait( &db_Y,  COMP_BUFF );
	double_buffering_wait( &db_B,  IO_BUFF );
#if VAMPIR
    VT_USER_END("WAIT_ALL");
#endif

    /* Clean-up */
    free( M );

    free( XL );
    free( B_t  );
    free( V_tl );
    free( tmpBs );
    free( tmpVs );

	double_buffering_destroy( &db_XR );
	double_buffering_destroy( &db_Y  );
	double_buffering_destroy( &db_B  );

    return 0;
}

#endif // WITH_GPU

