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
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

#include <omp.h>
#include <aio.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "GWAS.h"
#include "wrappers.h"
#include "timing.h"
#include "double_buffering.h"
#include "ooc_BLAS.h"
#include "GWAS.h"
#include "utils.h"
#include "databel.h"
#include "fgls_eigen.h"

#if VAMPIR
  #include "vt_user.h"
#endif

#include <stdint.h>

static void ooc_loops( eigen_loops_t *loops_t );

/*
 * Eigen-based computation of the Feasible Generalized Least-Squares problem
 */
int fgls_eigen( FGLS_config_t *cf )
{
    // Loops computation
    eigen_loops_t loops_t;

#if VAMPIR
    VT_USER_START("PRELOOP");
#endif
    // Compute the pre-loop operations
    set_multi_threaded_BLAS( cf->num_threads );

    // Premultiply Z' XR
    ooc_gemm( cf->n, cf->m * cf->wXR, cf->ooc_b, cf->Z, cf->XR_data_path, cf->ZtXR_path,
			cf->threshold, "SNP", &cf->XR_fvi->fvi_data[cf->n*NAMELENGTH], NAMELENGTH, cf->num_threads );
    
#if VAMPIR
    VT_USER_END("PRELOOP");
#endif

    // Single threaded BLAS for the loops
    set_single_threaded_BLAS();

#if VAMPIR
    VT_USER_START("LOOPS");
#endif

    cf->ZtXR = fgls_fopen( cf->ZtXR_path, "rb" );

	eigen_loops_init( &loops_t, 0, cf ); // id -> 0
	ooc_loops( &loops_t );

#if VAMPIR
    VT_USER_END("LOOPS");
#endif

	eigen_loops_destroy( &loops_t );

    return 0;
}

void ooc_loops( eigen_loops_t *loops_t ) 
{
	/*eigen_loops_t *loops_t = ( eigen_loops_t* ) in;*/
    FGLS_config_t *cf = loops_t->cf;
  
    /* Dimensions of the problem */
    size_t m = cf->m,
           n = cf->n,
           p = cf->p,
           t = cf->t,
           wXL = cf->wXL,
           wXR = cf->wXR,
           x_b = cf->x_b,
           y_b = cf->y_b;
    int nths = cf->num_threads;

	size_t size_one_b_record = cf->p + (cf->p*(cf->p+1))/2;
  
    // Pointers to simplify code and improve readibility
    double *Winv    = loops_t->Winv;
    double *XR_copy = loops_t->XR_copy;
    double *XL_copy = loops_t->XL_copy;
    double *B_t     = loops_t->B_t;
    double *V_tl    = loops_t->V_tl;
  
    // OOC data
    double *XR_comp, *Y_comp, *B_comp;
    double_buffering *db_XR = &loops_t->db_XR,
                     *db_Y  = &loops_t->db_Y,
                     *db_B  = &loops_t->db_B;
  
    /* BLAS / LAPACK constants */
    double ZERO = 0.0;
    double ONE  = 1.0;
    int   iONE  = 1;
    /* LAPACK error value */
    int info;
    // Integer version of some operands for BLAS
    int int_n   = (int)n,
        int_p   = (int)p,
        int_wXL = (int)wXL,
        int_wXR = (int)wXR;
  
    /* iterators and auxiliar vars */
    size_t ib, jb, i, j, k, l, ll;
    size_t x_inc, y_inc;

	double *XR_copy_all = malloc( nths * n * wXR * sizeof(double) );
	int id;

    // Read first block of XR's
    double_buffering_read_XR( db_XR, IO_BUFF, 0, (size_t)MIN( x_b, m ) - 1 );
    double_buffering_swap( db_XR );

    // Read first Y
	y_inc = y_b >= t ? t : y_b;
    double_buffering_read_Y( db_Y, IO_BUFF, 0, y_inc - 1 );
    double_buffering_swap( db_Y );

	for ( jb = 0; jb < t; jb += y_b )
    {
		y_inc = (jb + y_b) >= t ? t - jb : y_b;
#if VAMPIR
        VT_USER_START("READ_Y");
#endif
        // Read next Y
        size_t next_jb, next_y_inc;
		next_jb = (jb + y_b) >= t ? 0 : jb + y_b;
		next_y_inc = (jb + y_b) >= t ? 1 : MIN( y_b, t - (jb + y_b) );
        double_buffering_read_Y( db_Y, IO_BUFF, next_jb, next_jb + next_y_inc - 1 );
#if VAMPIR
        VT_USER_END("READ_Y");
#endif

        //
        // Probably all this computation could be done one by one
        // in the j loop. Would it affect performance? It would
        // reduce slightly memory consumption
        //

        /* Copy XL */
        for (ll = 0; ll < y_inc; ll++)
            memcpy( &XL_copy[ll * wXL * n], cf->ZtXL, wXL * n * sizeof(double) );
#if VAMPIR
        VT_USER_START("WAIT_Y"); // Can be delayed a bit more if necessary
#endif
        /* Wait until the current Y is available */
        double_buffering_wait( db_Y, COMP_BUFF );
#if VAMPIR
        VT_USER_END("WAIT_Y");
#endif
        Y_comp = double_buffering_get_comp_buffer( db_Y );
#if VAMPIR
        VT_USER_START("COMP_LOOP_Y_CODE");
#endif
        /* Set the scalars alpha and beta to compute: 
         *     Winv = alpha W + beta I */
        for (k = 0; k < y_inc; k++)
        {
            loops_t->alpha[k] = cf->res_sigma2[jb + k] * cf->h2[jb + k];
            loops_t->beta[k]  = cf->res_sigma2[jb + k] * (1 - cf->h2[jb + k]);
        }

        /* Winv := sqrt( inv( alpha W - beta I ) ) */
        // Possibly GER
        for (k = 0; k < y_inc; k++)
            for (l = 0; l < n; l++)
                Winv[k*n + l] = sqrt(1.0 / (loops_t->alpha[k] * cf->W[l] + loops_t->beta[k]));
		//checkNoNans(n*y_inc, Winv, "Nans in Winv\n");

        /* y := sqrt(Winv) * Z' * y */
        for (k = 0; k < y_inc; k++)
            for (l = 0; l < n; l++)
                Y_comp[k*n+l] *= Winv[k*n+l];

        /* XL := sqrt(Winv) * Z' * XL */
        for (ll = 0; ll < y_inc; ll++)
          for ( k = 0; k < wXL; k++ )
              for ( l = 0; l < n; l++ )
                  XL_copy[ ll * wXL * n + k * n + l ] *= Winv[ll * n + l];
		//checkNoNans(n*wXL, XL_copy, "Nans in XL\n");
          
        /* B_t := XL' * y */
        for (ll = 0; ll < y_inc; ll++)
          dgemv_("T", &int_n, &int_wXL, &ONE, &XL_copy[ll * wXL * n], &int_n, &Y_comp[ll * n], &iONE, &ZERO, &B_t[ll * wXL], &iONE);

        /* V_tl := Compute Top-left part of V */
        for (ll = 0; ll < y_inc; ll++)
        {
            dsyrk_("L", "T", // LOWER, NO_TRANS, 
                    &int_wXL, &int_n, // n, k
                    &ONE, &XL_copy[ll * wXL * n], &int_n, // KL KL' 
                    &ZERO, &V_tl[ll * wXL * wXL], &int_wXL); // V_TL
        }

#if VAMPIR
        VT_USER_END("COMP_LOOP_Y_CODE");
#endif
        for (ib = 0; ib < m; ib += x_b) 
        {
            x_inc = MIN(x_b, m - ib);
#if VAMPIR
            VT_USER_START("READ_X");
#endif
        /* Read the next X */
            size_t next_x_from = ((size_t)ib + x_b) >= m ?  0 : (size_t)ib + x_b;
            size_t next_x_to   = ((size_t)ib + x_b) >= m ? MIN( (size_t)x_b, (size_t)m ) - 1 : 
                                                           next_x_from + MIN( (size_t)x_b, (size_t)m - next_x_from ) - 1;
            double_buffering_read_XR( db_XR, IO_BUFF, next_x_from, next_x_to );
#if VAMPIR
            VT_USER_END("READ_X");
#endif

#if VAMPIR
            VT_USER_START("WAIT_X");
#endif
			/*read_clock( &start_read_x);*/
            /* Wait until the current X is available */
			double_buffering_wait( db_XR, COMP_BUFF );
#if VAMPIR
            VT_USER_END("WAIT_X");
#endif
            XR_comp = double_buffering_get_comp_buffer( db_XR );
			//checkNoNans(n*x_inc, XR_comp, "Nans in XR\n");

			double *Bij;
			int jt, it, xtile = cf->x_tile, ytile = cf->y_tile;
			B_comp = double_buffering_get_comp_buffer( db_B );

           	#pragma omp parallel private(Bij, i, j, it, jt, k, l, info, XR_copy, id) firstprivate(B_comp) num_threads(nths) 
			{
				id = omp_get_thread_num();
				XR_copy = &XR_copy_all[ id * wXR * n ];

				double *oneB = &loops_t->oneB[p*id];
				double *oneV = &loops_t->oneV[p*p*id];

				#pragma omp for schedule (guided) collapse(2)
				for ( jt = 0; jt < y_inc; jt+=ytile )
				{
					for (it = 0; it < x_inc; it+=xtile)
					{
						int jlim = MIN(jt + ytile, y_inc);
						int ilim = MIN(it + xtile, x_inc);
						for ( j = jt; j < jlim ; j++ )
						{
							for (i = it; i < ilim; i++)
							{
								/*for ( k = 0; k < wXR; k++ )*/
									for ( l = 0; l < n; l++ )
										XR_copy[ l ] = XR_comp[ i * n + l ] * Winv[j*n + l];
							
								Bij = &B_comp[ j*size_one_b_record*x_inc   + i*size_one_b_record];
			  
								/* Building B */
								// Copy B_T
								memcpy(oneB, &B_t[j * wXL], wXL * sizeof(double));
								// B_B := XR' * y
								dgemv_("T", 
									   &int_n, &int_wXR, 
									   &ONE, XR_copy, &int_n, &Y_comp[j* n], &iONE, 
									   &ZERO, &oneB[wXL], &iONE);

								//checkNoNans(p, oneB, "Nans in oneB\n");
			  
								  /* Building V */
								  // Copy V_TL
								  for( k = 0; k < wXL; k++ )
									  dcopy_(&int_wXL, &V_tl[j*wXL*wXL + k*wXL], &iONE, &oneV[k*p], &iONE); // TL
								  // V_BL := XR' * XL
								  dgemm_("T", "N",
										 &int_wXR, &int_wXL, &int_n, // m, n, k
										 &ONE, XR_copy, &int_n, &XL_copy[j*wXL*n], &int_n, // KR KL'
										 &ZERO, &oneV[wXL], &int_p); // BL
								  // V_BR := XR' * XR
								  dsyrk_("L", "T", 
										 &int_wXR, &int_n, 
										 &ONE, XR_copy, &int_n, 
										 &ZERO, &oneV[wXL * p + wXL], &int_p);

								//checkNoNans(p*p, oneV, "Nans in oneV\n");
			  
								  /* B := inv(V) * y 
									int ii, jj;
									for (ii = 0; ii < p; ii++)
									{
										for (jj = 0; jj < p; jj++)
											printf("%.16e ", oneV[jj*p + ii]);
										printf("\n");
									}
									printf("\n");*/
								  dpotrf_(LOWER, &int_p, oneV, &int_p, &info);
								  if (info != 0)
								  {
										for ( k = 0; k < (p+ p*(p+1)/2); k++ )
											Bij[k] = 0.0/0.0; //nan("char-sequence");
										continue;
								  }
								  dtrsv_(LOWER, NO_TRANS, NON_UNIT, &int_p, oneV, &int_p, oneB, &iONE);
								  dtrsv_(LOWER,    TRANS, NON_UNIT, &int_p, oneV, &int_p, oneB, &iONE);
			  
								  // V = res_sigma * inv( X' inv(M) X)
								  dpotri_(LOWER, &int_p, oneV, &int_p, &info);
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
						}
					}
				}
            }

#if VAMPIR
      VT_USER_START("WRITE_B");
#endif
      /* Write current B */
	      double_buffering_write_B( db_B, COMP_BUFF, ib, ib+x_inc - 1, jb, jb+y_inc - 1 );
#if VAMPIR
      VT_USER_END("WRITE_B");
#endif

#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
      if ( jb > 0 || ib > 0 )
		double_buffering_wait( db_B, IO_BUFF );
#if VAMPIR
      VT_USER_END("WAIT_B");
#endif

      /* Swap buffers */
	  double_buffering_swap( db_XR );
	  double_buffering_swap( db_B  );
    }
    /* Swap buffers */
    double_buffering_swap( db_Y );
  }
    /* Wait for the remaining IO operations issued */
#if VAMPIR
      VT_USER_START("WAIT_X");
#endif
    double_buffering_wait( db_XR, COMP_BUFF );
#if VAMPIR
      VT_USER_END("WAIT_X");
#endif
#if VAMPIR
      VT_USER_START("WAIT_Y");
#endif
    double_buffering_wait( db_Y,  COMP_BUFF );
#if VAMPIR
      VT_USER_END("WAIT_Y");
#endif
#if VAMPIR
      VT_USER_START("WAIT_B");
#endif
    double_buffering_wait( db_B,  IO_BUFF );

#if VAMPIR
      VT_USER_END("WAIT_B");
#endif

	  free(XR_copy_all);
}

// 
// EIGEN LOOPS
//
void eigen_loops_init( eigen_loops_t *loops, int id, FGLS_config_t *cf )
{
	size_t size_one_b_record = cf->p + (cf->p*(cf->p+1))/2;

	loops->alpha = ( double * ) fgls_malloc ( cf->y_b * sizeof(double) );
    loops->beta  = ( double * ) fgls_malloc ( cf->y_b * sizeof(double) );
    loops->Winv  = ( double * ) fgls_malloc ( cf->y_b * cf->n * sizeof(double) );

    loops->XR_copy = ( double * ) fgls_malloc ( cf->wXR * cf->n * cf->x_b * sizeof(double) );
    loops->XL_copy = ( double * ) fgls_malloc ( cf->wXL * cf->n * cf->y_b * sizeof(double) );
    loops->B_t  = ( double * ) fgls_malloc ( cf->wXL * cf->y_b * sizeof(double) );
    loops->V_tl = ( double * ) fgls_malloc ( cf->wXL * cf->wXL * cf->y_b * sizeof(double) );

    loops->oneB = ( double * ) fgls_malloc ( cf->p * cf->num_threads * sizeof(double) );
    loops->oneV = ( double * ) fgls_malloc ( cf->p * cf->p * cf->num_threads * sizeof(double) );

	double_buffering_init( &loops->db_XR, (size_t)cf->n * cf->wXR * cf->x_b * sizeof(double),
			                cf->ZtXR, cf ); // _fp
	double_buffering_init( &loops->db_Y,  (size_t)cf->n * cf->y_b * sizeof(double),
			                cf->ZtY,  cf );
	double_buffering_init( &loops->db_B,  (size_t)size_one_b_record * cf->x_b * cf->y_b * sizeof(double),
			                cf->B,  cf );

	loops->cf = cf;
	loops->id = id;
}

void eigen_loops_destroy( eigen_loops_t *loops )
{
	free( loops->alpha );
	free( loops->beta );
	free( loops->Winv );
	free( loops->XR_copy );
	free( loops->XL_copy );
	free( loops->B_t );
	free( loops->V_tl );
	free( loops->oneB );
	free( loops->oneV );

	double_buffering_destroy( &loops->db_XR );
	double_buffering_destroy( &loops->db_Y );
	double_buffering_destroy( &loops->db_B );
}
