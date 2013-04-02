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

#include "options.h"
#include "blas.h"
#include "lapack.h"
#include "ooc_BLAS.h"
#include "wrappers.h"
#include "utils.h"
#include "databel.h"
#include "statistics.h"
#include "optimization.h"
#include "REML.h"

// Auxiliar functions
static void eigen_mr3(int n, double *Phi, double *Z, double *W);
static void premultiplyXL( int n, int wXL, double *Z, double *XL, double *ZtXL );

/*
 * Eigendecomposition-based REML
 */
void estimates_eig( FGLS_config_t *cf )
{
	estimates_eig_t data;
	double eps; // machine epsilon
	int j, k;
	
	// Dimensions
	data.n   = cf->n;
	data.wXL = cf->wXL;
	// Load Phi
	data.Phi = cf->Phi;
	// Compute eigenvectors and eigenvalues in data.Z and data.W
	cf->Z = data.Z = fgls_malloc( data.n * data.n * sizeof(double) );
	cf->W = data.W = fgls_malloc( data.n * sizeof(double) );
	set_multi_threaded_BLAS(cf->num_threads); // Probably can be refined for the iterations below
	eigen_mr3(data.n, data.Phi, data.Z, data.W);
	// Phi should be SPD. We set the negative eigenvalues to 
	//   "epsilon * max( positive eigvalue )"
	// Eigenvalues are sorted in ascendent order:
	//   - min eigenvalue in W[0]
	//   - max eigenvalue in W[n-1]
	eps = get_epsilon();
	j = 0;
	while ( cf->W[j] < 0.0 && j < cf->n )
	{
		cf->W[j] = eps * cf->W[cf->n-1];
		j++;
	}
	// Phi is not used anymore
	// It is now that chol also uses this method to estimate h2 and res_sigma2
	/*free( cf->Phi ); cf->Phi = NULL;*/

	// Load XL
	data.X = cf->XL;
	cf->ZtXL = data.ZtX = fgls_malloc( data.n * data.wXL * sizeof(double) );
	// Premultiply ZtXL
	premultiplyXL( data.n, data.wXL, data.Z, data.X, data.ZtX );
	// Premultiply ZtY
	ooc_gemm( data.n, cf->t, cf->ooc_b, data.Z, cf->Y_data_path, cf->ZtY_path,
			cf->threshold, "TRAIT", &cf->Y_fvi->fvi_data[cf->n*NAMELENGTH], NAMELENGTH, cf->num_threads);

	// Iterate over phenotypes (Y)
	data.Y = fgls_malloc( data.n * sizeof(double) );
	data.ZtY = fgls_malloc( data.n * sizeof(double) );
	data.beta = fgls_malloc( cf->wXL * sizeof(double) );
	cf->ZtY = fgls_fopen( cf->ZtY_path, "rb" );
	for ( j = 0; j < cf->t; j++ )
	{
		sync_read( data.Y, cf->Y, cf->n, cf->n * j );
		// Sanity
		average( data.Y, cf->n, 1, cf->threshold, "TRAIT", &cf->Y_fvi->fvi_data[(cf->n+j)*NAMELENGTH], NAMELENGTH, 0, cf->num_threads );

		sync_read( data.ZtY, cf->ZtY, cf->n, cf->n * j );
		data.sigma = variance(data.Y, data.n);

		// range [a,b] = [0,0.99], tolerance = 10^{-8}
		minimize( 0.0, 0.99, 1e-8, &polygenic_REML_logLik_eig_wrapper, (void*)&data );
#if 0 //DEBUG
		printf("h: %.15e - sigma: %.15e - res_sigma: %.15e\n", data.h, data.sigma, data.res_sigma);
#endif
		cf->h2[j] = data.h;
		cf->sigma2[j] = data.sigma;
		cf->res_sigma2[j] = data.res_sigma;
		for ( k = 0; k < cf->wXL; k++ )
			cf->beta_ests[cf->wXL*j+k] = data.beta[k];
	}

	// Encapsulate (maybe reuse)
	free( data.Y );
	free( data.beta );
	free( data.ZtY );
}

double polygenic_REML_logLik_eig (
        int n, int widthXL,
		double *X, double *Y,
        double *Z, double *W, double *ZtX, double *ZtY,
		double sigma, double h2,
        double *loglik, double *res_sigma, double *beta
)
{
    double alpha, gamma, // to scale W
           *ZtY_upd, *ZtX_upd, *YmXB, *V, *D, det, // temporary
           ZERO = 0.0,
           ONE = 1.0,
           MINUS_ONE = -1.0;

    int    wXL = widthXL,
           iONE = 1, info,
           i, j;

    V       = (double *) fgls_malloc ( wXL * wXL * sizeof(double) );
	ZtY_upd = (double *) fgls_malloc ( n * sizeof(double) );
	ZtX_upd = (double *) fgls_malloc ( n * wXL * sizeof(double) );
	D = (double *) fgls_malloc ( n * sizeof(double) );
	YmXB = (double *) fgls_malloc ( n * sizeof(double) );

	/* 
	 * W = sqrt( inv ( W ) ) 
	 * det(M) = sum ( log( eigvalues(M) ) ) 
	 */
    alpha = h2; // * sigma2;
    gamma = (1 - h2); // * sigma2;
	det = 0.0;
	for ( i = 0; i < n; i++ )
	{
		D[i] = alpha * W[i] + gamma;
		det += log( D[i] );
		D[i] = sqrt( 1.0 / D[i] );
	}

	/* W * ZtX */
	/* W * ZtY */
	for (i = 0; i < n; i++ )
	{
		for (j = 0; j < wXL; j++ )
			ZtX_upd[j*n + i] = ZtX[j*n + i] * D[i];
		ZtY_upd[i] = ZtY[i] * D[i];
	}

    /* 5) beta := X' * y */
    dgemv_(TRANS, &n, &widthXL, &ONE, ZtX_upd, &n, ZtY_upd, &iONE, &ZERO, beta, &iONE);

    /* 6) V := X' * X */
	dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, ZtX_upd, &n, &ZERO, V, &wXL);

	/*dposv_(LOWER, &wXL, &iONE, V, &wXL, beta, &wXL, &info);*/
	dpotrf_(LOWER, &wXL,  V, &wXL, &info);
    if (info != 0)
    {
        fprintf(stderr, __FILE__ ": [ERROR] inv(V) y failed - info: %d\n", info);
        exit(-1);
    }
    /* beta := inv( V ) beta */
	dtrsv_(LOWER, NO_TRANS, NON_UNIT, &wXL, V, &wXL, beta, &iONE);
	dtrsv_(LOWER,    TRANS, NON_UNIT, &wXL, V, &wXL, beta, &iONE);

    // y - X beta
	memcpy( YmXB, Y, n * sizeof(double) );
    dgemv_( NO_TRANS, 
            &n, &wXL, 
            &MINUS_ONE, X, &n, beta, &iONE,
            &ONE, YmXB, &iONE );

    // residual sigma and loglik
	*res_sigma = variance( YmXB, n );
	// loglik = a + b
	//  a -> log(det(M))
	//  b -> YmXB' inv(M) YmXB
    dgemv_(TRANS, 
            &n, &n, 
            &ONE, Z, &n, YmXB, &iONE,
            &ZERO, ZtY_upd, &iONE);
	// YmXB' * inv( M ) * YmXB
	*loglik = 0.0;
	for (i = 0; i < n; i++ )
		*loglik += ZtY_upd[i] * D[i] * D[i] * ZtY_upd[i];
	*loglik = (*loglik) / (*res_sigma); // <- b
	// a + b -> det + loglik above
	//   det incomplete, we must take into account the sigma scaling M -> "+ n * log(*res_sigma)"
    *loglik = det + n * log(*res_sigma) + (*loglik) ;

    // Clean up
	free( V );
	free( ZtY_upd );
	free( ZtX_upd );
	free( YmXB );
	free( D );

	return *loglik;
}


double polygenic_REML_logLik_eig_wrapper ( double h, void *data )
{
	estimates_eig_t *eig_t = (estimates_eig_t*)data;
	size_t n   = eig_t->n, 
	       wXL = eig_t->wXL;
	double *X = eig_t->X,
		   *Y = eig_t->Y,
	       *Z   = eig_t->Z,
		   *W   = eig_t->W,
	       *ZtX = eig_t->ZtX,
		   *ZtY = eig_t->ZtY;
	double sigma = eig_t->sigma;
	double loglik;

	polygenic_REML_logLik_eig (
	        n, wXL,
	        X, Y, Z, W, ZtX, ZtY,
            sigma, h,
	        &loglik, &eig_t->res_sigma, eig_t->beta
	);

	eig_t->h = h;

#if 0
	printf("Loglik: %.15e\n", loglik);
	int i;
	for ( i = 0; i < wXL; i++ )
		printf("Beta[%d]: %.15e\n", i, eig_t->beta[i]);
#endif
	return loglik;
}

/*
 * MR3 eigendecomposition
 */
void eigen_mr3(int n, double *Phi, double *Z, double *W)
{
    int nb = 192;
    int idummy, nCompPairs, *isuppz, *iwork, info,
        lwork = n * (nb + 6),
        liwork = 10 * n;
    double ddummy = -1.0, *work;

    work = (double *) fgls_malloc ( lwork * sizeof(double) );
    iwork   = (int *) fgls_malloc ( liwork * sizeof(int) );
    isuppz  = (int *) fgls_malloc ( 2 * n * sizeof(int) );

    dsyevr_("V", "A", "L", &n, Phi, &n, 
			&ddummy, &ddummy, &idummy, &idummy, &ddummy, 
			&nCompPairs, W, Z, &n, isuppz, 
			work, &lwork, iwork, &liwork, &info);
    if (info != 0)
    {
        fprintf(stderr, "Error factoring Phi\n");
        exit(-1);
    }

	free( work );
	free( iwork );
	free( isuppz );
}

// Assumes leading dimension = n for every matrix
void premultiplyXL( 
		int n, int wXL,
		double *Z, // Eigenvectors
		double *XL,
		double *ZtXL
)
{
	int ldZ, ldXL, ldZtXL;
	
	double ZERO = 0.0,
	       ONE  = 1.0;

	ldZ = ldXL = ldZtXL = n;

	// Z' * X
	dgemm_(
			TRANS, NO_TRANS, 
			&n, &wXL, &n, 
			&ONE, Z, &ldZ, XL, &ldXL, 
			&ZERO, ZtXL, &ldZtXL
	);
}
