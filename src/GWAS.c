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

#include <sys/time.h>

#include <assert.h>

#include "utils.h"
#include "wrappers.h"
#include "GWAS.h"

static void free_and_nullify( double **p );
static void close_and_nullify( FILE **f );
static void load_databel_info( FGLS_config_t *cf );

// 
// GWAS config + data
//
void initialize_config(
        FGLS_config_t *cf, 
		char *cov_base, char *phi_base, char *snp_base, char *pheno_base, char *out_base//,
//		char *var, int num_threads, int xtile, int ytile, int xb, int yb, int write_output
)
{
	load_databel_info( cf );
    // Problem dimensions
	cf->n = cf->XL_fvi->fvi_header.numObservations;
	cf->p = cf->XL_fvi->fvi_header.numVariables + 1; // Intercept included
	cf->m = cf->XR_fvi->fvi_header.numVariables;
	cf->t = cf->Y_fvi->fvi_header.numVariables;
	// Assuming wXR = 1
	cf->wXL = cf->p - 1;
	cf->wXR = 1;
    // Algorithm parameters
	/*cf->x_b = x_b;*/
	/*cf->y_b = y_b;*/
	/*cf->num_threads = num_threads;*/
	get_main_memory_size( &cf->totalMem, &cf->availMem );
	estimate_block_sizes( cf, cf->var, cf->x_b == -1 || cf->y_b == -1 ); // if any is -1, estimate them
/*	if ( !(cf->x_b == -1 || cf->y_b == -1) )
	{
		cf->x_b = xb;
		cf->y_b = yb;
	}
	cf->x_tile = xtile;
	cf->y_tile = ytile;
*/
	/*cf->x_b = 10;*/
	/*cf->y_b = 10;*/
    // In/Out Files
  /*  sprintf( cf->Phi_data_path, "%s.fvd",  phi_base );
    sprintf( cf->Phi_info_path, "%s.fvi",  phi_base );
    sprintf( cf->XL_data_path,  "%s.fvd",  cov_base );
    sprintf( cf->XL_info_path,  "%s.fvi",  cov_base );
    sprintf( cf->XR_data_path,  "%s.fvd",  snp_base );
    sprintf( cf->XR_info_path,  "%s.fvi",  snp_base );
    sprintf( cf->Y_data_path,   "%s.fvd",  pheno_base );
    sprintf( cf->Y_info_path,   "%s.fvi",  pheno_base );
	if ( write_output )
	{
		sprintf( cf->B_data_path, "%s.out",  out_base );
		sprintf( cf->B_info_path, "%s.iout", out_base );
	}
	else
	{
		sprintf( cf->B_data_path,      "/dev/null" );
		sprintf( cf->B_info_path,      "/dev/null" );
	}
*/	
    // Temporary files
    snprintf( cf->ZtXL_path, STR_BUFFER_SIZE, "%s.tmp", cov_base );
    snprintf( cf->ZtXR_path, STR_BUFFER_SIZE, "%s.tmp", snp_base );
    snprintf( cf->ZtY_path,  STR_BUFFER_SIZE, "%s.tmp", pheno_base );

	// Allocate memory and load in-core data
	// NOTE: only data which is input to GWAS
	//    
	// Intermediate data as ZtXL will be handled on the fly
	// by the corresponding routines:
	//   * ZtXL, only available after REML_eigen (if so)
	//   * Z and W will be allocated and filled up during
	//     the first eigen-decomposition of Phi, (if so)
	
	// ests needed in both cases chol and eigen
	cf->ests = (double *) fgls_malloc ( (3 + cf->wXL) * cf->t * sizeof(double) );
	cf->h2         =  cf->ests;
	cf->sigma2     = &cf->ests[ cf->t ];
	cf->res_sigma2 = &cf->ests[ 2 * cf->t ];
	cf->beta_ests  = &cf->ests[ 3 * cf->t ];

	cf->Phi = (double *) fgls_malloc ( cf->n * cf->n * sizeof(double) );
	load_file ( cf->Phi_data_path,  "rb", cf->Phi,  cf->n * cf->n );
	// Sanity check (Phi)
	checkNoNans(cf->n * cf->n, cf->Phi, "[ERROR] NaNs not allowed in Phi\n");
	cf->XL   = (double *) fgls_malloc ( cf->n * cf->wXL * sizeof(double) );
	load_file ( cf->XL_data_path,   "rb", cf->XL,   cf->n * cf->wXL );
	// Sanity check (XL)
	average( cf->XL, cf->n, cf->wXL, cf->threshold, "Covariate", 
			&cf->XL_fvi->fvi_data[cf->n*NAMELENGTH], NAMELENGTH, 1, cf->num_threads );
	cf->ZtXL = NULL;
	cf->Z    = NULL;
	cf->W    = NULL;

	// Open files for out-of-core data
	// Again, only for data input/output to GWAS. Intermediate, on the fly
	cf->XR = fgls_fopen( cf->XR_data_path, "rb" );
	cf->ZtXR = NULL;
	cf->Y  = fgls_fopen( cf->Y_data_path, "rb" );
	cf->ZtY = NULL;
	cf->B = fgls_fopen( cf->B_data_path, "wb" );

    return;
}

void finalize_config( FGLS_config_t *cf )
{
	// Free databel headers (fvi)
	free_databel_fvi( &cf->Phi_fvi );
	free_databel_fvi( &cf->XL_fvi );
	free_databel_fvi( &cf->XR_fvi );
	free_databel_fvi( &cf->Y_fvi );
	// Free in-core data buffers
	free_and_nullify( &cf->ests );
	cf->h2 = NULL;
	cf->sigma2 = NULL;
	cf->res_sigma2 = NULL;
	cf->beta_ests = NULL;
	free_and_nullify( &cf->Phi );
	free_and_nullify( &cf->XL );
	free_and_nullify( &cf->ZtXL );
	free_and_nullify( &cf->Z );
	free_and_nullify( &cf->W );
    // Close out-of-core data files
	close_and_nullify( &cf->XR );
	close_and_nullify( &cf->ZtXR );
	close_and_nullify( &cf->Y );
	close_and_nullify( &cf->ZtY );
	close_and_nullify( &cf->B );
//	close_and_nullify( &cf->V );
}

static void free_and_nullify( double **p )
{
	if ( *p )
	{
		free( *p );
		*p = NULL;
	}
}

static void close_and_nullify( FILE **f )
{
	if ( *f )
	{
		fclose( *f );
		*f = NULL;
	}
}


void load_databel_info( FGLS_config_t *cf )
{
	cf->Phi_fvi = load_databel_fvi( cf->Phi_info_path );
	cf->XL_fvi = load_databel_fvi( cf->XL_info_path );
	cf->XR_fvi = load_databel_fvi( cf->XR_info_path );
	cf->Y_fvi = load_databel_fvi( cf->Y_info_path );
	// Check integrity, i.e., dimension n matches
	assert( cf->XL_fvi->fvi_header.numObservations == cf->XR_fvi->fvi_header.numObservations );
	assert( cf->XL_fvi->fvi_header.numObservations == cf->Y_fvi->fvi_header.numObservations );
	assert( cf->XL_fvi->fvi_header.numObservations == cf->Phi_fvi->fvi_header.numObservations );
	assert( cf->XL_fvi->fvi_header.numObservations == cf->Phi_fvi->fvi_header.numVariables );
	// Assert type is double
}
