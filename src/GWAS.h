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

#ifndef GWAS_H
#define GWAS_H

#define DEBUG 0
#define TIMING 0
#define VAMPIR 0
#define CHOL_MIX_PARALLELISM 1

#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define STR_BUFFER_SIZE 256

#include <stdio.h>
#include "databel.h"


// GWAS configuration and data carried throughout the computaiton
typedef struct
{
	// Paths to each data files
	// Zt* are temporary, no databel info + data
    char Phi_info_path[STR_BUFFER_SIZE];
    char Phi_data_path[STR_BUFFER_SIZE];
    char XL_info_path[STR_BUFFER_SIZE];
    char XL_data_path[STR_BUFFER_SIZE];
    char XR_info_path[STR_BUFFER_SIZE];
    char XR_data_path[STR_BUFFER_SIZE];
    char Y_info_path[STR_BUFFER_SIZE];
    char Y_data_path[STR_BUFFER_SIZE];
    char B_info_path[STR_BUFFER_SIZE];
    char B_data_path[STR_BUFFER_SIZE];
	// Temporary data
    char ZtXL_path[STR_BUFFER_SIZE];
    char ZtXR_path[STR_BUFFER_SIZE];
    char ZtY_path[STR_BUFFER_SIZE];

	// Headers
	databel_fvi *Phi_fvi;
	databel_fvi *XL_fvi;
	databel_fvi *XR_fvi;
	databel_fvi *Y_fvi;
	// Data ready for runtime computations
	//   in-core data     -> double *
	//   out-of-core data -> FILE *
	double *ests;
	double *h2; 
	double *sigma2;
	double *res_sigma2; // est_sigma2
	double *beta_ests;
	double *Phi;
	double *XL;
	double *ZtXL;
	FILE *XR;
	FILE *ZtXR;
	FILE *Y;
	FILE *ZtY;
	FILE *B;
	// Temporary data
	double *Z; // Eigenvectors of Phi
	double *W; // Eigenvalues of Phi

	// Problem dimensions
    size_t n;
    size_t p;
    size_t m;
    size_t t;
    size_t wXL; // width of XL
    size_t wXR; // width of XR
    size_t x_b; // Block-size along X
    size_t y_b; // Block-size along Y
    size_t ooc_b; // Block-size for ooc gemms

	// TO-DO
	// For now, the size of the openmp block
	// Should be the size of the ooc tile, and xb yb the size of the omp block
    size_t x_tile;
    size_t y_tile;

	size_t totalMem;
	size_t availMem;
    int num_threads;
	char var[STR_BUFFER_SIZE];
	int threshold;
} FGLS_config_t;

void initialize_config(
        FGLS_config_t *cf, 
		char *cov_base, char *phi_base, char *snp_base, char *pheno_base, char *out_base//,
//		char *var, int num_threads, int xtile, int ytile, int xb, int yb, int write_output
);

void finalize_config( FGLS_config_t *cf );

#endif // GWAS_H
