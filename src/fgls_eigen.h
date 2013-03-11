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

#ifndef FGLS_EIGEN_H
#define FGLS_EIGEN_H

#include "GWAS.h"
#include "double_buffering.h"

int fgls_eigen( FGLS_config_t *cf );

// Data structure for the computation of the loops in the
// GWAS_eigen algorithm
//
// Extends the basic configuration with additional temporary variables, and
// double buffering handlers
//
// One of this structure is needed per thread
typedef struct {
	// Temporary
    double *alpha; // sigma2 * h
    double *beta;  // sigma2 * ( 1 - h )
    double *Winv;  // Inverse of the shifted eigenvalues

	double *XR_copy; // Blocks of XR are reused across the a block of Ys
    double *XL_copy; // Same for XL
    double *B_t;  // Fixed top part of B
    double *V_tl; // Fixed top-left part of V
	
	double *oneB;
	double *oneV;

	// Double buffering
	double_buffering db_XR;
	double_buffering db_Y;
	double_buffering db_B;
//	double_buffering db_V;

    // Configuration
    FGLS_config_t *cf;

	// Thread id
    int id;
} eigen_loops_t;

void eigen_loops_init( eigen_loops_t *loops, int id, FGLS_config_t *cf );
void eigen_loops_destroy( eigen_loops_t *loops );


#endif // FGLS_EIGEN_H
