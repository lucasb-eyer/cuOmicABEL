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

#ifndef REML_H
#define REML_H

#include "GWAS.h"

/*
 * REML based on the Eigendecomposition approach
 *
 * Data structure, Eigen REML, and wrapper
 */
typedef struct
{
	// Problem dimensions
	size_t n;
	size_t wXL;

	// Data to fit
	double *Phi;
	double *X;
	double *Y;
	double sigma;

	// Precomputed values
	//
	// Z are the eigenvectors of Phi
	// W are the eigenvalues of Phi
	double *Z;
	double *W;
	double *ZtX;
	double *ZtY;

	// computed estimates
	double h;
	double res_sigma;
	double *beta;
} estimates_eig_t;

void estimates_eig( FGLS_config_t *cf );
double polygenic_REML_logLik_eig_wrapper ( double h, void *data );
double polygenic_REML_logLik_eig(
        int n, int widthXL,
        double *X, double *Y, 
		double *Z, double *W, double *ZtX, double *ZtY,
        double sigma2, double h2,
        double *loglik, double *res_sigma, double *beta
);


#endif // REML_H
