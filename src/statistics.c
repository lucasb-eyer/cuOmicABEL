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

#include <assert.h>

#include "statistics.h"

double mean( double *v, int n )
{
	int i;
	double mean = 0.0;

	assert( n > 0 );

	for ( i = 0; i < n; i++ )
		mean += v[i];
	mean /= n;

	return mean;
}

double variance( double *v, int n )
{
	int i;
	double m, var = 0.0;

	assert( n > 1 );

	m = mean( v, n );
	for ( i = 0; i < n; i++ )
		var += (m - v[i]) * (m - v[i]);
	var /= (n - 1);

	return var;
}
