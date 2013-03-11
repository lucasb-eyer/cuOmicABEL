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

#include <stdlib.h>

#include "databel.h"
#include "wrappers.h"

struct databel_fvi * load_databel_fvi( char *path )
{
	FILE *f;
	databel_fvi *fvi;
	size_t data_size;

	f = fgls_fopen( path, "r" );
	
	fvi = (databel_fvi*) fgls_malloc( sizeof(databel_fvi) );
	// Header
	fread( &fvi->fvi_header, sizeof(databel_fvi_header), 1, f );
	// Labels
	data_size = (fvi->fvi_header.numVariables +fvi->fvi_header.numObservations ) * 
		         fvi->fvi_header.namelength * sizeof(char);
	fvi->fvi_data = (char *) fgls_malloc ( data_size );
	// Load labels
	fread( fvi->fvi_data, 1, data_size, f );

	fclose( f );

	return fvi;
}

void free_databel_fvi( struct databel_fvi **fvi )
{
	free ((*fvi)->fvi_data);
	free (*fvi);
	*fvi = NULL;
}

