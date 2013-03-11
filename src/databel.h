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

#ifndef DATABEL_H
#define DATABEL_H

#define UNSIGNED_SHORT_INT_TYPE
#define SHORT_INT_TYPE
#define UNSIGNED_INT_TYPE
#define INT_TYPE
#define FLOAT_TYPE
#define DOUBLE_TYPE
#define SIGNED_CHAR_TYPE
#define UNSIGNED_CHAR_TYPE

#define NAMELENGTH 32
#define RESERVEDSPACE 5

typedef struct databel_fvi_header
{
	unsigned short int type;
	unsigned int nelements;
	unsigned int numObservations;
	unsigned int numVariables;
	unsigned int bytesPerRecord;
	unsigned int bitsPerRecord;
	unsigned int namelength;
	unsigned int reserved[RESERVEDSPACE];
} databel_fvi_header;

typedef struct databel_fvi
{
	databel_fvi_header  fvi_header;
	char               *fvi_data;
} databel_fvi;

struct databel_fvi * load_databel_fvi( char *path );
void free_databel_fvi( struct databel_fvi **fvi );

#endif // DATABEL_H
