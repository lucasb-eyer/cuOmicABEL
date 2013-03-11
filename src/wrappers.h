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

#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <stdio.h>

// Wrappers for:
//   - fwrite
//   - fread
//   - fopen
void sync_write(double *buffer, FILE* file, size_t size, size_t start);
void sync_read(double *buffer, FILE* file, size_t size, size_t start);
FILE * fgls_fopen( const char *path, const char *mode );

// Loads an entire file
void load_file( const char *path, const char *mode, double *buff, size_t size );

// Guess :P
#define fgls_malloc(size) fgls_malloc_impl(__FILE__, __LINE__, size)
void * fgls_malloc_impl( const char* file, long line, size_t size );

// Error handling - To implement properly
void error_msg(char *msg, int abort);

#endif // WRAPPERS_H
