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

#include <sys/time.h>

#if defined OSX
  #include <sys/sysctl.h>
#endif

#include <omp.h>
#if defined MKL
  #include <mkl.h>
#endif

#ifdef WITH_GPU
  #include <cuda_runtime_api.h>
  #include <cublas_v2.h>
#endif

#include "GWAS.h"
#include "wrappers.h"
#include "utils.h"

void set_multi_threaded_BLAS( int nths )
{
    char nths_str[STR_BUFFER_SIZE];

    snprintf( nths_str, STR_BUFFER_SIZE, "%d", nths );
#if defined GOTO
    setenv("GOTO_NUM_THREADS", nths_str, 1);
    setenv("OMP_NUM_THREADS", "1", 1);
#elif defined MKL
	mkl_set_num_threads(nths); // Set MKL to use nths for multithreaded BLAS
	omp_set_num_threads(nths); // Set OMP to use nths for the openmp parallel directives
#else
    setenv("OMP_NUM_THREADS", nths_str, 1);
#endif
}

void set_single_threaded_BLAS( void )
{
#if defined GOTO
    setenv("GOTO_NUM_THREADS", "1", 1);
#elif defined MKL
	mkl_set_num_threads(1);
	omp_set_num_threads(1);
#else
    setenv("OMP_NUM_THREADS",  "1", 1);
#endif
}

double get_epsilon( void )
{
	double eps = 1.0;

	while ( (1.0 + eps/2.0) > 1.0 )
		eps = eps / 2.0;

	return eps;
}

void get_main_memory_size( size_t *totalMem, size_t *availMem )
{
    FILE *fp;
    int found = 0;
    size_t size = 0;
    char buff[STR_BUFFER_SIZE], field[STR_BUFFER_SIZE], unit[STR_BUFFER_SIZE];

#if defined LINUX
    fp = fgls_fopen( "/proc/meminfo", "rb" );
    fgets( buff, STR_BUFFER_SIZE, fp );
    while ( found < 2 && !feof(fp) ) // Until 2 pieces of info found: Total Mem and Available Mem
    {
        sscanf( buff, "%s %zu %s", field, &size, unit);
        if ( !strncmp( field, "MemFree:", STR_BUFFER_SIZE ) )
        {
            if ( !strncmp( unit, "kB", STR_BUFFER_SIZE ) )
                *availMem = size * 1024;
			else if ( !strncmp( unit, "mB", STR_BUFFER_SIZE ) )
                *availMem = size * 1024 * 1024;
            found++;
        }
		if ( !strncmp( field, "MemTotal:", STR_BUFFER_SIZE ) )
        {
            if ( !strncmp( unit, "kB", STR_BUFFER_SIZE ) )
                *totalMem = size * 1024;
			else if ( !strncmp( unit, "mB", STR_BUFFER_SIZE ) )
                *totalMem = size * 1024 * 1024;
            found++;
        }
        fgets( buff, STR_BUFFER_SIZE, fp );
    }
    fclose( fp );
#elif defined OSX
    int mib[2];
    u_int namelen;
    size_t len;

    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    namelen = 2;
    len = sizeof(*totalMem);
    sysctl(mib, namelen, totalMem, &len, NULL, 0);

    mib[0] = CTL_HW;
    mib[1] = HW_USERMEM;
    namelen = 2;
    len = sizeof(*availMem);
    sysctl(mib, namelen, availMem, &len, NULL, 0);
#else
   // Something
#endif
}

void estimate_block_sizes( FGLS_config_t *cf, const char* var, int estimate_inc )
{
    size_t avail_mem = cf->availMem;
    size_t y_b = 100, ratio = 50, x_b = ratio * y_b;
	size_t ooc_b;
	int nths = cf->num_threads;

    size_t mem_usage_per_th;
    int converged = 0;
    if ( !strcmp( var, "eigen" ) )
    {
        avail_mem = avail_mem * .95; // "Safety"
        avail_mem = avail_mem - ( cf->n * cf->n ); // Z
        avail_mem = avail_mem - 100 * cf->n; // Extra little vars
		/*avail_mem = avail_mem / cf->num_threads; // In eigen, blocks indicate memory per thread*/

		x_b = y_b = 25600;
        // Estimate (x_b, y_b), favoring rectangular shape (column panel)
		if (estimate_inc)
		{
			while ( ! converged )
			{
				// Check XR_copy
				// Add oneB and oneV * nths
				mem_usage_per_th = ( y_b + // alpha
									 y_b + // beta
									 y_b * cf->n + // Winv
									 /*cf->wXR * cf->n * x_b + // XR_copy*/
									 cf->wXR * cf->n * x_b + // XR_copy
									 cf->wXL * cf->n * y_b + // XL_copy
									 cf->wXL * y_b +         // B_t
									 cf->wXL * cf->wXL * y_b + // V_tl
									 // Double buffering
									 2 * x_b * cf->n * cf->wXR +  // XR
									 2 * y_b * cf->n +            // Y
									 2 * x_b * y_b * (cf->p + cf->p*(cf->p+1)/2)     // B
									 /*2 * x_b * y_b * cf->p +      // B*/
									 /*2 * x_b * y_b *cf->p * cf->p // V*/
								   ) * sizeof(double);
				if ( mem_usage_per_th < avail_mem )
				{
					cf->x_b = MIN( cf->m, x_b );
					cf->y_b = MIN( cf->t, y_b );
					converged = 1;
				}
				else if ( y_b == 1 )
				{
					fprintf( stderr, "[ERROR] Not enought memory (y_b). Please try a system with larger memory.\n" );
					/*fprintf( stderr, "        For testing purposes, please reduce the value of (sample size).\n" );*/
					exit(EXIT_FAILURE);
				}
				else 
				{
					x_b = y_b = y_b / 2;
					/*if ( y_b < 64 ) // 32 on*/
					/*ratio = 20;*/
					/*else */
					if ( y_b < (160*sqrt(nths))  ) // 32 on
					{
						fprintf( stderr, "[WARNING] small value for (x_b, y_b), this might affect performance.\n" );
						fprintf( stderr, "          We suggest to run this test in a system with a larger amount of memory.\n" );
						/*fprintf( stderr, "          For testing purposes, if you observe delays, try reducing n (sample size).\n" );*/
					}
					/*x_b = y_b * ratio;*/
				}
			}
		}
		// Estimate the block size for the ooc gemms
		ooc_b = 1024 * nths; // 1024 columns per thread
		avail_mem = avail_mem / 2; // double buffering for ooc
		converged = 0;
		while ( !converged )
		{
			if ( ooc_b * cf->n < avail_mem ) // 2 buffers, n columns (Z' Y or Z' XR)
			{
				cf->ooc_b = ooc_b; // * cf->num_threads;
				converged = 1;
			}
			else if ( ooc_b == 1 )
			{
				fprintf( stderr, "[ERROR] Not enought memory (ooc_b). Please try a system with larger memory.\n" );
				/*fprintf( stderr, "        For testing purposes, please reduce the value of (sample size).\n" );*/
				exit(EXIT_FAILURE);
			}
			else
			{
				ooc_b = ooc_b / 2;
				if ( ooc_b < 512 )
				{
					fprintf( stderr, "[WARNING] ooc_b below 512, this might affect performance.\n" );
					fprintf( stderr, "          We suggest to run this test in a system with a larger amount of memory.\n" );
					/*fprintf( stderr, "          For testing purposes, if you observe delays, try reducing n (sample size).\n" );*/
				}
			}
		}
    }
	else if ( !strcmp( var, "chol" ) )
	{
		avail_mem = cf->availMem;
        avail_mem = avail_mem * .95; // "Safety"
        avail_mem = avail_mem - 2 * ( cf->n * cf->n ); // Phi and M
        avail_mem = avail_mem - 100 * cf->n; // Extra little vars

		y_b = 1, x_b = 1024 * cf->num_threads;
		/*ooc_b = 0;*/
		// Estimate the block size for the ooc gemms
		ooc_b = 1024 * nths; // 1024 columns per thread
		avail_mem = avail_mem / 2; // double buffering for ooc
		converged = 0;
		while ( !converged )
		{
			/*printf("ooc_b: %d\n", ooc_b);*/
			if ( ooc_b * cf->n < avail_mem ) // 2 buffers, n columns (Z' Y or Z' XR)
			{
				cf->ooc_b = ooc_b; // * cf->num_threads;
				converged = 1;
			}
			else if ( ooc_b == 1 )
			{
				fprintf( stderr, "[ERROR] Not enought memory (ooc_b). Please try a system with larger memory.\n" );
				/*fprintf( stderr, "        For testing purposes, please reduce the value of (sample size).\n" );*/
				exit(EXIT_FAILURE);
			}
			else
			{
				ooc_b = ooc_b / 2;
				if ( ooc_b < 512 )
				{
					fprintf( stderr, "[WARNING] ooc_b below 512, this might affect performance.\n" );
					fprintf( stderr, "          We suggest to run this test in a system with a larger amount of memory.\n" );
					/*fprintf( stderr, "          For testing purposes, if you observe delays, try reducing n (sample size).\n" );*/
				}
			}
		}

		converged = 0;
		size_t mem_usage;
        // Estimate x_b, y_b = 1 in chol
		if (estimate_inc)
		{
			while ( ! converged )
			{
				// Double buffering
				mem_usage = 2 * ( 
						// Check little vars (negligible, just accuracy)
								 x_b * cf->n * cf->wXR +  // XR
								 y_b * cf->n +            // Y
								 x_b * y_b * (cf->p + cf->p*(cf->p+1)/2) // B
								 /*x_b * y_b * cf->p +      // B*/
								 /*x_b * y_b *cf->p * cf->p // V*/
							) * sizeof(double);
				if ( mem_usage < avail_mem )
				{
					cf->x_b = MIN( cf->m, x_b );
					cf->y_b = y_b;
					cf->ooc_b = ooc_b;
					converged = 1;
				}
				else if ( x_b == 1 )
				{
					fprintf( stderr, "[ERROR] Not enought memory (x_b). Please try a system with larger memory.\n" );
					/*fprintf( stderr, "        For testing purposes, please reduce the value of (sample size).\n" );*/
					exit(EXIT_FAILURE);
				}
				else 
				{
					x_b = x_b / 2;
					if ( x_b < 512 * cf->num_threads )
					{
						fprintf( stderr, "[WARNING] x_b below 512 per thread, this might affect performance.\n" );
						fprintf( stderr, "          We suggest to run this test in a system with a larger amount of memory.\n" );
						/*fprintf( stderr, "          For testing purposes, if you observe delays, try reducing n (sample size).\n" );*/
					}
				}
			}
		}
	}
#ifdef WITH_GPU
	else if ( !strcmp( var, "chol_gpu" ) ) {
		// TODO(lucasb): This variant needs three buffers in main memory and two on the GPU.
		estimate_block_sizes( cf, "chol", estimate_inc );

		// See how many GPUs we have and what's the least memory available we can use.
		int ngpus = 0;
		if(cudaGetDeviceCount(&ngpus) != cudaSuccess) {
			fprintf( stderr, "[ERROR] Cannot determine the number of GPUs.\n");
			exit(EXIT_FAILURE);
		}

		size_t maxmem = 0;
		int igpu = 0;
		for(igpu = 0 ; igpu < ngpus ; ++igpu) {
			size_t free, total;
			cudaSetDevice(igpu);
			if(cudaMemGetInfo(&free, &total) != cudaSuccess) {
				fprintf(stderr, "[ERROR] Cannot get free memory amount of GPU %d.\n", igpu);
				exit(EXIT_FAILURE);
			}
			maxmem = free > maxmem ? free : maxmem;
		}

		// How much of that memory could we use for a block?
		// There are only two blocks and the triangular L matrix on each device.
		if(maxmem < cf->n*cf->n) {
			fprintf(stderr, "[ERROR] Need at the very least %ld MB on the GPU but it can only hold %ld MB!\n", cf->n*cf->n*sizeof(double)/1024/1024, maxmem/1024/1024);
			exit(EXIT_FAILURE);
		}
		size_t maxblock = (maxmem - cf->n*cf->n)/2;

		// Each GPU can handle such a block.
		cf->x_b = MIN(cf->x_b, maxblock*ngpus);

		// The problem with this is that IT STILL DOESN'T MEAN IT WILL WORK!
		// Because of memory fragmentation, it might not be possible to allocate
		// blocks of this size.
	}
#endif
}

void print_timestamp( void )
{
	struct timeval t;

	gettimeofday( &t, NULL );
	printf( "%lld:%lld seconds\n", (long long)t.tv_sec, (long long)t.tv_usec );
}

void average( double *data, int n, int ncols, int threshold, const char *obj_type, char *obj_name, int namelength, int verbose )
{
	int i, j;
	double sum, avg;
	int nans, infs;

	for ( j = 0; j < ncols; j++ )
	{
		sum = 0.0;
		nans = 0;
		infs = 0;
		for ( i = 0; i < n; i++ )
		{
			if ( isnan(data[j*n + i]) )
				nans++;
			else if ( isinf(data[j*n + i]) )
			{
				fprintf(stderr, "[WARNING] Infinity found in %s %s\n", obj_type, &obj_name[j*namelength]);
				infs++;
			}
			else
				sum += data[j*n + i];
		}
		avg = sum / (n-nans-infs);
		for ( i = 0; i < n; i++ )
		{
			if ( isnan(data[j*n + i]) || isinf(data[j*n + i]) )
				data[j*n + i] = avg;
		}
		if ( ((float)nans / n) > (100 - threshold)/(float)100 && verbose )
			fprintf(stderr, "\n[WARNING] %.1f%% of NaNs in %s %s exceeds threshold (%d%%)\n", 
					100*(float)nans/n, obj_type, &obj_name[j*namelength], 100-threshold);
	}
}

void checkNoNans( size_t n, double *buff, const char* err_msg)
{
	size_t i;
	for ( i = 0; i < n; i++ )
		if ( isnan( buff[i] ) || isinf( buff[i] ) )
		{
			fprintf( stderr, err_msg );
			exit( EXIT_FAILURE );
		}
}
