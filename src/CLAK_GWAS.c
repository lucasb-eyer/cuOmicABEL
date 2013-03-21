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

#include <sys/stat.h>

#include "wrappers.h"
#include "utils.h"
#include "timing.h"
#include "REML.h"
#include "fgls_chol.h"
#include "fgls_chol_gpu.h"
#include "fgls_eigen.h"

void usage( void );
int parse_input( int argc, char *argv[], char *var, 
		         char *cov_base, char *phi_base, char *snp_base, char *pheno_base, char *out_base,
				 int *nths, int *thres, int *xtile, int *ytile, int *xb, int *yb, int *write_output );
void check_input_integrity( FGLS_config_t *cf, char *var,
				         char *cov_base, char *phi_base, char *snp_base, char *pheno_base, char *out_base,
						 int nths, int thres, int xtile, int ytile, int xb, int yb, int write_output );
void print_info( FGLS_config_t *cf );
void append_estimates( FGLS_config_t *cf );
void write_output_info_file( FGLS_config_t *cf );
void cleanup( FGLS_config_t *cf );

int main( int argc, char *argv[] ) 
{
    char cov_base[STR_BUFFER_SIZE] = "", 
		 phi_base[STR_BUFFER_SIZE] = "", 
		 snp_base[STR_BUFFER_SIZE] = "", 
		 pheno_base[STR_BUFFER_SIZE] = "", 
		 out_base[STR_BUFFER_SIZE] = "", 
		 var[STR_BUFFER_SIZE] = "chol";
	int nths = 1, write_output = 1;
	int xb = -1, yb = -1, xtile = 160, ytile = 160;
	int thres = 95;
	int status;

    FGLS_config_t cf;
    struct timeval start, end;

	status = parse_input( argc, argv, var, 
			              cov_base, phi_base, snp_base, pheno_base, out_base,
						  &nths, &thres, &xtile, &ytile, &xb, &yb, &write_output );
    if ( status == 0 ) 
        exit(EXIT_FAILURE);

	check_input_integrity( &cf, var,
			              cov_base, phi_base, snp_base, pheno_base, out_base,
						  nths, thres, xtile, ytile, xb, yb, write_output );

	initialize_config(
			&cf, 
			cov_base, phi_base, snp_base, pheno_base, out_base //,
		//    var, nths, xtile, ytile, xb, yb, write_output
	);

	/*if ( !strcmp( var, "chol" ) )*/
	/*{*/
	/*printf("Chol variant is broken. Please use -var eigen\n");*/
	/*exit(EXIT_FAILURE);*/
	/*}*/
	
	/*printf("\n**************************************************************************\n");*/
	/*printf("*** This is an alpha version under development. Use it at your own risk!!!\n");*/
	/*printf("**************************************************************************\n\n");*/
	print_info( &cf );

	// Compute h and sigma
	//
	// Perform GWAS
	if ( !strcmp( var, "eigen" ) )
	{
		printf( "Estimating GWAS parameters: heritability and variance... " );
		fflush( stdout );
		gettimeofday(&start, NULL);
		estimates_eig( &cf ); // Computation
    	gettimeofday(&end, NULL);
		printf( "Done ( took %.3f secs )\n", get_diff_ms(&start, &end)/1000. );
		fflush( stdout );

		printf( "Performing the study... " );
		fflush( stdout );
		gettimeofday(&start, NULL);
		fgls_eigen( &cf ); // Computation
    	gettimeofday(&end, NULL);
		printf( "Done ( took %.3f secs )\n", get_diff_ms(&start, &end)/1000. );
		fflush( stdout );
	}
	else if ( !strcmp( var, "chol" ) || !strcmp( var, "chol_gpu" ) )
	{
		printf( "Estimating GWAS parameters: heritability and variance... " );
		fflush( stdout );
		gettimeofday(&start, NULL);
		/*estimates_chol( &cf ); // Computation*/
		estimates_eig( &cf ); // Computation
    	gettimeofday(&end, NULL);
		printf( "Done ( took %.3f secs )\n", get_diff_ms(&start, &end)/1000. );
		fflush( stdout );

		printf( "Performing the study... " );
		fflush( stdout );
		gettimeofday(&start, NULL);
#ifdef WITH_GPU
		if ( !strcmp( var, "chol_gpu" ) )
			fgls_chol_gpu( cf );
		else
#endif
			fgls_chol( cf ); // Computation
    	gettimeofday(&end, NULL);
		printf( "Done ( took %.3f secs )\n", get_diff_ms(&start, &end)/1000. );
		fflush( stdout );
	}

	append_estimates( &cf );
	write_output_info_file( &cf );
	if ( !strcmp( var, "eigen" ) )
		cleanup( &cf );

	finalize_config( &cf );

    return 0;
}

void usage( void ) 
{
    fprintf(stderr, "\nUsage: CLAK-GWAS <arguments> [options]\n\n");
    fprintf(stderr, "Following arguments are mandatory:\n\n");
    fprintf(stderr, "  -cov <path>    base path to the file containing the covariates\n");
    fprintf(stderr, "  -phi <path>    base path to the file containing the relationship matrix\n");
    fprintf(stderr, "  -snp <path>    base path to the file containing the SNPs\n");
    fprintf(stderr, "  -pheno <path>  base path to the file containing the phenotypes\n");
    fprintf(stderr, "  -out <path>    base path to the file where the output will be stored\n\n");
    fprintf(stderr, "Following options might be given:\n\n");
#ifdef WITH_GPU
    fprintf(stderr, "  -var [chol | chol_gpu | eigen] Default is chol.\n");
#else
    fprintf(stderr, "  -var [chol | chol_gpu | eigen] Default is chol.\n");
#endif
    fprintf(stderr, "  -nths <num>         Default is 1 thread.\n");
    fprintf(stderr, "  -thres <num>        Default is 95(%%).\n");
	/*fprintf(stderr, "  -no-output               Default is to write the output files.\n");*/
    fprintf(stderr, "  -h                  Show this help and exit\n\n");
}

int parse_input( int argc, char *argv[], char *var, 
		         char *cov_base, char *phi_base, char *snp_base, char *pheno_base, char *out_base,
				 int *nths, int *thres, int *xtile, int *ytile, int *xb, int *yb, int *write_output )
{
	int arg = 1;
	int ok = 1;
	char *err = NULL;

	while (arg < argc && ok)
    {
      if ( !strcmp(argv[arg], "-h") )  
	  {
          usage();
          ok = 0;
      }
      else if ( ! strcmp(argv[arg], "-var") )  
	  {
          if ( ++arg < argc ) 
		  {
			  if ( !strcmp( argv[arg], "chol" ) )
			  {
				  strncpy( var, argv[arg], STR_BUFFER_SIZE );
				  /*fprintf( stderr, "Chol variant is currently broken\n" );*/
				  /*exit( EXIT_FAILURE );*/
			  }
#ifdef WITH_GPU
			  else if ( !strcmp( argv[arg], "chol_gpu" ) )
				  strncpy( var, argv[arg], STR_BUFFER_SIZE );
#endif
			  else if ( !strcmp( argv[arg], "eigen" ) )
				  strncpy( var, argv[arg], STR_BUFFER_SIZE );
			  else
			  {
                 fprintf( stderr, "Unrecognized variant: %s\n\n", argv[arg] );
				 ok = 0;
			  }
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-cov") )  
	  {
          if ( ++arg < argc ) 
		  {
			  /*dir = argv[arg];*/
		      strncpy( cov_base, argv[arg], STR_BUFFER_SIZE );
/*			  if ( ! ( stat( dir, &st ) == 0 && S_ISDIR( st.st_mode ) ) )
			  {
				  fprintf( stderr, "%s does not exist or it is not a directory\n\n", argv[arg] );
				  ok = 0;
			  }*/
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-phi") )  
	  {
          if ( ++arg < argc ) 
		  {
		      strncpy( phi_base, argv[arg], STR_BUFFER_SIZE );
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-snp") )  
	  {
          if ( ++arg < argc ) 
		  {
		      strncpy( snp_base, argv[arg], STR_BUFFER_SIZE );
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-pheno") )  
	  {
          if ( ++arg < argc ) 
		  {
		      strncpy( pheno_base, argv[arg], STR_BUFFER_SIZE );
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-out") )  
	  {
          if ( ++arg < argc ) 
		  {
		      strncpy( out_base, argv[arg], STR_BUFFER_SIZE );
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-nths") )  
	  {
          if ( ++arg < argc ) 
		  {
			  *nths = atoi(argv[arg]);
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-thres") )  
	  {
          if ( ++arg < argc ) 
		  {
			  *thres = atoi(argv[arg]);
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-xtile") )  
	  {
          if ( ++arg < argc ) 
		  {
			  *xtile = atoi(argv[arg]);
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-ytile") )  
	  {
          if ( ++arg < argc ) 
		  {
			  *ytile = atoi(argv[arg]);
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-xb") )  
	  {
          if ( ++arg < argc ) 
		  {
			  *xb = atoi(argv[arg]);
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-yb") )  
	  {
          if ( ++arg < argc ) 
		  {
			  *yb = atoi(argv[arg]);
			  arg++;
		  }
		  else
		  {
			  ok = 0;
			  err = argv[arg-1];
		  }
	  }
      else if ( ! strcmp(argv[arg], "-no-output") )  
	  {
		  *write_output = 0;
	       arg++;
	  }
      else 
	  {
          fprintf(stderr, "Unrecognized input: %s\n\n", argv[arg]);
		  ok = 0;
          arg++;
      }
  
      if ( err ) 
	  {
          fprintf(stderr, "Error: Missing value for argument '%s'\n\n", err);
          err = NULL;
      }
  }
    return ok;
}

void check_input_integrity( FGLS_config_t *cf, char *var,
				         char *cov_base, char *phi_base, char *snp_base, char *pheno_base, char *out_base,
						 int nths, int thres, int xtile, int ytile, int xb, int yb, int write_output )
{
	struct stat st;

	// Variant: chol or eigen
	strncpy( cf->var, var, STR_BUFFER_SIZE );
	
	// Paths provided?
	if ( strcmp( phi_base, "" ) == 0 || strcmp( cov_base, "" ) == 0 || 
			 strcmp( snp_base, "" ) == 0 || strcmp( pheno_base, "" ) == 0 ||
			(strcmp( out_base, "" ) == 0 && write_output) )
	{
		fprintf( stderr, "You must provide all mandatory arguments (see -h option for details)\n\n" );
		exit( EXIT_FAILURE );
	}
//	if ( strcmp( out_base, "" ) == 0 && write_output )
//	{
//		fprintf( stderr, "Output base path required (see -h option for details)\n\n" );
//		exit( EXIT_FAILURE );
//	}
			
	// Phi data
    snprintf( cf->Phi_data_path, STR_BUFFER_SIZE, "%s.fvd", phi_base );
	if ( ! ( stat( cf->Phi_data_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->Phi_data_path );
		exit( EXIT_FAILURE );
	}
	// Phi info
    snprintf( cf->Phi_info_path, STR_BUFFER_SIZE, "%s.fvi", phi_base );
	if ( ! ( stat( cf->Phi_info_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->Phi_info_path );
		exit( EXIT_FAILURE );
	}
	// Covariates (XL) data
    snprintf( cf->XL_data_path, STR_BUFFER_SIZE, "%s.fvd", cov_base );
	if ( ! ( stat( cf->XL_data_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->XL_data_path );
		exit( EXIT_FAILURE );
	}
	// Covariates (XL) info
    snprintf( cf->XL_info_path, STR_BUFFER_SIZE, "%s.fvi", cov_base );
	if ( ! ( stat( cf->XL_info_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->XL_info_path );
		exit( EXIT_FAILURE );
	}
	// SNPs (XR) data
    snprintf( cf->XR_data_path, STR_BUFFER_SIZE, "%s.fvd", snp_base );
	if ( ! ( stat( cf->XR_data_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->XR_data_path );
		exit( EXIT_FAILURE );
	}
	// SNPs (XR) info
    snprintf( cf->XR_info_path, STR_BUFFER_SIZE, "%s.fvi", snp_base );
	if ( ! ( stat( cf->XR_info_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->XR_info_path );
		exit( EXIT_FAILURE );
	}
	// Phenotypes (Y) data
    snprintf( cf->Y_data_path, STR_BUFFER_SIZE, "%s.fvd", pheno_base );
	if ( ! ( stat( cf->Y_data_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->Y_data_path );
		exit( EXIT_FAILURE );
	}
	// SNPs (XR) info
    snprintf( cf->Y_info_path, STR_BUFFER_SIZE, "%s.fvi", pheno_base );
	if ( ! ( stat( cf->Y_info_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
	{
		fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->Y_info_path );
		exit( EXIT_FAILURE );
	}
	// Output (B) data
	if ( write_output )
	{
		snprintf( cf->B_data_path, STR_BUFFER_SIZE, "%s.out",  out_base );
		snprintf( cf->B_info_path, STR_BUFFER_SIZE, "%s.iout", out_base );
		// What to do if exists?
/*		if ( ! ( stat( cf->Y_data_path, &st ) == 0 && S_ISREG( st.st_mode ) ) )
		{
			fprintf( stderr, "%s does not exist or it is not a regular file\n\n", cf->Y_data_path );
			exit( EXIT_FAILURE );
		}
*/	}
	else
	{
		sprintf( cf->B_data_path, "/dev/null" );
		sprintf( cf->B_info_path, "/dev/null" );
	}

	// Threads
	cf->num_threads = nths;
	if ( cf->num_threads < 1 )
	{
		fprintf( stderr, "Number of threads (%d) must be >= 1\n\n", cf->num_threads );
		exit( EXIT_FAILURE );
	}

	// Threshold
	cf->threshold = thres;
	if ( cf->threshold < 1 || cf->threshold > 100 )
	{
		fprintf( stderr, "Threshold (%d) must be an integer value between 1 and 100\n\n", cf->threshold );
		exit( EXIT_FAILURE );
	}

	// Tiling and Blocking (hidden options for performance testing)
	cf->x_b = xb;
	cf->y_b = yb;
	cf->x_tile = xtile;
	cf->y_tile = ytile;
}

void print_info( FGLS_config_t *cf )
{
	printf( "\nRunning a Genome-Wide Association Study of the following size:\n" );
	printf( "  sample size:      %zu\n", cf->n);
	printf( "  # of covariates:  %zu\n", cf->p-2);
	printf( "  # of SNPs:        %zu\n", cf->m);
	printf( "  # of phenotypes:  %zu\n", cf->t);
	fflush( stdout );

	if ( !strcmp( cf->var, "chol" ) )
		printf( "\nWill use CLAK-Chol with the following parameters:\n" );
#ifdef WITH_GPU
	else if ( !strcmp( cf->var, "chol_gpu" ) )
		printf( "\nWill use CLAK-Chol on GPU(s) with the following parameters:\n" );
#endif
	else
		printf( "\nWill use CLAK-Eig with the following parameters:\n" );
	printf( "  x_b: %zu\n", cf->x_b );
	printf( "  y_b: %zu\n", cf->y_b );
	printf( "  o_b: %zu\n", cf->ooc_b );
	printf( "  x_tile: %zu\n", cf->x_tile );
	printf( "  y_tile: %zu\n", cf->y_tile );
	printf( "  # of threads: %d\n",  cf->num_threads );
	printf( "  Available memory %f GBs (out of %f GBs)\n\n", ((double)cf->availMem) / (1L<<30), ((double)cf->totalMem) / (1L<<30) );
	fflush( stdout );
}

void write_output_info_file( FGLS_config_t *cf )
{
	FILE *f;
	int type = 6; //DOUBLE;
	int nbytes = sizeof(double); // bytes per record
	/*int size_one_output_record = cf->p + cf->p*(cf->p+1)/2;*/
	int namelength = 32;
	char buff[STR_BUFFER_SIZE];
	int i, j;

	f = fgls_fopen( cf->B_info_path, "wb" );

	fwrite( &type, sizeof(int), 1, f );
	fwrite( &nbytes, sizeof(int), 1, f );
	fwrite( &cf->p, sizeof(int), 1, f );
	/*fwrite( &size_one_output_record, sizeof(int), 1, f );*/
	fwrite( &cf->m, sizeof(int), 1, f );
	fwrite( &cf->t, sizeof(int), 1, f );
	fwrite( &cf->x_b, sizeof(int), 1, f );
	fwrite( &cf->y_b, sizeof(int), 1, f );
	fwrite( &namelength, sizeof(int), 1, f );
	// Output labels
	// beta
	for ( i = 0; i < cf->wXL; i++ )
	{
		size_t cov_pos = (cf->XL_fvi->fvi_header.numObservations + i)*namelength;
		snprintf( buff, STR_BUFFER_SIZE, "beta_%s", &cf->XL_fvi->fvi_data[cov_pos] );
		fwrite( buff, namelength, 1, f );
	}
	snprintf( buff, STR_BUFFER_SIZE, "beta_SNP" );
	fwrite( buff, namelength, 1, f );
	// se
	for ( i = 0; i < cf->wXL; i++ )
	{
		size_t cov_pos = (cf->XL_fvi->fvi_header.numObservations + i)*namelength;
		snprintf( buff, STR_BUFFER_SIZE, "se_%s", &cf->XL_fvi->fvi_data[cov_pos] );
		fwrite( buff, namelength, 1, f );
	}
	snprintf( buff, STR_BUFFER_SIZE, "se_SNP" );
	fwrite( buff, namelength, 1, f );

	// cov
	for ( j = 0; j < cf->wXL; j++ )
	{
		size_t cov_pos_j = (cf->XL_fvi->fvi_header.numObservations + j)*namelength;
		for ( i = j+1; i < cf->wXL; i++ )
		{
			size_t cov_pos_i = (cf->XL_fvi->fvi_header.numObservations + i)*namelength;
			snprintf( buff, STR_BUFFER_SIZE, "cov_%s_%s", 
					&cf->XL_fvi->fvi_data[cov_pos_i],
					&cf->XL_fvi->fvi_data[cov_pos_j] );
			fwrite( buff, namelength, 1, f );
		}
		snprintf( buff, STR_BUFFER_SIZE, "cov_SNP_%s", 
				&cf->XL_fvi->fvi_data[cov_pos_j] );
		fwrite( buff, namelength, 1, f );
	}

	size_t vars_pos = cf->XR_fvi->fvi_header.numObservations*namelength;
	fwrite( &cf->XR_fvi->fvi_data[vars_pos], namelength, cf->m, f );
	vars_pos = cf->Y_fvi->fvi_header.numObservations*namelength;
	fwrite( &cf->Y_fvi->fvi_data[vars_pos], namelength, cf->t, f );

	fclose( f );
}

void append_estimates( FGLS_config_t *cf )
{
	int size_of_one_record = cf->p + cf->p*(cf->p+1)/2;
	sync_write( cf->ests, cf->B, (3+cf->wXL)*cf->t, size_of_one_record * cf->m * cf->t );
}

void cleanup( FGLS_config_t *cf )
{
	if ( remove( cf->ZtXR_path ) == -1 )
		perror("[ERROR] Couldn't delete temporary file for XR");
	if ( remove( cf->ZtY_path ) == -1 )
		perror("[ERROR] Couldn't delete temporary file for Y");
}
