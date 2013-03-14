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
 *   Lucas Beyer (lucas-b.eyer.be)
 */

#ifdef WITH_GPU
// Do nothing if there is no GPU support.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <sys/time.h>
#include <time.h>

#include <aio.h>

#include <cuda_runtime_api.h>
#include <cublas_v2.h>

#include <omp.h>

#include "blas.h"
#include "lapack.h"
#include "options.h"
#include "GWAS.h"
#include "wrappers.h"
#include "timing.h"
#include "double_buffering.h"
#include "utils.h"
#include "fgls_chol_gpu.h"

#if VAMPIR
    #include "vt_user.h"
#endif

#include <stdarg.h>
static void sync_gpus(int ngpus);
static void start_section(struct timeval* t_start, int ngpus, const char* vt_id, const char* text, ...);
static void end_section(struct timeval* t_start, int ngpus, const char* vt_id);
#if 1
#define START_SECTION(VT_ID, MSG) start_section(&t_start, ngpus, VT_ID, MSG "... ")
#define START_SECTION2(VT_ID, MSG, ...) start_section(&t_start, ngpus, VT_ID, MSG "... ", __VA_ARGS__)
#define END_SECTION(VT_ID) end_section(&t_start, ngpus, VT_ID)
#else
#define START_SECTION(VT_ID, MSG) (void*)0
#define START_SECTION2(VT_ID, MSG, ...) (void*)0
#define END_SECTION(VT_ID) (void*)0
#endif

/*
 * Builds Phi as an SPD matrix, after the eigenvalues were fixed 
 * during the REML estimation
 */
void build_SPD_Phi( int n, double *eigVecs, double *eigVals, double *Phi );

/*
 * Cholesky-based GPU solution of the 
 *  sequence of Feasible Generalized Least-Squares problem
 *  in the context of GWAS:
 */
int fgls_chol_gpu( FGLS_config_t cf )
{
	int n = cf.n,
		   m = cf.m,
		   p = cf.p,
		   t = cf.t,
		   x_b = cf.x_b,
		   /*y_b = cf.y_b,*/
		   wXL = cf.wXL,
		   wXR = cf.wXR;
    /* In-core operands */
    double *Phi;
    double *M;
    double *ests;
    double *h2;
	double *res_sigma;
    double alpha;
    double beta;

    /* Out-of-core operands */
    double *Bij; // Auxiliary variables

    /* Reusable data thanks to constant XL */
    double *XL;
    double *XL_orig; // XL and a copy (XL is overwritten at every iteration of j)
    double *B_t;  // Top part of b ( in inv(S) b )
    double *V_tl; // Top-Left part of V

    /* BLAS / LAPACK constants */
    double ZERO = 0.0;
    double ONE = 1.0;
    int iONE = 1;
    /* LAPACK error value */
    int info;

    /* iterators and auxiliar vars */
    int i, j, k, l; // size_t
    int nn = cf.n * cf.n; // size_t
	size_t size_one_b_record = p + (p*(p+1))/2;

	// Threading
	int id;
	double *tmpBs, *tmpVs; // Buffer with one B and one V per thread
	double *oneB, *oneV;   // Each thread pointer to its B and V

    // For measuring times.
    struct timeval t_start;

    if ( cf.y_b != 1 )
	{
        fprintf(stderr, "\n[Warning] y_b not used (set to 1)\n");
		cf.y_b = 1;
	}

    if ( cf.t > 1 ) {
        // Note/TODO: to make it work, you need to fix the "fetching next Xr
        // from disk" code part to restart from offset zero at some point.
        // That's why the loop over t is still left in.
        fprintf(stderr, "\n[ERROR] The chol_gpu variant doesn't support multiple phenotypes. Use the chol or eigen variants.\n");
        exit(EXIT_FAILURE);
    }

    ///////
    // GPU: Initializing the GPU(s)
    int ngpus = 0, igpu = 0;
    cublasHandle_t cu_handle;
    cublasStatus_t cu_status;
    if((cu_status = cublasCreate(&cu_handle)) != CUBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "\n[ERROR] cublasCreate() failed (info: %d)\n", cu_status);
        exit(EXIT_FAILURE);
    }

    cudaError_t cu_error;
    if((cu_error = cudaGetDeviceCount(&ngpus)) != cudaSuccess) {
        fprintf(stderr, "\n[ERROR] Can't get the cuda device count. Are there any? (info: %d)", cu_error);
        exit(EXIT_FAILURE);
    }

    //fprintf(stdout, "\n[Info] Using %d GPUs\n", ngpus);

    // Create two streams for each GPU: one computation stream and one data
    // transfer stream, so those two can work in parallel.
    cudaStream_t* cu_trans_streams = fgls_malloc(ngpus*sizeof(cudaStream_t));
    cudaStream_t* cu_comp_streams = fgls_malloc(ngpus*sizeof(cudaStream_t));

    for(igpu = 0 ; igpu < ngpus ; ++igpu) {
        cudaSetDevice(igpu);
        if((cu_error = cudaStreamCreate(cu_trans_streams + igpu)) != cudaSuccess) {
            fprintf(stderr, "\n[ERROR] Something is wrong with your GPUs: couldn't create a transfer stream on GPU %d (info: %d)\n", igpu, cu_error);
            exit(EXIT_FAILURE);
        }
        if((cu_error = cudaStreamCreate(cu_comp_streams + igpu)) != cudaSuccess) {
            fprintf(stderr, "\n[ERROR] Something is wrong with your GPUs: couldn't create a computation stream on GPU %d (info: %d)\n", igpu, cu_error);
            exit(EXIT_FAILURE);
        }
    }

    // We can already allocate GPU space for L here.
    double** L_gpus = fgls_malloc(ngpus*sizeof(void*));
    size_t L_gpu_bytes = (size_t)cf.n*cf.n * sizeof(double);

    // Aswell as for the streamed Xrs.
    double** Xr_alpha_gpus = fgls_malloc(ngpus*sizeof(double*));
    double** Xr_beta_gpus = fgls_malloc(ngpus*sizeof(double*));
    size_t Xr_gpu_bytes = (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double);

    for(igpu = 0 ; igpu < ngpus ; igpu++) {
        cudaSetDevice(igpu);
        if((cu_error = cudaMalloc((void**)&L_gpus[igpu], L_gpu_bytes)) != cudaSuccess) {
            fprintf(stderr, "\n[ERROR] Not enough memory to allocate %ld bytes for L on GPU %d (info: %d)\n", L_gpu_bytes, igpu, cu_error);
            exit(EXIT_FAILURE);
        }
        if((cu_error = cudaMalloc((void**)&Xr_alpha_gpus[igpu], Xr_gpu_bytes)) != cudaSuccess) {
            fprintf(stderr, "\n[ERROR] Not enough memory to allocate %ld bytes for Xr on GPU %d (info: %d)\n", Xr_gpu_bytes, igpu, cu_error);
            exit(EXIT_FAILURE);
        }
        if((cu_error = cudaMalloc((void**)&Xr_beta_gpus[igpu], Xr_gpu_bytes)) != cudaSuccess) {
            fprintf(stderr, "\n[ERROR] Not enough memory to allocate %ld bytes for Xr on GPU %d (info: %d)\n", Xr_gpu_bytes, igpu, cu_error);
            exit(EXIT_FAILURE);
        }
    }
    // /GPU
    ///////

    /* Memory allocation */
    // In-core
	build_SPD_Phi( cf.n, cf.Z, cf.W, cf.Phi );
	Phi   = cf.Phi;
    M     = ( double * ) fgls_malloc ( (size_t)cf.n * cf.n * sizeof(double) );
    ests  = cf.ests;

	h2 = ests;
	res_sigma = &ests[2*cf.t];

    XL_orig = cf.XL;
    XL      = ( double * ) fgls_malloc ( cf.wXL * cf.n * sizeof(double) );
    B_t  = ( double * ) fgls_malloc ( cf.wXL * sizeof(double) );
    V_tl = ( double * ) fgls_malloc ( cf.wXL * cf.wXL * sizeof(double) );

	// Temporary storage prior to copying in db_B
    tmpBs = ( double * ) fgls_malloc ( cf.p * cf.num_threads * sizeof(double) );
    tmpVs = ( double * ) fgls_malloc ( cf.p * cf.p * cf.num_threads * sizeof(double) );

    /* Files and pointers for out-of-core */
    double *XR_comp, *Y_comp, *B_comp;

    /* Asynchronous IO data structures */
	double_buffering db_Y, db_B;
	double_buffering_init( &db_Y, (size_t)cf.n * cf.y_b * sizeof(double),
			                cf.Y,  &cf );
	double_buffering_init( &db_B, (size_t)size_one_b_record * cf.x_b * cf.y_b * sizeof(double),
			                cf.B,  &cf );

    ///////
    // GPU
    // For the triple buffering:
    double* Xr[3];
    size_t A = 0, B = 1, C = 2; // Indices for buffer rotation.
    cudaError_t cu_error1 = cudaHostAlloc((void**)&Xr[A], (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cudaHostAllocPortable);
    cudaError_t cu_error2 = cudaHostAlloc((void**)&Xr[B], (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cudaHostAllocPortable);
    cudaError_t cu_error3 = cudaHostAlloc((void**)&Xr[C], (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cudaHostAllocPortable);
    if(cu_error != cudaSuccess || cu_error2 != cudaSuccess || cu_error3 != cudaSuccess) {
        fprintf(stderr, "\n[ERROR] Not enough memory to allocate page-locked triple buffers (3*%ld MB, info: %d)\n", (size_t)cf.x_b * cf.wXR * cf.n * sizeof(double), cu_error | cu_error2 | cu_error3);
        exit(EXIT_FAILURE);
    }
    struct aiocb aio_Xr; // Always only have one async read (in A).
    // /GPU
    ///////

// #if VAMPIR
//     VT_USER_START("READ_X");
// #endif
//     /* Read first block of XR's */
// 	double_buffering_read_XR( &db_XR, IO_BUFF, 0, (size_t)MIN( cf.x_b, cf.m ) - 1 );
// 	double_buffering_swap( &db_XR );
// #if VAMPIR
//     VT_USER_END("READ_X");
// #endif
#if VAMPIR
    VT_USER_START("READ_Y");
#endif
    /* Read first Y */
	double_buffering_read_Y( &db_Y, IO_BUFF, 0, 0 );
	double_buffering_swap( &db_Y );
#if VAMPIR
    VT_USER_END("READ_Y");
#endif

    int iter = 0;
    for ( j = 0; j < t; j++ )
    {
        /* Set the number of threads for the multi-threaded BLAS */
		set_multi_threaded_BLAS( cf.num_threads );

#if VAMPIR
        VT_USER_START("READ_Y");
#endif
        /* Read next Y */
		size_t next_j = (j+1) >= t ? 0 : j+1;
		double_buffering_read_Y( &db_Y, IO_BUFF, next_j, next_j );
#if VAMPIR
        VT_USER_END("READ_Y");
#endif

#if VAMPIR
        VT_USER_START("COMP_J");
#endif
        /* M := sigma * ( h^2 Phi - (1 - h^2) I ) */
        memcpy( M, Phi, (size_t)n * n * sizeof(double) );
		alpha = res_sigma[j] * h2[j];
        beta  = res_sigma[j] * (1 - h2[j]);
        dscal_(&nn, &alpha, M, &iONE);
        for ( i = 0; i < n; i++ )
            M[i*n + i] = M[i*n + i] + beta;

        /* L * L' = M */
        dpotrf_(LOWER, &n, M, &n, &info);
        if (info != 0)
        {
            char err[STR_BUFFER_SIZE];
            snprintf(err, STR_BUFFER_SIZE, "dpotrf(M) failed (info: %d)", info);
            error_msg(err, 1);
        }

        /* XL := inv(L) * XL */
        memcpy( XL, XL_orig, wXL * n * sizeof(double) );
        dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &wXL, &ONE, M, &n, XL, &n);

#if VAMPIR
        VT_USER_START("WAIT_Y");
#endif
        // Wait until current Y is available for computation
		double_buffering_wait( &db_Y, COMP_BUFF );
#if VAMPIR
        VT_USER_END("WAIT_Y");
#endif

        /* y := inv(L) * y */
		Y_comp = double_buffering_get_comp_buffer( &db_Y );
		// Sanity check
		average( Y_comp, n, 1, cf.threshold, "TRAIT", 
				&cf.Y_fvi->fvi_data[n*NAMELENGTH], NAMELENGTH, 0 );
        dtrsv_(LOWER, NO_TRANS, NON_UNIT, &n, M, &n, Y_comp, &iONE);

        /* B_t := XL' * y */
        dgemv_(TRANS, &n, &wXL, &ONE, XL, &n, Y_comp, &iONE, &ZERO, B_t, &iONE);

        /* V_tl := XL' * XL */
        dsyrk_(LOWER, TRANS, &wXL, &n, &ONE, XL, &n, &ZERO, V_tl, &wXL);
#if VAMPIR
        VT_USER_END("COMP_J");
#endif
		/* Solve for x_b X's at once */

        // Please refer to Lucas Beyer's paper and/or thesis to understand the loop.
        // The dependency colors refer to the "timeline-perspective" figure.
        // Basically, the idea is the following (read from top to bottom, then left to right):
        //
        //       b = -1             b = 0             b = 1         b = 2..last-2       b = last-1          b = last         b = last+1
        //
        //        ---                ---            wait beta         wait alpha        wait beta          wait alpha           ---
        //        ---                ---               ---            wait beta         wait alpha         wait beta         wait alpha
        //        ---                ---            trsm beta         trsm alpha        trsm beta          trsm alpha           ---
        // read X[b+2] -> A   read X[b+2] -> B  read X[b+2] -> C  read X[b+2] -> A         ---                ---               ---
        //        ---                ---               ---        recv_s beta -> B  recv_s alpha -> C  recv_s beta -> A  recv_s alpha -> B
        //        ---         wait X[b+1] -> A  wait X[b+1] -> B  wait X[b+1] -> C  wait X[b+1] -> A          ---               ---
        //        ---          send A -> beta    send B -> alpha   send C -> beta    send A -> alpha          ---               ---
        //        ---                ---               ---         S-loop B -> B^    S-loop C -> C^     S-loop A -> A^    S-loop B -> B^
        //        ---                ---               ---             write B^          write C^           write A^          write B^
        //
        // So we just need the following rotations at the end of each iteration:
        //    A->B->C->A
        //   alpha<->beta
        int blockcount = m % x_b == 0 ? m/x_b : m/x_b+1;
        int iblock = 0;
        for (iblock = -1 ; iblock <= blockcount+1 ; ++iblock)
        {
            // cu_trsm_wait alpha
            // Wait for the previous GPU computation to be done...
            // (The first GPU computation happens at i == 1 thus we first wait at i == 2)
            // (bordeaux dependency, also implied by single lifeline of GPU.)
            if(2 <= iblock) {
                START_SECTION2("GPU_trsm", "%d: Waiting for the trsms (%d) on the GPUs to be done", iblock, iblock-1);
//                 for(igpu = 0 ; igpu < ngpus ; ++igpu) {
//                     // cudaSetDevice(igpu); // Unnecessary according to "cuda_webinar_multi_gpu.pdf", p. 6
//                     cudaStreamSynchronize(cu_comp_streams[igpu]);
//                 }
                END_SECTION("GPU_trsm");
            }

            // cu_send_wait C -> beta
            // wait for the sending of the previous data-block to the GPU to be done.
            // (olive dependency)
            if(1 <= iblock && iblock <= blockcount) {
                START_SECTION2("GPU_send_Xr", "%d: Waiting for the sending of Xr (%d) parts to the GPUs to be done", iblock, iblock);
//                 for(igpu = 0 ; igpu < ngpus ; ++igpu) {
//                     cudaStreamSynchronize(cu_trans_streams[igpu]);
//                 }
                END_SECTION("GPU_send_Xr");
            }

            // alpha <- cu-trsm_async L_gpu, alpha
            // dispatch the TRSM on the GPU for the block which was just sent there.
            if(1 <= iblock && iblock <= blockcount) {
                START_SECTION2("GPU_trsm", "%d: Starting the trsms (%d) on the GPUs", iblock, iblock);
                int curr_Xr_block_length = MIN(x_b, m - (iblock-1)*x_b);
                int rhss  = wXR * curr_Xr_block_length;

                // TODO(lucasb): sanity check! (call to average, replaces NaNs by avg.)

                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
//                     cudaSetDevice(igpu);
//                     cublasSetStream(cu_handle, cu_comp_streams[igpu]);
//                     if((cu_status = cublasDtrsm(cu_handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, rhss/ngpus, &ONE, L_gpus[igpu], n, Xr_beta_gpus[igpu], n)) != CUBLAS_STATUS_SUCCESS) {
//                         fprintf(stderr, "\n[ERROR] Computing the trsm on the GPU %d failed (info: %d)\n", igpu, cu_status);
//                         exit(EXIT_FAILURE);
//                     }
                }
                END_SECTION("GPU_trsm");
            }

            // aio_read Xr[b+2] -> A
            // Read the second-next block from HDD to main memory.
            if(iblock <= blockcount-2) {
                START_SECTION2("READ_X", "%d: Starting the read of the Xr block (%d)", iblock, iblock+2);
//                 fgls_aio_read( &aiocb_x,
//                                fileno( XR_fp ), X[hdd_cpu_buff],
//                                MIN((size_t)x_b, m - ((size_t)x_b*iblock)) * wXR * n * sizeof(double),
//                                ((off_t)x_b*iblock) * wXR * n * sizeof(double) );
                END_SECTION("READ_X");
            }

            // cu_recv B <- beta
            // synchronuously get the results of the previous TRSM, sync since
            // we need them for further computation anyway.
            if(2 <= iblock) {
                START_SECTION2("GPU_recv_LXr", "%d: Getting inv(L) * Xr (%d) back to the CPU", iblock, iblock-1);
                //size_t prev_Xr_elems_per_device = prev_Xr_block_length*wXR*n/ngpus;
                for(igpu = 0 ; igpu < ngpus ; ++igpu) {
                    // TODO(lucasb): sync. Async makes almost no sense and isn't what the paper says.
//                     cudaSetDevice(igpu);
//                     if((cu_status = cublasGetVectorAsync(prev_Xr_elems_per_device, sizeof(double), Xr_alpha_gpus[igpu], 1, X[cpu_gpu_buff] + igpu*prev_Xr_elems_per_device, 1, cu_trans_streams[igpu])) != CUBLAS_STATUS_SUCCESS) {
//                         fprintf(stderr, "\n[ERROR]sending GPU %d's part of inv(L)*Xr to the CPU failed (info: %d)\n", igpu, cu_status);
//                         exit(EXIT_FAILURE);
//                     }
                }
                END_SECTION("GPU_recv_LXr");
            }

            // aio_wait Xr[b+1] -> C
            // cu_send_async C -> beta
            // waits for the next X-block to be loaded and then sends it to
            // the GPU so that it can be TRSMed in the next iteration.
            // (blue dependency)
            if(0 <= iblock && iblock <= blockcount-1) {
                START_SECTION2("WAIT_X", "%d: Waiting for the Xr block's (%d) arrival in main memory with departure from HDD. Weight is %.2f MB", iblock, iblock+1, 0.0/0.0);//, aiocb_x.aio_nbytes/1024.0/1024.0);
//                 fgls_aio_suspend( &aiocb_x, 1, NULL );
                END_SECTION("WAIT_X");

                // acu_send_async C -> beta
                // send the next X-block we just waited for to the GPU.
                START_SECTION2("GPU_send_Xr", "%d: Distributing the Xr block (%d) equally amongst the %d GPUs", iblock, iblock+1, ngpus);
//                 size_t Xr_elems_per_device = n*(wXR*x_b)/ngpus;
//                 for(igpu = 0 ; igpu < ngpus ; ++igpu) {
//                     cudaSetDevice(igpu);
//                     cublasSetStream(cu_handle, cu_trans_streams[igpu]);
//                     if((cu_status = cublasSetVectorAsync(Xr_elems_per_device, sizeof(double), X[hdd_cpu_buff]+igpu*Xr_elems_per_device, 1, Xr_alpha_gpus[igpu], 1, cu_trans_streams[igpu])) != CUBLAS_STATUS_SUCCESS) {
//                         char err[STR_BUFFER_SIZE];
//                         snprintf(err, STR_BUFFER_SIZE, "sending part of Xr to the GPU %d failed (info: %d)", igpu, cu_status);
//                         error_msg(err, 1);
//                     }
//                 }
                END_SECTION("GPU_send_Xr");
            }

// #if VAMPIR
//             VT_USER_START("READ_X");
// #endif
//             /* Read next block of XR's */
// 			size_t next_x_from = ((size_t)ib + x_b) >= m ?  0 : (size_t)ib + x_b;
// 			size_t next_x_to   = ((size_t)ib + x_b) >= m ? MIN( (size_t)x_b, (size_t)m ) - 1 : 
// 				                                           next_x_from + MIN( (size_t)x_b, (size_t)m - next_x_from ) - 1;
// 			double_buffering_read_XR( &db_XR, IO_BUFF, next_x_from, next_x_to );
// #if VAMPIR
//             VT_USER_END("READ_X");
// #endif

// #if VAMPIR
//             VT_USER_START("WAIT_X");
// #endif
//             /* Wait until current block of XR's is available for computation */
// 			double_buffering_wait( &db_XR, COMP_BUFF );
// #if VAMPIR
//             VT_USER_END("WAIT_X");
// #endif

#if 0
            /* Set the number of threads for the multi-threaded BLAS */
			set_multi_threaded_BLAS( cf.num_threads );

#if VAMPIR
            VT_USER_START("COMP_IB");
#endif
            /* XR := inv(L) XR */
			XR_comp = double_buffering_get_comp_buffer( &db_XR );
			// Auxiliar variables
            int x_inc = MIN(x_b, m - ib);
            int rhss  = wXR * x_inc;
			// Sanity check
			average( XR_comp, n, x_inc, cf.threshold, "SNP", 
					&cf.XR_fvi->fvi_data[(n+ib)*NAMELENGTH], NAMELENGTH, 1 );
			// Computation
            dtrsm_(LEFT, LOWER, NO_TRANS, NON_UNIT, &n, &rhss, &ONE, M, &n, XR_comp, &n);

#if VAMPIR
            VT_USER_END("COMP_IB");
#endif

#endif // 0

            // S-loop B -> B^
            // aio_wait r[b-2]
            // aio_write r[b-2]
            // Perform the "remaining" computations using the previous TRSM's
            // result on the CPU and then schedule their writing to HDD.
            if(2 <= iblock) {
                START_SECTION2("COMP_I", "%d: Computing the S-loop of block %d", iblock, iblock-1);
//                sloop(out_r, Xr_block, Xl, y, Stl);
                END_SECTION("COMP_I");
            }

#if 0

#if CHOL_MIX_PARALLELISM
            /* Set the number of threads for the multi-threaded BLAS to 1.
             * The innermost loop is parallelized using OPENMP */
			set_single_threaded_BLAS();
#endif
#if VAMPIR
            VT_USER_START("COMP_I");
#endif
            B_comp = double_buffering_get_comp_buffer( &db_B );
#if CHOL_MIX_PARALLELISM
            #pragma omp parallel for private(Bij, oneB, oneV, i, k, info, id) schedule(static) num_threads(cf.num_threads)
#endif
            for (i = 0; i < x_inc; i++)
            {
				id = omp_get_thread_num();
				oneB = &tmpBs[ id * p ];
				oneV = &tmpVs[ id * p * p ];
                Bij = &B_comp[i * size_one_b_record];

                // Building B
                // Copy B_T
                memcpy(oneB, B_t, wXL * sizeof(double));
                // B_B := XR' * y
                dgemv_("T", 
                        &n, &wXR, 
                        &ONE, &XR_comp[i * wXR * n], &n, Y_comp, &iONE, 
                        &ZERO, &oneB[wXL], &iONE);

                // Building V
                // Copy V_TL
                for( k = 0; k < wXL; k++ )
                    dcopy_(&wXL, &V_tl[k*wXL], &iONE, &oneV[k*p], &iONE); // V_TL
                // V_BL := XR' * XL
                dgemm_("T", "N",
                        &wXR, &wXL, &n,
                        &ONE, &XR_comp[i * wXR * n], &n, XL, &n,
                        &ZERO, &oneV[wXL], &p); // V_BL
                // V_BR := XR' * XR
                dsyrk_("L", "T", 
                        &wXR, &n, 
                        &ONE, &XR_comp[i * wXR * n], &n, 
                        &ZERO, &oneV[wXL * p + wXL], &p); // V_BR

                // B := inv(V) * B
                dpotrf_(LOWER, &p, oneV, &p, &info);
                if (info != 0)
                {
					for ( k = 0; k < size_one_b_record; k++ )
						Bij[k] = 0.0/0.0; //nan("char-sequence");
					continue;
                }
                dtrsv_(LOWER, NO_TRANS, NON_UNIT, &p, oneV, &p, oneB, &iONE);
                dtrsv_(LOWER,    TRANS, NON_UNIT, &p, oneV, &p, oneB, &iONE);

                /* V := res_sigma * inv( X' inv(M) X) */
                dpotri_(LOWER, &p, oneV, &p, &info);
                if (info != 0)
                {
                    char err[STR_BUFFER_SIZE];
                    snprintf(err, STR_BUFFER_SIZE, "dpotri failed (info: %d)", info);
                    error_msg(err, 1);
                }

				// Copy output
				for ( k = 0; k < p; k++ )
					Bij[k] = (float) oneB[k];
                for ( k = 0; k < p; k++ )
                    Bij[p+k] = (float)sqrt(oneV[k*p+k]);
				int idx = 0;
				for ( k = 0; k < p-1; k++ ) // Cols of V
					for ( l = k+1; l < p; l++ ) // Rows of V
					{
						Bij[p+p+idx] = (float)oneV[k*p+l];
						idx++;
					}
#if 0
			  printf("Chi square: %.6f\n", ( (oneB[p-1] / Bij[p+p-1]) * (oneB[p-1] / Bij[p+p-1]) ) );
#endif
            }
#if VAMPIR
            VT_USER_END("COMP_I");
#endif

#if VAMPIR
            VT_USER_START("WAIT_BV");
#endif
            /* Wait until the previous blocks of B's and V's are written */
            if ( iter > 0)
                double_buffering_wait( &db_B, IO_BUFF );
#if VAMPIR
            VT_USER_END("WAIT_BV");
#endif
            /* Write current blocks of B's and V's */
#if VAMPIR
            VT_USER_START("WRITE_BV");
#endif
			double_buffering_write_B( &db_B, COMP_BUFF, ib, ib+x_inc - 1, j, j );
#if VAMPIR
            VT_USER_END("WRITE_BV");
#endif

#endif

            /* Swap buffers */
//			double_buffering_swap( &db_XR );
			double_buffering_swap( &db_B  );
            iter++; // TODO: don't need.
        }
        /* Swap buffers */
		double_buffering_swap( &db_Y );
    }

#if VAMPIR
    VT_USER_START("WAIT_ALL");
#endif
    /* Wait for the remaining IO operations issued */
// 	double_buffering_wait( &db_XR, COMP_BUFF );
	double_buffering_wait( &db_Y,  COMP_BUFF );
	double_buffering_wait( &db_B,  IO_BUFF );
#if VAMPIR
    VT_USER_END("WAIT_ALL");
#endif

    // GPU: clean-up
    for(igpu = 0 ; igpu < ngpus ; ++igpu) {
        cudaSetDevice(igpu);
        cudaFree(L_gpus[igpu]);
        cudaFree(Xr_alpha_gpus[igpu]);
        cudaFree(Xr_beta_gpus[igpu]);
    }
    cudaFreeHost(Xr[A]);
    cudaFreeHost(Xr[B]);
    cudaFreeHost(Xr[C]);
    cublasDestroy(cu_handle);

    free(L_gpus);
    free(Xr_alpha_gpus);
    free(Xr_beta_gpus);

    /* Clean-up */
    free( M );

    free( XL );
    free( B_t  );
    free( V_tl );
    free( tmpBs );
    free( tmpVs );

// 	double_buffering_destroy( &db_XR );
	double_buffering_destroy( &db_Y  );
	double_buffering_destroy( &db_B  );

    return 0;
}

static void sync_gpus(int ngpus)
{
    int igpu = 0;
    for(igpu = 0 ; igpu < ngpus ; ++igpu) {
        cudaSetDevice(igpu);
        cudaStreamSynchronize(0);
    }
}

static void start_section(struct timeval* t_start,
                          int ngpus,
                          const char* vt_id,
                          const char* text, ...)
{
#ifdef FGLS_GPU_SERIAL
    sync_gpus(ngpus);
#endif
#ifdef TIMING
    read_clock(t_start);
#endif
#if defined(DEBUG) || defined(TIMING)
    va_list argp;
    va_start(argp, text);
    vprintf(text, argp);
    fflush(stdout);
#endif

#ifdef VTRACE
    VT_USER_START(vt_id);
#endif
}

static void end_section(struct timeval* t_start,
                        int ngpus,
                        const char* vt_id)
{
#ifdef FGLS_GPU_SERIAL
    sync_gpus(ngpus);
#endif
#ifdef VTRACE
    VT_USER_END(vt_id);
#endif
#ifdef TIMING
    struct timeval t_end;
    read_clock(&t_end);
    int dt = elapsed_time(t_start, &t_end);
    printf("done in %fs (%fms, %dus)\n", dt*0.000001, dt*0.001, dt);
    fflush(stdout);
#elif defined (DEBUG)
    printf("done\n");
    fflush(stdout);
#endif
}

#endif // WITH_GPU

