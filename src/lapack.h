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

#ifndef __LAPACK_H
#define __LAPACK_H

/* Cholesky */
extern void dpotrf_(char *uplo, 
                    int *n, double *A, int *ldA, 
                    int *info);

extern void dpotri_(char *uplo, 
                    int *n, double *A, int *ldA, 
                    int *info);

/* Triangular inverse */
extern void dtrtri_(char *uplo, char *diag, 
                    int *n, double *A, int *ldA, 
                    int *info);
extern void dtrti2_(char *uplo, char *diag, 
                    int *n, double *A, int *ldA, 
                    int *info);

/* 
 * One/Frobenius/Infinity norm or
 * Largest Absolute Value
 */
extern double dlange_(char *norm, 
                      int *m, int *n, double *A, int *ldA,
                      double *work);

/*
 * Eigenvalue decomposition.
 *
 * Solution of a real generalized symmetric-definite eigenproblem, of the form
 * A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x,
 * where symm(A) and SPD(B)
 */
extern void dsygvd_( int *itype, char *jobz, char *uplo, 
                     int *n, double *A, int *ldA, double *B, int *ldB,
                     double *W, double *work, int *lwork, double *iwork, int *liwork, 
                     int *info );
/*
 * Reduction of a real symmetric-definite generalized eigenproblem
 * to standard form
 */
extern void dsygst_(int *itype, char *uplo,
                    int *n, double *A, int *ldA, double *B, int *ldB,
                    int *info);

extern void dsygs2_(int *itype, char *uplo,
                    int *n, double *A, int *ldA, double *B, int *ldB,
                    int *info);
/*
 * Eigendecomposition of a real symmetric matrix
 */
extern void dsyevd_( char *jobz, char *uplo, 
                     int *n, double *A, int *lda, 
                     double *W, double *work, int *lwork, int *iwork, int *liwork, 
                     int *info );
/*
 * Eigendecomposition of a real symmetric matrix
 */
extern void dsyevr_(char *jobz, char *range, char *uplo, 
                    int *n, double *A, int *ldA, 
                    double *vl, double *vu, int *il, int *iu,
                    double *abstol, int *m, double *W, double *Z, int *ldZ, int *isuppz, 
                    double *work, int *lwork, int *iwork, int *liwork,
                    int *info);
/*
 * Solution to a real system of linear equations A * X = B,
 * where symm(A)
 */
extern void dsysv_(char *uplo, 
                   int *n, int *nrhs, double *A, int *ldA, int *ipiv, 
                   double *B, int *ldB, double *work, int *lwork, 
                   int *info);

extern void dsytrs_( char *uplo, 
		             int *n, int *nrhs, double *A, int *lda, int *ipiv, 
					 double *B, int *ldb, int *info );

/*
 * Solution to a real SPD system of linear equations A * X = B ( SPD(A) )
 */
extern void dposv_(char *uplo, 
                   int *n, int *nrhs, double *A, int *ldA,
                   double *B, int *ldB,
                   int *info);

/*
 * Generalized least squares problem
 */
extern void dgels_(char *trans, int *m, int *n,
                   int *nrhs, double *A, int *lda,
                              double *B, int *ldb,
                   double *work, int *lwork, int *info);

extern void dggglm_(int *n, int *m, int *p,
                    double *A, int *lda,
                    double *B, int *ldb,
                    double *D, double *X, double *Y,
                    double *work, int *lwork, int *info);

extern void dlaset_(char *uplo, int *m, int *n,
                    double *alpha, double *beta, double *A, int *lda);

/* Solution of a general system of equations */
extern void dgetrf_( int *m, int *n, double *A, int *lda,
                     int *ipiv, int *info );

extern void dgetrs_( char *trans, int *n, int *nrhs,
                     double *A, int *lda, int *ipiv,
                     double *B, int *ldb, int *info );

extern void dlaswp_( int *n, double *A, int *lda,
                     int *k1, int *k2, int *ipiv, int *incx );

/* QR */

extern void dgeqrf_( int *m, int *n, double *A, int *lda,
                     double *tau, double *work, int *lwork, int *info );

extern void dormqr_( char *side, char *trans,
                     int *m, int *n, int *k,
                     double *A, int *lda,
                     double *tau, double *C, int *ldc,
                     double *work, int *lwork, int *info );
    

#endif // __LAPACK_H
