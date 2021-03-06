//
//  dsyevx_example.cpp
//  misdp
//
//  Created by Qi Zhang on 3/9/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#include "dsyevx_example.h"
//< Table Of Contents Intel® Math Kernel Library LAPACK Examples

/*******************************************************************************
 *  Copyright (C) 2009-2014 Intel Corporation. All Rights Reserved.
 *  The information and material ("Material") provided below is owned by Intel
 *  Corporation or its suppliers or licensors, and title to such Material remains
 *  with Intel Corporation or its suppliers or licensors. The Material contains
 *  proprietary information of Intel or its suppliers and licensors. The Material
 *  is protected by worldwide copyright laws and treaty provisions. No part of
 *  the Material may be copied, reproduced, published, uploaded, posted,
 *  transmitted, or distributed in any way without Intel's prior express written
 *  permission. No license under any patent, copyright or other intellectual
 *  property rights in the Material is granted to or conferred upon you, either
 *  expressly, by implication, inducement, estoppel or otherwise. Any license
 *  under such intellectual property rights must be express and approved by Intel
 *  in writing.
 *
 ********************************************************************************
 */
/*
 DSYEVX Example.
 ==============
 
 Program computes the smallest eigenvalues and the corresponding
 eigenvectors of a real symmetric matrix A:
 
 6.29  -0.39   0.61   1.18  -0.08
 -0.39   7.19   0.81   1.19  -0.08
 0.61   0.81   5.48  -3.13   0.22
 1.18   1.19  -3.13   3.79  -0.26
 -0.08  -0.08   0.22  -0.26   0.83
 
 Description.
 ============
 
 The routine computes selected eigenvalues and, optionally, eigenvectors of
 an n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
 
 A*v(j) = lambda(j)*v(j)
 
 where lambda(j) is its eigenvalue. The computed eigenvectors are
 orthonormal.
 Eigenvalues and eigenvectors can be selected by specifying either a range
 of values or a range of indices for the desired eigenvalues.
 
 Example Program Results.
 ========================
 
 DSYEVX Example Program Results
 
 The total number of eigenvalues found: 3
 
 Selected eigenvalues
 0.71   0.82   6.58
 
 Selected eigenvectors (stored columnwise)
 0.22   0.09  -0.95
 0.21   0.08  -0.04
 -0.52  -0.22  -0.29
 -0.73  -0.21  -0.09
 -0.32   0.94   0.01
 */
#include <stdlib.h>
#include <stdio.h>

/* DSYEVX prototype */
extern void dsyevx( char* jobz, char* range, char* uplo, int* n, double* a,
                   int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                   int* m, double* w, double* z, int* ldz, double* work, int* lwork,
                   int* iwork, int* ifail, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, double* a, int lda );

/* Parameters */
#define N 5
#define NSELECT 3
#define LDA N
#define LDZ N

/* Main program */
int main() {
    /* Locals */
    int n = N, il, iu, m, lda = LDA, ldz = LDZ, info, lwork;
    double abstol, vl, vu;
    double wkopt;
    double* work;
    /* Local arrays */
    /* iwork dimension should be at least 5*n */
    int iwork[5*N], ifail[N];
    double w[N], z[LDZ*NSELECT];
    double a[LDA*N] = {
        6.29,  0.00,  0.00,  0.00,  0.00,
        -0.39,  7.19,  0.00,  0.00,  0.00,
        0.61,  0.81,  5.48,  0.00,  0.00,
        1.18,  1.19, -3.13,  3.79,  0.00,
        -0.08, -0.08,  0.22, -0.26,  0.83
    };
    /* Executable statements */
    printf( " DSYEVX Example Program Results\n" );
    /* Negative abstol means using the default value */
    abstol = -1.0;
    /* Set il, iu to compute NSELECT smallest eigenvalues */
    il = 1;
    iu = NSELECT;
    /* Query and allocate the optimal workspace */
    lwork = -1;
    dsyevx( "Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
           &abstol, &m, w, z, &ldz, &wkopt, &lwork, iwork, ifail, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    /* Solve eigenproblem */
    dsyevx( "Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
           &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info );
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
    /* Print the number of eigenvalues found */
    printf( "\n The total number of eigenvalues found:%2i\n", m );
    /* Print eigenvalues */
    print_matrix( "Selected eigenvalues", 1, m, w, 1 );
    /* Print eigenvectors */
    print_matrix( "Selected eigenvectors (stored columnwise)", n, m, z, ldz );
    /* Free workspace */
    free( (void*)work );
    exit( 0 );
} /* End of DSYEVX Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
    int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
        printf( "\n" );
    }
}