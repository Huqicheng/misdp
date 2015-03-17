//
//  Rectangle.cpp
//  misdp
//
//  Created by Qi Zhang on 3/8/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#include "Rectangle.h"
#include <iostream>

Rectangle::Rectangle () {
    width = 5;
    height = 5;
}

Rectangle::Rectangle (int a, int b) {
    width = a;
    height = b;
}

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    
    //Computes the inverse of a triangular matrix.
    void dtrtri_(char *uplo, char*diag, int *n, double *a, int *lda, int *info);
    
    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
    //    void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
    
    
    //    void dpotrs_(char* UPLOp, long* Np, long* NRHSp, double* A, long* LDAp, double* B, long* LDBp, long* infop);
    void dpotrs_(char* UPLOp, int* Np, int* NRHSp, double* A, int* LDAp, double* B, int* LDBp, int* infop);
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
    void dpotri_(char *uplo, int *n, double *a, int *lda, int *info);
    
}


extern "C" {

/* DSYEVX prototype */
extern void dsyevx( char* jobz, char* range, char* uplo, int* n, double* a,
                   int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                   int* m, double* w, double* z, int* ldz, double* work, int* lwork,
                   int* iwork, int* ifail, int* info );
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, double* a, int lda );
}
/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
    int i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ )
            cout << a[i+j*lda] << ", " ;
//            printf( " %6.2f", a[i+j*lda] );
        printf( "\n" );
    }
}

void Rectangle::test_LAPACK(){
#define N 5
#define NSELECT 2
#define LDA N
#define LDZ N
//    int N = 5;
//    int NSELECT = 3;
//    int LDA, LDZ;
//    LDA = N; LDZ = N;
//    
    
    /* Main program */
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
//        exit( 0 );

}