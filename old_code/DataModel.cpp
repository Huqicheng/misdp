//
//  DataModel.cpp
//  SparseRegression
//
//  Created by Qi on 5/27/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#include "DataModel.h"
//#define debugMode_sigma 0
//#define debugMode_pause 1

/*
#define debug_logDet_IplusDiagzW 1
#define debug_gradient 1
*/
/*
#define debugMode_eigen 1
#define debugMode_eigen2 1
#define debugMode 1
*/


//extern "C"{
//    void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
////    int dpotrf_(const char* UPLO, const int& N, double* A, const int& LDA, int& INFO);
////    int dpotrf_(char* UPLO, int& N, double* A, int& LDA, int& INFO);
////    int dpotrf_(const char* UPLO, const integer& N, double* A, const integer& LDA, int& INFO);
////    void dsptrf_( char *uplo, int *n, double *ap, int *ipiv, int *info );
////    void dsptri_( char *uplo, int *n, double *ap, int *ipiv, double *work, int *info );
//
//}
////
////extern "C" int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
////extern "C" int dpotri_(char *uplo, int *n, double *a, int *lda, int *info);
////
//#define debug_logDet_IplusDiagzW2 1
//#define debug_logDet_IplusDiagzW3 1
//#define debug_logDet_IplusDiagzW4 1



// take out on oct 10
//#define SIGN_output_info 1
//
//#define debug_print1 0
//#define debug_uMinusLGradient 0
//
//#define debug_20140911_1pm_findSegmentationFault 0
//#define print_LU 0
//
//
//#define debug_star 1
//
//#define function_title 1


#define SIGN_maxEig1 0

#define SIGN_testTempScale 0
#define SIGN_invSum 0


#define debug_newton 0
#define debug_newton_level2 0
#define debug_newton_l3 0
#define debug_newton_level4 0

#define betaTest 1


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
    void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);
}




/*
 file dtrtri.f  dtrtri.f plus dependencies
 prec double
 for  Computes the inverse of a triangular matrix.
 gams d2a3
 
 */


//extern "C" {
//    /* DSYEV prototype */
//    void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
//}

/* DSYEV prototype */
//  extern "C" void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );


//#define D 4



using namespace std;

void DataModel :: initializerStar()
{
    vector< vector< double > > U ;
    U.resize(n);
    for (int i = 0; i < n; i++)
        U[i].resize(n);
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            U[iRow][iCol] = UMaxEigen1[iRow][iCol]    ;
        }
    }
    for (int iRow = 0; iRow < n; iRow++){
        U[iRow][iRow] -= 1;
    }
    
    vector< double > UplusUsquare(n*n) ;
    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++) {
            for (int k = 0; k < n; k++) {
                UplusUsquare[i*n+j] += U[i][k] * U[k][j];
            }
            if (i != j){
                UplusUsquare[j*n+i] = UplusUsquare[i*n+j];
            }
        }
    }
    inverse( &UplusUsquare[0], n);
    
    vector< double > Ucopy(n*n) ;
    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++) {
            Ucopy[i*n+j] = U[i][j];
            if (i != j){
                Ucopy[j*n+i] = Ucopy[i*n+j];
            }
        }
    }
    inverse( &Ucopy[0], n);
    
    
    vector< double > UinvEta(n) ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            UinvEta[i] += Ucopy[i*n+j] * eta[j];
        }
    }
    
    // r_star 分子
    double numerator = 0;
    for (int i = 0; i < n; i++){
        numerator += UplusUsquare[i*n+i]*UinvEta[i]*UinvEta[i];
        for (int j = i+1; j < n; j++) {
            numerator += 2* UplusUsquare[i*n+j]*UinvEta[i]*UinvEta[j];
        }
    }
    double denominator = 0;
    for (int i = 0; i < n; i++){
        denominator += UinvEta[i] * eta[i];
    }
    denominator = sigma_yy - denominator;
    
    
    rStar = 1 + numerator / denominator;
    
#if debug_star
    cout << "rStar = " << rStar << endl;
    {
        int tempp;
        cin >> tempp;
    }
#endif
    
    Ucopy.end();
    UplusUsquare.end();
    UinvEta.end();
    
};

void DataModel :: initializeSigma()
{
    // Sigma = [X'X   X'y;
    //          y'X   y'y ]
    Sigma.resize(n+1);
    for (int i = 0; i < n+1; i++)
        Sigma[i].resize(n+1);
    
    for (int i = 0; i < n; i++) {
        if ((i % 100) == 0)
            cout << "Sigma[i] rows done! " << i << endl;
        for (int j = i; j < n; j++) {
            double currentSpot = 0;
            for (int k = 0; k < m; k++) {
                currentSpot += X[k][i] * X[k][j];
            }
            Sigma[i][j] = currentSpot / (m);
            if (i != j)
                Sigma[j][i] = currentSpot / (m);
            
            /// force Sigma full rank !!! !!! !!! !!! !!! !!! !!!
//            if (i == j)
//                Sigma[i][i] *= 1.1;
            
        }
    }
    
 
    for (int i = 0; i < n; i++) {
        double currentSpot = 0;

        for (int k = 0; k < m; k++) {
            currentSpot += X[k][i] * y[k];
        }

        Sigma[i][n] = currentSpot / (m);
        Sigma[n][i] = currentSpot / (m);
    }
    

    {
        double currentSpot = 0;
        for (int k = 0; k < m; k++) {
            currentSpot += y[k] * y[k];
        }
        Sigma[n][n] = currentSpot / (m);


    }
    
    flag_if_sigma_initialized = true;
    
#if 0
    cout << "\n ** \n ** \n in initializeSigma " << " Sigma initialiezd" << " (tag fqo4fjq2o)" << endl;
    cout << " * * * Matrix Sigma = \n";
    for (int nr = 0; nr <= n; nr++){
        for (int nc = 0; nc <= n; nc++) {
            cout << Sigma[nr][nc] << "  " ;
        }
        cout << "\n";
    }
    cout << "\n\n\n\n\n\n// ********* ********* ********* ********* ********* ********* // \n\n\n\n\n\n";

//    cout << " * * * Matrix X = \n";
//    for (int nr = 0; nr < m; nr++){
//        for (int nc = 0; nc < n; nc++) {
//            cout << X[nr][nc] << "\t" ;
//        }
//        cout << "\n";
//    }
    
#if 0 && 1 // ! ! !
    {
        cout << "please check with Sigma (in DataModel::initializeSigma )\n input something:" << endl;
        int temp;
        cin >> temp;
    }
#endif
#endif
    
    eta.resize(n);
    sigma_yy = Sigma[n][n];
    for (int i = 0; i< n; i++) {
        eta[i]  = Sigma[n][i];
    };

};

vector< double > DataModel :: eig ( const vector< vector < double > > & matrix)
{
#if debug_20140911_1pm_findSegmentationFault
    cout << "DataModel :: eig tag (fq4f4q3f4) " ;
#endif
    double* a;
    a = (double*)malloc( n*n*sizeof(double) );
    
//    for (int iRow = 0; iRow < n; iRow++){
//        for (int iCol = 0; iCol < n; iCol++){
//            cout << matrix[iRow][iCol] << "\t";
//        }
//        cout << "\n";
//    }
    
    for (int iRow = 0, indexMatrix = 0; iRow < n; iRow ++)
        for (int iCol = 0; iCol < n; iCol ++)
            a[indexMatrix++] = matrix[iRow][iCol];
    
    double* eigens;
    eigens = (double*)malloc( n * sizeof(double) );
    double wkopt;
    int lwork = -1, info;
    char jos = 'N';
    char uplo = 'U';
    int dim = n;
    
    dsyev_( &jos, &uplo, &dim, a, &dim, eigens, &wkopt, &lwork, &info );
    //    dsyev_( &jos, &uplo, &n, a, &n, eigens, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    double* work;
    
    work = (double*)malloc( lwork*sizeof(double) );
    
    /* Solve eigenproblem */
    dsyev_( &jos, &uplo, &dim, a, &dim, eigens, work, &lwork, &info );
    //    dsyev( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );
    
    
    // Sigma = [X'X   X'y;
    //          y'X   y'y ]
    maxLambda = eigens[n - 1];
    minEigenvalue = eigens[0];
//    cout << "minEigenvalue = " << minEigenvalue ;
//    cout << ", maxEigenvalue = " << eigens[n - 1] << " (tag fj092jf2039j)"<< endl;

    vector< double > result(n);
    for (int i = 0; i < n; i++) {
        result[i] = eigens[i];
    }
    
//    for (int i = 0; i < n; i++) {
//        cout << result[i] << "\t " ;
//    }
//    cout << "\n";
    sort(result.begin(), result.end());
//    for (int i = 0; i < n; i++) {
//        cout << result[i] << "\t " ;
//    }
//    cout << "\n";
    return  result;
};

void DataModel :: scalizeSigma()
{
#if debug_20140911_1pm_findSegmentationFault
    cout << "DataModel :: scalizeSigma tag (fniewaufnaw) " ;
#endif
    SigmaOriginal.resize(n+1);
    for (int i = 0; i < n+1; i++) {
        SigmaOriginal[i].resize(n+1);
    }
    for (int iRow = 0; iRow < n+1; iRow++)
        for (int iCol = 0; iCol < n+1; iCol++)
            SigmaOriginal[iRow][iCol] = Sigma[iRow][iCol];
    
    OriginalU.resize(n);
    for (int i = 0; i < n; i++) {
        OriginalU[i].resize(n);
    }
    for (int iRow = 0; iRow < n; iRow++)
        for (int iCol = 0; iCol < n; iCol++)
            OriginalU[iRow][iCol] = Sigma[iRow][iCol];

    OriginalL.resize(n);
    for (int i = 0; i < n; i++)
        OriginalL[i].resize(n);
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            OriginalL[iRow][iCol] = SigmaOriginal[iRow][iCol] - SigmaOriginal[iRow][n]*SigmaOriginal[n][iCol] / SigmaOriginal[n][n] ;
        }
    }

    
    double* a;
    a = (double*)malloc( n*n*sizeof(double) );
    
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            cout << Sigma[iRow][iCol] << "\t";
        }
        cout << "\n";
    }


    
    for (int iRow = 0, indexMatrix = 0; iRow < n; iRow ++)
        for (int iCol = 0; iCol < n; iCol ++)
            a[indexMatrix++] = Sigma[iRow][iCol];
    
    
    double* eigens;
    eigens = (double*)malloc( n*sizeof(double) );
    double wkopt;
    int lwork = -1, info;
    char jos = 'N';
    char uplo = 'U';
    int dim = n;

    dsyev_( &jos, &uplo, &dim, a, &dim, eigens, &wkopt, &lwork, &info );
//    dsyev_( &jos, &uplo, &n, a, &n, eigens, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    double* work;

    work = (double*)malloc( lwork*sizeof(double) );
    
    /* Solve eigenproblem */
    dsyev_( &jos, &uplo, &dim, a, &dim, eigens, work, &lwork, &info );
//    dsyev( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );

    
    // Sigma = [X'X   X'y;
    //          y'X   y'y ]
    maxLambda = eigens[n - 1];
    minEigenvalue = eigens[0];
    cout << "minEigenvalue = " << minEigenvalue ;
    cout << ", maxEigenvalue = " << eigens[n - 1] << endl;
    {
#if 0 && 1 // ! ! !
        cout << " input something here: "  << endl;
        int tempp32r23;
        cin >> tempp32r23;
#endif
    }
#if debugMode_eigen2
    cout << "\nSigma before tuned = \n" ;
    for (int iRow = 0; iRow < n+1; iRow++){
        for (int iCol = 0; iCol < n+1; iCol++){
            cout << Sigma[iRow][iCol] << "\t";
        }
        cout << "\n";
    }
#endif

    
//    for (int iRow = 0; iRow < n; iRow++)
//        for (int iCol = 0; iCol < n; iCol++)
//            Sigma[iRow][iCol] /= maxLambda;
//    double sqrtMaxLambda = sqrt(maxLambda);
//    for (int i = 0; i < n; i++){
//        Sigma[n][i] /= sqrtMaxLambda;
//        Sigma[i][n] /= sqrtMaxLambda;
//    }
    for (int iRow = 0; iRow < n; iRow++)
        for (int iCol = 0; iCol < n; iCol++)
            Sigma[iRow][iCol] /= minEigenvalue;
    double sqrtMinEigenvalue = sqrt(minEigenvalue);
    for (int i = 0; i < n; i++){
        Sigma[n][i] /= sqrtMinEigenvalue;
        Sigma[i][n] /= sqrtMinEigenvalue;
    }

    
#if debugMode_eigen2
    cout << "maxLambda = " << maxLambda << endl;
    cout << "sqrtMaxLambda = " << sqrtMaxLambda << endl;
#if debugMode_pause
    {
        cout << " input sth here: " << endl;
        int temp;
        cin >> temp;
    }
#endif

#endif

    // to here, Sigma settled !
    
    
    
//  void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* INFO);

#if 0
    cout <<" wkopt = " << wkopt << "\n" ;
    cout <<"\neigens = \n";
    for (int i = 0; i < n; i++) {
        cout << "\t" << eigens[i] ;
    }
    cout << "\nSigma after tuned = \n" ;
    for (int iRow = 0; iRow < n+1; iRow++){
        for (int iCol = 0; iCol < n+1; iCol++){
            cout << Sigma[iRow][iCol] << "\t";
        }
        cout << "\n";
    }

    
    
    cout << " We are at the [DataModel :: scalizeSigma]. check the eigenvalues " << endl;
#if debugMode_pause
    {
        cout << " input sth here: " << endl;
        int temp;
        cin >> temp;
    }
#endif

#endif

//    void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

    free( (void*)work );
    free( (void*)eigens );
    free( (void*)a );
    
    
    

    //////////////////////// //////////////////////// //////////////////////// //////////////////////// ////////////////////////
    
    
//                    ######## ####  ######   ######## ##    ##
//                    ##        ##  ##    ##  ##       ###   ##
//                    ##        ##  ##        ##       ####  ##
//                    ######    ##  ##   #### ######   ## ## ##
//                    ##        ##  ##    ##  ##       ##  ####
//                    ##        ##  ##    ##  ##       ##   ###
//                    ######## ####  ######   ######## ##    ##
    
    U.resize(n);
    for (int i = 0; i < n; i++)
        U[i].resize(n);
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            U[iRow][iCol] = Sigma[iRow][iCol] ;
        }
    }
    
    L.resize(n);
    for (int i = 0; i < n; i++)
        L[i].resize(n);
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            L[iRow][iCol] = Sigma[iRow][iCol] - Sigma[iRow][n]*Sigma[n][iCol] / Sigma[n][n] ;
        }
    }

//#if debug_print1
//    cout << "\nL before tuned = \n" ;
//    for (int iRow = 0; iRow < n; iRow++){
//        for (int iCol = 0; iCol < n; iCol++){
//            cout << L[iRow][iCol] << "\t";
//        }
//        cout << "\n";
//    }
//    cout << "\nU before tuned = \n" ;
//    for (int iRow = 0; iRow < n; iRow++){
//        for (int iCol = 0; iCol < n; iCol++){
//            cout << U[iRow][iCol] << "\t";
//        }
//        cout << "\n";
//    }
//#endif
    vector < double > eigenL(eig(L));
#if debug_print1
    cout << " * eigen value = ";
    for (int i = 0; i < n; i++) {
        cout << eigenL[i] << "\t" << endl;
    }
#endif
//    double minEigen = eigenL[0] * n ;
//    double minEigen = eigenL[0] * n * n;
    double minEigen = eigenL[0] ;
    
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            U[iRow][iCol] /= minEigen ;
            L[iRow][iCol] /= minEigen ;
        }
    }
#if print_LU
    cout << "original eigenL[0] = " << eigenL[0] << endl;
    cout << "minEigen = " << minEigen << endl;

    cout << "\nL after tuned = \n" ;
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            cout << L[iRow][iCol] << "\t";
        }
        cout << "\n";
    }
    cout << "\nU after tuned = \n" ;
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            cout << U[iRow][iCol] << "\t";
        }
        cout << "\n";
    }
    cout << "eigenL[*] = \n" ;
    for (int i = 0; i < n; i++) {
        cout << eigenL[i] << "\t";
    }
    cout << "\n";

    cout << "original eigenL[0] = " << eigenL[0] << endl;
    cout << "minEigen = " << minEigen << endl;

    #if 0 && 1 // ! ! !
        {
            cout << " input sth here: " << endl;
            int temp;
            cin >> temp;
        }
    #endif
#endif
    
//                ######## ####  ######   ######## ##    ##    ##     ##    ###    ##     ##       ##
//                ##        ##  ##    ##  ##       ###   ##    ###   ###   ## ##    ##   ##      ####
//                ##        ##  ##        ##       ####  ##    #### ####  ##   ##    ## ##         ##
//                ######    ##  ##   #### ######   ## ## ##    ## ### ## ##     ##    ###          ##
//                ##        ##  ##    ##  ##       ##  ####    ##     ## #########   ## ##         ##
//                ##        ##  ##    ##  ##       ##   ###    ##     ## ##     ##  ##   ##        ##
//                ######## ####  ######   ######## ##    ##    ##     ## ##     ## ##     ##     ######

    
    

    UMaxEigen1.resize(n);
    for (int i = 0; i < n; i++)
        UMaxEigen1[i].resize(n);
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            UMaxEigen1[iRow][iCol] = U[iRow][iCol]  ;
        }
    }

    vector < double > eigenU(eig(UMaxEigen1));

    maxEigenU_ieSigma = eigenU[n-1];

    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            UMaxEigen1[iRow][iCol] /= maxEigenU_ieSigma  ;
        }
    }
    
    LMaxEigen1.resize(n);
    for (int i = 0; i < n; i++)
        LMaxEigen1[i].resize(n);
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            LMaxEigen1[iRow][iCol] = ( Sigma[iRow][iCol] - Sigma[iRow][n]*Sigma[n][iCol] / Sigma[n][n] ) / maxEigenU_ieSigma;
        }
    }
#if maxEigenU_ieSigma // ! ! !
    {
        cout << " input sth here: (tag fwqfaw4ves5nrs5mndr6fq34g34q)" << endl;
        cout << " maxEigenU_ieSigma =  " << maxEigenU_ieSigma << endl;
        int temp;
        cin >> temp;
    }
#endif
#if 1
    cout << " eigenU (tag fjoawef03nj98)  \n";
    for (int i=0; i<n; i++) {
        cout << eigenU[i] << "\t";
    }
    cout << "  \n";

#endif
//
//    for (int iRow = 0; iRow < n; iRow++){
//        for (int iCol = 0; iCol < n; iCol++){
//            UMaxEigen1[iRow][iCol] /= maxEigenL ;
//            LMaxEigen1[iRow][iCol] /= maxEigenL ;
//        }
//    }

    
};



// added 10/09/2014
void DataModel :: initializeOriginalInvOfInvSum(const TuneParameters &tune)
{
#if SIGN_invSum
    cout << "DataModel :: initializeOriginalInvOfInvSum tag (faw4f4ta4gtbyrt54w) " ;
#endif
    
    vector < double > tempU(n*n); // = OriginalU;
    vector < double > tempL(n*n); // = OriginalL;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            tempU[i*n+j] = OriginalU[i][j];
            tempL[i*n+j] = OriginalL[i][j];
        }
    }
#if SIGN_invSum
    cout << " \n\n\n\n tag (4be577aw3v4) " << endl;
    cout << " tempU = " ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << tempU[i*n+j] << "\t" ;
        cout << endl;
    }
    cout << " tempL = " << endl ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << tempL[i*n+j] << "\t" ;
        cout << endl;
    }
#endif
    
    inverse(&tempU[0], n);
    inverse(&tempL[0], n);
    
    
    OriginalUinverse.resize(n);
    for (int i =0; i<n; i++) {
        OriginalUinverse[i].resize(n);
    }
    OriginalLinverse.resize(n);
    for (int i =0; i<n; i++) {
        OriginalLinverse[i].resize(n);
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            OriginalUinverse[i][j] = tempU[i*n+j];
            OriginalLinverse[i][j] = tempL[i*n+j];
        }
    }


#if SIGN_invSum
    cout << " \n\n\n\n tag (z4bdtfyrst4wagstr4), after inverse "  << endl;
    cout << " tempU = "  << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << tempU[i*n+j] << "\t" ;
        cout << endl;
    }
    cout << " tempL = " << endl ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << tempL[i*n+j] << "\t" ;
        cout << endl;
    }
#endif

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            tempU[i*n+j] += tempL[i*n+j];
        }
    }
#if SIGN_invSum
    cout << " \n\n\n\n tag (aegbrst4wgt), sum up " << endl ;
    cout << " tempU = " << endl ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << tempU[i*n+j] << "\t" ;
        cout << endl;
    }
#endif

    inverse(&tempU[0], n);
#if SIGN_invSum
    cout << " \n\n\n\n tag (aestr ydrse4ut), inverse " << endl ;
    cout << " tempU = " << endl ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << tempU[i*n+j] << "\t" ;
        cout << endl;
    }
#endif

//    for (int i = 0; i < n; i++){
//        for (int j = 0; j < n; j++){
//            tempU[i*n+j] *= 2;
//        }
//    }
    
    
    OriginalInvOfInvSum.resize(n);
    for (int i = 0; i < n; i++) {
        OriginalInvOfInvSum[i].resize(n);
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            OriginalInvOfInvSum[i][j] = tempU[i*n+j] * 2;
        }
    }
#if SIGN_invSum
    cout << " \n\n\n\n tag (f43g5rturg5th67) " << endl ;
    cout << " OriginalInvOfInvSum = " << endl ;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << OriginalInvOfInvSum[i][j] << "\t" ;
        cout << endl;
    }
#endif

};


vector < vector < double > > DataModel :: scaleBothSide(const vector < vector < double > > &A, const vector < double > &s)
{
    vector < vector < double > > B = A;
    for (int i = 0; i < B.size(); i++) {
        for(int j = 0; j < B[i].size(); j++){
            B[i][j] *= s[i]*s[j];
        }
    }
    return B;
};







// added 10/09/2014
void DataModel :: initializeUandLatNSD(const TuneParameters &tune)
{
    
    /*
     // in this function, we initialize the following
     vector< vector < double > > UatNegSemiDef ;
     vector< vector < double > > LatNegSemiDef ;
     vector< vector < double > > scalerUatNSD ;
     vector< vector < double > > scalerLatNSD ;
     */
    
#if function_title
    cout << " at the begining of DataModel :: initializeUandLatNSD (tag foqw8f230f23q0u)" << endl;
#endif
    
    LatNegSemiDef.resize(n);
    for (int i = 0; i < n; i++)
        LatNegSemiDef[i].resize(n);
    UatNegSemiDef.resize(n);
    for (int i = 0; i < n; i++)
        UatNegSemiDef[i].resize(n);
    scalerUatNSD.resize(n);
    scalerLatNSD.resize(n);
    
    // L first
    for (int iRow = 0; iRow < n; iRow++)
        for (int iCol = 0; iCol < n; iCol++)
            LatNegSemiDef[iRow][iCol] = SigmaOriginal[iRow][iCol] - SigmaOriginal[iRow][n]*SigmaOriginal[n][iCol] / SigmaOriginal[n][n] ;
    
    
    
#if betaTest
    printf(" SigmaOriginal tag g534wgw45 \n");
    for (int i=0; i<n+1; i++) {
        for(int j=0; j<n+1; j++)
            printf(" %f\t", SigmaOriginal[i][j]);
        printf(" \n");
    }
#endif

    
    
#if SIGN_output_info
    printf(" LatNegSemiDef tag 43fqw34g45whgw45 \n");
    for (int i=0; i<n; i++) {
        for(int j=0; j<n; j++)
            printf(" %f\t", LatNegSemiDef[i][j]);
        printf(" \n");
    }
#endif
    for (int i = 0; i < n; i++) {
        scalerLatNSD[i] = sqrt(LatNegSemiDef[i][i]);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            LatNegSemiDef[i][j] /= scalerLatNSD[i];
            LatNegSemiDef[j][i] /= scalerLatNSD[i];
        }
    }    // to here, diag(LatNegSemiDef) == 1 1 1 1 1 1 1 1
    
    
    
    
    
#if SIGN_output_info
    printf(" LatNegSemiDef tag wf4fq234f2q  \n");
    for (int i=0; i<n; i++) {
        for(int j=0; j<n; j++)
            printf(" %f\t", LatNegSemiDef[i][j]);
        printf(" \n");
    }
#endif
#if SIGN_output_info || 1
    printf(" scalerLatNSD tag q34fq43gq34 \n");
    for (int i=0; i<n; i++) {
        printf("%f\t", scalerLatNSD[i]);
    }
    printf("\n");
#endif

    vector < double > eigenTempL(eig(LatNegSemiDef));
#if SIGN_output_info
    printf(" eigenTempL tag fw34f34 \n");
    for (int i=0; i<n; i++) {
        printf("%f\t", eigenTempL[i]);
    }
    printf("\n");
#endif
    
    double maxEigValueTempL = eigenTempL[n-1];
    for (int i = 0; i < n; i++) {
        scalerLatNSD[i] *= sqrt(maxEigValueTempL);
        for (int j = 0; j < n; j++){
            LatNegSemiDef[i][j] /= maxEigValueTempL;
        }
    }    // to here, max eig(LatNegSemiDef) == 1, diagonal the same (try my best to maintain good condition number
    for (int i = 0; i < n; i++)
        scalerLatNSD[i] = 1.0 / scalerLatNSD[i];
#if SIGN_output_info
    printf(" LatNegSemiDef tag tafaw4awe final L  \n");
    for (int i=0; i<n; i++) {
        for(int j=0; j<n; j++)
            printf(" %f\t", LatNegSemiDef[i][j]);
        printf(" \n");
    }
#endif

#if SIGN_output_info
    printf(" scalerLatNSD  tag fq34fq34f43 final scaler of L\n");
    for (int i=0; i<n; i++) {
        printf("%f\t", scalerLatNSD[i]);
    }
    printf("\n");
#endif
    

    
    
    // U here
    int N = n;
    vector< double> sqrtU(N * N);
//    sqrtU.resize(N * N);
    
    int indexMatrix = 0;
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            sqrtU[indexMatrix] = SigmaOriginal[iRow][iCol];
            indexMatrix++;
        }
    }
    
    
    
    // ! ! ! take out the 4.4
    for (int i = 0; i < n; i++) {
        scalerUatNSD[i] = 1.0/sqrt(sqrtU[i*N+i]);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            sqrtU[i*N+j] *= scalerUatNSD[i];
            sqrtU[j*N+i] *= scalerUatNSD[i];
        }
    }    // to here, diag(LatNegSemiDef) == 1 1 1 1 1 1 1 1

    
    
    
    
    FILE *out;
    int LDA, INFO;
    LDA = N;
//    int M, N, LDA, INFO;
    out = stderr;
#if SIGN_output_info
    fprintf(out, "N = %d\n", N);
    fprintf(out, "A before call, original U or Sigma :\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
//            fprintf(out, "%f\t", sqrtU[i+N*j]);
            fprintf(out, "%f\t", sqrtU[j+N*i]);
        fprintf(out, "\n");
    }
#endif
    char uplo = 'L';
#if SIGN_output_info
    cout << " before dpotrf_ (tag jfoq48jfq2340):" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    //void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
    //    dpotrf_(&uplo, &N, A, &LDA, &INFO);
    dpotrf_(&uplo, &N, &sqrtU[0], &LDA, &INFO);
    //    dgetrf_(&M, &N, A, &LDA, &IPIV[0], &INFO);
#if SIGN_output_info
    cout << " after dpotrf_ :" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    
#if SIGN_output_info
    fprintf(out, "INFO: %d\n", INFO);
    cout << " N = " << N <<  endl;
    fprintf(out, "A after call (&lower tri&):\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", sqrtU[j+N*i]);
        fprintf(out, "\n");
    }
#endif
    
    
//    dtrtri_
//    void dtrtri_(char *uplo, char*diag, int *n, double *a, int *lda, int *info);
    char diag = 'N';
    dtrtri_(&uplo, &diag, &N, &sqrtU[0], &LDA, &INFO);

#if SIGN_output_info
    fprintf(out, "INFO: %d\n", INFO);
    cout << " N = " << N <<  endl;
    fprintf(out, "A after call (be the inverse &lower tri&), i.e., sqrtU, chol(Sigma) :\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", sqrtU[j+N*i]);
        fprintf(out, "\n");
    }
#endif

    vector< vector < double > > tempU ; // L'^-.5 * LatNSD * L^-.5
    tempU.resize(n);
    for (int i = 0; i < n; i++)
        tempU[i].resize(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            tempU[i][j] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // assign tempU[i][j]
            for (int k = 0; k <= j; k++) {
//                tempU[i][j] += LatNegSemiDef[i][k]*sqrtU[k+N*j]; // ik * kj
                tempU[i][j] += LatNegSemiDef[i][k]*sqrtU[j+N*k]; // ik * kj
            }
        }
    }
#if SIGN_output_info
    fprintf(out, "tempU = JSJ * BB \n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", tempU[i][j]);
        fprintf(out, "\n");
    }
#endif
    vector< vector < double > > tempU2 ; // L'^-.5 * LatNSD * L^-.5
    tempU2.resize(n);
    for (int i = 0; i < n; i++)
        tempU2[i].resize(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            tempU2[i][j] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // assign tempU[i][j]
            for (int k = 0; k <= i; k++) {
//                tempU2[i][j] += sqrtU[k+N*i] * tempU[k][j]; // ik * kj
                tempU2[i][j] += sqrtU[i+N*k] * tempU[k][j]; // ik * kj
            }
        }
    }
    
#if SIGN_output_info
    fprintf(out, " sigma^-.5 * (J Sigma - etaeta' J) sigma^-.5 , i.e., tempU2\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", tempU2[i][j]);
        fprintf(out, "\n");
    }
#endif

    
    
    
    
    
    vector < double > eigenTempU(eig(tempU2));
    double tempScalerUatNSD = sqrt( eigenTempU[0] );
    for (int i = 0; i < n; i++) {
        scalerUatNSD[i] *= tempScalerUatNSD;
    }
    


    
    
#if SIGN_output_info
    printf(" eigenTempU tag fawf232q3f \n");
    for (int i=0; i<n; i++) {
        printf("%f\t", eigenTempU[i]);
    }
    printf("\n");
#endif

    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
//            UatNegSemiDef[iRow][iCol] = SigmaOriginal[iRow][iCol] * scalerUatNSD[iRow] ;
//            UatNegSemiDef[iCol][iRow] = SigmaOriginal[iCol][iRow] * scalerUatNSD[iRow] ;
            UatNegSemiDef[iRow][iCol] = SigmaOriginal[iRow][iCol] * scalerUatNSD[iRow] * scalerUatNSD[iCol]  ;
        }
    }
#if SIGN_output_info
    fprintf(out, " UatNegSemiDef (tag fq928fj3498)\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", UatNegSemiDef[i][j]);
        fprintf(out, "\n");
    }
#endif
    
#if SIGN_testTempScale
    
    for (int i = 0; i < n; i++) {
        scalerUatNSD[i] *= tune.testNumber1 ;
    }
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            UatNegSemiDef[iRow][iCol] = SigmaOriginal[iRow][iCol] * scalerUatNSD[iRow] * scalerUatNSD[iCol]  ;
        }
    }
    fprintf(out, " UatNegSemiDef after scaled (tag aw4vgrs4wvt) \n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", UatNegSemiDef[i][j]);
        fprintf(out, "\n");
    }
#endif
    
#if SIGN_output_info
    fprintf(out, " LatNegSemiDef (tag 4320q9fjq3409j)\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", LatNegSemiDef[i][j]);
        fprintf(out, "\n");
    }
#endif
    
#if SIGN_output_info
    cout << " tune.testNumber1 = " << tune.testNumber1 << endl;
    
    fprintf(out, " scalerLatNSD (tag ra3wv4rydc4w)\n");
    for (int i=0; i<N; i++) {
        fprintf(out, "%f\t", scalerLatNSD[i]);
    }
    fprintf(out, "\n");
    fprintf(out, " scalerUatNSD (tag tdu65eb yf4w)\n");
    for (int i=0; i<N; i++) {
        fprintf(out, "%f\t", scalerUatNSD[i]);
    }
#endif
    
    if ( n > 10)
        return;
#if betaTest
    fprintf(out, " SigmaOriginal (tag fawv4 dyfd4c td)\n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout <<SigmaOriginal[i][j] << "\t" ;
        }
        cout << "\n";
    }

#endif
#if SIGN_output_info
    fprintf(out, "\n");
    fprintf(out, " maxEigenU_ieSigma \n");
    fprintf(out, "%f\t", maxEigenU_ieSigma);
    fprintf(out, "\n");
    fprintf(out, " SigmaOriginal (tag fawv4 dyfd4c td)\n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout <<SigmaOriginal[i][j] << "\t" ;
        }
        cout << "\n";
    }
    fprintf(out, "\n");
    fprintf(out, " Sigma \n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout << Sigma[i][j] << "\t" ;
        }
        cout << "\n";
    }
    fprintf(out, "\n");


    fprintf(out, " \n L \n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout << L[i][j] << "\t" ;
        }
        cout << "\n";
    }
    fprintf(out, "\n");
    
    
    fprintf(out, " \n U \n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout << U[i][j] << "\t" ;
        }
        cout << "\n";
    }
    fprintf(out, "\n");
    
    
    fprintf(out, " \n LatNegSemiDef (tag 3vwrve4 5e4)\n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout << LatNegSemiDef[i][j] << "\t" ;
        }
        cout << "\n";
    }
    fprintf(out, "\n");


    fprintf(out, " \n UatNegSemiDef (tag 234 wt65vwe r) \n");
    for (int i=0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            cout << UatNegSemiDef[i][j] << "\t" ;
        }
        cout << "\n";
    }
    fprintf(out, "\n");


    {
        int temp;
        cin >> temp;
    }

#endif

    
};







double DataModel :: logDetDelta( const bitset< _CONST_INT_MAX_FACTOR_SIZE > & i, const bitset< _CONST_INT_MAX_FACTOR_SIZE > & S )
{
#if debug_20140911_1pm_findSegmentationFault
    cout << " tag 0f92j30923 " << endl;
#endif
    // delta_i(S) = f(S\cup i)-f(S\i)
//    std::cout << " S = " << S << ", i = " << i << ", (S|i) = " << (S|i) << ", (S&(~i)) = " << (S&(~i)) << endl;
//    std::cout << " logDetSubmatrix(S|i) = " << logDetSubmatrix(S|i) << ", logDetSubmatrix(S&(~i)) = " << logDetSubmatrix(S&(~i)) << endl;
    return (  logDetSubmatrix(S|i) - logDetSubmatrix(S&(~i)) );
};

double DataModel :: logDetDelta( const int i, bitset< _CONST_INT_MAX_FACTOR_SIZE > S )
{
    // delta_i(S) = f(S\cup i)-f(S\i)
    return (  logDetSubmatrix(S.set(i)) - logDetSubmatrix(S.reset(i)) );
};




void DataModel :: push_back(const vector < double > &a){
    Sigma.push_back(a);
};

void DataModel :: push_back_y( const double a ){
    y.push_back(a);
};

void DataModel :: push_back_X( const vector < double > &a ){
    X.push_back(a);
};


void DataModel :: setCols(int & a){
    n = a;
};

void DataModel :: setRows(int & a){
    m = a;
    Sigma.reserve(m);
};

int DataModel :: getRows() const{
    return(m);
};

int DataModel :: getCols() const{
    return(n);
};


double DataModel :: logDetSubmatrix_with_y (const bitset<_CONST_INT_MAX_FACTOR_SIZE> & setsIndicators)
{
//    cout << "  DataModel :: logDetSubmatrix_with_y tag (g432b 24) " << endl;
    bitset<_CONST_INT_MAX_FACTOR_SIZE> s = setsIndicators;
    s.set(n);
    return logDetSubmatrix(s);
};


double DataModel :: logDetSubmatrix (const bitset< _CONST_INT_MAX_FACTOR_SIZE > & setsIndicators)
{
#if debug_20140911_1pm_findSegmentationFault
    cout << "  DataModel :: logDetSubmatrix tag (f3be45n6) " << endl;
#endif
    
    
#if debugMode
    cout << "\n ** \n ** \n in logDetSubmatrix " << " setsIndicators = " << setsIndicators << " (tag ymu67ju75e6h)" << endl;
#endif
    if (setsIndicators.none())
        return 0.0;

#if debug_20140911_1pm_findSegmentationFault
    cout << "   tag (gwn5w345nw45) " << endl;
#endif
    unordered_map< bitset<_CONST_INT_MAX_FACTOR_SIZE>, double >::const_iterator got = calculatedValues.find(setsIndicators);
    if ( got != calculatedValues.end() )
        return got->second;
#if debug_20140911_1pm_findSegmentationFault
    cout << "   tag (he46nmhw45n) " << endl;
#endif

    double result = 0;
    FILE *out;
//    int M, N, LDA, IPIV[ _CONST_INT_MAX_FACTOR_SIZE ], INFO;
    int M, N, LDA, INFO;
    N = (int) setsIndicators.count();
    M = N;
    LDA = M;
#if debug_20140911_1pm_findSegmentationFault
    cout << "   tag (g34wbw34n) " << endl;
#endif


#if debugMode
    cout << "\n ** \n ** \n in logDetSubmatrix " << " (tag f3bwhe57nr68)" << endl;
#endif
    vector<int> seletedIndex;
    seletedIndex.reserve(N);
#if debugMode
    cout << "\n ** \n ** \n in logDetSubmatrix " << " (tag h43y6h4hrg3w)" << endl;
#endif
    int n = getCols();

    for (int i = 0; i < n+1; i++)
        if (setsIndicators[i])
            seletedIndex.push_back(i);
#if debugMode
    cout << "\n ** \n ** \n in logDetSubmatrix " << " (tag g5nj765bhw4)" << endl;
    cout << " _CONST_INT_MAX_FACTOR_SIZE = " << _CONST_INT_MAX_FACTOR_SIZE << endl;
#endif
    //    double  A[ _CONST_INT_MAX_FACTOR_SIZE * _CONST_INT_MAX_FACTOR_SIZE ];
    vector< double> A(N * N);
//    A.resize(N * N);
    
    int indexMatrix = 0;
    for (vector<int>::iterator itRow = seletedIndex.begin(); itRow != seletedIndex.end(); itRow ++ ){
        for (vector<int>::iterator itCol = seletedIndex.begin(); itCol != seletedIndex.end(); itCol ++ ){
#if 0
            cout << "\n ** \n ** \n in logDetSubmatrix " << " (tag f3vg56n685)" << endl;
            cout << " seletedIndex.size() = " << seletedIndex.size() << endl;
            for (int tempp = 0; tempp < seletedIndex.size(); tempp++)
                cout << seletedIndex.at(tempp) << "\t";
            cout << endl;
            cout << " [*itRow] = " << *itRow << ", [*itCol] = " << *itCol << endl;
#endif
            A[indexMatrix] = Sigma[*itRow][*itCol];
            indexMatrix++;
        }
    }
#if debug_20140911_1pm_findSegmentationFault
    cout << "   tag (fwbw45ne56) " << endl;
#endif
#if debugMode
    cout << "A = \n" << endl;
    int indexMatrix_temp = 0;
    for (vector<int>::iterator itRow = seletedIndex.begin(); itRow != seletedIndex.end(); itRow ++ ){
        for (vector<int>::iterator itCol = seletedIndex.begin(); itCol != seletedIndex.end(); itCol ++ ){
            cout << A[indexMatrix_temp++] << "\t";
        }
        cout<<endl;
    }
#endif

    
#if debugMode
    cout << "\n ** \n ** \n in logDetSubmatrix (Cholesky decomposition, to do logdet) " << " (tag ikr5vwerfq4wfv)" << endl;
#endif

/////////////////////
    
    
    out = stderr;
#if SIGN_output_info
    fprintf(out, "A before call:\n");
    fprintf(out, "N = %d\n", N);
#endif
    
#if SIGN_output_info
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
#if debug_20140911_1pm_findSegmentationFault
    cout << "   tag (g34qnw3n) " << endl;
#endif

    char uplo = 'L';

#if SIGN_output_info
    cout << " before dpotrf_ (tag 2h873q9q34h):" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    //void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
//    dpotrf_(&uplo, &N, A, &LDA, &INFO);
    dpotrf_(&uplo, &N, &A[0], &LDA, &INFO);
//    dgetrf_(&M, &N, A, &LDA, &IPIV[0], &INFO);
#if SIGN_output_info
    cout << " after dpotrf_ :" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
#if debug_20140911_1pm_findSegmentationFault
    cout << "   tag (fq34bg3q4n) " << endl;
#endif

#if SIGN_output_info
    fprintf(out, "INFO: %d\n", INFO);
    cout << " N = " << N <<  endl;
    fprintf(out, "A after call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif

#if SIGN_output_info
    fprintf(out, "VECR:\n");
#endif
    double piii=1;
    for (int i=0; i<N; i++) {
#if SIGN_output_info
        fprintf(out, "%f\t", A[i+N*i]);
#endif
        piii = piii*A[i+N*i];
        result += log(A[i*N+i]);
        /* !!! why negative diagonal */
    }
    
    // cholesky give's L'*L = A, we need to double up
    result *= 2;
    
    
#if 0 // SIGN_output_info
    fprintf(out, ", result = %f\t", result);

    fprintf(out, ", piii = %f\n",piii);
    fprintf(out, "P:\n");
    for(int j=0; j<N; j++)
        fprintf(out, "%d\t", IPIV[j]);
    fprintf(out, "\n");
    
    cout << endl << "det = " << piii << ", logdet = " << result << endl<< endl<< endl;
#endif
    
    
    
    
#if debugMode
    cout << "\n ** \n ** \n in logDetSubmatrix (end) " << " (tag f3q4vgw45bhe56b)" << endl;
    cout << " result = " << result << endl;
#if debugMode_pause
    {
        int temp;
        cin >> temp;
    }
#endif

#endif


    calculatedValues.emplace(setsIndicators, result);
    
    return result;
};

void DataModel :: inverse(double* A, int N)
{
#if debug_newton
    cout << " start DataModel :: inverse (tag foaiwf890430)" << endl;
#endif
    //    int *IPIV = new int[N+1];
    
    vector<int> IPIV(N+1);
    int LWORK = N*N;
    vector<double> WORK(LWORK);
    //    double *WORK = new double[LWORK];
    int INFO;
    dgetrf_(&N,&N,A,&N,&IPIV[0],&INFO);
    dgetri_(&N,A,&N,&IPIV[0],&WORK[0],&LWORK,&INFO);
    
    //    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    //    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
#if debug_newton
    cout << " finished dgetri_ " << endl;
#endif
    
    //    delete IPIV;
#if debug_newton
    cout << " sdfawrgearger sefwaegaewrgare " << endl;
#endif
    
    //    delete WORK;
#if debug_newton
    cout << " finished DataModel :: inverse" << endl;
#endif
    
};

void DataModel :: inverse( vector<double> & A, int N)
{
#if debug_newton
    cout << " start DataModel :: 3" << endl;
#endif
    vector<int> IPIV(N+1);
    int LWORK = N*N;
    vector<double> WORK(LWORK);
    int INFO;
    dgetrf_(&N,&N,&A[0],&N,&IPIV[0],&INFO);
    dgetri_(&N,&A[0],&N,&IPIV[0],&WORK[0],&LWORK,&INFO);
    
    
#if debug_newton
    
    cout << " (tag fsefaerwgarew)\n A = " << endl;
    int n = N;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << A[i*n+j] << ",  " ;
        cout << endl;
    }
    
    
    cout << " finished DataModel :: inverse dfawfaerg"  << endl;
#endif
    
};

void DataModel :: inversePSD( vector<double> & A, int n)
{
#if debug_newton
    cout << " start DataModel :: inversePSD" << endl;
#endif
//    vector<int> IPIV(n);
//    int LWORK = N*N;
//    vector<double> WORK(LWORK);
    int INFO;
    int LDA = n;
    
    char uplo = 'L';
//    void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
//    void dpotri_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info)
#if debug_newton
    
    cout << " (tag fsefaerwgarew)\n before inverse A = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << A[i*n+j] << ",  " ;
        cout << endl;
    }
    
    
    cout << " finished DataModel :: inversePSD 4th56nj6r"  << endl;
#endif

    
    dpotrf_(&uplo, &n, &A[0], &LDA, &INFO);
    dpotri_(&uplo, &n, &A[0], &LDA, &INFO);
    
//    dpotrf_(&N,&N,&A[0],&N,&IPIV[0],&INFO);
//    dpotri_(&N,&A[0],&N,&IPIV[0],&WORK[0],&LWORK,&INFO);
    
#if debug_newton
    
    cout << " (tag fsefaerwgarew)\n after inverse A = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << A[i*n+j] << ",  " ;
        cout << endl;
    }
#endif

    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            A[i*n+j] = A[j*n+i];

    
#if debug_newton
    
    cout << " (tag fsefaerwgarew)\n after inverse A = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << A[i*n+j] << ",  " ;
        cout << endl;
    }
    
    
    cout << " finished DataModel :: inversePSD 4th56nj6r"  << endl;
#endif
    
    
    
    
    
};




double DataModel :: logDet_IplusDiagzW( const IloNumArray &valz )
{
    return logDet_IplusDiagzW(valz, Sigma);
};


double DataModel ::  logDet_IplusDiagzW( const vector < double > & valz , const vector< vector < double > > & Sigma )
{
    double result = 0;
    //    FILE *out;
    vector< double> A( n*n );
    
    vector< double> valzSQRT (n);
    for (int i = 0; i < n; i++) {
        if (valz[i] <= 0)
            valzSQRT[i] = 0;
        else
            valzSQRT[i] = sqrt(valz[i]);
    }
    
    int N, LDA, IPIV[ _CONST_INT_MAX_FACTOR_SIZE ], INFO;
    N = n;
    LDA = N;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                A[i*n+j] = 1 + ( Sigma[i][j] - 1 ) * valzSQRT[i] * valzSQRT[j];
            } else {
                A[i*n+j] = Sigma[i][j] * valzSQRT[i] * valzSQRT[j];
            }
        }
    };
    
    // here A is symmetric
    
#if debug_logDet_IplusDiagzW
    cout << " > > in {DataModel::logDet_IplusDiagzW} " << endl;
    cout << " * diaganol of the matrix = [ " << endl;
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            cout << A[i*N+j] << ", " ;
        }
        cout <<endl;
    }
    cout << " ]"<<  endl;
#endif
    
    
    
    //    out = stderr;
#if debug_logDet_IplusDiagzW2
    cout <<" ( tag fqo4f89q34980q ) " << endl;
    printf("A before call:\n");
    printf("N = %d\n", N);
#endif
#if debug_logDet_IplusDiagzW2
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            printf( "%f\t", A[i+N*j]);
        printf("\n");
    }
#endif
    
    dgetrf_(&N, &N, &A[0], &LDA, IPIV, &INFO);
    
    
    // use dpbtrf
    
    
#if debug_logDet_IplusDiagzW2
    printf("INFO: %d\n", INFO);
    cout << " N = " << N << endl;
    printf("A after call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            printf( "%f\t", A[i+N*j]);
        printf( "\n");
    }
#endif
    
#if debug_logDet_IplusDiagzW2
    printf( "VECR:\n");
#endif
    double piii=1;
    for (int i=0; i<N; i++) {
#if debug_logDet_IplusDiagzW2
        printf( "%f\t", A[i+N*i]);
#endif
        piii=piii*A[i+N*i];
        result += log(A[i*N+i]);
    }
#if debug_logDet_IplusDiagzW2
    printf( "%f\t", result);
    
    printf( "%f\n",piii);
    printf( "P:\n");
    for(int j=0; j<N; j++)
        printf( "%d\t", IPIV[j]);
    printf( "P(end)\n");
    
#endif
#if debug_logDet_IplusDiagzW3
    printf( "* * * Adiag:\n");
    for(int i=0; i<N; i++)
        printf( "%f\t", A[i*N+i]);
    printf( "Adiag(end)\n");
    
#endif
    
#if debugMode_pause
    {
        int temp;
        cin >> temp;
    }
#endif
    
    
#if debug_logDet_IplusDiagzW
    cout << " > > in {DataModel::logDet_IplusDiagzW} " << endl;
    cout << " * diaganol of the matrix = [ " ;
    for (int i = 0; i < N; i++) {
        cout << A[i*N+i] << ", " ;
    }
    cout << " ]" <<  endl;
    cout << " * valzSQRT = [ " ;
    for (int i = 0; i < N; i++) {
        cout << valzSQRT[i] << ", " ;
    }
    cout << " ]" <<  endl;
    
    
#endif
    
    
    return result;
};

double DataModel ::  logDet_IplusDiagzW( const IloNumArray &valz , const vector< vector < double > > & Sigma )
{
    double result = 0;
    //    FILE *out;
    vector< double> A( n*n );
    
    vector< double> valzSQRT (n);
    for (int i = 0; i < n; i++) {
        if (valz[i] <= 0)
            valzSQRT[i] = 0;
        else
            valzSQRT[i] = sqrt(valz[i]);
    }
    
//    int N, LDA, IPIV[ _CONST_INT_MAX_FACTOR_SIZE ], INFO;
    int N, LDA, INFO;
    int M;
    M = n;
    N = n;
    LDA = N;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                A[i*n+j] = 1 + ( Sigma[i][j] - 1 ) * valzSQRT[i] * valzSQRT[j];
            } else {
                A[i*n+j] = Sigma[i][j] * valzSQRT[i] * valzSQRT[j];
            }
        }
    };
    
    // here A is symmetric
    
    //
    //
    //
    //
    //
    
    FILE *out;
    
    
    out = stderr;
#if SIGN_output_info
    fprintf(out, "A before call:\n");
    fprintf(out, "N = %d\n", N);
#endif
    
#if SIGN_output_info
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
#if debug_logDet_IplusDiagzW4
    fprintf(out, "A before call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
    
    char uplo = 'L';
    
#if SIGN_output_info
    cout << " before dpotrf_ :" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    //void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
    dpotrf_(&uplo, &N, &A[0], &LDA, &INFO);
    //    dgetrf_(&M, &N, A, &LDA, &IPIV[0], &INFO);
#if SIGN_output_info
    cout << " after dpotrf_ :" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    
#if SIGN_output_info
    fprintf(out, "INFO: %d\n", INFO);
    cout << " N = " << N <<  endl;
    fprintf(out, "A after call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
    
#if SIGN_output_info
    fprintf(out, "VECR:\n");
#endif
    double piii=1;
    for (int i = 0; i < N; i++) {
#if SIGN_output_info
        fprintf(out, "%f\t", A[i+N*i]);
#endif
        piii = piii*A[i+N*i];
        result += log(A[i*N+i]);
        /* !!! why negative diagonal */
    }
    
    // cholesky give's L'*L = A, we need to double up
    result *= 2;
    
    
#if debug_logDet_IplusDiagzW4 
    printf( "* * * Adiag:\n");
    for(int i=0; i<N; i++)
        printf( "%f\t", A[i*N+i]);
    printf( "Adiag(end)\n");
    cout << "result (log det) = " << result << endl;
#endif
    
    
 
    
    
    return result;
    
};





double DataModel ::  logDet_IplusDiagzW_old_LU( const IloNumArray &valz , const vector< vector < double > > & Sigma )
{
    double result = 0;
//    FILE *out;
    vector< double> A( n*n );
    
    vector< double> valzSQRT (n);
    for (int i = 0; i < n; i++) {
        if (valz[i] <= 0)
            valzSQRT[i] = 0;
        else
            valzSQRT[i] = sqrt(valz[i]);
    }
    
    int N, LDA, IPIV[ _CONST_INT_MAX_FACTOR_SIZE ], INFO;
    N = n;
    LDA = N;
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                A[i*n+j] = 1 + ( Sigma[i][j] - 1 ) * valzSQRT[i] * valzSQRT[j];
            } else {
                A[i*n+j] = Sigma[i][j] * valzSQRT[i] * valzSQRT[j];
            }
        }
    };
    
    // here A is symmetric
    
#if debug_logDet_IplusDiagzW
    cout << " > > in {DataModel::logDet_IplusDiagzW} " << endl;
    cout << " * diaganol of the matrix = [ " << endl;
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            cout << A[i*N+j] << ", " ;
        }
        cout <<endl;
    }
    cout << " ]"<<  endl;
#endif

    
    
//    out = stderr;
#if debug_logDet_IplusDiagzW2
    cout <<" ( tag fqo4f89q34980q ) " << endl;
    printf("A before call:\n");
    printf("N = %d\n", N);
#endif
#if debug_logDet_IplusDiagzW2
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            printf( "%f\t", A[i+N*j]);
        printf("\n");
    }
#endif

    dgetrf_(&N, &N, &A[0], &LDA, IPIV, &INFO);

    
    // use dpbtrf
    
    
#if debug_logDet_IplusDiagzW2
    printf("INFO: %d\n", INFO);
    cout << " N = " << N << endl;
    printf("A after call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            printf( "%f\t", A[i+N*j]);
        printf( "\n");
    }
#endif
    
#if debug_logDet_IplusDiagzW2
    printf( "VECR:\n");
#endif
    double piii=1;
    for (int i=0; i<N; i++) {
#if debug_logDet_IplusDiagzW2
        printf( "%f\t", A[i+N*i]);
#endif
        piii=piii*A[i+N*i];
        result += log(A[i*N+i]);
    }
#if debug_logDet_IplusDiagzW2
    printf( "%f\t", result);
    
    printf( "%f\n",piii);
    printf( "P:\n");
    for(int j=0; j<N; j++)
        printf( "%d\t", IPIV[j]);
    printf( "P(end)\n");
    
#endif
#if debug_logDet_IplusDiagzW3
    printf( "* * * Adiag:\n");
    for(int i=0; i<N; i++)
        printf( "%f\t", A[i*N+i]);
    printf( "Adiag(end)\n");
    
#endif
    
#if debugMode_pause
    {
        int temp;
        cin >> temp;
    }
#endif
    
    
#if debug_logDet_IplusDiagzW 
    cout << " > > in {DataModel::logDet_IplusDiagzW} " << endl;
    cout << " * diaganol of the matrix = [ " ;
    for (int i = 0; i < N; i++) {
        cout << A[i*N+i] << ", " ;
    }
    cout << " ]" <<  endl;
    cout << " * valzSQRT = [ " ;
    for (int i = 0; i < N; i++) {
        cout << valzSQRT[i] << ", " ;
    }
    cout << " ]" <<  endl;
    
    
#endif
    
    
    return result;
    
};




vector< double> DataModel :: R2uMinusLGradient( const IloNumArray &valz, double r, const vector<vector< double>> &L, const vector<vector< double>> &U )
{
    vector< double> grad(n);
    
    vector< double> gradL(n);
    vector< double> gradU(n);
    vector< double> IplusLZ(n*n);
    vector< double> IplusUZ(n*n);
    
    /*   comput for L */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                IplusLZ[i*n+j] = 1 + ( L[i][j] - 1 ) * valz[i];
            } else {
                IplusLZ[i*n+j] = L[i][j] * valz[i];
            }
        }
    };
#if debug_uMinusLGradient
    cout << "IplusLZ = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusLZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    inverse( &IplusLZ[0], n);
#if debug_uMinusLGradient
    cout << "IplusLZ_inverse = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusLZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    
    for (int i = 0; i < n; i++){
        double currentDiag = 0;
        for (int j = 0; j < n; j++) {
            //            currentDiag += L[i][j] * IplusLZ[j*n+i];
            currentDiag += L[j][i] * IplusLZ[j*n+i];
        }
        currentDiag -= IplusLZ[i*n+i];
        gradL[i] = currentDiag;
    };
    
    
    
#if debug_uMinusLGradient
    cout << "U = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << U[i][j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    
    /*   comput for U */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                IplusUZ[i*n+j] = 1 + ( U[i][j] - 1 ) * valz[i];
            } else {
                IplusUZ[i*n+j] = U[i][j] * valz[i];
            }
        }
    };
#if debug_uMinusLGradient
    cout << "IplusUZ = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusUZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    inverse( &IplusUZ[0], n);
#if debug_uMinusLGradient
    cout << "IplusUZ_inverse = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusUZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    for (int i = 0; i < n; i++){
        double currentDiag = 0;
        for (int j = 0; j < n; j++) {
            currentDiag += U[j][i] * IplusUZ[j*n+i];
        }
        currentDiag -= IplusUZ[i*n+i];
        gradU[i] = currentDiag;
    };
    
    
    
    
    
    
    double r2 = r*r;
    
    /*   comput for L - U */
    for (int i = 0; i < n; i++)
        grad[i] = r2 * gradL[i] - gradU[i];
    
#if debug_uMinusLGradient
    // pinrt out something
    cout << "valz = [ ";
    for (int i = 0; i < n; i++)
        cout << valz[i] << ", ";
    cout << " ]\n";
    cout << "gradL = [ ";
    for (int i = 0; i < n; i++)
        cout << gradL[i] << ", ";
    cout << " ]\n";
    cout << "gradU = [ ";
    for (int i = 0; i < n; i++)
        cout << gradU[i] << ", ";
    cout << " ]\n";
    cout << "grad = [ ";
    for (int i = 0; i < n; i++)
        cout << grad[i] << ", ";
    cout << " ]\n";
#endif
    
    
    
    
    gradL.end();
    gradU.end();
    IplusLZ.end();
    IplusUZ.end();
    
    
    
    
    return grad;
};

vector< double> DataModel :: uMinusLGradient( const IloNumArray &valz)
{
    return uMinusLGradient(valz, L, U);
};

vector< double> DataModel :: uMinusLGradient( const IloNumArray &valz, const vector< vector < double > > &L, const vector< vector < double > > &U)
{
    vector< double> grad(n);
    
    vector< double> gradL(n);
    vector< double> gradU(n);
    vector< double> IplusLZ(n*n);
    vector< double> IplusUZ(n*n);
    
    /*   comput for L */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                IplusLZ[i*n+j] = 1 + ( L[i][j] - 1 ) * valz[i];
            } else {
                IplusLZ[i*n+j] = L[i][j] * valz[i];
            }
        }
    };
#if debug_uMinusLGradient
    cout << "IplusLZ = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusLZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    inverse( &IplusLZ[0], n);
#if debug_uMinusLGradient
    cout << "IplusLZ_inverse = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusLZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    
    for (int i = 0; i < n; i++){
        double currentDiag = 0;
        for (int j = 0; j < n; j++) {
            //            currentDiag += L[i][j] * IplusLZ[j*n+i];
            currentDiag += L[j][i] * IplusLZ[j*n+i];
        }
        currentDiag -= IplusLZ[i*n+i];
        gradL[i] = currentDiag;
    };
    
    
    
#if debug_uMinusLGradient
    cout << "U = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << U[i][j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    
    /*   comput for U */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++){
            if (j == i) {
                IplusUZ[i*n+j] = 1 + ( U[i][j] - 1 ) * valz[i];
            } else {
                IplusUZ[i*n+j] = U[i][j] * valz[i];
            }
        }
    };
#if debug_uMinusLGradient
    cout << "IplusUZ = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusUZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    inverse( &IplusUZ[0], n);
#if debug_uMinusLGradient
    cout << "IplusUZ_inverse = \n[ ";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            cout << IplusUZ[i*n+j] << ", " ;
        }
        cout << endl;
    }
    cout << " ]\n";
#endif
    for (int i = 0; i < n; i++){
        double currentDiag = 0;
        for (int j = 0; j < n; j++) {
            currentDiag += U[j][i] * IplusUZ[j*n+i];
        }
        currentDiag -= IplusUZ[i*n+i];
        gradU[i] = currentDiag;
    };
    
    
    
    
    
    
    
    /*   comput for L - U */
    for (int i = 0; i < n; i++)
        grad[i] = gradL[i] - gradU[i];
    
#if debug_uMinusLGradient
    // pinrt out something
    cout << "valz = [ ";
    for (int i = 0; i < n; i++)
        cout << valz[i] << ", ";
    cout << " ]\n";
    cout << "gradL = [ ";
    for (int i = 0; i < n; i++)
        cout << gradL[i] << ", ";
    cout << " ]\n";
    cout << "gradU = [ ";
    for (int i = 0; i < n; i++)
        cout << gradU[i] << ", ";
    cout << " ]\n";
    cout << "grad = [ ";
    for (int i = 0; i < n; i++)
        cout << grad[i] << ", ";
    cout << " ]\n";
#endif
    
    
    
    
    gradL.end();
    gradU.end();
    IplusLZ.end();
    IplusUZ.end();
    
    
    
    
    return grad;
};


vector < double > DataModel :: gradient( const IloNumArray &valz )
{
    return gradient(valz, Sigma);
};



vector < double > DataModel :: gradient( const IloNumArray &valz , const vector< vector < double > > & Sigma)
{
    vector< double> grad(n);
//    vector< vector< double> > IplusDiagxW(n, vector< double>(n));
    vector< double> IplusDiagxW(n*n);

    
    
#if SIGN_output_info || debug_gradient
    cout << endl << "Sigma = (tag fq3bw45bgw45)" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << Sigma[i][j] << "\t";
        cout << endl;
    };
    cout << endl << "valz = (tag 45wbgw45bw45b )" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << valz[i] << "\t";
    };
    cout << endl;
#endif

    
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            if (j == i) {
                IplusDiagxW[i*n+j] = 1 + ( Sigma[i][j] - 1 ) * valz[i];
            } else {
                IplusDiagxW[i*n+j] = Sigma[i][j] * valz[i];
            }
        }
    };
#if SIGN_output_info || debug_gradient
    cout << endl << "before computing is: (IplusDiagxW) (tag 345vg64bheh57n)" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << IplusDiagxW[i*n+j] << "\t";
        cout << endl;
    };
    cout << endl;
#endif

    inverse( &IplusDiagxW[0], n);
    
    
#if SIGN_output_info || debug_gradient
    cout << endl << "The computed inverse is: (IplusDiagxW) (tag 0243f034mjc83)" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << IplusDiagxW[i*n+j] << "\t";
        cout << endl;
    };
    cout << endl;
#endif
    for (int i = 0; i < n; i++){
        double currentDiag = 0;
        for (int j = 0; j < n; j++) {
            currentDiag += Sigma[i][j] * IplusDiagxW[j*n+i];
        }
        currentDiag -= IplusDiagxW[i*n+i];
        grad[i] = currentDiag;
    };

#if debugMode_pause || debug_gradient
    {
        int temp;
        cin >> temp;
    }
#endif
    
    return grad;
};


vector < double > DataModel :: gradientLDYplusM( const vector < double > &Y , const vector < vector < double > > &M)
{
    vector< double> grad(n);
    vector< double> matrix(n*n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            if (i == j)
                matrix[i*n+j] = M[i][j] + Y[i];
            else
                matrix[i*n+j] = M[i][j];
        }
    }
    inverse( &matrix[0], n);
    
    for (int i = 0; i < n; i++){
        grad[i] = matrix[i*n+i];
    };
    return grad;
};



vector < vector < double > > DataModel :: KUKplusJLJ( const vector < double > &K , const vector < double > &J , const vector < vector < double > > &U, const vector < vector < double > > &L )
{
    unsigned long n = U.size();
    vector < vector < double > > result(n,vector<double>(n));
    for (int i = 0; i < n; i++) {
        for(int j = i; j < n; j++){
            result[i][j] = U[i][j]*K[i]*K[j]+L[i][j]*J[i]*J[j];
            if (j != i)
                result[j][i] = result[i][j];
        }
    }
    return result;
};
vector < vector < double > > DataModel :: KUKminusJLJ( const vector < double > &K , const vector < double > &J , const vector < vector < double > > &U, const vector < vector < double > > &L )
{
    unsigned long n = U.size();
    vector < vector < double > > result(n,vector<double>(n));
    for (int i = 0; i < n; i++) {
        for(int j = i; j < n; j++){
            result[i][j] = U[i][j]*K[i]*K[j]-L[i][j]*J[i]*J[j];
            if (j != i)
                result[j][i] = result[i][j];
        }
    }
    return result;
};

vector < vector < double > > DataModel :: KUK( const vector < double > &K ,const vector < vector < double > > &U)
{
    unsigned long n = U.size();
    vector < vector < double > > result(n,vector<double>(n));
    for (int i = 0; i < n; i++) {
        for(int j = i; j < n; j++){
            result[i][j] = U[i][j]*K[i]*K[j];
            if (j != i)
                result[j][i] = result[i][j];
        }
    }
    return result;
};


vector < vector < double > > DataModel :: sumMatrices( const vector < vector < double > > &A ,const vector < vector < double > > &B)
{
    vector < vector < double > > result(A);
    for (int i = 0; i < A.size(); i++) {
        for(int j = i; j < A.size(); j++){
            result[i][j] += B[i][j];
            if (j != i)
                result[j][i] = result[i][j];
        }
    }
    return A;
};









double DataModel :: logDet_AplusX( const vector < double > &x , const vector< vector < double > > &matrix )
{
    double result = 0;
    vector< double> A( n*n );
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++){
            if (j == i) {
                A[i*n+j] = matrix[i][j]+ x[i];
            } else {
                A[i*n+j] = matrix[i][j];
            }
        }
    
    //    int N, LDA, IPIV[ _CONST_INT_MAX_FACTOR_SIZE ], INFO;
    int N, LDA, INFO;
    int M;
    M = n;
    N = n;
    LDA = N;
    
    
    // here A is symmetric
    
    //
    //
    //
    //
    //
    
    FILE *out;
    
    
    out = stderr;
#if SIGN_output_info
    fprintf(out, "A before call:\n");
    fprintf(out, "N = %d\n", N);
#endif
    
#if SIGN_output_info
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
#if debug_logDet_IplusDiagzW4
    fprintf(out, "A before call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
    
    char uplo = 'L';
    
#if SIGN_output_info
    cout << " before dpotrf_ :" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    //void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
    dpotrf_(&uplo, &N, &A[0], &LDA, &INFO);
    //    dgetrf_(&M, &N, A, &LDA, &IPIV[0], &INFO);
#if SIGN_output_info
    cout << " after dpotrf_ :" << endl;
    cout << "uplo = " << uplo << endl;
    cout << "N = " << N << endl;
    cout << "lda = " << LDA << endl;
    cout << "info = " << INFO   << endl;
#endif
    
#if SIGN_output_info
    fprintf(out, "INFO: %d\n", INFO);
    cout << " N = " << N <<  endl;
    fprintf(out, "A after call:\n");
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++)
            fprintf(out, "%f\t", A[i+N*j]);
        fprintf(out, "\n");
    }
#endif
    
#if SIGN_output_info
    fprintf(out, "VECR:\n");
#endif
    double piii=1;
    for (int i = 0; i < N; i++) {
#if SIGN_output_info
        fprintf(out, "%f\t", A[i+N*i]);
#endif
        piii = piii*A[i+N*i];
        result += log(A[i*N+i]);
        /* !!! why negative diagonal */
    }
    
    // cholesky give's L'*L = A, we need to double up
    result *= 2;
    
    
#if debug_logDet_IplusDiagzW4
    printf( "* * * Adiag:\n");
    for(int i=0; i<N; i++)
        printf( "%f\t", A[i*N+i]);
    printf( "Adiag(end)\n");
    cout << "result (log det) = " << result << endl;
#endif
    
    
    
    
    
    return result;

};








// added 10/27/2014
double DataModel :: minEigAsqrtinvBAsqrtinv(const vector < vector < double > > &A ,const vector < vector < double > > &B)
{
    // input A , B
    // aA <= B
    // find a = A^-.5 * B * A^-.5
    
#if function_title
    cout << " at the begining of DataModel :: minEigAsqrtinvBAsqrtinv (tag fw4vbtdfva4e5svwv4)" << endl;
#endif
//    
//    LatNegSemiDef.resize(n);
//    for (int i = 0; i < n; i++)
//        LatNegSemiDef[i].resize(n);
//    UatNegSemiDef.resize(n);
//    for (int i = 0; i < n; i++)
//        UatNegSemiDef[i].resize(n);
//    scalerUatNSD.resize(n);
//    scalerLatNSD.resize(n);
//    
//    // L first
//    for (int iRow = 0; iRow < n; iRow++)
//        for (int iCol = 0; iCol < n; iCol++)
//            LatNegSemiDef[iRow][iCol] = SigmaOriginal[iRow][iCol] - SigmaOriginal[iRow][n]*SigmaOriginal[n][iCol] / SigmaOriginal[n][n] ;
//    
//    
//    
//    for (int i = 0; i < n; i++) {
//        scalerLatNSD[i] = sqrt(LatNegSemiDef[i][i]);
//    }
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++){
//            LatNegSemiDef[i][j] /= scalerLatNSD[i];
//            LatNegSemiDef[j][i] /= scalerLatNSD[i];
//        }
//    }    // to here, diag(LatNegSemiDef) == 1 1 1 1 1 1 1 1
//    
//    
//    
//    
//    vector < double > eigenTempL(eig(LatNegSemiDef));
//    
//    double maxEigValueTempL = eigenTempL[n-1];
//    for (int i = 0; i < n; i++) {
//        scalerLatNSD[i] *= sqrt(maxEigValueTempL);
//        for (int j = 0; j < n; j++){
//            LatNegSemiDef[i][j] /= maxEigValueTempL;
//        }
//    }    // to here, max eig(LatNegSemiDef) == 1, diagonal the same (try my best to maintain good condition number
//    for (int i = 0; i < n; i++)
//        scalerLatNSD[i] = 1.0 / scalerLatNSD[i];
//
    
    
    // U here
    int N = n;
    vector< double> sqrtA(N * N);
    int indexMatrix = 0;
    for (int iRow = 0; iRow < n; iRow++){
        for (int iCol = 0; iCol < n; iCol++){
            sqrtA[indexMatrix] = A[iRow][iCol];
            indexMatrix++;
        }
    }
    
    FILE *out;
    int LDA, INFO;
    LDA = N;
    //    int M, N, LDA, INFO;
    out = stderr;
    char uplo = 'L';
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    //void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
    //    dpotrf_(&uplo, &N, A, &LDA, &INFO);
    dpotrf_(&uplo, &N, &sqrtA[0], &LDA, &INFO);
    //    dgetrf_(&M, &N, A, &LDA, &IPIV[0], &INFO);
    
    //    dtrtri_
    //    void dtrtri_(char *uplo, char*diag, int *n, double *a, int *lda, int *info);
    char diag = 'N';
    dtrtri_(&uplo, &diag, &N, &sqrtA[0], &LDA, &INFO);
    
    vector< vector < double > > RHS1(n, vector < double >(n)) ;
    // A'^-.5 * B * A^-.5
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // assign tempU[i][j]
            for (int k = 0; k <= j; k++) {
                // tempU[i][j] += LatNegSemiDef[i][k]*sqrtU[k+N*j]; // ik * kj
                RHS1[i][j] += B[i][k]*sqrtA[j+N*k]; // ik * kj
            }
        }
    }
    vector< vector < double > > RHS2(n, vector < double >(n)) ;
//    vector< vector < double > > tempU2 ; // L'^-.5 * LatNSD * L^-.5
//    tempU2.resize(n);
//    for (int i = 0; i < n; i++)
//        tempU2[i].resize(n);
//    for (int i = 0; i < n; i++)
//        for (int j = 0; j < n; j++)
//            tempU2[i][j] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) { // assign tempU[i][j]
            for (int k = 0; k <= i; k++) {
                // tempU2[i][j] += sqrtU[k+N*i] * tempU[k][j]; // ik * kj
                RHS2[i][j] += sqrtA[i+N*k] * RHS1[k][j]; // ik * kj
            }
        }
    }
    
    vector < double > eigenRHS(eig(RHS2));
    double tempScaler = sqrt( eigenRHS[0] );
    return tempScaler;
    
};






// added 10/29/2014
void DataModel :: newtonDirectionPSDSeparationProblem( vector<double> &dir,const vector<double> &x,const vector<vector<double>> &A,const vector<vector<double>> &B,const vector<vector<double>> &C, const double mu)
{
    vector < double > inverseXplusA(n*n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                inverseXplusA[i*n+i] = A[i][i]+x[i];
            else
                inverseXplusA[i*n+j] = A[i][j];
#if debug_newton
    {
    cout << " (tag faw98fajw)\ninverseXplusA = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseXplusA[i*n+j] << "  " ;
        cout << endl;
    }
    }
#endif
    inversePSD(inverseXplusA, n);
#if debug_newton
    {
        cout << " (tag hsetge4fwae)\ninverseXplusA = " << endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << inverseXplusA[i*n+j] << "  " ;
            cout << endl;
        }
    }
#endif

    vector < double > inverseXplusB(n*n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                inverseXplusB[i*n+i] = B[i][i]+x[i];
            else
                inverseXplusB[i*n+j] = B[i][j];
#if debug_newton
    cout << " (tag wa4geshe5rh)\nnverseXplusB = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseXplusB[i*n+j] << "  " ;
        cout << endl;
    }
#endif
    inversePSD(inverseXplusB, n);
#if debug_newton
    cout << " (tag hestherswatfwa)\nnverseXplusB = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseXplusB[i*n+j] << "  " ;
        cout << endl;
    }
#endif

    vector < double > inverseCminusX(n*n);
    
#if debug_newton
    cout << " (tag bxdgbserhryd)\ninverseCminusX = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseCminusX[i*n+j] << "  " ;
        cout << endl;
    }
    cout << " (tag awfwa4)\n x = " << endl;
    for (int i = 0; i < n; i++){
        cout << x[i] << "  " ;
    }
    cout << endl;

#endif
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i == j)
                inverseCminusX[i*n+i] = C[i][i]-x[i];
            else
                inverseCminusX[i*n+j] = C[i][j];
#if debug_newton
    cout << " (tag bxdgbserhryd)\ninverseCminusX = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseCminusX[i*n+j] << "  " ;
        cout << endl;
    }
#endif
    inversePSD(inverseCminusX, n);
#if debug_newton
    cout << " (tag nrtherses)\ninverseCminusX = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseCminusX[i*n+j] << "  " ;
        cout << endl;
    }
#endif

    
    
    vector < double > gradient(n);
    vector < double > hessian(n*n);
    
    
    
    
    
    
#if debug_newton
    cout << " (tag fwgege5h56u78)\nnverseXplusA = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseXplusA[i*n+j] << "  " ;
        cout << endl;
    }
    cout << " (tag kutikyftbhr)\nnverseXplusB = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << inverseXplusB[i*n+j] << "  " ;
        cout << endl;
    }
#endif

    
    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++){
            if (i == j){
                int ini = i*n+i;
                gradient[i] = inverseXplusA[ini] - inverseXplusB[ini] + mu * ( -inverseCminusX[ini]+1.0/x[i] );
                hessian[ini] = - inverseXplusA[ini]*inverseXplusA[ini] +inverseXplusB[ini]*inverseXplusB[ini] -mu*( inverseCminusX[ini]*inverseCminusX[ini]+1.0/x[i]/x[j] );

            } else {
                int inj = i*n+j;
                hessian[inj] = - inverseXplusA[inj]*inverseXplusA[inj] +inverseXplusB[inj]*inverseXplusB[inj] -mu*( inverseCminusX[inj]*inverseCminusX[inj] );
                hessian[j*n+i] = hessian[inj];
            }
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            hessian[i*n+j] = - hessian[i*n+j];
#if debug_newton
    cout << " (tag sgesrhrshsr)\ngradient = " << endl;
    for (int i = 0; i < n; i++){
        cout << gradient[i] << "  " ;
    }
    cout << endl;
    
    cout << " hessian = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << hessian[i*n+j] << "  " ;
        cout << endl;
    }
    
    cout << " 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! " << endl;
    {
        int temp;
        cin >> temp;
        
    }
#endif

#if debug_newton_level2
    vector<vector<double>> hessianVV(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            hessianVV[i][j] = hessian[i*n+j];
        }
    }
    cout << " eig(hessianVV) = " << endl;
    for (int i = 0; i < n; i++) {
        cout << eig(hessianVV)[i] << ", " ;
    }
    cout << endl;
    if (eig(hessianVV)[0] < 1e-8){
        cout << " 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! 检查一下这里 ! ! " << endl;
        cout << " eig(hessianVV)[0] = " << eig(hessianVV)[0] << endl;
        {
            int temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            cin >> temp;
            
        }

    }
#endif
    
            
            
    int LDA, LDB, INFO;
    int N = n;
    int NRHS = 1;
    LDA = N;
    LDB = N;
    
    //    int M, N, LDA, INFO;
    char uplo = 'L';
    //Computes the Cholesky factorization of a symmetric positive definite matrix.
    //    dpotrf_(&uplo, &N,  A   , &LDA, &INFO);
#if debug_newton_level2
    cout << " (tag fawefaw4gaw4) before before partition \n gradient = " << endl;
    for (int i = 0; i < n; i++){
        cout << gradient[i] << "  " ;
    }
    cout << endl;
    
    cout << " hessian = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << hessian[i*n+j] << "  " ;
        cout << endl;
    }
#endif
    dpotrf_(&uplo, &N, &hessian[0], &LDA, &INFO);
    
#if debug_newton_level4
//    if (1) { // (INFO == 0){
//        cout << "INFO = " << INFO <<"   ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== "<< endl;
////        int temp;
////        cin >> temp;
//        system("pause");
//    }
#endif
    
#if debug_newton_level2
    cout << " (tag hdrtyhrd6ns5e) during problem  problem \n gradient = " << endl;
    for (int i = 0; i < n; i++){
        cout << gradient[i] << "  " ;
    }
    cout << endl;
    
    cout << " hessian = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << hessian[i*n+j] << "  " ;
        cout << endl;
    }
#endif
    

    //    dpotrs_(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, &info);
    dpotrs_(&uplo, &N, &NRHS, &hessian[0], &LDA, &gradient[0], &LDB, &INFO);
    // ! ! ! ! ! ! ! ! 问题可能出在这里 ! ! ! ! ! ! ! ! hessian 可能不正定
    
    

    
#if debug_newton_level2
    cout << " (tag fe5bd6rndr6) after solve the problem \n solution = " << endl;
    for (int i = 0; i < n; i++){
        cout << gradient[i] << "  " ;
    }
    cout << endl;
    
    cout << " hessian = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
            cout << hessian[i*n+j] << "  " ;
        cout << endl;
    }
#endif

#if debug_newton
    cout << " (tag zsrger)\n x = " << endl;
    for (int i = 0; i < n; i++){
            cout << x[i] << "  " ;
    }        cout << endl;

#endif

    
//    for (int i = 0; i < n; i++)
//        gradient[i] = - gradient[i];

    
//    vector < double > direction(n);

//    void inverse(double* A, int N);

    
#if debug_newton_level2
    cout << " (tag f4q3t34bs5ehn) finishing newtonDirectionPSDSeparationProblem  @ @ @ @ @ @ @ @ @ @ @ " << endl;
    {
        int temp;
        cin >> temp;
    }
#endif

    
    
    dir.resize(n);
    for (int i = 0; i < n; i++)
        dir[i] = gradient[i];
//    
//    vector<double> result(gradient);
//    return result;
};
































