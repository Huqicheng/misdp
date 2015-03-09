//
//  DataModel.h
//  SparseRegression
//
//  Created by Qi on 5/27/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __SparseRegression__DataModel__
#define __SparseRegression__DataModel__

#include <iostream>
#include <bitset>
#include <vector>
#include <cmath>
#include <algorithm>    // std::min

#include <unordered_map>
#include "SystemParameters.h"
#include <ilcplex/ilocplex.h>
#include "TuneParameters.h"



using namespace std;

class DataModel{
public:
    
    void push_back( const vector < double > &a );
    void push_back_X( const vector < double > &a );
    void push_back_y( const double a );
    void setCols(int & a);
    void setRows(int & a);
    int getRows() const;
    int getCols() const;
    
    void initializeSigma();
    void scalizeSigma();
    
//    bitset<_CONST_INT_MAX_FACTOR_SIZE> at(const int &i) const;

    double logDetSubmatrix( const bitset<_CONST_INT_MAX_FACTOR_SIZE> & setsIndicators ); // log det (X_SS)
//    double logDetSubmatrix_real( const bitset<_CONST_INT_MAX_FACTOR_SIZE> & setsIndicators ); // log det (X_SS)
    double logDetSubmatrix_with_y( const bitset<_CONST_INT_MAX_FACTOR_SIZE> & setsIndicators ); // log det (X_SS)
    double logDetDelta( const bitset<_CONST_INT_MAX_FACTOR_SIZE> & i, const bitset<_CONST_INT_MAX_FACTOR_SIZE> & S ); // delta_i(S) = f(S\cup i)-f(S\i)
    double logDetDelta( const int i, bitset<_CONST_INT_MAX_FACTOR_SIZE> S ); // delta_i(S) = f(S\cup i)-f(S\i)
//    double * gradient();
    
    vector < double > gradient( const IloNumArray &valz ) ;
    vector < double > gradient( const IloNumArray &valz , const vector < vector < double > > & Sigma) ;
    vector < double > gradientLDYplusM( const vector < double > &Y , const vector < vector < double > > &M ) ;


    
    vector < double > uMinusLGradient( const IloNumArray &valz ) ;
    vector < double > uMinusLGradient( const IloNumArray &valz, const vector<vector<double>> &L, const vector<vector<double>> &U ) ;
    vector < double > R2uMinusLGradient( const IloNumArray &valz, double r, const vector <vector < double>> &L, const vector <vector < double>> &U ) ;
    void inverse(double* A, int N);
    void inverse( vector<double> & A, int N);
    void inversePSD( vector<double> & A, int N);
    double logDet_IplusDiagzW( const IloNumArray &valz , const vector < vector < double > > & Sigma ) ;
    double logDet_IplusDiagzW( const vector < double > & valz , const vector < vector < double > > & Sigma ) ;
    double logDet_IplusDiagzW_old_LU( const IloNumArray &valz , const vector < vector < double > > & Sigma ) ;
    double logDet_IplusDiagzW( const IloNumArray &valz ) ;
    
    
    vector < double> eig ( const vector < vector < double > > & matrix) ;
    
// private:
//    string dataType;
    int n; // #cols, number of features,
    int m; // #rows, number of observations
    unordered_map< bitset<_CONST_INT_MAX_FACTOR_SIZE>, double > calculatedValues;
    
    
    
    vector < vector < double > > SigmaOriginal ; // Sigma = [X y]'[X y]
    vector < vector < double > > Sigma ; // Sigma = [X y]'[X y]
    vector < vector < double > > U ; // Sigma = X'X
    vector < vector < double > > L ; // Sigma = X'X - 1/yy *
    vector < double > eta ;
    double sigma_yy ;
    
    
    double maxLambda ;
    double maxEigenvalue, minEigenvalue;
    vector < vector < double > > X ;
    vector < double > y ;
    
    
    
    vector < vector < double > > UMaxEigen1 ; // U = X'X
    vector < vector < double > > LMaxEigen1 ; // L = X'X - 1/yy * ...
    double maxEigenU_ieSigma; // compare with U
    double rStar;
    void initializerStar();
    
    
    
    vector < vector < double > > OriginalU ; // Sigma = [X y]'[X y] -> X'X    // *finished*
    vector < vector < double > > OriginalL ; // Sigma = [X y]'[X y] -> X'X - 1/yy *    // *finished*
    vector < vector < double > > OriginalUinverse ; // Sigma = [X y]'[X y] -> X'X    // *finished*
    vector < vector < double > > OriginalLinverse ; // Sigma = [X y]'[X y] -> X'X - 1/yy *    // *finished*
    vector < vector < double > > OriginalInvOfInvSum ; // Sigma = [X y]'[X y] -> X'X - 1/yy *    // *finished*
    void initializeOriginalInvOfInvSum(const TuneParameters &tune);
    vector < vector < double > > scaleBothSide(const vector < vector < double > > &A, const vector < double > &s);
    
    
    void initializeUandLatNSD(const TuneParameters &tune);

//    vector < bitset<_CONST_INT_MAX_FACTOR_SIZE> > bitData;
    
    
    
    // obj = LD(I + (%L%-I) X) - 2sum log(sL)X  -  LD(I + (%U%-I) X) + 2sum log(sU)X
    vector < vector < double > > UatNegSemiDef ;
    vector < vector < double > > LatNegSemiDef ;
    vector < double > scalerUatNSD ;
    vector < double > scalerLatNSD ;

    
    
    
    bool flag_if_sigma_initialized;
    
    
    
    // some basic matrix function
    vector < vector < double > > KUKplusJLJ( const vector < double > &K , const vector < double > &J , const vector < vector < double > > &U, const vector < vector < double > > &L ) ;
    vector < vector < double > > KUKminusJLJ( const vector < double > &K , const vector < double > &J , const vector < vector < double > > &U, const vector < vector < double > > &L ) ;
    vector < vector < double > > KUK( const vector < double > &K ,const vector < vector < double > > &U);
    vector < vector < double > > sumMatrices( const vector < vector < double > > &A ,const vector < vector < double > > &B);
    
    
    double logDet_AplusX( const vector < double > &x , const vector< vector < double > > &A );

    double minEigAsqrtinvBAsqrtinv(const vector < vector < double > > &A ,const vector < vector < double > > &B);
    
    
    
    
    // added 10/29/2014
    void newtonDirectionPSDSeparationProblem(vector<double> &dir,const vector<double> &x,const vector<vector<double>> &A,const vector<vector<double>> &B,const vector<vector<double>> &C, const double mu);

};



#endif /* defined(__SparseRegression__DataModel__) */














































