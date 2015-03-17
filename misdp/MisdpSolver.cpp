//
//  MisdpSolver.cpp
//  misdp
//
//  Created by Qi Zhang on 3/7/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//
#define title_everywhere 1
#define info_level1 1
#define test_input 1
#define test_building 1
#define test_build_initCons 0
#define test_lazy1 1
#define test_lazy2 1

#include "MisdpSolver.h"



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
    
    /* DSYEVX prototype */
    extern void dsyevx( char* jobz, char* range, char* uplo, int* n, double* a,
                       int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
                       int* m, double* w, double* z, int* ldz, double* work, int* lwork,
                       int* iwork, int* ifail, int* info );
    /* Auxiliary routines prototypes */
    extern void print_matrix( char* desc, int m, int n, double* a, int lda );

}


/*
     min  <C_0, X>
     s.t. <A_i, X> \geq b_i, i = 1...m
          X >= 0
 */


ILOLAZYCONSTRAINTCALLBACK3(myLazy, IloNumVarArray,x, MisdpModel &,misdpmodel, MisdpTuneParameters &, tune)
{
    //    if (getNnodes64() > 10)
    //        return;
#if title_everywhere
    cout << " ILOLAZYCONSTRAINTCALLBACK3 tag (fq34345wgw345) " << endl;
#endif
    
    if (getNnodes64() == 0){
        tune.rootBestValue = getBestObjValue();
    };
    
    // want to find negative eigen-value (vector) of X = [x]
    
    IloEnv env = getEnv();
    IloNumArray valx_ilo(env);
    try{
#if test_lazy2
        double currObjValue = getObjValue();
        cout << " currObjValue = " << currObjValue << endl;
#endif
        getValues( valx_ilo, x );
        int n = misdpmodel.n;

        
        vector<double> valx(n*n);
        for (int i = 0; i < n*n; i++) {
            valx[i] = valx_ilo[i];
        }

        //    #define N 5
        //    #define NSELECT 2
        //    #define LDA N
        //    #define LDZ N
        //    int N = 5;
        //    int NSELECT = 3;
        //    int LDA, LDZ;
        //    LDA = N; LDZ = N;
        //
        
        /* Main program */
        /* Locals */
        int nSelect = n;
        
        //    int n = N, il, iu, m, lda = LDA, ldz = LDZ, info, lwork;
        int il, iu, m, lda = n, ldz = n, info, lwork;
        double abstol, vl, vu;
        double wkopt;
        double* work;
        /* Local arrays */
        /* iwork dimension should be at least 5*n */
        //    int iwork[5*N], ifail[N];
        vector<int> iwork(5*n, 0);
        vector<int> ifail(n, 0);
        //    double w[N], z[LDZ*NSELECT];
        vector<double> w(n,0);
        vector<double> z(ldz*nSelect, 0);
        //    double a[LDA*N] = {
        //        6.29,  0.00,  0.00,  0.00,  0.00,
        //        -0.39,  7.19,  0.00,  0.00,  0.00,
        //        0.61,  0.81,  5.48,  0.00,  0.00,
        //        1.18,  1.19, -3.13,  3.79,  0.00,
        //        -0.08, -0.08,  0.22, -0.26,  0.83
        //    };
        //        vector<double> a(n*n, 0);
        //        a  --> & valx[0]
        
        /* Executable statements */
        if (n < 10)
            printf( " DSYEVX Example Program Results\n" );
        
        /* Print original matrix */
        if (n < 10){
            cout << valx_ilo << endl;
            print_matrix( "original matrix (before)", n, n, &valx[0], n );
        }

        
        /* Negative abstol means using the default value */
        abstol = -1.0;
        /* Set il, iu to compute NSELECT smallest eigenvalues */
        il = 1;
        iu = nSelect;
        /* Query and allocate the optimal workspace */
        lwork = -1;
        dsyevx( "Vectors", "Indices", "Upper", &n, &valx[0], &lda, &vl, &vu, &il, &iu,
               &abstol, &m, &w[0], &z[0], &ldz, &wkopt, &lwork, &iwork[0], &ifail[0], &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Solve eigenproblem */
        dsyevx( "Vectors", "Indices", "Upper", &n, &valx[0], &lda, &vl, &vu, &il, &iu,
               &abstol, &m, &w[0], &z[0], &ldz, &work[0], &lwork, &iwork[0], &ifail[0], &info );
        /* Check for convergence */
        if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
        }
        /* Print the number of eigenvalues found */
        printf( "\n The total number of eigenvalues found:%2i\n", m );
        /* Print original matrix */
        if (n < 10){
            cout << valx_ilo << endl;
            print_matrix( "original matrix (after)", n, n, &valx[0], n );
        }
        /* Print eigenvalues */
        print_matrix( "Selected eigenvalues", 1, m, &w[0], 1 );
        /* Print eigenvectors */
        if (n < 10)
            print_matrix( "Selected eigenvectors (stored columnwise)", n, m, &z[0], ldz );
        /* Free workspace */
//        free( (void*)work );
        
        
        
#if test_lazy1
        cout << " ready to build cuts2BeAdded (tag of9wj90)" << endl;
#endif
        
        IloConstraintArray cuts2BeAdded(env);
        nSelect = 1;
        for (int iEigen = 0; iEigen < nSelect; iEigen++){
            if (w[iEigen] < -1e-2) {
#if test_lazy2
                cout << " iEigen = " << iEigen;
                cout << ", w[iEigen] = " << w[iEigen] ;
                cout << ". z[iEigen] = " << endl << "[";
                for (int i = 0; i < n; i++) {
                    cout << z[i+iEigen*n] << ", " ;
                }
                cout  << "]" << endl;
#endif
                
                
                IloExpr lhs(env);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
#if 0 & test_lazy1
//                        cout << " lhs = " << lhs << endl;
                        IloExpr temp(env);
                        temp = z[j+iEigen*n] * z[i+iEigen*n] * x[i*n+j];
                        cout << " temp = " << temp << endl;
#endif
                        lhs += n  *  z[j+iEigen*n] * z[i+iEigen*n] * x[i*n+j];
                    }
                }
//                IloExpr lhs(env);
//                for (int i = 0; i < n; i++) {
//                    for (int j = 0; j < n; j++) {
//#if 0 & test_lazy1
//                        //                        cout << " lhs = " << lhs << endl;
//                        IloExpr temp(env);
//                        temp = z[j+iEigen*n] * z[i+iEigen*n] * x[i*n+j];
//                        cout << " temp = " << temp << endl;
//#endif
//                        if ( i == j )
//                            lhs += n * z[j+iEigen*n] * z[i+iEigen*n] * x[i*n+j];
//                        else if (i > j)
//                            lhs += 2 * n * z[j+iEigen*n] * z[i+iEigen*n] * x[i*n+j];
//                    }
//                }

                IloRange cut;
                cut = ( lhs >= 1e-4 );
#if test_lazy1
                if (n < 100)
                    cout << " cut = " << cut << endl;
#endif
                cuts2BeAdded.add( cut );
                tune.num_tangentcut ++ ;
            }
        }
        
#if test_lazy1
        cout << " ready to add cuts. " << "cuts2BeAdded.getSize() = " << cuts2BeAdded.getSize() << " (tag g3qw5gfq34fq432)" << endl;
#endif
        for (int i = 0; i < cuts2BeAdded.getSize(); i++){
            add( cuts2BeAdded[i] ).end();
        }
#if test_lazy1
        cout << " finished adding cuts (tag g3qw5gfq34fq432)" << endl;
#endif

        
        valx_ilo.end();
        cuts2BeAdded.end();

#if test_lazy2
        cout << " currObjValue = " << currObjValue << endl;
#endif
        
#if title_everywhere
        cout << " end of user cut callback " << endl;
        cout << " - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
        cout << "\n\n\n" << endl;
        {
            int temp;
            cin >> temp;
        }
#endif

        //        exit( 0 );
        
    } catch (...) {
        valx_ilo.end();
        //        for (int i = 0; i < cuts2BeAddedOnU.getSize(); i++)
        //            cuts2BeAddedOnU[i].end();
        //        for (int i = 0; i < cuts2BeAddedOnL.getSize(); i++)
        //            cuts2BeAddedOnL[i].end();
        //        for (int i = 0; i < cuts2BeAddedOnLandU.getSize(); i++)
        //            cuts2BeAddedOnLandU[i].end();
        //
        //        cuts2BeAddedOnL.end();
        //        cuts2BeAddedOnLandU.end();
        //        cuts2BeAddedOnU.end();
        throw;
    }
    //#if debugMode_display_cut_info_summary
    //    {
    //        cout << "  end of user cut round jwaoejfaw9jfaw4" << endl ;
    //    }
    //#endif
    //
    //#if debugMode_display_pause || debugMode_usercut
    //    {
    //        cout << " - at the end of user cut callback" << endl ;
    //        cout << " - type something here: ( 45hn667rhe )" << endl;
    //        int temp;
    //        cin >> temp;
    //    }
    //#endif
    
    
    
}

void MisdpSolver::solve () {
    int n = misdpmodel.n;
    int m = misdpmodel.m;
#if info_level1
    cout << " problem info: (tag 23cc45ybr67mty8) " << endl;
    cout << " m = " << m << ", n = " << n  << endl;
#endif
    try{
        IloModel        model(env);
        IloCplex        cplex(env);

        IloNumVarArray x(env, n*n, -IloInfinity, IloInfinity);
//        IloNumVarArray x(env, n*n);
                // x[i][j] = x[j*n+i]
        IloBoolVar tempBoolVar(env);
        
        namingVar(model, x);
        
        build(model, x, tempBoolVar);
        
//        cplex.use( myUser( env, x, misdpmodel, tune) );
        cplex.use( myLazy( env, x, misdpmodel, tune) );
//        cplex.use( myIncumbent( env, x, misdpmodel, tune) );
//        cplex.use( mySimplex( env, x, misdpmodel, tune) );
//        cplex.use( myMipinfo( env, x, misdpmodel, tune) );
//        cplex.use( myCallback( env, x, misdpmodel, tune) );

        /*
         ########     ###    ########     ###    ##     ##
         ##     ##   ## ##   ##     ##   ## ##   ###   ###
         ##     ##  ##   ##  ##     ##  ##   ##  #### ####
         ########  ##     ## ########  ##     ## ## ### ##
         ##        ######### ##   ##   ######### ##     ##
         ##        ##     ## ##    ##  ##     ## ##     ##
         ##        ##     ## ##     ## ##     ## ##     ##
         */
        //        cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
//        cplex.setParam(IloCplex::MIPDisplay, 3);
//        cplex.setParam(IloCplex::MIPInterval, 1);
//        cplex.setParam(IloCplex::NumericalEmphasis, 1);
//        cplex.setParam(IloCplex::ParallelMode, 1);
//        cplex.setParam(IloCplex::Threads, 4);
//        cplex.setParam(IloCplex::MIPEmphasis, 4);
        //        cplex.setParam(IloCplex::EpGap, .1);   // Any number from 0.0 to 1.0; default: 1e-04.
        
        //        cplex.setParam(IloCplex::EpAGap, 1e-12);  // Any nonnegative number; default: 1e-06.
        //        cplex.setParam(IloCplex::EpGap, 1e-12);   // Any number from 0.0 to 1.0; default: 1e-04.
        //        cplex.setParam(IloCplex::EpLin, 1e-12);   // Any positive value greater than zero; default: 1e-3.
        //        cplex.setParam(IloCplex::EpMrk, 0.99999); // Any number from 0.0001 to 0.99999; default: 0.01.
        //        cplex.setParam(IloCplex::EpInt, 0 ); // Any number from 0.0 to 0.5; default: 1e-05.
        //        cplex.setParam(IloCplex::EpPer, 1e-8); // Any positive number greater than or equal to 1e-8; default: 1e-6.
//        cplex.setParam(IloCplex::EpRHS, 1e-9); // Any number from 1e-9 to 1e-1; default: 1e-06.
        //        cplex.setParam(IloCplex::EpRelax, 1e-12); // Any nonnegative value; default: 1e-6.
        //        cplex.setParam(IloCplex::EpOpt, 1e-09); // Any number from 1e-9 to 1e-1; default: 1e-06.
        //
        //        cplex.setParam(IloCplex::BtTol, 0); // Any number from 1e-9 to 1e-1; default: 1e-06.
        
        
        //        cplex.setParam(IloCplex::SIM, 1e-09); // Any number from 1e-9 to 1e-1; default: 1e-06.
        
        //        cplex.setParam(IloCplex::EpMrk, 1e-6);
        //        cplex.setParam(IloCplex::Ep, 1e-6);
        
                cplex.setParam(IloCplex::PreInd, 0);
        // cplex.setParam(IloCplex::MIPEmphasis, 1);
        // this is for cplex 12.6
        //        cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
        //        cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
        //        cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
        
        cplex.extract(model);
        cplex.exportModel ("model1.lp");
#if title_everywhere
        cout << " ready to solve tag (fawefawefaw) " << endl;
#endif
        
        
        
#if 1
        cout << " * * * * * * * * * * * * * * * * * * * * * * * * * * * *  " << endl;
        cout << "\n\n\n" << endl;
        cout << " ready to solve cplex " << endl;
        cout << "\n\n\n" << endl;
        cout << " * * * * * * * * * * * * * * * * * * * * * * * * * * * *  " << endl;
        {
            int temp;
            cin >> temp;
        }
#endif

        cplex.solve();
#if title_everywhere
        cout << " solved tag (faw34f43) " << endl;
#endif
        cplex.writeSolution("./testWriteSol.sol");
        
        output(cplex, x, tune);

        
        
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();
}


void MisdpSolver::output(IloCplex &cplex, IloNumVarArray &x, MisdpTuneParameters &tune )
{
#if title_everywhere
    cout << " MisdpSolver :: output tag (fawef) " << endl;
#endif
//    int n = misdpmodel.n;
    IloNumArray valx(env);
    cplex.getValues(valx, x);
    env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << "Solution value  = " << cplex.getObjValue() << endl;
    env.out() << "# nodes = " << cplex.getNnodes64() << endl;
//    env.out() << "val.x  = " << valx << endl;
    env.out() << "tune.num_tangentcut = " << tune.num_tangentcut << endl;
    
    int n = misdpmodel.n;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cout << valx[i * n + j] << ", ";
        }
        cout << endl;
    }


}

MisdpSolver::MisdpSolver () {
    if_modelInputed = false;
}


void MisdpSolver::readdata(string filename) {
    misdpmodel.readdata(filename);
    
}



void MisdpSolver :: namingVar (IloModel &model, IloNumVarArray &x)
{
#if title_everywhere
    cout << " MisdpSolver :: namingVar tag (f23f23gq2) \n\n\n " << endl;
#endif
    int n = misdpmodel.n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            string name = "x_";
            name += to_string( i+1 );
            name += "_";
            name += to_string( j+1 );
            x[i*n+j].setName(name.c_str());
#if test_building
            cout << name.c_str() << endl;
#endif
        }
    }
//    
//    for (int i = 0; i < n*n; i++) {
//        string name = "x";
//        name += to_string( i+1 );
//        x[i].setName(name.c_str());
//#if test_building
//        cout << x[i] << endl;
//        cout << name.c_str() << endl;
//#endif
//    };

};




void MisdpSolver::build (IloModel &model, IloNumVarArray &x, IloBoolVar &tempBoolVar) {
    build_obj(model, x, tempBoolVar);
    build_initCons(model, x);
}



void MisdpSolver::build_obj (IloModel &model, IloNumVarArray &x, IloBoolVar &tempBoolVar)
{
#if title_everywhere
    cout << " MisdpSolver::build_obj tag (f24wfwaf) " << endl;
#endif
    IloExpr exprObj(env);
#if test_building
    cout << "obj = " << " (size)=" << misdpmodel.obj.coefficients.size()<<endl;
#endif

    for (int i = 0; i < misdpmodel.obj.coefficients.size(); i++){
#if test_building
        cout << " (" << misdpmodel.obj.coefficients[i] << ") , "   ;
        cout << x[i] << ", "  ;
//        cout << misdpmodel.obj.coefficients[i]* x[i] << ", " << endl  ;
//        cout << exprObj << ", " << endl  ;
#endif
        
        exprObj += misdpmodel.obj.coefficients[i] * x[i];
#if test_building
        if ( misdpmodel.obj.coefficients.size() < 2 )
            cout << exprObj << ", " << endl ;
#endif
    }
#if test_building
    cout  <<endl;
#endif

    exprObj += tempBoolVar;
#if test_building
    if (x.getSize() < 100)
        cout  << " exprObj = " << exprObj << endl;;
#endif
    
    model.add(IloMinimize(env, exprObj));
#if test_building
    cout  << " \n \n \n * * * * * * * * * \n \n \n exprObj = " << exprObj << endl;;
    {
        int temp;
        cin >> temp;
    }
#endif

}

void MisdpSolver::build_initCons (IloModel &model, IloNumVarArray &x) {
#if title_everywhere
    cout << " MisdpSolver::build_initCons tag (fq234f34gw35) " << endl;
#endif
    int n = misdpmodel.n;
    // symmetric
    for (int i = 0; i < n; i++ ) {
        for (int j = 0; j < n; j++ ) {
            if (i >= j)
                model.add( x[i*n+j] == x[j*n+i] );
        }
    }
    
#if test_build_initCons
    cout << " MisdpSolver::build_initCons tag (fwrges5hw54) " << endl;
#endif
    // non-negative diagonal
    for (int i = 0; i < n; i++) {
        x[i*n+i].setBounds(0, IloInfinity);
    }
#if test_build_initCons
    cout << " MisdpSolver::build_initCons tag (h54h45whw35) " << endl;
#endif

    // initial cons in misdpmodel
    for (int iCons =0; iCons < misdpmodel.m; iCons++) {
        IloExpr lhs(env);
#if test_build_initCons
        cout << " MisdpSolver::build_initCons tag (he56j6e545w) " << endl;
        cout << " misdpmodel.m = " << misdpmodel.m << endl;
#endif
        for (int iCoeff = 0; iCoeff < n*n; iCoeff++) {
#if 0 &  test_build_initCons
            cout << "misdpmodel.cons[iCons].coefficients.size() = " << misdpmodel.cons[iCons].coefficients.size() << endl;
//            cout << "misdpmodel.cons[iCons].coefficients[iCoeff] = " << misdpmodel.cons[iCons].coefficients[iCoeff] << endl;
//            cout << "misdpmodel.cons[iCons].coefficients[1] = " << misdpmodel.cons[iCons].coefficients[1] << endl;
//            cout << "misdpmodel.cons[iCons].coefficients[2] = " << misdpmodel.cons[iCons].coefficients[2] << endl;
            cout << "x[iCoeff] = " << x[iCoeff] << endl;
//            cout << "x[1] = " << x[1] << endl;
//            cout << "x[2] = " << x[2] << endl;
#endif
            lhs += misdpmodel.cons[iCons].coefficients[iCoeff] * x[iCoeff];
#if test_build_initCons
            cout << "x[iCoeff] = " << x[iCoeff] << endl;
#endif
        }
        double rhs = misdpmodel.cons[iCons].rhs;
#if test_build_initCons
        cout << " MisdpSolver::build_initCons tag (h56he4hw45) " << endl;
#endif

        switch (misdpmodel.cons[iCons].sign) {
//            case '>':
            case 'g':
                model.add(lhs >= rhs);
                break;
//            case '<':
            case 'l':
                model.add(lhs <= rhs);
                break;
//            case '=':
            case 'e':
                model.add(lhs == rhs);
                break;
        }
#if test_build_initCons
        cout << " MisdpSolver::build_initCons tag (gw35g56jr57) " << endl;
#endif

    }
#if test_build_initCons
    cout << " MisdpSolver::build_initCons tag (k86r7je5f34f) " << endl;
#endif

    // 2-dim submatrix
    for (int i = 0; i < n; i++ ) {
        for (int j = i+1; j < n; j++ ) {
            if (i != j){
                model.add( x[i*n+i] + x[j*n+j] >= 2 * x[j*n+i] );
                model.add( x[i*n+i] + x[j*n+j] >= - 2 * x[j*n+i] );
            }
        }
    }
    
    
    
    
}
