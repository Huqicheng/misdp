//
//  MisdpSolver.cpp
//  misdp
//
//  Created by Qi Zhang on 3/7/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//
#define title_everywhere 1
#define info_level1 1

#include "MisdpSolver.h"
/*
     min  <C_0, X>
     s.t. <A_i, X> \geq b_i, i = 1...m
          X >= 0
 */
ILOUSERCUTCALLBACK3(myUser, IloNumVarArray,x, MisdpModel &,misdpmodel, MisdpTuneParameters &, tune)
{
    //    if (getNnodes64() > 10)
    //        return;
#if title_everywhere
    cout << " ILOUSERCUTCALLBACK3 tag (fq34345wgw345) " << endl;
#endif
    
    if (getNnodes64() <= 1){
        tune.rootBestValue = getBestObjValue();
    };
}


//
//    
//    
//    
//    tune.any_cut_added = false;
//    tune.num_cut_added_this_round = 0;
//    tune.cumulative_violated_amount = 0;
//    
//#if debugMode_positionTitle
//    cout << " \n * * \n * * user cut callback: upper gradient + lower facet (tag f234v57henu68r)  " << endl;
//#endif
//    
//#if debugMode_display_callback_info
//    cout << "  getNnodes64 = " << getNnodes64() << ", getNremainingNodes64 = " << getNremainingNodes64() << endl;
//#endif
//    IloEnv env = getEnv();
//    IloNumArray valz(env);
//    IloNum valf = 0;
//    IloNum valg = 0;
//    IloNum valmse = 0;
//    
//    
//    IloConstraintArray cuts2BeAddedOnU(env);
//    IloConstraintArray cuts2BeAddedOnL(env);
//    IloConstraintArray cuts2BeAddedOnLandU(env);
//    
//    
//    try {
//        // now we may have a integer feasible solution
//        valz  = IloNumArray(env);
//        valf  = getValue( f );
//        valg  = getValue( g );
//        getValues( valz, z );
//        valmse  = getValue( mse );
//        
//        int n = data.getCols();
//        
//        
//#if debugMode_display_callback_sol
//        cout << " valz = " << valz;
//        cout << ", valf = " << valf;
//        cout << ", valg = " << valg << endl;
//#endif
//        double distFromLastValz = 0;
//        for (int i = 0; i < n ; i++){
//            distFromLastValz += pow(lastValz[i] - valz[i], 2);
//        }
//        distFromLastValz = sqrt(distFromLastValz);
//#if debugMode_display_callback_sol
//        cout << " lastValz = [ " ;
//        for (int i = 0; i < n; i++)
//            cout << lastValz[i] << ", " ;
//        cout << "]  distFromLastValz = " << distFromLastValz << endl;
//#endif
//        if (distFromLastValz < .01*n)
//            return;
//        
//        
//        
//        
//        bitset<_CONST_INT_MAX_FACTOR_SIZE> currentSelection;
//        //        bitset<_CONST_INT_MAX_FACTOR_SIZE> currentNode;
//        for (int j = 0; j < n; j++)
//            if ( valz[j] > 0.5 )
//                currentSelection.set(j);
//#if debugMode_display_sol // debugMode_display_sol
//        cout << "    sol information (tag fq4293jfq3400): " << endl;
//        cout << ",   currentSelection = " << currentSelection << endl;
//        cout << "    valf = " << valf << ", data.logDetSubmatrix(currentSelection) = " << data.logDetSubmatrix(currentSelection) << endl;
//        cout << "    valg = " << valg << ", data.logDetSubmatrix_with_y(currentSelection) = " << data.logDetSubmatrix_with_y(currentSelection) << endl;
//#endif
//        
//        
//        
//#if debugMode_usercut && 0
//        {
//            cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl ;
//            cout << " - - - - - a new label to test when to add user gradient cut" << endl<< endl<< endl ;
//            
//            cout << " f.getUB() = " << f.getUB() ;
//            cout << ",  f.getLB() = " << f.getLB() << endl;
//            
//            double c = data.logDet_IplusDiagzW( valz ) ;
//            cout << " c = " << c << ", f = " << valf << endl << endl;
//            
//            int temp;
//            cin >> temp;
//        }
//#endif
//        
//        // test
//        cout << " valz = " <<  valz << endl;
//        //
//        if  ( tune.cutTypeFlag_User_adaptive_NSD ) {
//            cout << " apply generateCuts_total_gradient_adaptive_nsd (user, tag fq4brste4wafers) \n";
//            SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
//            cout << " end generateCuts_total_gradient_adaptive_nsd (user, tag fa3w4egvrv545dt) \n";
//        }
//        
//        if  ( tune.cutTypeFlag_User_NSD ) {
//            cout << " apply generateCuts_total_gradient_nsd (user, tag fawfeawfaewfaw ) \n";
//            SparseRegressionProblem :: generateCuts_total_gradient_nsd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
//        }
//        
//        
//        
//        if  ( tune.cutTypeFlag_User_adaptive_PSD ) {
//            cout << " apply generateCuts_total_gradient_adaptive_psd (user, tag fawrbs5enhdr6n) \n";
//            SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
//            cout << " end generateCuts_total_gradient_adaptive_psd (user, tag fzsebfdr6m7) \n";
//        }
//        
//        
//        if  ( tune.cutTypeFlag_User_PSD ) {
//            cout << " apply generateCuts_total_gradient_psd (user, tag fawefawfawgae) \n";
//            SparseRegressionProblem :: generateCuts_total_gradient_psd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
//        }
//        
//        
//        
//        //        double c = data.logDet_IplusDiagzW( valz ) ;
//        
//        
//        if (tune.cutTypeFlag_User_uppergradient_eig1) {
//            cout << " apply generateCuts_uppercuts_maxEig1_gradient (user, tag 2039f203j) \n";
//            SparseRegressionProblem :: generateCuts_uppercuts_maxEig1_gradient(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf);
//        }
//        
//        
//        
//        // upper gradient
//        //        if ( ( tune.cutTypeFlag_User_uppergradient )
//        //                && ( valf > c + tune.DOUBLE_EPSILON_ADD_USER_CUT_THRESHOLD ) ) {
//        if (  tune.cutTypeFlag_User_uppergradient  ) {
//            cout << " apply generateCuts_uppercuts_gradient (user, tag 2309r23jjf2) \n";
//            SparseRegressionProblem :: generateCuts_uppercuts_gradient(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf);
//        }
//        
//        
//        
//        
//        // lower facet
//        if ( (tune.cutTypeFlag_User_lowerfacet) ){
//            cout << " apply generateCuts_lowercuts (user) \n";
//            SparseRegressionProblem :: generateCuts_lowercuts(env, cuts2BeAddedOnL, data, tune, z, g, valz, valg);
//        }
//        
//        // exp grad
//        if (1 && (tune.kSparse >= 0)){
//            cout << " apply generateCuts_exp_gradient (user) \n";
//            SparseRegressionProblem :: generateCuts_exp_gradient(env, cuts2BeAddedOnLandU, tune, mse, valmse, f, valf, g, valg);
//            cout << " end generateCuts_exp_gradient (user) \n";
//        }
//        
//        if  ( 0 ) {
//            cout << " apply generateCuts_partial_r_gradient (user) \n";
//            SparseRegressionProblem :: generateCuts_partial_r_gradient(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
//        }
//        
//        
//        
//        
//#if debugMode_display_cut_info_summary
//        cout << " number of cuts added in this round of USER on f = " << cuts2BeAddedOnU.getSize() << endl;
//        cout << " number of cuts added in this round of USER on g = " << cuts2BeAddedOnL.getSize() << endl;
//        cout << " number of cuts added in this round of USER on f-g = " << cuts2BeAddedOnLandU.getSize() << endl;
//#endif
//        
//        for (int i = 0; i < cuts2BeAddedOnU.getSize(); i++){
//#if debugMode_display_cut_info_summary
//            cout << " in cuts2BeAddedOnU loop, i = " << i << endl;
//#endif
//            if (n < 5)
//                cout << cuts2BeAddedOnU[i] << endl;
//            
//            add( cuts2BeAddedOnU[i] ).end();
//            tune.num_user_cut++;
//            tune.num_user_cut_on_f++;
//        }
//        
//        for (int i = 0; i < cuts2BeAddedOnL.getSize(); i++){
//#if debugMode_display_cut_info_summary
//            cout << " in cuts2BeAddedOnL loop, i = " << i << endl;
//#endif
//            if (n < 5)
//                cout << cuts2BeAddedOnL[i] << endl;
//            
//            add( cuts2BeAddedOnL[i] ).end();
//            tune.num_user_cut++;
//            tune.num_user_cut_on_g++;
//        }
//        for (int i = 0; i < cuts2BeAddedOnLandU.getSize(); i++){
//#if debugMode_display_cut_info_summary
//            cout << " in cuts2BeAddedOnLandU loop, i = " << i << endl;
//#endif
//            if (n < 5)
//                cout << cuts2BeAddedOnLandU[i] << endl;
//            
//            add( cuts2BeAddedOnLandU[i] ).end();
//            tune.num_user_cut++;
//            tune.num_user_cut_on_fandg++;
//        }
//        
//#if debugMode_display_pause || debugMode_usercut
//        cout << " valz = " << valz;
//        cout << ", valf = " << valf;
//        cout << ", valg = " << valg << endl;
//#endif
//        
//        
//        
//        for (int i = 0; i < n ; i++)
//            lastValz[i] = valz[i];
//        
//        valz.end();
//        cuts2BeAddedOnL.end();
//        cuts2BeAddedOnLandU.end();
//        cuts2BeAddedOnU.end();
//        
//    } catch (...) {
//        valz.end();
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
//        throw;
//        
//    }
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
//    
//    
//};





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
        IloNumVarArray x(env=env, n=n*n);                    // x[i][j] = x[j*n+i]
        
        cplex.use( myUser( env, x, misdpmodel, tune) );

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
        
        //        cplex.setParam(IloCplex::PreInd, 0);
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
    env.out() << "val.x  = " << valx << endl;


}

MisdpSolver::MisdpSolver () {
    if_modelInputed = false;
}


void MisdpSolver::readdata(string filename) {
    misdpmodel.readdata(filename);
    
}


void MisdpSolver::build (IloModel &model, IloBoolVarArray &x) {
    build_obj(model, x);
    build_initCons(model, x);
}



void MisdpSolver::build_obj (IloModel &model, IloBoolVarArray &x) {
#if title_everywhere
    cout << " MisdpSolver::build_obj tag (fq34gw3q) " << endl;
#endif
    IloNumExpr exprObj(env);
    for (int i = 0; i < misdpmodel.obj.coefficients.size(); i++)
        exprObj += misdpmodel.obj.coefficients[i] * x[i];
    
    model.add(IloMinimize(env, exprObj));
}

void MisdpSolver::build_initCons (IloModel &model, IloBoolVarArray &x) {
#if title_everywhere
    cout << " MisdpSolver::build_initCons tag (fq234f34gw35) " << endl;
#endif
    int n = misdpmodel.n;
    // symmetric
    for (int i = 0; i < n; i++ ) {
        for (int j = 0; j < n; j++ ) {
            if (i != j)
                model.add( x[i*n+j] == x[j*n+i] );
        }
    }
    
    // 2-dim submatrix
    for (int i = 0; i < n; i++ ) {
        for (int j = 0; j < n; j++ ) {
            if (i != j){
                model.add( x[i*n+i] >= x[j*n+i] );
                model.add( x[i*n+i] >= -x[j*n+i] );
            }
        }
    }
    
}
