//
//  SparseRegressionProblem.h
//  SparseRegression
//
//  Created by Qi on 5/27/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __SparseRegression__SparseRegressionProblem__
#define __SparseRegression__SparseRegressionProblem__

#include <iostream>
#include <fstream>
#include <string>
#include <stack>          // std::stack


#include "SystemParameters.h"
#include "TuneParameters.h"
#include "DataModel.h"
#include "Qi_IO.h"
#include "Qi_Algorithm.h"
#include "Qi_Stat.h"
// #include "TuneParameters.h"
#include <stdio.h>




#include <ilcplex/ilocplex.h>

using namespace std;


class SparseRegressionProblem{
public:
    void readData(const string & filename);
    void core(bool flagCPX=true);
    
    void buildModel();
    void build_SR_VarNames (IloModel &model, IloBoolVarArray &z, IloNumVar &f, IloNumVar &g , IloNumVar &mse ) ;
    void build_SR_Objective (IloModel &model, IloBoolVarArray &z, IloNumVar &f, IloNumVar &g , IloNumVar &mse ) ;
    void build_SR_ULboth_Constraints (IloModel &model, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l ) ;
    void build_SR_Exp_Constraints (IloModel &model, IloNumVar &mse, IloNumVar &u, IloNumVar &l) ;
    void build_SR_UpperConstraints (IloModel &model, IloBoolVarArray &z, IloNumVar &u) ;
    void build_SR_LowerConstraints (IloModel &model, IloBoolVarArray &z, IloNumVar &g) ;
    void build_SR_Cardinality_Constraints (IloModel &model, IloBoolVarArray &z) ;
//    void build_SR_PenaltyConstraints (IloModel &model, IloBoolVarArray &x, IloNumVarArray &p) ;
//    void build_SR_AcyclicConstraints (IloModel &model, IloBoolVarArray &x) ;
    
//    void feed_SR_Solutions ( IloCplex &cplex, IloBoolVarArray &x, IloNumVarArray &f, IloNumVarArray &g, IloNumVarArray &p ) ;
    void output_SR_Solutions ( IloCplex &cplex, IloBoolVarArray &z, IloNumVar &f, IloNumVar &g, IloNumVar &mse, TuneParameters &tune  ) ;
    
    
    
    void static generateCuts_lowercuts( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &g, IloNumArray &valz, IloNum valg);


    void static generateCuts_uppercuts_submodular1( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf, bitset<_CONST_INT_MAX_FACTOR_SIZE> &currentSelection);
    void static generateCuts_uppercuts_submodular2( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf, bitset<_CONST_INT_MAX_FACTOR_SIZE> &currentSelection);

    void static generateCuts_uppercuts_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf);
    void static generateCuts_uppercuts_maxEig1_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf);
    
    void static generateCuts_exp_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, TuneParameters &tune, IloNumVar &mse, IloNum valmse, IloNumVar &f, IloNum valf, IloNumVar &g, IloNum valg);
    
    void static generateCuts_total_gradient_psd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_adaptive_psd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_adaptive_psd_copy20141031_0138pm( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);

    void static generateCuts_total_gradient_adaptive_psd_copy2( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_nsd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_adaptive_nsd_old( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_adaptive_nsd_old_negneg( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_adaptive_nsd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    void static generateCuts_total_gradient_adaptive_nsd_copy2( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);
    
    void static generateCuts_partial_r_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz);

    
    void setTuneParameters();
    void test();
    //    void writeSolution();
    
    
    int n, m;

    // last
    vector<double> lastValz;

    
//private:
    TuneParameters tune;
    DataModel data;
    IloEnv env;
    
    
};





#endif /* defined(__SparseRegression__SparseRegressionProblem__) */
