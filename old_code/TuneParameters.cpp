//
//  TuneParameters.cpp
//  BNSL
//
//  Created by Qi on 4/26/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#include "TuneParameters.h"
using namespace std;


TuneParameters :: TuneParameters(void)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: setTuneParameters tag (gw345nw453) " << endl;
#endif
    sigma2Noise = 1.0;
    
    how_many_usercut_in_one_incumbent = 4;
    DOUBLE_EPSILON_USER_CUT_LOWER_THRESHOLD = 0 ; // 1e-4;
    DOUBLE_EPSILON_ADD_USER_CUT_THRESHOLD = 0 ; // 1e-4;
    DOUBLE_EPSILON_ON_CUT = 0 ; // 1e-4;
    DOUBLE_EPSILON = 0; // 1e-4;
    
    DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD = 1e-12; // 1e-8;
    
    flag_ifcount_num_of_all_cuts = true;
    
    num_lazy_cut = 0;
    num_user_cut = 0;
    
    num_lazy_cut_on_f = 0;
    num_lazy_cut_on_g = 0;
    num_lazy_cut_on_fandg = 0;
    num_user_cut_on_f = 0;
    num_user_cut_on_fandg = 0;
    num_user_cut_on_g = 0;
    num_sbm_lower_cut_on_g = 0;
    num_sbm_upper_cut_on_f = 0;
    num_gradient_upper_cut_on_f = 0;
    num_gradient_upper_cut_on_LMinusU_psd = 0;
    num_gradient_upper_cut_on_LMinusU_nsd = 0;
    num_gradient_upper_cut_on_LMinusU_partial = 0;
    num_gradient_cut_on_LMinusU_adaptive_psd = 0;
    num_gradient_cut_on_LMinusU_adaptive_nsd = 0;
    num_gradient_exp = 0;
    
    
    rootGap = 0;
    rootBestValue = -1e10;
    incumbentBestValue = 1e10;
    rootFlag = true;
    
    
    
    cutTypeFlag_Lazy_upper1 = true;
    cutTypeFlag_Lazy_upper2 = true;
    cutTypeFlag_Lazy_lowerfacet = true;
    cutTypeFlag_Lazy_uppergradient = true;
    cutTypeFlag_Lazy_NSD = false;
    cutTypeFlag_Lazy_PSD = false;
    cutTypeFlag_Lazy_uppergradient_eig1 = false;
    
    
    cutTypeFlag_User_lowerfacet = false;
    cutTypeFlag_User_uppergradient = false;
    cutTypeFlag_User_uppergradient_eig1 = true;
    cutTypeFlag_User_NSD = true;
    cutTypeFlag_User_PSD = true;
    cutTypeFlag_User_adaptive_PSD = true;
    cutTypeFlag_User_adaptive_NSD = true;
    
    initialConstraintFlag_totalGrad_PSD = true;
    initialConstraintFlag_upper = true;
    initialConstraintFlag_lower = true;

    limitNumber_Nodes = 1000;

    // temp
    testNumber1 = 1;

    
};


void TuneParameters :: infoPrinter()
{
    cout << "cutTypeFlag_Lazy_upper1 = " << cutTypeFlag_Lazy_upper1 << endl;
    cout << "cutTypeFlag_Lazy_upper2 = " << cutTypeFlag_Lazy_upper2 << endl;
    cout << "cutTypeFlag_Lazy_lowerfacet = " << cutTypeFlag_Lazy_lowerfacet << endl;
    cout << "cutTypeFlag_Lazy_uppergradient = " << cutTypeFlag_Lazy_uppergradient << endl;
    cout << "cutTypeFlag_Lazy_NSD = " << cutTypeFlag_Lazy_NSD << endl;
    cout << "cutTypeFlag_Lazy_PSD = " << cutTypeFlag_Lazy_PSD << endl;
    cout << "cutTypeFlag_User_lowerfacet = " << cutTypeFlag_User_lowerfacet << endl;
    cout << "cutTypeFlag_User_uppergradient = " << cutTypeFlag_User_uppergradient << endl;
    cout << "cutTypeFlag_User_NSD = " << cutTypeFlag_User_NSD << endl;
    cout << "cutTypeFlag_User_PSD = " << cutTypeFlag_User_PSD << endl;
    cout << "cutTypeFlag_User_adaptive_PSD = " << cutTypeFlag_User_adaptive_PSD << endl;
    cout << "cutTypeFlag_User_adaptive_NSD = " << cutTypeFlag_User_adaptive_NSD << endl;
    
    cout << "cutTypeFlag_User_uppergradient_eig1 = " << cutTypeFlag_User_uppergradient_eig1 << endl;
    
    
    cout << "cutTypeFlag_Lazy_upper1 = " << cutTypeFlag_Lazy_upper1 << endl;
    cout << "cutTypeFlag_Lazy_upper1 = " << cutTypeFlag_Lazy_upper1 << endl;
    cout << "cutTypeFlag_Lazy_upper1 = " << cutTypeFlag_Lazy_upper1 << endl;
    cout << "cutTypeFlag_Lazy_upper1 = " << cutTypeFlag_Lazy_upper1 << endl;
    
    
    cout << "initialConstraintFlag_totalGrad_PSD = " << initialConstraintFlag_totalGrad_PSD << endl;
    cout << "initialConstraintFlag_upper = " << initialConstraintFlag_upper << endl;
    cout << "initialConstraintFlag_lower = " << initialConstraintFlag_lower << endl;

    
    cout << "tune.limitNumber_Nodes = " << limitNumber_Nodes << endl; // atoi(argv[++i]); // the number of variables

};
