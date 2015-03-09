//
//  TuneParameters.h
//  BNSL
//
//  Created by Qi on 4/26/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __SR__TuneParameters__
#define __SR__TuneParameters__

#include <iostream>
#include <string>

class TuneParameters{
public:
    TuneParameters();
    void infoPrinter();
    
    /* model parameters */
    int kSparse; // k in k-sparsity, if 0, then BIC
    double lambdaPenalty; //  || Y-X beta ||^2 + \lambda * cardinality
    double sigma2Noise; //  || Y-X beta ||^2 + log(m)/2 * sigma2 * cardinality
    std::string infilename;
    std::string outfilename;
    std::string constraintType;

    
    /* tune parameters*/
    double DOUBLE_EPSILON;
    double DOUBLE_EPSILON_USER_CUT_LOWER_THRESHOLD;
    double DOUBLE_EPSILON_ADD_USER_CUT_THRESHOLD;
    double DOUBLE_EPSILON_ON_CUT;
    double DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD;

    
    int how_many_usercut_in_one_incumbent;
    
    bool flag_ifcount_num_of_all_cuts;
    int num_lazy_cut_on_f;
    int num_lazy_cut_on_g;
    int num_user_cut_on_f;
    int num_lazy_cut_on_fandg;
    int num_user_cut_on_g;
    int num_user_cut_on_fandg;
    
    int num_user_cut;
    int num_lazy_cut;

    int num_sbm_lower_cut_on_g;
    int num_sbm_upper_cut_on_f;
    int num_gradient_upper_cut_on_f;
    int num_gradient_upper_cut_on_LMinusU_psd;
    int num_gradient_upper_cut_on_LMinusU_nsd;
    int num_gradient_upper_cut_on_LMinusU_partial;
    int num_gradient_cut_on_LMinusU_adaptive_psd;
    int num_gradient_cut_on_LMinusU_adaptive_nsd;
    
    
    int num_gradient_exp;
//    int num_cut_on_p;
//    int num_lazycut_on_p;
//    int num_cut_on_cycle;
//    int num_heuristic_provided;
    
    
    
    
    
    
    //temp in a round,
    bool any_cut_added ;
    int num_cut_added_this_round ;
    double cumulative_violated_amount;
    
    

    double rootGap;
    double rootBestValue;
    double incumbentBestValue;
    bool rootFlag;
    
    
    // type of cut
    bool cutTypeFlag_Lazy_upper1;
    bool cutTypeFlag_Lazy_upper2;
    bool cutTypeFlag_Lazy_lowerfacet;
    bool cutTypeFlag_Lazy_uppergradient;
    bool cutTypeFlag_Lazy_NSD;
    bool cutTypeFlag_Lazy_PSD;
    bool cutTypeFlag_Lazy_uppergradient_eig1;
    
    bool cutTypeFlag_User_lowerfacet;
    bool cutTypeFlag_User_uppergradient;
    bool cutTypeFlag_User_uppergradient_eig1;
    bool cutTypeFlag_User_NSD;
    bool cutTypeFlag_User_PSD;
    bool cutTypeFlag_User_adaptive_PSD;
    bool cutTypeFlag_User_adaptive_NSD;
    
    bool initialConstraintFlag_totalGrad_PSD;
    bool initialConstraintFlag_upper;
    bool initialConstraintFlag_lower;
    
    
    long limitNumber_Nodes;

    
    // temp
    double testNumber1;
};



#endif /* defined(__SR__TuneParameters__) */
