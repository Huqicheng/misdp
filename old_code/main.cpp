//
//  main.cpp
//  SparseRegression
//
//  Created by Qi on 5/27/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

// #include <iostream>
#include "main.h"





int main(int argc, const char * argv[])
{
    std::cout << "Ready to kick off!\n";
    std::cout << "_CONST_INT_MAX_FACTOR_SIZE = " << _CONST_INT_MAX_FACTOR_SIZE << endl;
    
    SparseRegressionProblem sparseRegressionInstance;

    /* initialize the tune */
    sparseRegressionInstance.tune.kSparse = 0;
//    sparseRegressionInstance.tune.infilename = "/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/prostateCancer/prostate.txt";
    sparseRegressionInstance.tune.infilename = "/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/GaussianData/test/test2d3M100full.txt";
    sparseRegressionInstance.tune.outfilename = "./out.dat";
    sparseRegressionInstance.tune.constraintType = "L0";

    std::cout << " default file:" << sparseRegressionInstance.tune.infilename << endl ;

    
    /* deal with the input */
    for (int i = 1; i < argc-1; i++) { // command line options
        if (!strncmp(argv[i],  "-cutUserPsdAdaptive", 19)) {
            sparseRegressionInstance.tune.cutTypeFlag_User_adaptive_PSD = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutUserNsdAdaptive", 19)) {
            sparseRegressionInstance.tune.cutTypeFlag_User_adaptive_NSD = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutUserPsd", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_User_PSD = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutUserNsd", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_User_NSD = atoi(argv[++i]); // the number of variables
        };

        if (!strncmp(argv[i],  "-cutUserUpEigMax1", 17)) {
            sparseRegressionInstance.tune.cutTypeFlag_User_uppergradient_eig1 = atoi(argv[++i]); // the number of variables
        };
        
        if (!strncmp(argv[i],  "-cutLazyUpEigMax1", 17)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_uppergradient_eig1 = atoi(argv[++i]); // the number of variables
        };
        
        
        if (!strncmp(argv[i],  "-cutLazyNsd", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_NSD = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutLazyPsd", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_PSD = atoi(argv[++i]); // the number of variables
        };
        
        if (!strncmp(argv[i],  "-cutLazyUp1", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_upper1 = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutLazyUp2", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_upper2 = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutLazyUpG", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_uppergradient = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-cutLazyLow", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_Lazy_lowerfacet = atoi(argv[++i]); // the number of variables
        };

        if (!strncmp(argv[i],  "-initConTGP", 11)) {
            sparseRegressionInstance.tune.initialConstraintFlag_totalGrad_PSD = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-initConUpp", 11)) {
            sparseRegressionInstance.tune.initialConstraintFlag_upper = atoi(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-initConLow", 11)) {
            sparseRegressionInstance.tune.initialConstraintFlag_lower = atoi(argv[++i]); // the number of variables
        };

        
        
        
        if (!strncmp(argv[i],  "-cutUserLow", 11)) {
            sparseRegressionInstance.tune.cutTypeFlag_User_lowerfacet = atoi(argv[++i]); // the number of variables
        };

        
        if (!strncmp(argv[i],  "-limitNodeNumber", 16)) {
            sparseRegressionInstance.tune.limitNumber_Nodes = atoi(argv[++i]); // the number of variables
        };

        if (!strncmp(argv[i],  "-test", 5)) {
            sparseRegressionInstance.tune.testNumber1 = atof(argv[++i]); // the number of variables
        };
        
        
        
        

        if (!strncmp(argv[i],  "-i", 2)) {
            sparseRegressionInstance.tune.infilename.assign(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-o", 2)) {
            sparseRegressionInstance.tune.outfilename.assign(argv[++i]); // the number of variables
        };
        if (!strncmp(argv[i],  "-k", 2)) {
            sparseRegressionInstance.tune.kSparse = atoi(argv[++i]);
        };
        if (!strncmp(argv[i],  "-lambda", 7)) {
            sparseRegressionInstance.tune.lambdaPenalty = atof(argv[++i]);
        };
        if (!strncmp(argv[i],  "-sigma2", 7)) {
            sparseRegressionInstance.tune.sigma2Noise = atof(argv[++i]);
        };
        


        if (!strncmp(argv[i],  "-c", 2)) {
            sparseRegressionInstance.tune.constraintType.assign(argv[++i]); // the number of variables
        };

    
};
    

//    sparseRegressionInstance.readData("/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/prostateCancer/prostate.txt");
//    sparseRegressionInstance.readData("/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/temp_2d/data2d.txt");
//    sparseRegressionInstance.readData("/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/YearPredictionMSD/YearPredictinoMSD_r1000_c21.txt");
    
    
//    sparseRegressionInstance.readData("/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/communitiesAndCrime/SigmaCommunitiesCrimeData.txt");
    
//    sparseRegressionInstance.readData("/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/RelativeLocationOfCT/dataSliceLocation.txt");

    
//    !!!!!
    std::cout << " ready to read file:" << sparseRegressionInstance.tune.infilename << endl ;

    sparseRegressionInstance.readData(sparseRegressionInstance.tune.infilename);

    std::cout << " ready to print info:" << endl;
    sparseRegressionInstance.tune.infoPrinter();
    
    std::cout << " ready to run > sparseRegressionInstance.core(); \n";
    int tempp;
    cin >> tempp;

    sparseRegressionInstance.core();

    std::cout << "All set!\n";
    return 0;
}






















































































































































































