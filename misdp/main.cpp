//
//  main.cpp
//  misdp
//
//  Created by Qi Zhang on 3/8/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#include "main.h"






int main(int argc, const char * argv[])
{
    std::cout << "Ready to kick off!\n";
    
    
//    Rectangle rect (3,4);
//    Rectangle rectb;
//    cout << "rect area: " << rect.area() << endl;
//    cout << "rectb area: " << rectb.area() << endl;
//    return 0;
    
    
    //    MisdpSolver misdpInstance;
    MisdpSolver misdpInstance;
    /* initialize the tune */
    //    misdpInstance.tune.kSparse = 0;
    //    sparseRegressionInstance.tune.infilename = "/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/prostateCancer/prostate.txt";
    //    misdpInstance.tune.infilename = "/Users/Qi/Dropbox/research/my work/causal inference + structure learning + sparse regression/data/GaussianData/test/test2d3M100full.txt";
    //    misdpInstance.tune.outfilename = "./out.dat";
    //    misdpInstance.tune.constraintType = "L0";

    misdpInstance.readdata();
    misdpInstance.solve();

    std::cout << "All set!\n";
    return 0;
}






//std::cout << " default file:" << sparseRegressionInstance.tune.infilename << endl ;
//
//
///* deal with the input */
//for (int i = 1; i < argc-1; i++) { // command line options
//    if (!strncmp(argv[i],  "-cutUserPsdAdaptive", 19)) {
//        sparseRegressionInstance.tune.cutTypeFlag_User_adaptive_PSD = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutUserNsdAdaptive", 19)) {
//        sparseRegressionInstance.tune.cutTypeFlag_User_adaptive_NSD = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutUserPsd", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_User_PSD = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutUserNsd", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_User_NSD = atoi(argv[++i]); // the number of variables
//    };
//
//    if (!strncmp(argv[i],  "-cutUserUpEigMax1", 17)) {
//        sparseRegressionInstance.tune.cutTypeFlag_User_uppergradient_eig1 = atoi(argv[++i]); // the number of variables
//    };
//
//    if (!strncmp(argv[i],  "-cutLazyUpEigMax1", 17)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_uppergradient_eig1 = atoi(argv[++i]); // the number of variables
//    };
//
//
//    if (!strncmp(argv[i],  "-cutLazyNsd", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_NSD = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutLazyPsd", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_PSD = atoi(argv[++i]); // the number of variables
//    };
//
//    if (!strncmp(argv[i],  "-cutLazyUp1", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_upper1 = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutLazyUp2", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_upper2 = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutLazyUpG", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_uppergradient = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-cutLazyLow", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_Lazy_lowerfacet = atoi(argv[++i]); // the number of variables
//    };
//
//    if (!strncmp(argv[i],  "-initConTGP", 11)) {
//        sparseRegressionInstance.tune.initialConstraintFlag_totalGrad_PSD = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-initConUpp", 11)) {
//        sparseRegressionInstance.tune.initialConstraintFlag_upper = atoi(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-initConLow", 11)) {
//        sparseRegressionInstance.tune.initialConstraintFlag_lower = atoi(argv[++i]); // the number of variables
//    };
//
//
//
//
//    if (!strncmp(argv[i],  "-cutUserLow", 11)) {
//        sparseRegressionInstance.tune.cutTypeFlag_User_lowerfacet = atoi(argv[++i]); // the number of variables
//    };
//
//
//    if (!strncmp(argv[i],  "-limitNodeNumber", 16)) {
//        sparseRegressionInstance.tune.limitNumber_Nodes = atoi(argv[++i]); // the number of variables
//    };
//
//    if (!strncmp(argv[i],  "-test", 5)) {
//        sparseRegressionInstance.tune.testNumber1 = atof(argv[++i]); // the number of variables
//    };
//
//
//
//
//
//    if (!strncmp(argv[i],  "-i", 2)) {
//        sparseRegressionInstance.tune.infilename.assign(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-o", 2)) {
//        sparseRegressionInstance.tune.outfilename.assign(argv[++i]); // the number of variables
//    };
//    if (!strncmp(argv[i],  "-k", 2)) {
//        sparseRegressionInstance.tune.kSparse = atoi(argv[++i]);
//    };
//    if (!strncmp(argv[i],  "-lambda", 7)) {
//        sparseRegressionInstance.tune.lambdaPenalty = atof(argv[++i]);
//    };
//    if (!strncmp(argv[i],  "-sigma2", 7)) {
//        sparseRegressionInstance.tune.sigma2Noise = atof(argv[++i]);
//    };
//
//
//
//    if (!strncmp(argv[i],  "-c", 2)) {
//        sparseRegressionInstance.tune.constraintType.assign(argv[++i]); // the number of variables
//    };
//
//
//};
//
////    !!!!!
//std::cout << " ready to read file:" << sparseRegressionInstance.tune.infilename << endl ;
//
//sparseRegressionInstance.readData(sparseRegressionInstance.tune.infilename);
//
//std::cout << " ready to print info:" << endl;
//sparseRegressionInstance.tune.infoPrinter();
//
//std::cout << " ready to run > sparseRegressionInstance.core(); \n";
//int tempp;
//cin >> tempp;