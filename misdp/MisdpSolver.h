//
//  MisdpSolver.h
//  misdp
//
//  Created by Qi Zhang on 3/7/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#ifndef __misdp__MisdpSolver__
#define __misdp__MisdpSolver__

#include <stdio.h>

//#include <iostream>
//#include <fstream>
//#include <string>
//#include <stack>          // std::stack


#include "MisdpTuneParameters.h"
#include "MisdpModel.h"
//#include "Qi_IO.h"
//#include "Qi_Algorithm.h"
//#include "Qi_Stat.h"


#include <ilcplex/ilocplex.h>

//using namespace std;




class MisdpSolver {
    MisdpModel misdpmodel;
    MisdpTuneParameters tune;
    bool if_modelInputed;
    
    IloEnv env;
public:
    MisdpSolver();
    
    void solve();
    void readdata(string filename="/Users/Qi/Dropbox/research/my work/MISDP-solver/misdp/src/data/data.txt");

};




#endif /* defined(__misdp__MisdpSolver__) */
