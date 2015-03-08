//
//  MisdpSolver.cpp
//  misdp
//
//  Created by Qi Zhang on 3/7/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#include "MisdpSolver.h"
/*
     min  <C_0, X>
     s.t. <A_i, X> \geq b_i, i = 1...m
          X >= 0
 */

MisdpSolver::MisdpSolver () {
    if_modelInputed = false;
}



void MisdpSolver::readdata(string filename) {
    misdpmodel.readdata(filename);
}

void MisdpSolver::solve () {
    ;
}




