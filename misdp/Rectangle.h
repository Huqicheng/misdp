//
//  Rectangle.h
//  misdp
//
//  Created by Qi Zhang on 3/8/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#ifndef __misdp__Rectangle__
#define __misdp__Rectangle__

#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include "MisdpTuneParameters.h"
#include "MisdpModel.h"


class Rectangle {
    int width, height;
    MisdpModel misdpmodel;
    MisdpTuneParameters tune;
    bool if_modelInputed;
    
    IloEnv env;
    

    
public:
    Rectangle ();
    Rectangle (int,int);
    int area (void) {return (width*height);}

    void test_LAPACK();
};

#endif /* defined(__misdp__Rectangle__) */
