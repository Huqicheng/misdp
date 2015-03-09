//
//  Qi_Stat.h
//  BNSL
//
//  Created by Qi on 4/23/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __BNSL__Qi_Stat__
#define __BNSL__Qi_Stat__

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "SystemParameters.h"



class Qi_Stat{
public:
    
    double square(double);
    static double generateRandomUniform01 (double left=0.0, double right=1.0);
    double generateRandomNormal01 ();
    //    double sum(double *);
    
};



#endif /* defined(__BNSL__Qi_Stat__) */
