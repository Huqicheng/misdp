//
//  Qi_Math.h
//  BNSL
//
//  Created by Qi on 4/23/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __BNSL__Qi_Math__
#define __BNSL__Qi_Math__

#include <iostream>



#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include <ilcplex/ilocplex.h>

#include "SystemParameters.h"


using namespace std;

class Qi_Math{
public:

    double square(double);
    long double square(long double);
//    double sum(double *);
//    template <typename T> vector<size_t> sort_indexes(const vector<T> &v) ;                       // move to Qi_Algorithm
//    vector<int> sort_IloNumArray_sub_indexes(const IloNumArray &v, const int from, int to);       // move to Qi_Algorithm
};






#endif /* defined(__BNSL__Qi_Math__) */
