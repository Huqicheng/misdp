//
//  Qi_Algorithm.h
//  BNSL
//
//  Created by Qi on 4/26/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __BNSL__Qi_Algorithm__
#define __BNSL__Qi_Algorithm__

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>    // std::sort

#include <ilcplex/ilocplex.h>

#include "SystemParameters.h"


using namespace std;

class Qi_Algorithm{
public:
    
    template <typename T> vector<size_t> sort_indexes(const vector<T> &v) ;
    static vector<int> sort_IloNumArray_sub_indexes_decreasing(const IloNumArray &v, const int from, int to);
};





#endif /* defined(__BNSL__Qi_Algorithm__) */
