//
//  Qi_IO.h
//  BNSL
//
//  Created by Qi on 4/23/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#ifndef __BNSL__Qi_IO__
#define __BNSL__Qi_IO__







#include <iostream>
//#include <stdio.h>      /* printf, scanf, puts, NULL */
//#include <stdlib.h>     /* srand, rand */
//#include <time.h>       /* time */
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>


#include "SystemParameters.h"
#include "DataModel.h"



using namespace std;


class Qi_IO_input{
public:
    static void readData(const string & filename, DataModel & data);

    static void writeData(const string & filename, const string & content, bool add =true);
    static void writeData(const string & filename, const stringstream & content, bool add =true);

    /*
     template <class T>
     T GetMax (T a, T b) {
     T result;
     result = (a>b)? a : b;
     return (result);
     // usage:
     //        int i=5, j=6, k;
     //        long l=10, m=5, n;
     //        k=GetMax<int>(i,j);
     //        n=GetMax<long>(l,m);
     }
     */
};


class ReadData{
public:
    // static bool readTableToLasso(string, Model *, int /* number of var*/);
    // static bool readTableToLassoV2(string, Model * /* number of var*/);
};
class WriteData{
public:
    static bool writeStringToFile(const string &, const string &, bool);
    static bool writeStringToFile(const string &, const string &);
};






#endif /* defined(__BNSL__Qi_IO__) */
