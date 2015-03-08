//
//  MisdpModel.h
//  misdp
//
//  Created by Qi Zhang on 3/7/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#ifndef __misdp__MisdpModel__
#define __misdp__MisdpModel__

#include <stdio.h>

#include <iostream>
#include <bitset>
#include <vector>
#include <cmath>
#include <algorithm>    // std::min



//#include <unordered_map>
//#include "TuneParameters.h"
#include <ilcplex/ilocplex.h>




using namespace std;

class Expression{
    vector<double> coefficients;
public:
//    Expression(vector<double>);
//    Expression( );
    void setCoefficients(const vector<double> &);
};


class Constraint : public Expression {
    double rhs;
    char sign; // 'g', 'l', 'e' ~ greater, less, equal to
public:
//    Constraint();
//    Constraint(vector<double>, char sign = 'g', double rhs = 0);
    void setRhs(double);
    void setSign(char);
};


class Objective : public Expression {
    string sense = "min";
public:
    Objective();
    Objective(vector<double>, string sense = "min");
};


class MisdpModel{
    Objective obj;
    vector<Constraint> cons;
    int m;
    int n;

public:
    MisdpModel();
    void readdata(string filename="/Users/Qi/Dropbox/research/my work/MISDP-solver/misdp/src/data/data.txt");

};




#endif /* defined(__misdp__MisdpModel__) */
