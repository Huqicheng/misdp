//
//  MisdpModel.cpp
//  misdp
//
//  Created by Qi Zhang on 3/7/15.
//  Copyright (c) 2015 Qi Zhang. All rights reserved.
//

#define test_input 1
#include "MisdpModel.h"

void MisdpModel :: readdata( string filename )
{
    /*
     // data example:
     obj 5
     -0.7590 0.4124 0.7353 0.0959 -0.5028 0.4124 0.1498 0.0129 -0.3608 0.2109 0.7353 0.0129 -0.9140 -0.4120 0.1673 0.0959 -0.3608 -0.4120 0.7516 0.1118 -0.5028 0.2109 0.1673 0.1118 -0.2621
     cons 1
     1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1
     < rhs
     */
    
    ifstream f(filename);
    string line;
    
    getline(f, line); // expect "obj, n"
    {
        istringstream iss(line);
        string first;
        iss >> first;
        if (! (first == "obj")) {
            cout << "iss = "  << endl;
            cout << "iss = " << iss.str() << endl;
            cout << "Error data file. First = " << first << "\n";
            cout << "filename = " << filename << endl;
            cout << "some action needed ! (tag q34f3w5g46) ";
            int temp;
            cin >> temp;
        }
        iss >> n;
#if test_input
        cout << "n = " << n << endl;
#endif
    }
    
    getline(f, line); // expect "coef of obj", n*n numbers
    {
        istringstream iss(line);
        vector< double > objCoef;
        double currCoef;
#if test_input
        cout << "coef = " << endl;
#endif
        for (int i = 0; i < n*n; i++) {
            iss >> currCoef;
            objCoef.push_back(currCoef);
#if test_input
            cout << currCoef << "\t"  ;
#endif
        }
        obj.setCoefficients(objCoef);
    }
#if test_input
    cout  << endl;
#endif

    getline(f, line); // expect "cons, m"
    {
        istringstream iss(line);
        string first;
        iss >> first;
        if (! (first == "cons")) {
            cout << "Error data file. First = " << first << "\n";
            cout << "filename = " << filename;
            cout << "some action needed ! (tag f34f34g34q) ";
            int temp;
            cin >> temp;
        }
        iss >> m;
#if test_input
        cout << "m = " << m << endl;
#endif
    }
    for (int iCons = 0; iCons < m; iCons++) {
        Constraint con;
        
        getline(f, line); // one line of constraints coef, expect 1 0 0 0 0 0 1 0 0 ...
        vector< double > conCoef;
        double currCoef;
        istringstream iss(line);
#if test_input
        cout << "coef = " << endl;
#endif
        for (int i = 0; i < n*n; i++) {
            iss >> currCoef;
            conCoef.push_back(currCoef);
#if test_input
            cout << currCoef << ", "  ;
#endif
        }
        con.setCoefficients(conCoef);
#if test_input
        cout  << endl;
#endif

        getline(f, line); // one line of constraints sign, rhs, expect "< rhs"
        istringstream iss2(line);
        char sign;
        double rhs;
        iss2 >> sign >> rhs;
        con.setRhs(rhs);
        con.setSign(sign);
        
        cons.push_back(con);
    }
};
































void Expression::setCoefficients(const vector<double> & coef)
{
    coefficients = coef;
#if test_input
    cout << "\nin Expression::setCoefficients. (tag fawoiefjoaw) coefficients = " << endl;
#endif
    for (int i = 0; i < coefficients.size(); i++) {
#if test_input
        cout << coefficients[i] << ", "  ;
#endif
    }

};

void Constraint::setSign(char s)
{
    switch (s) {
        case '>':
        case 'g':
            sign = 'g';
            break;
        case '<':
        case 'l':
            sign = 'l';
            break;
        case '=':
        case 'e':
            sign = 'e';
            break;

        default:
            {
                cout << "Error sign. s = " << s << "\n";
                cout << "some action needed ! (tag fq23fq23fawef) ";
                int temp;
                cin >> temp;
            }
            break;
    }
};

void Constraint::setRhs(double r)
{
    rhs = r;
};



MisdpModel::MisdpModel(){
    ;
}
Objective::Objective(){
    ;
}