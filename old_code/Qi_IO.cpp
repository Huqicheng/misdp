//
//  Qi_IO.cpp
//  BNSL
//
//  Created by Qi on 4/23/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//
#define debugMode_display_pause 0


#include "Qi_IO.h"




void Qi_IO_input :: readData(const string &filename, DataModel & data)
{
#if 1
    cout << " Qi_IO_input :: readData (tag 3298ru239)" << endl;
    cout << " we are going to read data. the data should be \ny1 x11 x12 x13 ... \ny2 x21 x22 x23 ... \n... " << endl;
#endif
    int m = 0;
    int n = 0; // m = #rows, n = #cols
    ifstream f(filename);

    string line;

    while(getline(f, line))
    {
        // getline. If it works then processes then try again.
        istringstream iss(line);
        m++;
        double currentNum;
        iss >> currentNum;
        data.push_back_y(currentNum);
        vector< double > currentRow;
        while ( iss >> currentNum ){
            n++;
            currentRow.push_back(currentNum);
        };
        data.push_back_X(currentRow);
    };
    
    // here we didn't check the legity of the data.
    n /= m;
    
    data.setCols(n);
    data.setRows(m);
}



void Qi_IO_input :: writeData(const string & filename, const string & content, bool add)
{
    ofstream myfile;
    if (add)
        myfile.open (filename.c_str(), ios::app);
    else
        myfile.open (filename.c_str());
    
    myfile << content;
    myfile.close();
}

void Qi_IO_input :: writeData(const string & filename, const stringstream & content, bool add)
{
    ofstream myfile;
    if (add)
        myfile.open (filename.c_str(), ios::app);
    else
        myfile.open (filename.c_str());
    
    myfile << content.rdbuf();;
    myfile.close();
}

//
//void Qi_IO_input::readData(const string &filename, DataModel & data)
//{
//    int m, n; // m = #rows, n = #cols
//    ifstream f(filename);
//    f >> m >> n;
//    data.setCols(n);
//    data.setRows(m);
//    
//    for (int i = 0; i < m; i++) {
//        bitset<_CONST_INT_MAX_FACTOR_SIZE> currentRow;
//        for (int j = 0; j < n; j++) {
//            int currentBit;
//            f >> currentBit;
//            if (currentBit)
//                currentRow.set( j, currentBit );
//        }
//        data.push_back(currentRow);
//    }
//}
//
//



