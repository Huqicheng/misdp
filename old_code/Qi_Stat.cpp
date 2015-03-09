//
//  Qi_Stat.cpp
//  BNSL
//
//  Created by Qi on 4/23/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#include "Qi_Stat.h"

double Qi_Stat :: generateRandomUniform01 (double left, double right)
{
    double f = (double)rand() / RAND_MAX;
    return left + f * (right - left);
}
//
//double Qi_Stat :: generateRandomNormal01 ()
////{
////    std::default_random_engine generator;
////    std::normal_distribution<double> distribution(0.0,1.0);
////    return distribution(generator);
////    //    return std::normal_distribution<double> distribution(0.0,1.0);
////}
//{
//    return 0;
//}
