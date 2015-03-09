//
//  Qi_Math.cpp
//  BNSL
//
//  Created by Qi on 4/23/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

#include "Qi_Math.h"






double Qi_Math :: square (double x)
{
    return x*x;
}


long double Qi_Math :: square (long double x)
{
    return x*x;
}
//
//
//double Qi_Math :: min (vector<double> a)
//{
//    return x*x;
//}


//double Qi_Math :: sum(double a[])
//{
//	double b = 0.0;
//	for (int i = 0; i < sizeof(a) / sizeof(a[0]); i++)
//	{
//		b += a[i];
//	}
//	return b;
//}
//

//    // reference http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
//    template <typename T>
//    vector<size_t> Qi_Math :: sort_indexes(const vector<T> &v) {
//        // initialize original index locations
//        vector<size_t> idx(v.size());
//        for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
//        
//        // sort indexes based on comparing values in v
//        sort(idx.begin(), idx.end(),
//             [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
//        
//        return idx;
//    }
//    //Now you can use the returned index vector in iterations such as
//    //for (auto i: sort_indexes(v)) {
//    //    cout << v[i] << endl;
//    //}
//
//
//
//
//
//    vector<int> Qi_Math :: sort_IloNumArray_sub_indexes(const IloNumArray &v, const int from=0, int to=-1) {
//        // initialize original index locations
//        if (to == -1)
//            to = (int) v.getSize();
//        
//        vector<int> idx(to-from);
//        for (int i = 0; i != idx.size(); ++i)
//            idx[i] = i;
//        
//        // sort indexes based on comparing values in v
//        sort(idx.begin(), idx.end(),
//             [&v, from](int i1, int i2) {return v[from + i1] < v[from + i2];});
//        
//        return idx;
//    }


