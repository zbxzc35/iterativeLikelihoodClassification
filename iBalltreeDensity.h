/* 
 * File:   iBalltreeDensity.h
 * Author: kaiwu
 *
 * Created on October 24, 2013, 11:13 AM
 */

#ifndef IBALLTREEDENSITY_H
#define IBALLTREEDENSITY_H

#include "iBalltree.h"
#include <vector>
#include <cstdlib>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#define SQRT2PI 2.50662827463100024
using namespace std;

class iBalltreeDensity: public iBalltree {
public:
    iBalltreeDensity(const char* fileName,int nTerms,vector<int>&beta,double errorTorlerance);
    virtual ~iBalltreeDensity();
    double setErrorTolerance(double errorTolerance){
        double r=_errorTolerance;
        _errorTolerance=errorTolerance;
        setPreCalculateErrorTolerance(_errorTolerance);
        return r;
    };
    int setMinBallSize(int minBallSize){
       if(minBallSize<pow(double(_nTerms),_nDimensions))
           _nMinimumBallSize=pow(double(_nTerms),_nDimensions);
       else
           _nMinimumBallSize=minBallSize;
    };
    double directEvaluate(vector<float> p);
    double combinationEvaluate(float* p);
    void init(int type=0)       ;
private:
    double _errorTolerance      ;
    int _hasInit                ;
    vector<double> _coefficients;
    vector<double> _expPart     ;
};

#endif  /* IBALLTREEDENSITY_H */

