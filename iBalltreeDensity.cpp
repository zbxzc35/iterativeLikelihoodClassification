/* 
 * File:   iBalltreeDensity.cpp
 * Author: kaiwu
 * 
 * Created on October 24, 2013, 11:13 AM
 */
#include "iBalltreeDensity.h"
#include "TProfile.h"
iBalltreeDensity::iBalltreeDensity(const char* fileName,int nTerms,vector<int>&beta,double errorTorlerance):iBalltree(fileName),_errorTolerance(errorTorlerance),_hasInit(false){
    _nTerms=nTerms;
    _nDimensions=beta.size();
    _beta.clear();
    _beta.resize(_nDimensions);
    setPreCalculateErrorTolerance(_errorTolerance);
    for(int i1=0;i1<beta.size();i1++){
        _beta[i1]=beta[i1];
    }
    cout<<_nTerms<<endl;
    preInit(fileName);
}

void iBalltreeDensity::init(int type) {
    _nMinimumBallSize=log(_nPoints)+6;
    buildTree();
    cout<<"build Tree OK!"<<endl;
    switch(type) {
        case 0:
            _root->buildFGTCoefficient();
            break;
        case 1:
        case 2:
            break;
    }
    _hasInit=true;
}
iBalltreeDensity::~iBalltreeDensity() {

}
double iBalltreeDensity::directEvaluate(vector<float> p){
    double ret=0.;
    if(!_hasInit){
        cout<<"The structure has not been initialized!"<<endl;
        return -1;
    }
    if(_nDimensions!=p.size()){
        cout<<"Need the same dimensions!"<<endl;
        return 0.;
    }
    double Hy=1.;
    for(int i1=0;i1<_nPoints;i1++){
        double r=0.;
        double h=1.;
        Hy=1.      ;
        for(int i2=0;i2<_nDimensions;i2++){
            double s=(p[i2]-_dataPoints[i1][i2])*_bandwidth[i1][i2];
            r+=s*s;
            h*=(INVERSESQRT2PI*_bandwidth[i1][i2]);
            Hy*=hermitePolynomial(s/SQRT2,_beta[i2])*pow(-1*INVERSESQRT2*_bandwidth[i1][i2],_beta[i2]);
        }
        r*=-0.5;
        ret+=exp(r)*h*Hy;
    }
    
    return 1./_nPoints*ret;
}
double iBalltreeDensity::combinationEvaluate(float* p){
    double ret=0.       ;
    if(!_hasInit){
        cout<<"The structure has not been initialized!"<<endl;
        return -1;
    }
    ret=_root->combinationEvaluate(p,_errorTolerance);
    return ret;
}