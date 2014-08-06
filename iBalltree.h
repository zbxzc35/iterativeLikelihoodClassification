/* 
 * File:   iBalltree.h
 * Author: kaiwu
 *
 * Created on October 24, 2013, 11:12 AM
 */

#ifndef IBALLTREE_H
#define IBALLTREE_H
#include <vector>
#include <iostream>
#include "math.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TBox.h"
#include <fstream>
#include "TMatrixD.h"
#include "TNamed.h"
#include "TTree.h"
#include "TMarker.h"
#include "stdlib.h"
#define SQRT2PI 2.50662827463100024
#define INVERSESQRT2PI 0.398942280401432703
#define SQRT2   1.41421356237309515
#define INVERSESQRT2 0.707106781186547462
using namespace std;
//some precalculated consts

class iBalltree {
public:
    iBalltree(const char* databaseFileName);
    virtual ~iBalltree();
    class iBalltreeNode{
    public:
          iBalltreeNode(int nDimensions,int nTerms,iBalltree* rootTree,vector<int>&beta);
          ~iBalltreeNode(){}         ;
          void calculateStatistics() ;
          void updateStatistics()    ;
          void reset()               ;
          int drawBox(int color)     ;
          double combinationEvaluate(float* p,double errorTolerance);
          void  calculateTargetCoefficients(float* p,double* c);
          double calculateTargetError(float* p) ;
          double maxContributions(float* p) ;
          void minDistance(float* p,float* minDist);
          int nPoints(){return _rightPoints-_leftPoints+1;};
          void buildFGTCoefficient()    ;
          void buildFGTCoefficient2()   ;
          
          //statistics of this nodes
          vector<double> _mean  ;
          vector<double> _sigma ;
          vector<double> _min   ;
          vector<double> _max   ;
          vector<double> _minBandwidth  ;
          vector<double> _maxBandwidth  ;
          vector<double> _meanBandwidth ;
          vector<int> _beta             ;
          double _maxScale              ;
          double _logMaxScale           ;
          double _minScale              ;
          double _logMinScale           ; 
          double _scale                 ;
          double _logScale              ;
          double _inverseHBeta          ;
          double _weight     ;
          double _logWeight  ;
          double _volume     ;
          double _probability;
          int _maxSigmaDimension;
          int _nDimensions      ;
          
          int _leftPoints       ;
          int _rightPoints      ;
          
          //variables of Fast Gauss Transform
          bool _hasFGTCoefficient       ;
          int _nTerms                   ;
          int _totalNTerms              ;
          
          vector<double> _coefficients          ;
          vector<double> _coefficients2         ;
          
          double _error                         ;
                  
          iBalltreeNode* _leftChild     ;
          iBalltreeNode* _rightChild    ;
          iBalltreeNode* _parent        ;
          iBalltree*    _rootTree       ;
    };
    //reset the interior attributes
    void reset()                ;
    void buildTree()            ;
    void setData(double* x,double* y,double* hx,double* hy,int nP);
    double hermitePolynomial(double x,int n);
    
    void printInfo()                                ;
    bool preInit(const char* databaseFileName)         ;
    
    void setPreCalculateErrorTolerance(double errorTolerance){ 
        _preCalculateErrorTolerance=log(errorTolerance);
        _errorTolerance=errorTolerance;
    }
    vector<double* > _dataPoints;
    vector<double* > _bandwidth ;    
    iBalltreeNode* _root        ;
    vector<int>    _beta        ;
    double _totalWeights        ;
    //dimension of data points
    int _nDimensions            ;
    //number of data points
    int _nPoints                ;
    int _nMinimumBallSize       ;
    int _nTerms                 ;
private:
    //swap 2 data points
    void swap(int i,int j)                                      ;
    //re-arrange the data points array
    void select(int dimension,int position, int low,int high)   ;
    //build ball
    void buildBall(int low,int high,iBalltreeNode* root)        ;
    
    //parameters for the interpolation function
    double* _paramE1     ;
    double* _paramE2     ;
    double  _xH          ;
    double  _maxX        ;
    int     _minP        ;
    int     _maxP        ;
    int     _maxBeta     ;
    int     _nX                 ;
    string _databaseFileName    ;
    
    double _preCalculateErrorTolerance  ; 
    double _errorTolerance              ;
};
#endif  /* IBALLTREE_H */