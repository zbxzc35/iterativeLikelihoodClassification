/* 
 * File:   iBalltree.cpp
 * Author: kaiwu
 * 
 * Created on October 24, 2013, 11:12 AM
 */

#include "iBalltree.h"
#include <TBits.h>
#include <iomanip>
#include <TGraph.h>
double inverseN[100];
iBalltree::iBalltree(const char* databaseFileName) : _root(NULL) {
    _databaseFileName=databaseFileName;
    for(int i1=0;i1<100;i1++){
        inverseN[i1]=1./(double(i1)+1.);
    }
    reset();
}

bool iBalltree::preInit(const char* databaseFileName) {
    bool ret = true;
    TFile* fdb = TFile::Open(databaseFileName);
    if(fdb==NULL){
        cout<<"Can not find "<<databaseFileName<<endl;
        exit(-1);
    }
    TNamed* name = (TNamed*) fdb->Get("h");
    _xH = atof(name->GetTitle())        ;
    name = (TNamed*) fdb->Get("minP")   ;
    _minP = atoi(name->GetTitle())      ;
    name = (TNamed*) fdb->Get("maxP")   ;
    _maxP = atoi(name->GetTitle())      ;
    name = (TNamed*) fdb->Get("maxBeta");
    _maxBeta = atoi(name->GetTitle())   ;
    name = (TNamed*) fdb->Get("maxX")   ;
    _maxX = atoi(name->GetTitle())      ;
    _nX = _maxX / _xH + 1               ;
    
    double v;
    if(!_paramE1)
        delete [] _paramE1;
    if(!_paramE2)
        delete [] _paramE2;
    _paramE1 = new double[(_maxBeta + 1) * _nX];
    _paramE2 = new double[(_maxBeta + 1) * _nX];
    for (int iBeta = 0; iBeta <= _maxBeta; iBeta++) {
        cout<<"Loading ..."<<Form("pB%d", iBeta)<<endl;
        TTree* tree = (TTree*) fdb->Get(Form("pB%d", iBeta));
        tree->SetBranchAddress("v", &v);
        for (int iX = 0; iX < tree->GetEntries(); iX++) {
            tree->GetEntry(iX);
            _paramE1[iBeta * _nX + iX] = v;
        }
        tree = (TTree*) fdb->Get(Form("p%dB%d", _nTerms, iBeta));
        cout<<"Loading ..."<<Form("p%dB%d", _nTerms, iBeta)<<endl;
        tree->SetBranchAddress("v", &v);
        for (int iX = 0; iX < tree->GetEntries(); iX++) {
            tree->GetEntry(iX);
            _paramE2[iBeta * _nX + iX] = v;
        }
    }
    fdb->Close();
    return ret;
}

iBalltree::~iBalltree() {
    if (_root != NULL) {
        _root->reset()  ;
        delete _root    ;
    }
    if(_paramE1)
        delete[] _paramE1;
    if(_paramE2)
        delete[] _paramE2;
}
void iBalltree::setData(double* x,double* y,double* hx,double* hy,int nP){
    _dataPoints.push_back(x);
    _dataPoints.push_back(y);
    _bandwidth.push_back(hx);
    _bandwidth.push_back(hy);
    _nPoints=nP;
}
//reset the interior attributes

iBalltree::iBalltreeNode::iBalltreeNode(int nDimensions, int nTerms, iBalltree* rootTree,vector<int>&beta) {
    _nDimensions = nDimensions          ;
    _mean.clear();
    _mean.resize(_nDimensions, 0.)      ;
    _sigma.clear();
    _sigma.resize(_nDimensions, 0.)     ;
    _min.clear();
    _min.resize(_nDimensions, 1.e300)   ;
    _max.clear();
    _max.resize(_nDimensions, -1.e300)  ;
    _minBandwidth.clear();
    _minBandwidth.resize(_nDimensions, 1.e300)  ;
    _maxBandwidth.clear();
    _maxBandwidth.resize(_nDimensions, -1.e300) ;
    _meanBandwidth.clear();
    _meanBandwidth.resize(_nDimensions, 0.)     ;
    _leftPoints = _rightPoints = 0;
    _leftChild = _rightChild = _parent = NULL   ;
    _maxSigmaDimension = 0;
    _rootTree = rootTree;
    _nTerms = nTerms;
    _totalNTerms = pow(_nTerms, _nDimensions)   ;
    _coefficients.clear();
    _coefficients.resize(_totalNTerms, 0.)      ;
    _coefficients2.clear();
    _coefficients2.resize(_totalNTerms, 0.);
    _hasFGTCoefficient = false;
    _scale = 1.;
    _logScale = 0.;
    _beta.clear();
    _beta.resize(_nDimensions);
    for(int i1=0;i1<_nDimensions;i1++)  _beta[i1]=beta[i1];
}

void iBalltree::reset() {
    _dataPoints.clear()         ;
    _bandwidth.clear()          ;
    _nPoints = 0                ;
    //clear nodes' tree
    if (_root != NULL) {
        _root->reset()  ;
        delete _root    ;
    }
    _root = NULL;
};

void iBalltree::iBalltreeNode::calculateStatistics() {
    double maxSigma = 0;
    _scale = 1.;
    _minScale = 1.;
    _maxScale = 1.;
    for (int i1 = _leftPoints; i1 <= _rightPoints; i1++) {
        for (int i2 = 0; i2 < _nDimensions; i2++) {
            _mean[i2] += _rootTree->_dataPoints[i2][i1] ;
            _sigma[i2] += (_rootTree->_dataPoints[i2][i1] * _rootTree->_dataPoints[i2][i1]);
            if (_min[i2] > _rootTree->_dataPoints[i2][i1]) _min[i2] = _rootTree->_dataPoints[i2][i1];
            if (_max[i2] < _rootTree->_dataPoints[i2][i1]) _max[i2] = _rootTree->_dataPoints[i2][i1];

            if (_minBandwidth[i2] < _rootTree->_bandwidth[i2][i1]) _minBandwidth[i2] = _rootTree->_bandwidth[i2][i1];
            if (_maxBandwidth[i2] > _rootTree->_bandwidth[i2][i1]) _maxBandwidth[i2] = _rootTree->_bandwidth[i2][i1];
            _meanBandwidth[i2] += _rootTree->_bandwidth[i2][i1];
        }
    }
    double w=nPoints();
    for (int i1 = 0; i1 < _nDimensions; i1++) {
        _mean[i1] /= w;
        _meanBandwidth[i1] /= w;
        _sigma[i1] /= w;
        _sigma[i1] -= _mean[i1] * _mean[i1];
        _sigma[i1] = sqrt(_sigma[i1]);
        if (maxSigma < _sigma[i1]) {
            maxSigma = _sigma[i1];
            _maxSigmaDimension = i1;
        }
        _scale *= (_meanBandwidth[i1] * INVERSESQRT2PI);
        _minScale *= (_maxBandwidth[i1] * INVERSESQRT2PI);
        _maxScale *= (_minBandwidth[i1] * INVERSESQRT2PI);
    }
    _inverseHBeta=1.;
    for(int i1=0;i1<_nDimensions;i1++)
        _inverseHBeta*=pow(-1*_meanBandwidth[i1]*INVERSESQRT2,_beta[i1]);
    _logScale = log(_scale)+log(fabs(_inverseHBeta));
    _logMinScale = log(_minScale);
    _logMaxScale = log(_maxScale);
    _logWeight   = log(_weight);
////    cout<<"calculateStatistics Ok!"<<endl;
}

void iBalltree::iBalltreeNode::updateStatistics() {
    double maxSigma = 0;
    _weight = 0.;
    _scale = 1.;
    _minScale = 1.;
    _maxScale = 1.;
    if (_leftChild == NULL && _rightChild == NULL) {
        calculateStatistics();
    } else {
        _leftChild->updateStatistics();
        for (int i1 = 0; i1 < _nDimensions; i1++) {
            _minBandwidth[i1] = _leftChild->_minBandwidth[i1];
            _maxBandwidth[i1] = _leftChild->_maxBandwidth[i1];
            _meanBandwidth[i1] = _leftChild->_meanBandwidth[i1] * _leftChild->_weight;
            _min[i1] = _leftChild->_min[i1];
            _max[i1] = _leftChild->_max[i1];
            _mean[i1] = _leftChild->_mean[i1] * _leftChild->_weight;
        }
        _weight += _leftChild->_weight;
        _rightChild->updateStatistics();
        for (int i1 = 0; i1 < _nDimensions; i1++) {
            if (_rightChild->_minBandwidth[i1] > _minBandwidth[i1])
                _minBandwidth[i1] = _rightChild->_minBandwidth[i1];
            if (_rightChild->_maxBandwidth[i1] < _maxBandwidth[i1])
                _maxBandwidth[i1] = _rightChild->_maxBandwidth[i1];
            _meanBandwidth[i1] += _rightChild->_meanBandwidth[i1] * _rightChild->_weight;
            if (_rightChild->_min[i1] < _min[i1])
                _min[i1] = _rightChild->_min[i1];
            if (_rightChild->_max[i1] > _max[i1])
                _max[i1] = _rightChild->_max[i1];
            _mean[i1] = _rightChild->_mean[i1] * _rightChild->_weight;
        }
        _weight += _leftChild->_weight;
        for (int i1 = 0; i1 < _nDimensions; i1++) {
            _mean[i1]/=_weight;
            _meanBandwidth[i1] /= _weight;
            _sigma[i1]=_leftChild->_weight*(_leftChild->_sigma[i1]*_leftChild->_sigma[i1]+_leftChild->_mean[i1]*_leftChild->_mean[i1])+_rightChild->_weight*(_rightChild->_sigma[i1]*_rightChild->_sigma[i1]+_rightChild->_mean[i1]*_rightChild->_mean[i1]);
            _sigma[i1]=sqrt(_sigma[i1]/_weight-_mean[i1]*_mean[i1]);
            if(_sigma[i1]>maxSigma){
                _maxSigmaDimension=i1;
                maxSigma=_sigma[i1];
            }
            _scale *= (_meanBandwidth[i1] * INVERSESQRT2PI);
            _minScale *= (_maxBandwidth[i1] * INVERSESQRT2PI);
            _maxScale *= (_minBandwidth[i1] * INVERSESQRT2PI);
        }
        _inverseHBeta = 1.;
        for (int i1 = 0; i1 < _nDimensions; i1++)
            _inverseHBeta *= pow(-1*_meanBandwidth[i1]*INVERSESQRT2, _beta[i1]);
        
        _logScale = log(_scale)+log(fabs(_inverseHBeta));
        _logMinScale = log(_minScale);
        _logMaxScale = log(_maxScale);
    }  
}

int iBalltree::iBalltreeNode::drawBox(int color) {
    int r;
    int rl = 0, rr = 0;
    if (_nDimensions == 2) {
        if (_leftChild != NULL)
            rl = _leftChild->drawBox(color + 1);
        if (_rightChild != NULL)
            rr = _rightChild->drawBox(color + 1);
        if (_leftPoints != _rightPoints) {
            r = rl < rr ? rr : rl;
            r++;
            TBox* b = new TBox(_min[0], _min[1], _max[0], _max[1]);
            b->Draw();
            b->SetLineColor(color);
            b->SetFillStyle(0);
            b->SetLineWidth(r);
        }
    }
    return r;
}

void iBalltree::iBalltreeNode::reset() {
    if (_leftChild != NULL) {
        _leftChild->reset();
        //cout<<"delete "<<_leftChild->_leftPoints<<"->"<<_leftChild->_rightPoints<<endl;
        delete _leftChild;
    }
    if (_rightChild != NULL) {
        _rightChild->reset();
        delete _rightChild;
    }
}

void iBalltree::buildTree() {
    //calculateStatistics();
    //initialize the building...
    _root = new iBalltreeNode(_nDimensions, _nTerms, this,_beta);
    _root->_volume = 1.;  
    buildBall(0, _nPoints - 1, _root);
}

double iBalltree::iBalltreeNode::maxContributions(float* p){
    double ret=0.       ;
    float minDist[2];
    minDistance(p,minDist);
    for(int i1=0;i1<_nDimensions;i1++){
        ret+=minDist[i1]*minDist[i1]*_meanBandwidth[i1]*_meanBandwidth[i1];
    }
    ret*=(-0.5);
    return ret;
}
void iBalltree::iBalltreeNode::minDistance(float* p,float* minDist){
    for(int i1=0;i1<_nDimensions;i1++){
        double d;
        if((_min[i1]-p[i1])*(_max[i1]-p[i1])<=0.)
            d=0.;
        else
            d=fabs(_min[i1]-p[i1])<fabs(_max[i1]-p[i1])?fabs(_min[i1]-p[i1]):fabs(_max[i1]-p[i1]);
        minDist[i1]=d;
    }
}
double iBalltree::iBalltreeNode::combinationEvaluate(float* p, double errorTolerance) {
    double leftValue = 0., rightValue = 0. ;
    double ret = 0.                        ;
    double _targetCoefficients[900]        ;
        if (_leftChild == NULL && _rightChild == NULL) {
                if (_hasFGTCoefficient) {
                        double targetError = calculateTargetError(p);
                        if (targetError*_error < errorTolerance) { //Please note that here we need to count the error contribution of this node to overall error 
                                calculateTargetCoefficients(p,_targetCoefficients);
                                ret = 0.;
                                for (int i1 = 0; i1 < _totalNTerms; i1++){
                                        ret += _coefficients[i1] * _targetCoefficients[i1];
                                }
                                double s=_scale*_inverseHBeta;
                                ret*=s*targetError*targetError;
                                return ret;
                        }
                }
                for (int i1 = _leftPoints; i1 <= _rightPoints; i1++) {
                    double temp2 = 0.;
                    double temp3=1.;
                    for (int i2 = 0; i2 < _nDimensions; i2++) {
                        temp2 += (p[i2]-_rootTree->_dataPoints[i2][i1])*(p[i2]-_rootTree->_dataPoints[i2][i1])*_rootTree->_bandwidth[i2][i1] * _rootTree->_bandwidth[i2][i1];
                        temp3*=_rootTree->hermitePolynomial((p[i2]-_rootTree->_dataPoints[i2][i1])*INVERSESQRT2*_rootTree->_bandwidth[i2][i1],_beta[i2]);
                     }
                     temp2 *= (-0.5);
                     ret += exp(temp2)* temp3*_inverseHBeta*_scale;
                }
//              cout<<"directly calculate "<<nPoints()<<endl;
                return ret;
        }
        double maxValue = maxContributions(p);
        if(maxValue<_rootTree->_preCalculateErrorTolerance-_logScale)
                return 0.;
        if (_hasFGTCoefficient) {
            double targetError = calculateTargetError(p);
            if (targetError*_error  < errorTolerance) { 
                //Please note that here we need to count the error contribution of this node to overall error
                calculateTargetCoefficients(p,_targetCoefficients)      ;
                ret = 0.;
                for (int i1 = 0; i1 < _totalNTerms; i1++){
                    ret += _coefficients[i1] * _targetCoefficients[i1];
                }
                double s=_scale*_inverseHBeta   ;
                ret*=s*targetError*targetError  ;
//                cout<<"fast calculate "<<nPoints()<<endl;
                //serval tests here
//                double ret2=0.;
//                for (int i1 = _leftPoints; i1 <= _rightPoints; i1++) {
//                  double temp2 = 0.;
//                  double temp3=1.;
//                  for (int i2 = 0; i2 < _nDimensions; i2++) {
//                      temp2 += (p[i2]-_rootTree->_dataPoints[i2][i1])*(p[i2]-_rootTree->_dataPoints[i2][i1])*_rootTree->_bandwidth[i2][i1] * _rootTree->_bandwidth[i2][i1];
//                      temp3*=_rootTree->hermitePolynomial((p[i2]-_rootTree->_dataPoints[i2][i1])*INVERSESQRT2*_rootTree->_bandwidth[i2][i1],_beta[i2]);
//                 }
//                 temp2 *= (-0.5);
//                 ret2 += exp(temp2)* temp3*_inverseHBeta*_scale;
//              }
////                if(fabs(ret2-ret)>0.01){
//                cout<<"ball ("<<_leftPoints<<", "<<_rightPoints<<") rError= "<<fabs(ret2-ret)<<" vs "<<errorTolerance*_weight<<", "<<errorTolerance*_weight-fabs(ret2-ret)<<", Predicted Error= "<<targetError*_error<<"("<<targetError<<","<<_error<<")"<<endl;
//                s=0.;
//                if(fabs(ret2-ret)>1.){
//                    for (int i1 = 0; i1 < _totalNTerms; i1++){
//                        printf("%.6e*%.6e\n",_coefficients[i1] , _targetCoefficients[i1]);
//                        s+=_coefficients[i1] * _targetCoefficients[i1];
//                    }
//                }
//                cout<<"sum= "<<s<<","<<targetError*targetError<<",["<<ret<<", "<<ret2<<"]"<<(p[0] - _mean[0]) *INVERSESQRT2 * _meanBandwidth[0]<<", "<<_scale*_inverseHBeta<<endl;
//                cout<<"Ball Mean= "<<_mean[0]<<" ("<<_min[0]<<","<<_max[0]<<") mean bandwidth= "<<1./_meanBandwidth[0]<<" nPoints= "<<nPoints()<<endl;
//                }
                return ret;
            }
        }   
        if (_leftChild != NULL)
            leftValue = _leftChild->combinationEvaluate(p, errorTolerance);
        if (_rightChild != NULL)
            rightValue = _rightChild->combinationEvaluate(p, errorTolerance);
        ret = leftValue + rightValue;
        return ret;
}

double iBalltree::iBalltreeNode::calculateTargetError(float* p) {
    double ret = 1.;
    double r   = 0.;
    for (int i1 = 0; i1 < _nDimensions; i1++)
        r += ((p[i1] - _mean[i1])*(p[i1] - _mean[i1]) * _meanBandwidth[i1] * _meanBandwidth[i1]);
    ret =exp(-0.25*r)       ;
    //ret=-0.25*r         ;
    return ret          ;
}
double iBalltree::hermitePolynomial(double x, int n){
    double h0=1.;
    double h1=2*x;
    double ret=0.;
    if(n==0)
        return h0;
    if(n==1)
        return h1;
    if(n>1){
        for(int i1=2;i1<=n;i1++){
            ret=2*x*h1-2*(i1-1)*h0;
            h0=h1;
            h1=ret;
        }
    }
    return ret;
}
void iBalltree::iBalltreeNode::calculateTargetCoefficients(float* p,double* _targetCoefficients) {
    //cout<<(p[0] - _mean[0]) / (SQRT2 * _meanBandwidth[0])<<endl;
    //cout<<(p[1] - _mean[1]) / (SQRT2 * _meanBandwidth[1])<<endl;
    int c1 = 1, c2, c3 = 1;
    //calculate h_{Beta}
    double prefixZero[2]={0,0};
    double prefixOne[2]={0,0};
    _targetCoefficients[0] = 1.;
    for (int i1 = 0; i1 < _nDimensions; i1++) {
        int b = _beta[i1];
        double h0 = 1.;
        double r =(p[i1] - _mean[i1]) *INVERSESQRT2 * _meanBandwidth[i1];
        double h2 = 1.,h1=2*r;
        if (b == 0){
            prefixZero[i1]=0.;
            prefixOne[i1]=1.;
        }
        else if (b == 1){
            prefixZero[i1]=1.;
            prefixOne[i1]=h1;
        }
        else {
            for (int i2 = 2; i2 <= b; i2++) {
                h2 = 2.*r* h1 - 2.*(i2 - 1)*h0;
                h0 = h1;
                h1 = h2;
            }
            prefixZero[i1]=h0;
            prefixOne[i1]=h1;
        }
    }
    _targetCoefficients[0]=1.;
    for (int i1 = 0; i1 < _nDimensions; i1++) {
        double r = (p[i1] - _mean[i1]) *INVERSESQRT2 * _meanBandwidth[i1];
        c2 = c1;
        for (int i3 = c2 - c3; i3 < c2; i3++) {
            _targetCoefficients[c1] = _targetCoefficients[i3]*(2.*r*prefixOne[i1]-2.*_beta[i1]*prefixZero[i1]);
            _targetCoefficients[i3]*=prefixOne[i1];
            c1++;
        }
        for (int i2 = 2+_beta[i1]; i2 < _nTerms+_beta[i1]; i2++) {
            c2 = c1;
            for (int i3 = c2 - c3; i3 < c2; i3++) {
                _targetCoefficients[c1] = 2.* r * _targetCoefficients[i3] - 2.* (i2 - 1) * _targetCoefficients[i3 - c3];
                c1++;
            }
        }
        c3 *= _nTerms;
    }
//    for (int i1 = 0; i1 < _totalNTerms; i1++){
//        //cout<<_targetCoefficients[i1]<<", ";
//        _targetCoefficients[i1] *= targetExpPart;
//    }
    //exit(0);
}

//maximum contribution of this node to the final value
//bandwidth will be chosen to be the maximum
//distance will be chose to be minimum

void iBalltree::swap(int i, int j) {
    double temp;
    for (int i2 = 0; i2 < _nDimensions; i2++) {
        temp = _dataPoints[i2][i];
        _dataPoints[i2][i] = _dataPoints[i2][j];
        _dataPoints[i2][j] = temp;

        temp = _bandwidth[i2][i];
        _bandwidth[i2][i] = _bandwidth[i2][j];
        _bandwidth[i2][j] = temp;
    }
}

void iBalltree::select(int dimension, int position, int low, int high) {
    int m, r;
    while (low < high) {
        r = (low + high) / 2;
        swap(r, low);
        m = low;
        for (int i = low + 1; i <= high; i++) {
            if (_dataPoints[dimension][i] < _dataPoints[dimension][low]) {
                m++;
                swap(m, i);
            }
        }
        swap(low, m);
        if (m <= position) low = m + 1;
        if (m >= position) high = m - 1;
    }
}

void iBalltree::printInfo() {
    //currently only visualized the 2D case
    if (_nDimensions == 2) {
        TFile* f = TFile::Open("iBalltreePrintInfo.root", "recreate");
        TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
        TGraph* gr = new TGraph();
        for (int i1 = 0; i1 < _nPoints; i1++) {
            gr->SetPoint(i1, _dataPoints[0][i1], _dataPoints[1][i1]);
        }
        gr->Draw("AP*");
        gr->SetMarkerStyle(20);
        int color = 0;
        _root->drawBox(color + 1);
        c1->Write("c1");
        gr->Write("gr");
        f->Close();
    }
}

void iBalltree::iBalltreeNode::buildFGTCoefficient2() {
    double error1 = 1.;
    double error2 = 1.;
    _error = 0.;
    for (int iPoint = _leftPoints; iPoint <= _rightPoints; iPoint++) {
        error1 = error2 = 1.;
        for (int i1 = 0; i1 < _nDimensions; i1++) {
                double s = fabs(_rootTree->_dataPoints[i1][iPoint] - _mean[i1])*_meanBandwidth[i1];
                if (s >= _rootTree->_maxX||_beta[i1]>_rootTree->_maxBeta){
                        _hasFGTCoefficient = false;
                        break;
                }
                int rlow = s/_rootTree->_xH              ;
                int rindLow  =_rootTree->_nX * _beta[i1] + rlow    ;
                int rindhigh = rindLow + 1              ;
                double deltaS=s - _rootTree->_xH * rlow            ;
                
                error1 *= (_rootTree->_paramE1[rindLow]*(_rootTree->_xH-deltaS) +_rootTree-> _paramE1[rindhigh]*deltaS)/_rootTree->_xH;
                error2 *= (_rootTree->_paramE2[rindLow]*(_rootTree->_xH-deltaS) + _rootTree->_paramE2[rindhigh]*deltaS)/_rootTree->_xH;
        }
        _error += fabs(error1 - error2);
    }
    _error *= (_scale*fabs(_inverseHBeta));
    if (!_leftChild && !_rightChild) {
        //calculate the FGT coefficients
        double* tempCoefficients = new double[_totalNTerms];
        for (int i1 = 0; i1 < _totalNTerms; i1++)
            _coefficients2[i1] = 0.;
        for (int iPoints = _leftPoints; iPoints <= _rightPoints; iPoints++) {
            int c1 = 1, c2, c3 = 1;
            tempCoefficients[0] = 1.;
            _coefficients2[0] += (tempCoefficients[0] );
            for (int i1 = 0; i1 < _nDimensions; i1++) {
                double r = (_rootTree->_dataPoints[i1][iPoints] - _mean[i1])*INVERSESQRT2 * _meanBandwidth[i1];
                for (int i2 = 1; i2 < _nTerms; i2++) {
                    c2 = c1;
                    for (int i3 = c2 - c3; i3 < c2; i3++) {
                        tempCoefficients[c1] = r * (tempCoefficients[i3] *inverseN[i2]);
                        _coefficients2[c1] += (tempCoefficients[c1] );
                        c1++;
                    }
                }
                c3 *= _nTerms;
            }
        }
        delete[] tempCoefficients;
    } else {
        if (_leftChild)
            _leftChild->buildFGTCoefficient2();
        if (_rightChild)
            _rightChild->buildFGTCoefficient2();
        int complexityDirect = 2 * nPoints();
        int complexityRecursive = 9 + _totalNTerms;
        if (complexityDirect < complexityRecursive) {
            double* tempCoefficients = new double[_totalNTerms];
            for (int i1 = 0; i1 < _totalNTerms; i1++)
                _coefficients2[i1] = 0.;
            for (int iPoints = _leftPoints; iPoints <= _rightPoints; iPoints++) {
                int c1 = 1, c2, c3 = 1;
                tempCoefficients[0] = 1.;
                _coefficients2[0] += (tempCoefficients[0] );
                for (int i1 = 0; i1 < _nDimensions; i1++) {
                    double r = (_rootTree->_dataPoints[i1][iPoints] - _mean[i1]) *INVERSESQRT2 * _meanBandwidth[i1];
                    for (int i2 = 1; i2 < _nTerms; i2++) {
                        c2 = c1;
                        for (int i3 = c2 - c3; i3 < c2; i3++) {
                            tempCoefficients[c1] = r * (tempCoefficients[i3] *inverseN[i2]);
                            _coefficients2[c1] += (tempCoefficients[c1]);
                            c1++;
                        }
                    }
                    c3 *= _nTerms;
                }
            }
            delete[] tempCoefficients;
        }else {
            double* tempCoefficients1 = new double[_totalNTerms];
            double* tempCoefficients2 = new double[_totalNTerms];
            double* tempCoefficients3 = new double[_totalNTerms];
            double* tempCoefficients4 = new double[_totalNTerms];
            for (int i1 = 0; i1 < _totalNTerms; i1++)
                _coefficients2[i1] = 0.;
            int c1 = 1, c2, c3 = 1;
            tempCoefficients1[0] = 1.;
            tempCoefficients2[0] = 1.;
            tempCoefficients3[0] = 1.;
            tempCoefficients4[0] = 1.;
            for (int i1 = 0; i1 < _nDimensions; i1++) {
                tempCoefficients3[0] *= (_meanBandwidth[i1]/_leftChild->_meanBandwidth[i1]  );
                tempCoefficients4[0] *= (_meanBandwidth[i1]/_rightChild->_meanBandwidth[i1] );
            }
            for (int i1 = 0; i1 < _nDimensions; i1++) {
                double r1 = (_leftChild->_mean[i1] - _mean[i1]) *INVERSESQRT2 * _leftChild->_meanBandwidth[i1];
                double r2 = (_rightChild->_mean[i1] - _mean[i1]) *INVERSESQRT2 * _rightChild->_meanBandwidth[i1];
                double r3 = _meanBandwidth[i1]/_leftChild->_meanBandwidth[i1] ;
                double r4 = _meanBandwidth[i1]/_rightChild->_meanBandwidth[i1];
                for (int i2 = 1; i2 < _nTerms; i2++) {
                    c2 = c1;
                    for (int i3 = c2 - c3; i3 < c2; i3++) {
                        tempCoefficients1[c1] = r1 * (tempCoefficients1[i3] );
                        tempCoefficients2[c1] = r2 * (tempCoefficients2[i3] );
                        tempCoefficients3[c1] = r3 * tempCoefficients3[i3];
                        tempCoefficients4[c1] = r4 * tempCoefficients4[i3];
                        c1++;
                    }
                }
                c3 *= _nTerms;
            }
            for (int i1 = 0; i1 < _totalNTerms; i1++) {
                double lv = 0.;
                double rv = 0.;
                for (int i2 = 0; i2 <= i1; i2++) {
                    lv += tempCoefficients1[i1 - i2] * _leftChild->_coefficients2[i2];
                    rv += tempCoefficients2[i1 - i2] * _rightChild->_coefficients2[i2];
                }
                _coefficients2[i1] = lv * tempCoefficients3[i1] + rv * tempCoefficients4[i1];
            }
            delete [] tempCoefficients1;
            delete [] tempCoefficients2;
            delete [] tempCoefficients3;
            delete [] tempCoefficients4;
        }
    }
}

void iBalltree::iBalltreeNode::buildFGTCoefficient() {
    //cout<<"Ball("<<_leftPoints<<", "<<_rightPoints<<") minBallSize= "<<_rootTree->_nMinimumBallSize<<"_nDim= "<<_nDimensions<<", nTerm= "<<_nTerms<<", totNTerms= "<<_totalNTerms<<endl;
    int l;
    if (_leftChild) {
        l = _leftChild->nPoints() * _leftChild->nPoints();
        _leftChild->buildFGTCoefficient();
    }
    if (_rightChild) {
        l = _rightChild->nPoints() * _rightChild->nPoints();
        _rightChild->buildFGTCoefficient();
    }
    l = nPoints();
    if (l > 0) {
        if (!_hasFGTCoefficient) return;
        double* tempCoefficients = new double[_totalNTerms];
        for (int i1 = 0; i1 < _totalNTerms; i1++)
            _coefficients[i1] = 0.;
        for (int iPoints = _leftPoints; iPoints <= _rightPoints; iPoints++) {
            int c1 = 1, c2, c3 = 1  ;
            tempCoefficients[0] = 1.;
            _coefficients[0] += (tempCoefficients[0] );
            for (int i1 = 0; i1 < _nDimensions; i1++) {
                double r = (_rootTree->_dataPoints[i1][iPoints] - _mean[i1]) *INVERSESQRT2 * _meanBandwidth[i1];
                for (int i2 = 1; i2 < _nTerms; i2++) {
                    c2 = c1;
                    for (int i3 = c2 - c3; i3 < c2; i3++) {
                        tempCoefficients[c1] = r * (tempCoefficients[i3] *inverseN[i2-1]);
                        _coefficients[c1] += (tempCoefficients[c1] );
                        c1++;
                    }
                }
                c3 *= _nTerms;
            }
        }
        delete[] tempCoefficients;
    }
//    cout<<"Ball ("<<_leftPoints<<", "<<_rightPoints<<"), "<<" ("<<_min[0]<<", "<<_mean[0]<<", "<<_max[0]<<") h= "<<_meanBandwidth[0]<<" Error= "<<_error<<": "<<endl;
//    for(int i1=0;i1<_totalNTerms;i1++)
//        printf("%.5e ",_coefficients[i1]);
//    cout<<endl;
}

void iBalltree::buildBall(int low, int high, iBalltreeNode* root) {
//    cout<<"Build Ball ("<<low<<", "<<high<<")"<<endl; 
    root->_leftPoints = low     ;
    root->_rightPoints = high   ;
    root->calculateStatistics() ;
    //calculate Error
        double error1 = 1.;
        double error2 = 1.;
        root->_error  = 0.;
        root->_hasFGTCoefficient=true;
        for (int i1 = 0; i1 < _nDimensions; i1++) {
                double s1=fabs(root->_max[i1]-root->_mean[i1])* root->_meanBandwidth[i1];
                double s2=fabs(root->_min[i1]-root->_mean[i1])* root->_meanBandwidth[i1];
                if(s1>= _maxX||s2>= _maxX){
                        root->_hasFGTCoefficient=false;
                }
        }
        for (int iPoint = low; iPoint <= high; iPoint++) {
            if (!root->_hasFGTCoefficient)
                break;
            error1 = error2 = 1.;
            for (int i1 = 0; i1 < _nDimensions; i1++) {
                double s = fabs(_dataPoints[i1][iPoint] - root->_mean[i1])* root->_meanBandwidth[i1];
                if (s >= _maxX||_beta[i1]>_maxBeta){
                        root->_hasFGTCoefficient = false;
                        break   ;
                }
                int rlow = s / _xH                      ;
                int rindLow  =_nX * _beta[i1] + rlow    ;
                int rindhigh = rindLow + 1              ;
                double deltaS=s - _xH * rlow            ;
                error1 *= (_paramE1[rindLow]*(_xH-deltaS) +_paramE1[rindhigh]*deltaS)/_xH;
                error2 *= (_paramE2[rindLow]*(_xH-deltaS) +_paramE2[rindhigh]*deltaS)/_xH;
               
            }
//            if(error1-error2<0)
//                cout<<error1-error2<<endl;
            root->_error += (error1 - error2);
            if(root->_error<0){
                root->_hasFGTCoefficient=false;
                break;
            }
        }
        
    if(root->_hasFGTCoefficient)
        root->_error *= root->_scale*fabs(root->_inverseHBeta);
    //leaf node
    //cout<<root->_error<<", "<< _errorTolerance<<", nPnts= "<<high-low+1<<endl;
    if ((root->_rightPoints - root->_leftPoints <= _nMinimumBallSize) ||root->_hasFGTCoefficient&&(root->_error<_errorTolerance)){
        root->_leftChild = NULL ;
        root->_rightChild = NULL;
    } else {
        int middle = int(0.5*(low + high));
        select(root->_maxSigmaDimension, middle, low, high);
        iBalltreeNode* leftChild = new iBalltreeNode(_nDimensions, _nTerms, this,_beta) ;
        iBalltreeNode* rightChild = new iBalltreeNode(_nDimensions, _nTerms, this,_beta);
        //leftChild->_volume = root->_volume * (_dataPoints[root->_maxSigmaDimension][middle] - root->_min[root->_maxSigmaDimension]) / (root->_max[root->_maxSigmaDimension] - root->_min[root->_maxSigmaDimension]);
        //rightChild->_volume = root->_volume * (root->_max[root->_maxSigmaDimension] - _dataPoints[root->_maxSigmaDimension][middle]) / (root->_max[root->_maxSigmaDimension] - root->_min[root->_maxSigmaDimension]);
        //cout<<"low= "<<low<<", middle= "<<middle<<", high= "<<high<<endl;
        buildBall(low, middle, leftChild);
        root->_leftChild = leftChild;
        leftChild->_parent = root;
        buildBall(middle + 1, high, rightChild);
        root->_rightChild = rightChild;
        rightChild->_parent = root;
    }
}

