/* 
 * File:   kMvaClassification.cxx
 * Author: kaiwu
 * 
 * Created on May 04, 2013, 4:37 PM
 */
#include <TH2.h>
#include <TMatrixDSymEigen.h>
#include <TH1.h>
#include <TAxis.h>
#include <TObject.h>
#include <TFile.h>
#include <TRandom.h>
#include <TH2D.h>
#include <omp.h>
#include <time.h> 
#include "kMvaClassification.h"
//External procedure: Parse the classification expression 
//Input:
// first        :       Classification Expression
// second       :       Reverse Polish Notation
//Return        :       0---> Ok
//                      1---> False
extern "C" int parseExpression(const char* , char* );
/*Contruction function
 input variable: configFile ---> configuration file name
 */
kMvaClassification::kMvaClassification(string configFile){
        _configFileName=configFile;
        bool ret=readConfig()     ;
        if(!ret){
                cout<<"Configuration failed, please check your configuration file!"<<endl;
                exit(0);
        }
        printConfig();
        fx=NULL;
        fy=NULL;
        fhx=NULL;
        fhy=NULL;
        fx0=NULL;
        fy0=NULL;
        fz0=NULL;
        if(_iBT!=NULL)  delete _iBT;
        vector<int> beta(2,0);
        _iBT=new iBalltreeDensity(_dataBase.c_str(),_nTerms,beta,_epsilon);
        _iBT->preInit(_dataBase.c_str());
}
kMvaClassification::~kMvaClassification(){
        if(fx!=NULL) delete[] fx        ;
        if(fy!=NULL) delete[] fy        ;
        if(fhx!=NULL) delete[] fhx      ;
        if(fhy!=NULL) delete[] fhy      ;
        if(fx0!=NULL) delete[] fx0      ;
        if(fy0!=NULL) delete[] fy0      ;
        if(fz0!=NULL) delete[] fz0      ;

}
//Initialize Data Structure
void kMvaClassification::initData(){
        _inputSigFile =TFile::Open(_inputSigFileName.c_str());
        _inputBkgFile =TFile::Open(_inputBkgFileName.c_str());
        _outputFile=TFile::Open(_outputFileName.c_str(),"recreate");
        _outputFile->cd();
        _sigTree   =new TTree("sigTree_1","sigTree_1");
        _bkgTree   =new TTree("bkgTree_1","bkgTree_1");

        float* valTemp;
        valTemp=new float[_nVariable+_nSeperator+1];

        TTree* _sigTree1  =(TTree*)_inputSigFile->Get(_sigTreeName.c_str());
        TTree* _bkgTree1  =(TTree*)_inputBkgFile->Get(_bkgTreeName.c_str());
        if(!_sigTree||!_bkgTree){
                cout<<"Can not load "<<_sigTreeName<<", from "<<_inputSigFileName<<endl;
                cout<<"Can not load "<<_bkgTreeName<<", from "<<_inputBkgFileName<<endl;
                exit(0);
        }
        //Bind variables
        for(int i1=0;i1<_nVariable;i1++){
                _sigTree1->SetBranchAddress(_variableNames[i1].c_str(),&(valTemp[i1]));
                _bkgTree1->SetBranchAddress(_variableNames[i1].c_str(),&(valTemp[i1]));
                _sigTree->Branch(_variableNames[i1].c_str(),&(valTemp[i1]),Form("%s/F",_variableNames[i1].c_str()));
                _bkgTree->Branch(_variableNames[i1].c_str(),&(valTemp[i1]),Form("%s/F",_variableNames[i1].c_str()));
        }
        //Bind seperators
        for(int i1=0;i1<_nSeperator;i1++){
                _sigTree1->SetBranchAddress(_seperatorNames[i1].c_str(),&(valTemp[_nVariable+i1]));
                _bkgTree1->SetBranchAddress(_seperatorNames[i1].c_str(),&(valTemp[_nVariable+i1]));
                _sigTree->Branch(_variableNames[i1].c_str(),&(valTemp[i1]),Form("%s/F",_variableNames[i1].c_str()));
                _bkgTree->Branch(_variableNames[i1].c_str(),&(valTemp[i1]),Form("%s/F",_variableNames[i1].c_str()));
        }
        //Covariance matrix
        covSig=new TMatrixDSym(_nVariable)                              ;
        TMatrixD *covSigV=new TMatrixD(_nVariable,_nVariable)           ;
        mSig  =new TVectorD(_nVariable)                 ;
        deSigR=new TMatrixD(_nVariable,_nVariable)      ;
        for(int i1=0;i1<_nVariable;i1++){
                for(int i2=0;i2<_nVariable;i2++){
                        (*covSigV)(i1,i2)=0.;
                        (*deSigR)(i1,i2)=0.;
                }
                (*mSig)(i1)=0.;
        }

        VAR v;
        v.min=1.0e300   ;
        v.max=-1.0e300  ;
        v.mean=0.       ;
        v.sigma=0.      ;
        v.transformType=_transFormType;

        for(int i1=0;i1<_nVariable;i1++){
                varInfoSig[_variableNames[i1]]=v;
                varInfoBkg[_variableNames[i1]]=v;
        }

        int entriesSig=_sigTree1->GetEntries();
        int entriesBkg=_bkgTree1->GetEntries();
        //Load data ...
        for(int i1=0;i1<entriesSig;i1++){
                _sigTree1->GetEntry(i1);
                for(int i2=0;i2<_nVariable;i2++){
                        (*mSig)(i2)=(((*mSig)(i2))*i1)/(1.0+i1)+valTemp[i2]/(1.0+i1);
                        for(int i3=0;i3<_nVariable;i3++){
                                (*covSigV)(i2,i3)=(((*covSigV)(i2,i3))*i1)/(1.0+i1)+valTemp[i2]*valTemp[i3]/(1.0+i1);
                        }
                }
                _sigTree->Fill();
        }
        for(int i1=0;i1<_nVariable;i1++){
                for(int i2=0;i2<_nVariable;i2++){
                        (*covSigV)(i1,i2)-=(*mSig)(i1)*(*mSig)(i2);
                        (*covSig)(i1,i2)=(*covSigV)(i1,i2);
                }
        }
        for(int i1=0;i1<entriesBkg;i1++){
                _bkgTree1->GetEntry(i1);
                _bkgTree->Fill();
        }

        TBranch** brTempSig     ;
        TBranch** brTempBkg     ;
        float* valTemp2 ;
        brTempSig=new TBranch*[_nVariable+1];
        brTempBkg=new TBranch*[_nVariable+1];
        valTemp2 =new float[_nVariable+1]  ; 
        v.transformType=0               ;
        //binding TBranch of decorrelated variables 
        for(int i1=0;i1<_nVariable;i1++){
                brTempSig[i1]=_sigTree->Branch(Form("%sD",_variableNames[i1].c_str()),&(valTemp2[i1]),Form("%sD/F",_variableNames[i1].c_str()));
                brTempBkg[i1]=_bkgTree->Branch(Form("%sD",_variableNames[i1].c_str()),&(valTemp2[i1]),Form("%sD/F",_variableNames[i1].c_str()));
                varInfoSig[Form("%sD",_variableNames[i1].c_str())]=v;
                varInfoBkg[Form("%sD",_variableNames[i1].c_str())]=v;
        }

        //get decorrelation matrix
        //_transFormType ==0,  no tranformation performed
        //_transFormType ==1,  covarivance matrix based decorrelation
        //_transFormType ==2,  principle component based decorrelation

        if(_transFormType==0){
                for(int i1=0;i1<_nVariable;i1++){
                        for(int i2=0;i2<_nVariable;i2++){
                                (*deSigR)(i1,i2)=0.;
                                if(i1==i2)
                                        (*deSigR)(i1,i2)=1.0;
                        }
                }
        }

        if(_transFormType==1){
                TMatrixDSymEigen* eigen = new TMatrixDSymEigen(*covSig);
                // D = ST C S
                TMatrixD* si = new TMatrixD(eigen->GetEigenVectors());
                TMatrixD* s = new TMatrixD(*si); // copy
                si->Transpose(*si); // invert (= transpose)

                // diagonal matrices
                TMatrixD* d = new TMatrixD(_nVariable, _nVariable);
                d->Mult((*si), *covSig);
                (*d)*= (*s)         ;
                // make exactly diagonal
                for (int i = 0; i < _nVariable; i++) for (int j = 0; j < _nVariable; j++) if (j != i) (*d)(i, j) = 0;
                // compute the square-root C' of covariance matrix: C = C'*C'
                for (int i = 0; i < _nVariable; i++) (*d)(i, i) = TMath::Sqrt((*d)(i, i));

                deSigR->Mult((*s), (*d))        ;
                (*deSigR) *= (*si)      ;
                // invert square-root matrices
                deSigR->Invert();

                delete eigen        ;
                delete si           ;
                delete s            ;
                delete d            ;
        }
        if(_transFormType==2){
                TMatrixDSymEigen* eigen = new TMatrixDSymEigen(*covSig);
                TMatrixD* si = new TMatrixD(eigen->GetEigenVectors());
                for(int i1=0;i1<_nVariable;i1++)
                        for(int i2=0;i2<_nVariable;i2++)
                                (*deSigR)(i1,i2)=(*si)(i1,i2);
                delete eigen    ;
                delete si       ;
        }
        //Performan decorrelation
        cout<<"Decorrelation Matrix: "<<endl;
        deSigR->Print();
        cout<<"mean value: "<<endl;
        mSig->Print();
        //Performan Decorrelation of Signal Samples and collect information(mean/sigma,etc) of decorrelated varibales 
        for(int i1=0;i1<entriesSig;i1++){
                _sigTree1->GetEntry(i1);
                for(int i2=0;i2<_nVariable;i2++){
                        valTemp2[i2]=0.;
                        for(int i3=0;i3<_nVariable;i3++)
                                valTemp2[i2]+=(*deSigR)(i2,i3)*(valTemp[i3]-(*mSig)(i3));
                        brTempSig[i2]->Fill();
                        if(valTemp2[i2]<varInfoSig[Form("%sD",_variableNames[i2].c_str())].min) varInfoSig[Form("%sD",_variableNames[i2].c_str())].min=valTemp2[i2];
                        if(valTemp2[i2]>varInfoSig[Form("%sD",_variableNames[i2].c_str())].max) varInfoSig[Form("%sD",_variableNames[i2].c_str())].max=valTemp2[i2];
                        double v=varInfoSig[Form("%sD",_variableNames[i2].c_str())].mean;
                        varInfoSig[Form("%sD",_variableNames[i2].c_str())].mean=v*i1/(1.0+i1)+valTemp2[i2]/(1.0+i1);
                        v=varInfoSig[Form("%sD",_variableNames[i2].c_str())].sigma      ;
                        varInfoSig[Form("%sD",_variableNames[i2].c_str())].sigma=v*i1/(1.0+i1)+valTemp2[i2]*valTemp2[i2]/(1.0+i1);
                }
        }
        for(int i1=0;i1<_nVariable;i1++){
                varInfoSig[Form("%sD",_variableNames[i1].c_str())].sigma-=(varInfoSig[Form("%sD",_variableNames[i1].c_str())].mean*varInfoSig[Form("%sD",_variableNames[i1].c_str())].mean);
                varInfoSig[Form("%sD",_variableNames[i1].c_str())].sigma=sqrt(varInfoSig[Form("%sD",_variableNames[i1].c_str())].sigma);
        }
        //Performan Decorrelation of Background Samples and collect information(mean/sigma,etc) of decorrelated varibales 
        for(int i1=0;i1<entriesBkg;i1++){
                _bkgTree1->GetEntry(i1);
                for(int i2=0;i2<_nVariable;i2++){
                        valTemp2[i2]=0.;
                        for(int i3=0;i3<_nVariable;i3++)
                                valTemp2[i2]+=((*deSigR)(i2,i3)*(valTemp[i3]-(*mSig)(i3)));
                        brTempBkg[i2]->Fill();
                        if(valTemp2[i2]<varInfoBkg[Form("%sD",_variableNames[i2].c_str())].min) varInfoBkg[Form("%sD",_variableNames[i2].c_str())].min=valTemp2[i2];
                        if(valTemp2[i2]>varInfoBkg[Form("%sD",_variableNames[i2].c_str())].max) varInfoBkg[Form("%sD",_variableNames[i2].c_str())].max=valTemp2[i2];
                        double v=varInfoBkg[Form("%sD",_variableNames[i2].c_str())].mean;
                        varInfoBkg[Form("%sD",_variableNames[i2].c_str())].mean=v*i1/(1.0+i1)+valTemp2[i2]/(1.0+i1);
                        v=varInfoBkg[Form("%sD",_variableNames[i2].c_str())].sigma      ;
                        varInfoBkg[Form("%sD",_variableNames[i2].c_str())].sigma=v*i1/(1.0+i1)+valTemp2[i2]*valTemp2[i2]/(1.0+i1);
                }
        }
        for(int i1=0;i1<_nVariable;i1++){
                varInfoBkg[Form("%sD",_variableNames[i1].c_str())].sigma-=(varInfoBkg[Form("%sD",_variableNames[i1].c_str())].mean*varInfoBkg[Form("%sD",_variableNames[i1].c_str())].mean);
                varInfoBkg[Form("%sD",_variableNames[i1].c_str())].sigma=sqrt(varInfoBkg[Form("%sD",_variableNames[i1].c_str())].sigma);
        }
}
//Print out configuration
void kMvaClassification::printConfig(){
        cout<<"################ Configure #######################"<<endl;
        cout<<"inputSigFile "<<_inputSigFileName<<endl;
        cout<<"inputBkgFile "<<_inputBkgFileName<<endl;
        cout<<"sigTree   "<<_sigTreeName<<endl  ;
        cout<<"bkgTree   "<<_bkgTreeName<<endl  ;
        cout<<"Total "<<_nVariable<<" training variables: ";
        for(int i1=0;i1<_nVariable;i1++)
                cout<<_variableNames[i1]<<", ";
        cout<<endl;
        cout<<"Total "<<_nSeperator<<" sparator variables: ";  
        for(int i1=0;i1<_nSeperator;i1++)
                cout<<_seperatorNames[i1]<<", "; 
        cout<<endl;
        cout<<"transFormType "<<_transFormType<<endl;
        cout<<"densityNBins  "<<_densityNBins<<endl;
        cout<<"cachePdfNBins "<<_cachePdfNBins<<endl;
        cout<<"nThreads      "<<_nThreads<<endl         ;
        cout<<"plotROC       "<<_plotRoc<<endl          ;
        cout<<"expression    "<<_expression<<endl       ;
        //cout<<"internal format "<<_expressionP<<endl  ;
        cout<<"outputFileName "<<_outputFileName<<endl  ;
        cout<<"dataBase       "<<_dataBase<<endl;
        cout<<"nTerms         "<<_nTerms<<endl;
        cout<<"epsilon        "<<_epsilon<<endl;
        cout<<"auto           "<<_auto<<endl            ; 
        cout<<"##################################################"<<endl;
}
//Read Configuration file
//Return:       true ---> Configuration file is Ok
//                      false---> Configuration file has some error
bool kMvaClassification::readConfig(){
        bool ret=true;
        char readWords[1000];
        FILE* fin=fopen(_configFileName.c_str(),"r");
        _nVariable=0            ;
        _nSeperator=0           ;
        _transFormType=2        ;
        _densityNBins=100       ;
        _cachePdfNBins=100      ;
        _nThreads=8             ;
        _doThecondPass=0        ;
        _plotRoc=1              ;
        _auto=1                 ;

        if(!fin){
                cout<<"Can not open "<<_configFileName<<endl;
                ret=false;
        }
        else{
                fscanf(fin,"%s",readWords);
                while(!feof(fin)){
                        if(strcmp(readWords,"inputSigFile")==0){
                                fscanf(fin,"%s",readWords);
                                _inputSigFileName=readWords;
                        }
                        else if(strcmp(readWords,"inputBkgFile")==0){
                                fscanf(fin,"%s",readWords);
                                _inputBkgFileName=readWords;
                        }
                        else if(strcmp(readWords,"sigTree")==0){
                                fscanf(fin,"%s",readWords);
                                _sigTreeName=readWords;
                        }
                        else if(strcmp(readWords,"bkgTree")==0){
                                fscanf(fin,"%s",readWords);
                                _bkgTreeName=readWords;
                        }
                        else if(strcmp(readWords,"nVariable")==0){
                                fscanf(fin,"%s",readWords);
                                _nVariable=atoi(readWords);
                        }
                        else if(strcmp(readWords,"variableNames")==0){
                                for(int i1=0;i1<_nVariable;i1++){
                                        fscanf(fin,"%s",readWords);
                                        _variableNames.push_back(readWords);
                                }
                        }
                        else if(strcmp(readWords,"nSeperator")==0){
                                fscanf(fin,"%s",readWords);
                                _nSeperator=atoi(readWords);
                        }
                        else if(strcmp(readWords,"seperatorNames")==0){
                                for(int i1=0;i1<_nSeperator;i1++){
                                        fscanf(fin,"%s",readWords);
                                        _seperatorNames.push_back(readWords);
                                }
                        }
                        else if(strcmp(readWords,"transFormType")==0){
                                fscanf(fin,"%s",readWords);
                                _transFormType=atoi(readWords);
                        }
                        else if(strcmp(readWords,"densityNBins")==0){
                                fscanf(fin,"%s",readWords);
                                _densityNBins=atoi(readWords);
                        }
                        else if(strcmp(readWords,"cachePdfNBins")==0){
                                fscanf(fin,"%s",readWords);
                                _cachePdfNBins=atoi(readWords);
                        }
                        else if(strcmp(readWords,"nThreads")==0){
                                fscanf(fin,"%s",readWords);
                                _nThreads=atoi(readWords);
                        }
                        else if(strcmp(readWords,"expression")==0){
                                fscanf(fin,"%s",readWords);
                                _expression=readWords;
                                char temp[1000];
                                if(parseExpression(readWords,temp)!=0){
                                        cout<<"Classification expression syntax error!"<<endl;
                                        ret=false;
                                        break;
                                }
                                _expressionP=temp;
                        }
                        else if(strcmp(readWords,"outputFile")==0){
                                fscanf(fin,"%s",readWords);
                                _outputFileName=readWords;
                        }
                        else if(strcmp(readWords,"doThecondPass")==0){
                                fscanf(fin,"%s",readWords);
                                _doThecondPass=atoi(readWords);
                        }
                        else if(strcmp(readWords,"auto")==0){
                                fscanf(fin,"%s",readWords);
                                _auto=atoi(readWords);
                        }
                        else if(strcmp(readWords,"plotRoc")==0){
                                fscanf(fin,"%s",readWords);
                                _plotRoc=atoi(readWords);
                        }
                        else if(strcmp(readWords,"dataBase")==0){
                                fscanf(fin,"%s",readWords);
                                _dataBase=readWords;
                        }
                        else if(strcmp(readWords,"nTerms")==0){
                                fscanf(fin,"%s",readWords);
                                _nTerms=atoi(readWords);
                        }
                        else if(strcmp(readWords,"epsilon")==0){
                                fscanf(fin,"%s",readWords);
                                _epsilon=atof(readWords);
                        }
                        else{
                                cout<<"Can not find configuration keyword "<<readWords<<endl;
                                ret=false;
                                break;
                        }
                        fscanf(fin,"%s",readWords);
                }
        }
        fclose(fin);
        return ret;
}

//Calculate the bandwidth of kernel density estimator
//Input:
//      _nPoints        : Number of data points
//      _x/_y           : x/y coordinates of the data points
//      _sigmaX/Y       : variance of x/y
//      _histPDF        : raw histogram PDF as an approxmination of PDF of the data points
//Output:
// _hx/_hy              : Bandwidth of data points.
void kMvaClassification::calculateBindWidth(int _nPoints,double* _x,double* _y,double _sigmaX,double _sigmaY,TH2D* _histPdf,double* _hx,double* _hy){
        double nd =pow(_nPoints,-1./6.);

        double sigSum       = _sigmaX*_sigmaX + _sigmaY*_sigmaY;
        double sqrtSum      = sqrt( sigSum );
        double sigProd      = _sigmaX*_sigmaY;

        double h  =nd*sqrt(sigSum/sigProd);
        double xhmin=h*_sigmaX*sqrt(2.)/10;
        double yhmin=h*_sigmaY*sqrt(2.)/10;
        double hSigmaX =h*pow(_sigmaX/sqrtSum,1.5);
        double hSigmaY =h*pow(_sigmaY/sqrtSum,1.5);
        double minX,maxX,minY,maxY;
        minX=_histPdf->GetXaxis()->GetBinLowEdge(1);
        maxX=_histPdf->GetXaxis()->GetBinUpEdge(_densityNBins);
        minY=_histPdf->GetYaxis()->GetBinLowEdge(1);
        maxY=_histPdf->GetYaxis()->GetBinUpEdge(_densityNBins);
        for (Int_t i = 0; i < _nPoints; i++){
                double f=0.;
                if(_x[i]<=minX||_x[i]>=maxX||_y[i]<=minY||_y[i]>=maxY)
                        f=0.;
                else
                        f=_histPdf->Interpolate(_x[i],_y[i]);
                _hx[i]=1./(hSigmaX*pow(f,-0.25));
                _hy[i]=1./(hSigmaY*pow(f,-0.25));
                if(_hx[i]<xhmin) _hx[i]=1./xhmin;
                if(_hy[i]<yhmin) _hy[i]=1./yhmin;
        }
}

//Two dimensinal kernel density estimation
//Input:
//      nPoints : number of data points
//      x               : x coordinate of the data points
//      y               : y coordinate of the data points
//      hx              : bandwidth of x coordinate of the data points
//      hy              : bandwifth of y coordinate of the data points
//      xmin/max: minimum/maximum of x coordinate
//      ymin/max: minimum/maximum of y coordinate
//Output:
//      mEval   : number of cached points
//      x0/y0   : x/y coordinates of the cached points
//      z0              : density of the cached points
//      nthreads: number of parallel threads
//Return        : Time(s) used
float kMvaClassification::kde2d(int nPoints,double*x,double* y,double* hx,double* hy,double xmin,double xmax,double ymin,double ymax,int mEval,const double* x0,const double* y0,double* z0,int nthreads){
        int i1  ;
        int i2=0 ;
        int m1   ;
        double _f;
        double rx2,ry2,zx,zy            ;
        _iBT->reset()                   ;
        _iBT->setData(x,y,hx,hy,nPoints);
        cout<<"set OK!"<<endl           ;
        _iBT->init(0)                   ;
        cout<<"init OK!"<<endl        ;
        //cout<<"Set nthreads "<<nthreads<<endl;
        omp_set_num_threads(nthreads);
        //cout<<"Total number of thread= "<<omp_get_max_threads()<<endl;
        time_t startTime,endTime;
        time(&startTime);
        //Evaluate the density for each cached points
//        #pragma omp parallel private(i1,i2,_f,rx2,ry2,zx,zy)
        #pragma omp parallel
        {
                #pragma omp for nowait 
                for(m1=0;m1<mEval;m1++){
//                        _f=0.;
//                        for(i1=0;i1<nPoints;i1++){
//                                zx=0.;zy=0.;rx2=0.;ry2=0.;
//                                if(hx[i1]!=0.0) rx2=(x0[m1]-x[i1])*hx[i1];
//                                if(hy[i1]!=0.0) ry2=(y0[m1]-y[i1])*hy[i1];
//
//                                if(hx[i1]!=0.0) zx=exp(-0.5*rx2*rx2)*hx[i1];
//                                if(hy[i1]!=0.0) zy=exp(-0.5*ry2*ry2)*hy[i1];
//                                _f+=zy*zx;
//                        }
//                        i2++     ;
//                        z0[m1]=_f;
                         float data[2] ;
                        data[0]=x0[m1];
                        data[1]=y0[m1];
                        z0[m1]=_iBT->combinationEvaluate(data);
//                        cout<<"m1= "<<m1<<endl;
//                        cout<<"z0= "<<z0[m1]<<", fgt= "<<_iBT->combinationEvaluate(data)<<", direct= "<<_iBT->directEvaluate(data)<<endl;
                }
        }
        time(&endTime)                          ;
        return difftime(endTime,startTime)      ;
}

//Build two dimension pdf for variable var1 and var2
//Input:
//      _var1/2/_newvar : variables' name
//      _varInfo                : Information of the input variables
//      _treeIn                 : Data Points TTree
//Output:
//      _pdfs                   : PDF array
void kMvaClassification::buildPdf2Parallel(const char* _var1, const char* _var2, const char* _newvar,map<string,VAR>& _varInfo, TTree* _treeIn, map<string,TH2D*>& _pdfs){
        string newVar(_newvar)      ;
        double varXmin,varYmin,varXmax,varYmax,varXSigma,varYSigma,varXMean,varYMean;
        float varXVal, varYVal      ;
        double varXD,varYD          ;
        TBranch* varXB              ;
        TBranch* varYB              ;
        _treeIn->SetBranchAddress(_var1,&varXVal,&varXB);
        _treeIn->SetBranchAddress(_var2,&varYVal,&varYB);
        varXmin=_varInfo[_var1].min     ;
        varXmax=_varInfo[_var1].max     ;
        varXMean=_varInfo[_var1].mean   ;
        varXSigma=_varInfo[_var1].sigma ;
        varYmin=_varInfo[_var2].min     ;
        varYmax=_varInfo[_var2].max     ;
        varYMean=_varInfo[_var2].mean   ;
        varYSigma=_varInfo[_var2].sigma ;
        int v1TransformType=_varInfo[_var1].transformType;
        int v2TransformType=_varInfo[_var2].transformType;
        TH2D* pdfHist= _pdfs[newVar.c_str()]    ;
        TH2D* pdfTemp=new TH2D("pdfTemp","pdfTemp",_densityNBins,varXmin-1./_densityNBins*(varXmax-varXmin),varXmax+1./_densityNBins*(varXmax-varXmin),_densityNBins,varYmin-1./_densityNBins*(varYmax-varYmin),varYmax+1./_densityNBins*(varYmax+varYmin));
        int entries=_treeIn->GetEntries()       ;
        int validateEntries=0                                   ;
        //Collect validate data points to build the two dimension pdf.
        for(int i1=0;i1<entries;i1++){
                varXB->GetEntry(i1);
                varYB->GetEntry(i1);
                if(v1TransformType==-1&&(varXVal==0||varXVal==1))
                        continue;
                if(v2TransformType==-1&&(varYVal==0||varYVal==1))
                        continue;
                pdfTemp->Fill(varXVal,varYVal)  ;
                fx[validateEntries]=varXVal     ;
                fy[validateEntries]=varYVal     ;
                validateEntries++;
        }
        //Build a histogram PDF as an approximation 
        pdfTemp->Smooth(1);
        double tot=pdfTemp->Integral();
        for(int i3=1;i3<=pdfTemp->GetXaxis()->GetNbins();i3++){
                for(int i4=1;i4<=pdfTemp->GetYaxis()->GetNbins();i4++){
                        double vp=pdfTemp->GetBinContent(i3,i4)         ;
                        double wx=pdfTemp->GetXaxis()->GetBinWidth(i3)  ;
                        double wy=pdfTemp->GetYaxis()->GetBinWidth(i4)  ;
                        double s =wx*wy                                 ;
                        if(vp>0)
                                pdfTemp->SetBinContent(i3,i4,vp/tot*(1./s));
                        else
                                pdfTemp->SetBinContent(i3,i4,1/tot*(1./s));
                }
        }
        //Use better pdf from kde to evaluate the bindwidth 
        if(_doThecondPass==1){
                calculateBindWidth(validateEntries,fx,fy,varXSigma,varYSigma,pdfTemp,fhx,fhy);
                for(int i3=1;i3<=pdfTemp->GetXaxis()->GetNbins();i3++){
                        for(int i4=1;i4<=pdfTemp->GetYaxis()->GetNbins();i4++){
                                fx0[(i3-1)*_densityNBins+i4-1]=pdfTemp->GetXaxis()->GetBinCenter(i3) ;
                                fy0[(i3-1)*_densityNBins+i4-1]=pdfTemp->GetYaxis()->GetBinCenter(i4) ;
                        }
                }
                kde2d(validateEntries,fx,fy,fhx,fhy,varXmin,varXmax,varYmin,varYmax,_densityNBins*_densityNBins,fx0,fy0,fz0,_nThreads);
                tot=0.;
                for(int i3=1;i3<=pdfHist->GetXaxis()->GetNbins();i3++){
                        for(int i4=1;i4<=pdfHist->GetYaxis()->GetNbins();i4++){
                                tot+=fz0[(i3-1)*_densityNBins+i4-1];
                                pdfTemp->SetBinContent(i3,i4,fz0[(i3-1)*_densityNBins+i4-1]);
                        }
                }
                for(int i3=1;i3<=pdfTemp->GetXaxis()->GetNbins();i3++){
                        for(int i4=1;i4<=pdfTemp->GetYaxis()->GetNbins();i4++){
                                double vp=pdfTemp->GetBinContent(i3,i4)         ;
                                double wx=pdfTemp->GetXaxis()->GetBinWidth(i3)  ;
                                double wy=pdfTemp->GetYaxis()->GetBinWidth(i4)  ;
                                double s =wx*wy                                 ;
                                if(vp>0)
                                        pdfTemp->SetBinContent(i3,i4,vp/tot*(1./s));
                                else
                                        pdfTemp->SetBinContent(i3,i4,1/tot*(1./s));
                        }
                }
        }
        calculateBindWidth(validateEntries,fx,fy,varXSigma,varYSigma,pdfTemp,fhx,fhy);
        for(int i3=1;i3<=pdfHist->GetXaxis()->GetNbins();i3++){
                for(int i4=1;i4<=pdfHist->GetYaxis()->GetNbins();i4++){
                        fx0[(i3-1)*_cachePdfNBins+i4-1]=pdfHist->GetXaxis()->GetBinCenter(i3) ;
                        fy0[(i3-1)*_cachePdfNBins+i4-1]=pdfHist->GetYaxis()->GetBinCenter(i4) ;
                }
        }
        kde2d(validateEntries,fx,fy,fhx,fhy,varXmin,varXmax,varYmin,varYmax,_cachePdfNBins*_cachePdfNBins,fx0,fy0,fz0,_nThreads);
        tot=0.;
        double minBinContent=1.0e300;
        for(int i3=1;i3<=pdfHist->GetXaxis()->GetNbins();i3++){
                for(int i4=1;i4<=pdfHist->GetYaxis()->GetNbins();i4++){
                        if(minBinContent>fz0[(i3-1)*_cachePdfNBins+i4-1]) minBinContent=fz0[(i3-1)*_cachePdfNBins+i4-1];
                        tot+=fz0[(i3-1)*_cachePdfNBins+i4-1];
                        pdfHist->SetBinContent(i3,i4,fz0[(i3-1)*_cachePdfNBins+i4-1]);
                }
        }
        //cout<<"minBinContent= "<<minBinContent<<endl;
        minBinContent/=2.;
        //Performan Normalization
        for(int i3=1;i3<=pdfHist->GetXaxis()->GetNbins();i3++){
                for(int i4=1;i4<=pdfHist->GetYaxis()->GetNbins();i4++){
                        double vp=pdfHist->GetBinContent(i3,i4)         ;
                        double wx=pdfHist->GetXaxis()->GetBinWidth(i3)  ;
                        double wy=pdfHist->GetYaxis()->GetBinWidth(i4)  ;
                        double s =wx*wy                                 ;
                        if(vp>0)
                                pdfHist->SetBinContent(i3,i4,vp/tot*(1./s));
                        else
                                pdfHist->SetBinContent(i3,i4,minBinContent/tot*(1./s));
                }
        }
        delete pdfTemp;
}
//Automatically Build the Classification Expression
void kMvaClassification::startAutoBuilding(){
        initData()                      ;
        string newVarName               ;
        string workingVar1              ;
        string workingVar2              ;
        vector<string> origVarArray     ;
        VAR V               ;
        V.max=-1.e300       ;
        V.min=1.e300        ;
        V.mean=0.           ;
        V.sigma=0.          ;
        V.transformType=0   ;
        _outputFile->cd()          ;
        int entries=_sigTree->GetEntries();
        if(entries<_bkgTree->GetEntries())
                entries=_bkgTree->GetEntries();
        fx =new double[entries];
        fy =new double[entries];
        fhx=new double[entries];
        fhy=new double[entries];
        if(_densityNBins<_cachePdfNBins){
                fx0=new double[_cachePdfNBins*_cachePdfNBins];
                fy0=new double[_cachePdfNBins*_cachePdfNBins];
                fz0=new double[_cachePdfNBins*_cachePdfNBins];
        }
        else{
                fx0=new double[_densityNBins*_densityNBins];
                fy0=new double[_densityNBins*_densityNBins];
                fz0=new double[_densityNBins*_densityNBins];
        }
        vector<string> bestExpression;
        for(int i1=0;i1<_nVariable;i1++){
                origVarArray.push_back((_variableNames[i1]+"D"));
                bestExpression.push_back(_variableNames[i1]);
        }
        int mini1,mini2         ;
        double minLoss          ;
        double lossRet          ;
        //search First two best variables
        mini1=-1;mini2=-1;minLoss=1.0e300;
        for(int i1=0;i1<origVarArray.size();i1++){
                for(int i2=i1+1;i2<origVarArray.size();i2++){
                        workingVar1=origVarArray[i1]            ;
                        workingVar2=origVarArray[i2]            ;
                        newVarName=workingVar1+workingVar2      ;
                        V.transformType = -1       ;
                        varInfoSig[newVarName] = V;
                        varInfoBkg[newVarName] = V;
                        lossRet=buildNewVariable2(workingVar1.c_str(), workingVar2.c_str(), newVarName.c_str(), varInfoSig, varInfoBkg, _sigTree, _bkgTree, pdfSig, pdfBkg);
                        if(lossRet<minLoss){
                                minLoss=lossRet ;
                                mini1=i1        ;
                                mini2=i2        ;
                        }
                }
        }
        cout << "Ret: Min "<<origVarArray[mini1]<<" + "<<origVarArray[mini2]<<" with Loss= "<<minLoss<<endl;
        string curExpr=Form("(%s+%s)",bestExpression[mini1].c_str(),bestExpression[mini2].c_str());
        newVarName=origVarArray[mini1]+origVarArray[mini2];
        if(_plotRoc==1){
                roc[newVarName] = new TGraph();
                plotRocCurve(newVarName,_sigTree,_bkgTree,roc);
                roc[newVarName]->Write(Form("%s_Roc",newVarName.c_str()),TObject::kOverwrite);
        }
        pdfSig[newVarName]->Write("",TObject::kOverwrite);
        pdfBkg[newVarName]->Write("",TObject::kOverwrite);
        origVarArray.erase(origVarArray.begin()+mini1)  ;
        origVarArray.erase(origVarArray.begin()+mini2-1);
        bestExpression.erase(bestExpression.begin()+mini1);
        bestExpression.erase(bestExpression.begin()+mini2-1);
        bestExpression.push_back(curExpr)       ;
        workingVar1=newVarName;         

        //Iteratively Combine varibles
        while(origVarArray.size()>0){
                mini1=-1;minLoss=1.0e300;
                for(int i1=0;i1<origVarArray.size();i1++){
                        workingVar2=origVarArray[i1]            ;
                        newVarName=workingVar1+workingVar2      ;
                        V.transformType = -1       ;
                        varInfoSig[newVarName] = V;
                        varInfoBkg[newVarName] = V;
                        lossRet=buildNewVariable2(workingVar1.c_str(), workingVar2.c_str(), newVarName.c_str(), varInfoSig, varInfoBkg, _sigTree, _bkgTree, pdfSig, pdfBkg);
                        if(lossRet<minLoss){
                                minLoss=lossRet ;
                                mini1=i1        ;
                        }
                }
                cout << "Ret: Min "<<workingVar1<<" + "<<origVarArray[mini1]<<" with Loss= "<<minLoss<<endl;

                string curExpr=Form("(%s+%s)",bestExpression.back().c_str(),bestExpression[mini1].c_str());
                newVarName=workingVar1+origVarArray[mini1];
                if(_plotRoc==1){
                        roc[newVarName] = new TGraph();
                        plotRocCurve(newVarName,_sigTree,_bkgTree,roc);
                        roc[newVarName]->Write(Form("%s_Roc",newVarName.c_str()),TObject::kOverwrite);
                }
                pdfSig[newVarName]->Write("",TObject::kOverwrite);
                pdfBkg[newVarName]->Write("",TObject::kOverwrite);
                origVarArray.erase(origVarArray.begin()+mini1)  ;
                bestExpression.erase(bestExpression.begin()+mini1);
                bestExpression.push_back(curExpr)       ;
                workingVar1=newVarName                  ;
                cout<<curExpr<<endl;
        }
        _expression=bestExpression.back();
        writeVarInfo();
        cout<<"Training finished! Classifier has been written into "<<_outputFileName<<endl;
        cout<<"Best expression is "<<bestExpression.back()<<endl;
        _outputFile->Close();
        delete[] fx;
        delete[] fy;
        delete[]fhx;
        delete[]fhy;
        delete[]fx0;
        delete[]fy0;
        delete[]fz0;
}
/*Start to build the classifier*/ 
void kMvaClassification::start(){
        if(_auto==1){
                return startAutoBuilding();
        }
        initData()                                      ;
        string formula=_expressionP     ;
        int totLength=formula.length()  ;
        vector<string> stackNames       ;
        char varTempName[500]           ;
        char c                          ;
        int varTempNameLength = 0       ;
        int cur = 0                     ;
        string newVarName               ;
        string workingVar1              ;
        string workingVar2              ;
        VAR V               ;
        V.max=-1.e300       ;
        V.min=1.e300        ;
        V.mean=0.           ;
        V.sigma=0.          ;
        V.transformType=0   ;
        _outputFile->cd()          ;
        int entries=_sigTree->GetEntries();
        if(entries<_bkgTree->GetEntries())
                entries=_bkgTree->GetEntries();
        fx =new double[entries];
        fy =new double[entries];
        fhx=new double[entries];
        fhy=new double[entries]; 
        if(_densityNBins<_cachePdfNBins){
                fx0=new double[_cachePdfNBins*_cachePdfNBins];
                fy0=new double[_cachePdfNBins*_cachePdfNBins];
                fz0=new double[_cachePdfNBins*_cachePdfNBins];
        }
        else{
                fx0=new double[_densityNBins*_densityNBins];
                fy0=new double[_densityNBins*_densityNBins];
                fz0=new double[_densityNBins*_densityNBins];
        }
        //According to the classification expression to build two PDFs
        while(cur < totLength) {
                c = formula[cur] ;
                switch (c) {
                        case '|':
                                varTempName[varTempNameLength] = 0                  ;
                                stackNames.push_back(string(varTempName)+"D")           ;
                                varTempNameLength = 0                               ;
                                break;
                        case '+':
                                //Start to form new variable
                                newVarName = "";
                                {
                                workingVar1=stackNames.back()   ;
                                stackNames.pop_back()           ;
                                workingVar2=stackNames.back()   ;
                                stackNames.pop_back()           ;
                                newVarName = workingVar1+workingVar2;
                                V.transformType = -1       ;
                                varInfoSig[newVarName] = V;
                                varInfoBkg[newVarName] = V;
                                buildNewVariable2(workingVar1.c_str(), workingVar2.c_str(), newVarName.c_str(), varInfoSig, varInfoBkg, _sigTree, _bkgTree, pdfSig, pdfBkg);
                                if(_plotRoc==1){
                                        roc[newVarName] = new TGraph();
                                        plotRocCurve(newVarName,_sigTree,_bkgTree,roc);
                                        roc[newVarName]->Write(Form("%s_Roc",newVarName.c_str()),TObject::kOverwrite);
                                }
                                pdfSig[newVarName]->Write("", TObject::kOverwrite);
                                pdfBkg[newVarName]->Write("", TObject::kOverwrite);
                                }
                                stackNames.push_back(newVarName);
                                varTempNameLength = 0;
                                //read one more char, to skip |
                                cur++;
                                break;
                        default:
                                varTempName[varTempNameLength] = c;
                                varTempNameLength++;
                                break;
                }
            cur++;
        }
        writeVarInfo();
        cout<<"Training finished! Classifier has been written into "<<_outputFileName<<endl;
        _outputFile->Close();
        delete[] fx;
        delete[] fy;
        delete[]fhx;
        delete[]fhy;
        delete[]fx0;
        delete[]fy0;
        delete[]fz0;
}

//Create a new variable by combinine two varibales
//Input:
//      _var1/_var2             : input variables' name
//      _newvar                 : new variable name
//      _varInfoSig             : variables information array of Signal samples
//      _varInfoBkg             : variables information array of Background samples
//      _treeSig/Bkg    : Data points' TTree
//Output:
//      _pdfSig                 : Signal pdf array
//      _pdfBkg                 : Background pdf array
//Return:
//      likehood value of the new variable
double kMvaClassification::buildNewVariable2(const char* _var1,const char* _var2,const char* _newvar,map<string,VAR>&_varInfoSig,map<string,VAR>&_varInfoBkg,TTree* _treeSig,TTree* _treeBkg,map<string, TH2D*>& _pdfSig,map<string, TH2D*>& _pdfBkg){
    char newVarName[10000]                   ;
    double retLoss=0.                      ;
    sprintf(newVarName,"%s",_newvar)       ;
    double hwx=(_varInfoSig[_var1].max-_varInfoSig[_var1].min)/_cachePdfNBins;
    double hwy=(_varInfoSig[_var2].max-_varInfoSig[_var2].min)/_cachePdfNBins;

        _pdfSig[newVarName]=new TH2D(Form("%sS",newVarName),Form("Sig:%s", newVarName),_cachePdfNBins,_varInfoSig[_var1].min-0.5*hwx,_varInfoSig[_var1].max+0.5*hwx,_cachePdfNBins,_varInfoSig[_var2].min-0.5*hwy,_varInfoSig[_var2].max+0.5*hwy);
    hwx=(_varInfoBkg[_var1].max-_varInfoBkg[_var1].min)/_cachePdfNBins;
    hwy=(_varInfoBkg[_var2].max-_varInfoBkg[_var2].min)/_cachePdfNBins;
    
        _pdfBkg[newVarName]=new TH2D(Form("%sB",newVarName),Form("Bkg:%s", newVarName),_cachePdfNBins,_varInfoBkg[_var1].min-0.5*hwx,_varInfoBkg[_var1].max+0.5*hwx,_cachePdfNBins,_varInfoBkg[_var2].min-0.5*hwy,_varInfoBkg[_var2].max+0.5*hwy);
    buildPdf2Parallel(_var1,_var2,newVarName,_varInfoSig,_treeSig,_pdfSig);
    buildPdf2Parallel(_var1,_var2,newVarName,_varInfoBkg,_treeBkg,_pdfBkg);
    float var1Val,var2Val,newVarVal     ;
    double pSig,pBkg,pLik               ;
    TBranch* var1BranchSig      ;
    TBranch* var1BranchBkg      ;
    TBranch* newVarBranchSig    ;
    TBranch* var2BranchSig      ;
    TBranch* var2BranchBkg      ;
    TBranch* newVarBranchBkg    ;
    _treeSig->SetBranchAddress(_var1,&var1Val,&var1BranchSig);
    _treeSig->SetBranchAddress(_var2,&var2Val,&var2BranchSig);
    newVarBranchSig=_treeSig->Branch(newVarName,&newVarVal,Form("%s/F",newVarName));
    _treeBkg->SetBranchAddress(_var1,&var1Val,&var1BranchBkg);
    _treeBkg->SetBranchAddress(_var2,&var2Val,&var2BranchBkg);
    newVarBranchBkg=_treeBkg->Branch(newVarName,&newVarVal,Form("%s/F",newVarName));
    
    int v1TransfromType=_varInfoSig[_var1].transformType;
    int v2TransfromType=_varInfoSig[_var2].transformType;
    int entriesSig=_treeSig->GetEntries()       ;
    double validateEntriesSig=0.                        ;
    double v1,v2                        ;
    double sigMinX=_pdfSig[newVarName]->GetXaxis()->GetBinLowEdge(1);
    double sigMinY=_pdfSig[newVarName]->GetYaxis()->GetBinLowEdge(1);
    double sigMaxX=_pdfSig[newVarName]->GetXaxis()->GetBinLowEdge(_cachePdfNBins);
    double sigMaxY=_pdfSig[newVarName]->GetYaxis()->GetBinLowEdge(_cachePdfNBins);
    double bkgMinX=_pdfBkg[newVarName]->GetXaxis()->GetBinLowEdge(1);
    double bkgMinY=_pdfBkg[newVarName]->GetYaxis()->GetBinLowEdge(1);
    double bkgMaxX=_pdfBkg[newVarName]->GetXaxis()->GetBinUpEdge(_cachePdfNBins);
    double bkgMaxY=_pdfBkg[newVarName]->GetYaxis()->GetBinUpEdge(_cachePdfNBins);
        //Loop all Signal data points and calculate the new varible
    for(int i1=0;i1<entriesSig;i1++){
        var1BranchSig->GetEntry(i1)     ;
        var2BranchSig->GetEntry(i1)     ;
        v1=var1Val;
        v2=var2Val;
                if(v1<_varInfoSig[_var1].min||v2<_varInfoSig[_var2].min||v1>_varInfoSig[_var1].max||v2>_varInfoSig[_var2].max) 
                        pSig=0.;
                else
                        pSig=_pdfSig[newVarName]->Interpolate(v1,v2);
                if ((_varInfoSig[_var1].transformType==-1&&v1==1)||(_varInfoSig[_var2].transformType==-1&&v2==1)) {
                        pSig=1.;
                        pBkg=0.;
                }

                if(v1<_varInfoBkg[_var1].min||v2<_varInfoBkg[_var2].min||v1>_varInfoBkg[_var1].max||v2>_varInfoBkg[_var2].max) 
                        pBkg=0.;
                else
                        pBkg=_pdfBkg[newVarName]->Interpolate(v1,v2); 
                if ((_varInfoBkg[_var1].transformType==-1&&v1==0)||(_varInfoBkg[_var2].transformType==-1&&v2==0)) {
                        pBkg=1.;
                        pSig=0.;
                }
                if(pSig+pBkg>0.)
                        pLik=pSig/(pSig+pBkg)      ;
                else
                        pLik=0.5;
                if (pLik>0.)
                        retLoss+=-1./double(entriesSig)*log(pLik);
                else
                        retLoss+=-1./double(entriesSig)*log(1.0e-300);
                newVarVal=0.5*(1.+tanh(0.3*atanh((2*pLik-1.))));
                if(newVarVal!=1.&&newVarVal!=0.){
                        _varInfoSig[newVarName].mean =(_varInfoSig[newVarName].mean*validateEntriesSig)/(1.0+validateEntriesSig)+newVarVal/(1.0+validateEntriesSig);
                        _varInfoSig[newVarName].sigma=(_varInfoSig[newVarName].sigma*validateEntriesSig)/(1.0+validateEntriesSig)+newVarVal*newVarVal/(1.0+validateEntriesSig);
                        validateEntriesSig+=1.;
                }
                newVarBranchSig->Fill()                 ;
    }
    _varInfoSig[newVarName].min=0.;
    _varInfoSig[newVarName].max=1.;
    _varInfoSig[newVarName].sigma-=(_varInfoSig[newVarName].mean*_varInfoSig[newVarName].mean);
    _varInfoSig[newVarName].sigma=sqrt(_varInfoSig[newVarName].sigma); 
    int entriesBkg=_treeBkg->GetEntries()       ;
    double validateEntriesBkg=0.                        ;
        //Loop all Background data points and calculate the new varible
    for(int i1=0;i1<entriesBkg;i1++){
        var1BranchBkg->GetEntry(i1)     ;
        var2BranchBkg->GetEntry(i1)     ;
        v1=var1Val;
        v2=var2Val;
                if(v1<_varInfoSig[_var1].min||v2<_varInfoSig[_var2].min||v1>_varInfoSig[_var1].max||v2>_varInfoSig[_var2].max) 
                        pSig=0.                                         ;
                else
                        pSig=_pdfSig[newVarName]->Interpolate(v1,v2);
                if ((_varInfoSig[_var1].transformType==-1&&v1==1)||(_varInfoSig[_var2].transformType==-1&&v2==1)) {
                        pSig=1.         ;
                        pBkg=0.         ;
                }

                if(v1<_varInfoBkg[_var1].min||v2<_varInfoBkg[_var2].min||v1>_varInfoBkg[_var1].max||v2>_varInfoBkg[_var2].max) pBkg=0.;
                else
                        pBkg=_pdfBkg[newVarName]->Interpolate(v1,v2);
                if ((_varInfoBkg[_var1].transformType==-1&&v1==0)||(_varInfoBkg[_var2].transformType==-1&&v2==0)) {
                        pBkg=1.         ;
                        pSig=0.         ;
                }

                if(pSig+pBkg>0.)
                        pLik=pSig/(pSig+pBkg)      ;
                else
                        pLik=0.5                ;
                if(pLik<1.)
                        retLoss+=-1./double(entriesBkg)*log(1.-pLik);
                else
                        retLoss+=-1./double(entriesBkg)*log(1.0e-300);
                newVarVal=0.5*(1.+tanh(0.3*atanh((2*pLik-1.))));
                if(newVarVal!=0.&&newVarVal!=1.){
                        _varInfoBkg[newVarName].mean =(_varInfoBkg[newVarName].mean*validateEntriesBkg)/(1.0+validateEntriesBkg)+newVarVal/(1.0+validateEntriesBkg)                          ;
                        _varInfoBkg[newVarName].sigma=(_varInfoBkg[newVarName].sigma*validateEntriesBkg)/(1.0+validateEntriesBkg)+newVarVal*newVarVal/(1.0+validateEntriesBkg)       ;
                        validateEntriesBkg+=1.          ;
        }
                newVarBranchBkg->Fill()         ;
    }
    _varInfoBkg[newVarName].min=0.;
    _varInfoBkg[newVarName].max=1.;
    _varInfoBkg[newVarName].sigma-=(_varInfoBkg[newVarName].mean*_varInfoBkg[newVarName].mean);
    _varInfoBkg[newVarName].sigma=sqrt(_varInfoBkg[newVarName].sigma);
    
    return retLoss;
}
//Generate ROC curve
//Input:
//      _varName        : Variable Name
//      _treeSig        : Signal TTree Name
//      _treeBkg        : Background TTree Name
//Output:
// _roc                 : ROC curve array
void kMvaClassification::plotRocCurve(string _varName,TTree* _treeSig,TTree* _treeBkg,map<string,TGraph*>& _roc){
    TBranch* brSig;
    TBranch* brBkg;
        float val;
        int entriesSig,entriesBkg;
        double sigEff[101];
    double bkgEff[101];
    _treeSig->SetBranchAddress(_varName.c_str(),&val,&brSig);
    _treeBkg->SetBranchAddress(_varName.c_str(),&val,&brBkg);

    entriesSig=brSig->GetEntries()              ;
    entriesBkg=brBkg->GetEntries()              ;    
    vector<float> likSigV(entriesSig)   ;
    vector<float> likBkgV(entriesBkg)   ;
    //load data
    for(int i1=0;i1<entriesSig;i1++){
                brSig->GetEntry(i1);
                likSigV[i1]=val;
    }
    for(int i1=0;i1<entriesBkg;i1++){
                brBkg->GetEntry(i1);
                likBkgV[i1]=val;
        }
        //Sort the data points
    std::sort(likSigV.begin(),likSigV.end())            ;
    std::sort(likBkgV.begin(),likBkgV.end())            ;
    int k=0             ;
        //roc points
    for(int i1=0;i1<101;i1++){
        int n=int((0.01*i1)*entriesSig);
        if(n==entriesSig) n=entriesSig-1;
        double v=likSigV[n];
        sigEff[i1]=1.0-double(n)/double(entriesSig);
        for(;k<entriesBkg;k++){
                        if(likBkgV[k]>v) break;
        }
        bkgEff[i1]=1.0-double(k)/double(entriesBkg);
    }
    for (int i1 = 0; i1 < 101; i1++)
        _roc[_varName]->SetPoint(i1,sigEff[i1],1.0-bkgEff[i1])    ;
}
//Write classifier data base file
void kMvaClassification::writeVarInfo(){
    //output variable's information
    _outputFile->cd();
    TTree* varInfoSigTree=new TTree("varInfoSigTree","varInfoSigTree");
    TTree* varInfoBkgTree=new TTree("varInfoBkgTree","varInfoBkgTree");
    TTree* variableNamesTree=new TTree("variableNames","variableNames");
    double _min         ;
    double _max         ;
    double _mean        ;
    double _sigma       ;
    int _transformType  ;
    char varName[100]   ;
    varInfoSigTree->Branch("varName",(void*)varName,"varName/C")   ;
    varInfoSigTree->Branch("min",&_min,"min/D")                    ;
    varInfoSigTree->Branch("max",&_max,"max/D")                    ;
    varInfoSigTree->Branch("mean",&_mean,"mean/D")                 ;
    varInfoSigTree->Branch("sigma",&_sigma,"sigma/D")              ;
    varInfoSigTree->Branch("transformType",&_transformType,"transformType/I");
    
    varInfoBkgTree->Branch("varName",(void*)varName,"varName/C")   ;
    varInfoBkgTree->Branch("min",&_min,"min/D")                    ;
    varInfoBkgTree->Branch("max",&_max,"max/D")                    ;
    varInfoBkgTree->Branch("mean",&_mean,"mean/D")                 ;
    varInfoBkgTree->Branch("sigma",&_sigma,"sigma/D")              ;
    varInfoBkgTree->Branch("transformType",&_transformType,"transformType/I");

    variableNamesTree->Branch("varName",(void*)varName,"varName/C")   ;
    
    for(map<string,VAR>::iterator iter=varInfoSig.begin();iter!=varInfoSig.end();iter++){
        sprintf(varName,"%s",iter->first.c_str());
        _min=iter->second.min   ;
        _max=iter->second.max   ;
        _mean=iter->second.mean ;
        _sigma=iter->second.sigma;
        _transformType=iter->second.transformType;
        varInfoSigTree->Fill();
    }
    
    for(map<string,VAR>::iterator iter=varInfoBkg.begin();iter!=varInfoBkg.end();iter++){
        sprintf(varName,"%s",iter->first.c_str());
        _min=iter->second.min;
        _max=iter->second.max;
        _mean=iter->second.mean;
        _sigma=iter->second.sigma;
        _transformType=iter->second.transformType;
        varInfoBkgTree->Fill();
    }
        for(int i1=0;i1<_variableNames.size();i1++){
                sprintf(varName,"%s",_variableNames[i1].c_str())        ;
                variableNamesTree->Fill()                       ;
        }
    varInfoSigTree->Write("",TObject::kOverwrite)       ;
    varInfoBkgTree->Write("",TObject::kOverwrite)       ;
    deSigR->Write("deSigR",TObject::kOverwrite)         ;
    mSig->Write("mSig",TObject::kOverwrite)             ;
    _sigTree->Write("",TObject::kOverwrite)             ;
    _bkgTree->Write("",TObject::kOverwrite)             ;
    variableNamesTree->Write("",TObject::kOverwrite)    ;
    TNamed expression("expresssion",_expression.c_str());
    expression.Write("",TObject::kOverwrite)    ;
}

