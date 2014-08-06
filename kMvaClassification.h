/* 
 * File:   kMvaClassification.h
 * Author: kaiwu
 * Email : kaiwu@cern.ch
 * Created on May 04, 2013, 4:37 PM
 */
#ifndef CLASSIFICATION_H
#define CLASSIFICATION_H
#include "TTree.h"
#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include "TNamed.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMatrixDSymEigen.h"
#include "iBalltreeDensity.h"

using namespace std;
class kMvaClassification{
public:
        /*Contruction function
          input variable: configFile ---> configuration file name
        */
        kMvaClassification(string configFile);
        ~kMvaClassification()                ;
        /*Start to build the classifier*/ 
        void start()                                     ;
private:
        string _configFileName;
        string _inputSigFileName ;
        string _inputBkgFileName ;
        TFile* _inputSigFile;
        TFile* _inputBkgFile;
        string _sigTreeName;
        string _bkgTreeName;
        TTree* _sigTree  ;
        TTree* _bkgTree  ;
        int    _nVariable;
        vector<string> _variableNames;
        int    _nSeperator;
        vector<string> _seperatorNames;
        int     _transFormType;
        int     _densityNBins;
        int     _cachePdfNBins;
        int     _nThreads       ;
        int _doThecondPass      ;
        int _plotRoc            ;
        string _dataBase        ;
        int    _nTerms          ;
        double _epsilon         ;
        int _auto               ;
        string  _expression             ;
        string  _expressionP            ;
        string  _outputFileName         ;
        TFile*  _outputFile             ;
        iBalltreeDensity* _iBT           ;
                
        TMatrixDSym* covSig       ;
        TVectorD* mSig         ;
        TMatrixD* deSigR       ;
        //Varibale information structure
        typedef struct{
                double min              ;
                double max              ;
                double mean             ;
                double sigma            ;
                int transformType       ; 
        }VAR;
        //Variables' information
        map<string,VAR> varInfoSig  ;
        map<string,VAR> varInfoBkg  ;
        //Receiver of response curve
        map<string, TGraph*> roc    ;
        //PDF array
        map<string, TH2D*>  pdfSig  ;
        map<string, TH2D*>  pdfBkg  ;

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
        float kde2d(int nPoints,double*x,double* y,double* hx,double* hy,double xmin,double xmax,double ymin,double ymax,int mEval,const double* x0,const double* y0,double* z0,int nthreads=1);

        //Build two dimension pdf for variable var1 and var2
        //Input:
        //      _var1/2/_newvar : variables' name
        //      _varInfo                : Information of the input variables
        //      _treeIn                 : Data Points TTree
        //Output:
        //      _pdfs                   : PDF array
        void buildPdf2Parallel(const char* _var1, const char* _var2, const char* _newvar,map<string,VAR>& _varInfo, TTree* _treeIn, map<string,TH2D*>& _pdfs);
        //Calculate the bandwidth of kernel density estimator
        //Input:
        //      _nPoints        : Number of data points
        //      _x/_y           : x/y coordinates of the data points
        //      _sigmaX/Y       : variance of x/y
        //      _histPDF        : raw histogram PDF as an approxmination of PDF of the data points
        //Output:
        // _hx/_hy              : Bandwidth of data points.
        void calculateBindWidth(int _nPoints,double* _x,double* _y,double _sigmaX,double _sigmaY,TH2D* _histPdf,double* _hx,double* _hy);
        //Generate ROC curve
        //Input:
        //      _varName        : Variable Name
        //      _treeSig        : Signal TTree Name
        //      _treeBkg        : Background TTree Name
        //Output:
        // _roc                 : ROC curve array
        void plotRocCurve(string _varName,TTree* _treeSig,TTree* _treeBkg,map<string,TGraph*>& _roc)   ;
        //Read Configuration file
        //Return:       true ---> Configuration file is Ok
        //                      false---> Configuration file has some error
        bool readConfig()     ;
        //Print out configuration
        void printConfig()    ;
        //Initialize Data Structure
        void initData()       ;
        //Write classifier data base file
        void writeVarInfo()   ;
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
        double buildNewVariable2(const char* _var1,const char* _var2,const char* _newvar,map<string,VAR>&_varInfoSig,map<string,VAR>&_varInfoBkg,TTree* _treeSig,TTree* _treeBkg,map<string, TH2D*>& _pdfSig,map<string, TH2D*>& _pdfBkg);
        void startAutoBuilding()             ;
        //arrays for calculating kernel density estimation
        double* fx          ;
        double* fy          ;
        double* fhx         ;
        double* fhy         ;
        double* fx0         ;
        double* fy0         ;
        double* fz0         ;
};
#endif
