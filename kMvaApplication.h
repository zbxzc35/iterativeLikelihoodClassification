/* 
 * File:   kMvaApplication.h
 * Author: kaiwu
 *
 * Created on May 4, 2013, 4:46 PM
 */

#ifndef KMVAAPPLICATION_H
#define KMVAAPPLICATION_H
#include "TTree.h"
#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include "TGraph.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMatrixDSymEigen.h"
#include "map"

using namespace std;

class kMvaApplication {
public:
        //Construct function
        //Input:
        //      _databaseFileName: Database root file name
        kMvaApplication(string _databaseFileName);
        virtual ~kMvaApplication()              ;
        //Add variables' names used in the classification expression
        //Input:
        //      _var    : Variable name
        //Return        : Current number of input varibles
        int addVariable(string _var)    ;
        //Evaluate iterative likelihood, and store intermidiate likelihood value for tracing back purpose
        //Input:
        //      _varVals:       input variables' value array
        //Output:
        //  _intermediateVariables: intermediate likelihood values' array
        //Return:       Likelihood value(performaned tanh transform to strech the 0&1 regions)
        double eval(vector<float>&_varVals,vector<float>& _intermediateVariables);
        //Evaluate iterative likelihood
        //Input:
        //      _varVals:       input variables' value array
        //Return:       Likelihood value(performaned tanh transform to strech the 0&1 regions)
        double eval(vector<float>&_varVals);
private:
        vector<string> variables ;
        map<string,float> varVals;
        map<string,int>   varInputOrder;
        map<string,bool> varValidate;
        int     nVariables      ;
        string expression       ;
        string expressionP      ;
        string databaseFileName ;
        TFile* databaseFile     ;
        TVectorD* mSig          ;
        TMatrixD* deSigR        ;
        typedef struct{
                double min              ;
                double max              ;
                double mean             ;
                double sigma            ;
                int transformType       ;
        }VAR;
        //variable's information
        map<string,VAR> varInfoSig  ;
        map<string,VAR> varInfoBkg  ;
        //pdf array
        map<string, TH2D*>  pdfSig  ;
        map<string, TH2D*>  pdfBkg  ;
        bool init()                             ;
        bool hasInit                ;
};

#endif  /* KMVAAPPLICATION_H */


