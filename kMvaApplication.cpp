/* 
 * File:   kMvaApplication.cxx
 * Author: kaiwu
 * 
 * Created on May 4, 2013, 4:46 PM
 */

#include <TObjArray.h>
#include <algorithm>
#include <TNamed.h>
#include "kMvaApplication.h"
//External procedure: Parse the classification expression 
//Input:
// first        :       Classification Expression
// second       :       Reverse Polish Notation
//Return        :   0---> Ok
//                              1---> False
extern "C" int parseExpression(const char* , char* );
/*Contruction function
 input variable: configFile ---> configuration file name
 */
kMvaApplication::kMvaApplication(string _databaseFileName){
        databaseFileName=_databaseFileName      ;
        databaseFile=TFile::Open(databaseFileName.c_str());
        if(!databaseFile){
                cout<<"Can not open "<<databaseFileName<<"!"<<endl;
                exit(-1);
        }
        TNamed* _expression=(TNamed*) databaseFile->Get("expresssion");
        expression      =_expression->GetTitle()        ;
        char temp1[1000];
        char temp2[1000];
        sprintf(temp1,"%s",expression.c_str())  ;
        if(parseExpression(temp1,temp2)!=0){
                cout<<"Classification expression syntax error!"<<endl;
                exit(-1);
        }
        expressionP=temp2;
        hasInit=false;
}
//Initialize...
bool kMvaApplication::init(){
        bool ret=true                   ;
        nVariables=variables.size()     ;
        deSigR=(TMatrixD*) databaseFile->Get("deSigR")  ;
        mSig  =(TVectorD*) databaseFile->Get("mSig")    ;
        //Load variable information
        TTree* varInfoSigTree=(TTree*) databaseFile->Get("varInfoSigTree");
        TTree* varInfoBkgTree=(TTree*) databaseFile->Get("varInfoBkgTree");
        TTree* variableNamesTree=(TTree*) databaseFile->Get("variableNames");
        int _nVariables;
        VAR V;
        char varName[100]  ;
        varInfoSigTree->SetBranchAddress("varName",varName)   ;
        varInfoSigTree->SetBranchAddress("min",&(V.min))      ;
        varInfoSigTree->SetBranchAddress("max",&(V.max))      ;
        varInfoSigTree->SetBranchAddress("mean",&(V.mean))    ;
        varInfoSigTree->SetBranchAddress("sigma",&(V.sigma))  ;
        varInfoSigTree->SetBranchAddress("transformType",&(V.transformType))  ;


        varInfoBkgTree->SetBranchAddress("varName",varName)   ;
        varInfoBkgTree->SetBranchAddress("min",&(V.min))      ;
        varInfoBkgTree->SetBranchAddress("max",&(V.max))      ;
        varInfoBkgTree->SetBranchAddress("mean",&(V.mean))    ;
        varInfoBkgTree->SetBranchAddress("sigma",&(V.sigma))  ;
        varInfoBkgTree->SetBranchAddress("transformType",&(V.transformType))  ;

        variableNamesTree->SetBranchAddress("varName",varName);
        _nVariables=variableNamesTree->GetEntries();
        if(_nVariables!=nVariables){
                cout<<"The classification method use "<<_nVariables<<" variables, while user input "<<nVariables<<" variables!"<<endl;
                cout<<"Please Check!"<<endl;
                exit(0);
        }
        for(int i1=0;i1<_nVariables;i1++){
                variableNamesTree->GetEntry(i1);
                varInputOrder[varName]=i1;
        }
        //clear map array
        varInfoSig.clear()  ;
        varInfoBkg.clear()  ;
        int entries         ;
        entries=varInfoSigTree->GetEntries();
        for(int i1=0;i1<entries;i1++){
            varInfoSigTree->GetEntry(i1)    ;
            varInfoSig[varName]=V           ;
        }
        entries=varInfoBkgTree->GetEntries();
        for(int i1=0;i1<entries;i1++){
            varInfoBkgTree->GetEntry(i1);
            varInfoBkg[varName]=V;
        }
        {
        //according to the expression, load required pdfs
        int totLength=expressionP.length()  ;
        vector<string> stackNames       ;
        char varTempName[500]           ;
        char c                          ;
        int varTempNameLength = 0       ;
        int cur = 0                     ;
        string newVarName               ;
        string newVarName2              ;
        string workingVar1              ;
        string workingVar2              ;
        while (cur < totLength) {
                        c= expressionP[cur] ;
                        switch (c) {
                                case '|':
                                        varTempName[varTempNameLength] = 0              ;
                                        stackNames.push_back(string(varTempName)+"D")   ;
                                        varTempNameLength = 0;
                                        break;
                                case '+':
                                        //Get new variable's name and load corresponding PDF for it
                                        newVarName = "";
                                        newVarName2= "";
                                        workingVar1=stackNames.back()   ;
                                        stackNames.pop_back()           ;
                                        workingVar2=stackNames.back()   ;
                                        stackNames.pop_back()           ;
                                        newVarName = workingVar1+workingVar2;
                                        pdfSig[newVarName]=(TH2D*) databaseFile->Get((newVarName+"S").c_str());
                                        pdfBkg[newVarName]=(TH2D*) databaseFile->Get((newVarName+"B").c_str());
                                        if(!pdfSig[newVarName]||!pdfBkg[newVarName]){
                                                cout<<"Can not load pdf for variable "<<newVarName<<" from "<<databaseFileName<<endl;
                                                ret=false;
                                                return ret;
                                        }
                                        stackNames.push_back(newVarName)    ;
                                        varTempNameLength = 0               ;
                                        cur++;
                                        break;
                                default:
                                        varTempName[varTempNameLength] = c;
                                        varTempNameLength++;
                                        break;
                        }
                        cur++;
                }
        }
        if(ret)
                hasInit=true;
        return ret;
}

kMvaApplication::~kMvaApplication(){
        if(databaseFile){
                if(databaseFile->IsOpen())      databaseFile->Close();
        }
}
//Add variables' names used in the classification expression
//Input:
//      _var    : Variable name
//Return        : Current number of input varibles
int kMvaApplication::addVariable(string _var){
        variables.push_back(_var);
        varVals[_var]=0.         ;
        return variables.size()  ;
}
//Evaluate iterative likelihood, and store intermidiate likelihood value for tracing back purpose
//Input:
//      _varVals:       input variables' value array
//Output:
//  _intermediateVariables: intermediate likelihood values' array
//Return:       Likelihood value(performaned tanh transform to strech the 0&1 regions)
double kMvaApplication::eval(vector<float>&_varVals,vector<float>& _intermediateVariables){
        double ret;
        if(_varVals.size()!=variables.size()){
                cout<<"Number of variables' values should be equal to the number you variables you have added!"<<endl;
                return -1;
        }
        if(!hasInit){
                if(!init()){
                        cout<<"Can not initialize! Please check!"<<endl;
                        exit(0);
                }
        }
        //make sure the order of the variables is the same of the order of the training
        vector<int> varOrder(_varVals.size());
        for(int i1=0;i1<_varVals.size();i1++)   varOrder[i1]=varInputOrder[variables[i1]];
        for(int i1=0;i1<_varVals.size();i1++){
                varVals[variables[i1]]=0.       ;
                for(int i2=0;i2<_varVals.size();i2++)
                        varVals[variables[i1]]+=((*deSigR)(varOrder[i1],varOrder[i2])*(_varVals[i2]-(*mSig)(varOrder[i2])));
        }
        //Start to calcuate the likelihood ratio value
        {

                int totLength = expressionP.length();
                vector<string> stackNames       ;
                vector<double> stackVals        ;
                _intermediateVariables.clear()  ;
                char varTempName[500]   ;
                char c                  ;
                int varTempNameLength = 0       ;
                int cur = 0                     ;
                string newVarName               ;
                string workingVar1              ;
                string workingVar2              ;
                double newVarVal                ;
                double workingVar1Val   ;
                double workingVar2Val   ;
                double pSig, pBkg,pLik  ;
                int var1Validate                ;
                int var2Validate                ;
                while (cur < totLength) {
                        c = expressionP[cur];
                        switch (c) {
                                case '|':
                                        varTempName[varTempNameLength] = 0              ;
                                        stackVals.push_back(varVals[varTempName])       ;
                                        stackNames.push_back(string(varTempName)+"D")   ;
                                        varTempNameLength = 0                           ;
                                        break   ;
                                case '+':
                                        //Get new variable's name and calculate the value of it
                                        newVarName = ""                         ;
                                        workingVar1Val = stackVals.back()       ;
                                        stackVals.pop_back()                    ;
                                        workingVar2Val = stackVals.back()       ;
                                        stackVals.pop_back()                    ;
                                        workingVar1 = stackNames.back()         ;
                                        stackNames.pop_back()                   ;
                                        workingVar2 = stackNames.back()         ;
                                        stackNames.pop_back()                   ;
                                        newVarName = workingVar1 + workingVar2  ;
                                        float v1, v2;
                                        v1=workingVar1Val;v2=workingVar2Val     ;
                                        if(v1<varInfoSig[workingVar1].min||v2<varInfoSig[workingVar2].min||v1>varInfoSig[workingVar1].max||v2>varInfoSig[workingVar2].max) pSig=0.;
                                        else pSig=pdfSig[newVarName]->Interpolate(v1,v2);
                                        if((varInfoSig[workingVar1].transformType==-1&&v1==1)||(varInfoSig[workingVar2].transformType==-1&&v2==1)){
                                                pSig=1.;
                                                pBkg=0.;
                                        }
                                        if(v1<varInfoBkg[workingVar1].min||v2<varInfoBkg[workingVar2].min||v1>varInfoBkg[workingVar1].max||v2>varInfoBkg[workingVar2].max) pBkg=0.;
                                        else pBkg=pdfBkg[newVarName]->Interpolate(v1,v2);
                                        if((varInfoBkg[workingVar1].transformType==-1&&v1==0.)||(varInfoBkg[workingVar2].transformType==-1&&v2==0.)){
                                                pBkg=1.;
                                                pSig=0.;
                                        }
                                        if(pSig+pBkg>0.)
                                                pLik=pSig/(pSig+pBkg)      ;
                                        else
                                                pLik=0.5                   ;
                                        newVarVal=0.5*(1.+tanh(0.3*atanh((2*pLik-1.))));
                                        stackNames.push_back(newVarName)        ;
                                        stackVals.push_back(newVarVal)          ;
                                        _intermediateVariables.push_back(newVarVal);
                                        varTempNameLength = 0                   ;
                                        cur++   ;
                                        break   ;   
                                default:
                                        varTempName[varTempNameLength] = c;
                                        varTempNameLength++;
                                        break;  
                        }
                        cur++;
            }
            ret=stackVals.back();
            stackVals.pop_back();
        }
        return ret;
}
//Evaluate iterative likelihood
//Input:
//      _varVals:       input variables' value array
//Return:       Likelihood value(performaned tanh transform to strech the 0&1 regions)
double kMvaApplication::eval(vector<float>&_varVals){
        double ret;
        if(_varVals.size()!=variables.size()){
                cout<<"Number of variables' values should be equal to the number you variables you have added!"<<endl;
                return -1;
        }
        if(!hasInit){
                if(!init()){
                        cout<<"Can not initialize! Please check!"<<endl;
                        exit(0);
                }
        }
        //make sure the order of the variables is the same of the order of the training
        vector<int> varOrder(_varVals.size());
        for(int i1=0;i1<_varVals.size();i1++)   varOrder[i1]=varInputOrder[variables[i1]];
        for(int i1=0;i1<_varVals.size();i1++){
                varVals[variables[i1]]=0.       ;
                for(int i2=0;i2<_varVals.size();i2++)
                        varVals[variables[i1]]+=((*deSigR)(varOrder[i1],varOrder[i2])*(_varVals[i2]-(*mSig)(varOrder[i2])));
        }
        //Start to calcuate the likelihood ratio value
        {

                int totLength = expressionP.length();
                vector<string> stackNames       ;
                vector<double> stackVals        ;
            char varTempName[500]       ;
                char c                  ;
                int varTempNameLength = 0       ;
                int cur = 0                     ;
                string newVarName               ;
                string workingVar1              ;
                string workingVar2              ;
                double newVarVal                ;
                double workingVar1Val   ;
                double workingVar2Val   ;
                double pSig, pBkg,pLik  ;
            int var1Validate            ;
            int var2Validate            ;
                while (cur < totLength) {
                        c = expressionP[cur];
                        switch (c) {
                                case '|':
                                        varTempName[varTempNameLength] = 0              ;
                                        stackVals.push_back(varVals[varTempName])       ;
                                        stackNames.push_back(string(varTempName)+"D")   ;
                                        varTempNameLength = 0                           ;
                                        break   ;
                                case '+':
                                        //Get new variable's name and calculate the value of it
                                        newVarName = ""                         ;
                                        workingVar1Val = stackVals.back()       ;
                                        stackVals.pop_back()                    ;
                                        workingVar2Val = stackVals.back()       ;
                                        stackVals.pop_back()                    ;
                                        workingVar1 = stackNames.back()         ;
                                        stackNames.pop_back()                   ;
                                        workingVar2 = stackNames.back()         ;
                                        stackNames.pop_back()                   ;
                                        newVarName = workingVar1 + workingVar2  ;
                                        float v1, v2;
                                        v1=workingVar1Val;v2=workingVar2Val     ;
                                        if(v1<varInfoSig[workingVar1].min||v2<varInfoSig[workingVar2].min||v1>varInfoSig[workingVar1].max||v2>varInfoSig[workingVar2].max) pSig=0.;
                                        else pSig=pdfSig[newVarName]->Interpolate(v1,v2);
                                if((varInfoSig[workingVar1].transformType==-1&&v1==1)||(varInfoSig[workingVar2].transformType==-1&&v2==1)){
                                                pSig=1.;
                                                pBkg=0.;
                                        }
                                if(v1<varInfoBkg[workingVar1].min||v2<varInfoBkg[workingVar2].min||v1>varInfoBkg[workingVar1].max||v2>varInfoBkg[workingVar2].max) pBkg=0.;
                                        else pBkg=pdfBkg[newVarName]->Interpolate(v1,v2);
                                        if((varInfoBkg[workingVar1].transformType==-1&&v1==0.)||(varInfoBkg[workingVar2].transformType==-1&&v2==0.)){
                                                pBkg=1.;
                                                pSig=0.;
                                        }
                                if(pSig+pBkg>0.)
                                                pLik=pSig/(pSig+pBkg)      ;
                                else
                                                pLik=0.5                   ;
                                newVarVal=0.5*(1.+tanh(0.3*atanh((2*pLik-1.))));
                                        stackNames.push_back(newVarName)        ;
                                        stackVals.push_back(newVarVal)          ;
                                        varTempNameLength = 0                   ;
                                        cur++   ;
                                        break   ;   
                                default:
                                        varTempName[varTempNameLength] = c;
                                        varTempNameLength++;
                                        break;  
                        }
                        cur++;
            }
            ret=stackVals.back();
            stackVals.pop_back();
        }
        return ret;
}
