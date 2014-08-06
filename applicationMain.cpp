/* 
 * File:   applicationMain.cpp
 * Author: kaiwu
 *
 * Created on October 24, 2013, 10:17 AM
 */

#include <iostream>
#include <string>
#include <vector>
#include "TObject.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "kMvaApplication.h" 
using namespace std;
int main(int argc,char*argv[]){
        string database_file_name       ;
        string input_file_name          ;
        string input_tree_name          ;
        string output_file_name         ;
        string output_tree_name         ;
        int nVariables                          ;
        vector<string>  variables       ;
        database_file_name=argv[1]      ;
        input_file_name   =argv[2]      ;
        input_tree_name   =argv[3]      ;
        output_file_name  =argv[4]      ;
        output_tree_name  =argv[5]      ;
        nVariables        =atoi(argv[6]);
        for(int i1=0;i1<nVariables;i1++){
                variables.push_back(argv[i1+7]);
        }
        cout<<"database_file_name "<<database_file_name<<endl   ;
        cout<<"input_file_name    "<<input_file_name<<endl      ;
        cout<<"input_tree_name    "<<input_tree_name<<endl      ;
        cout<<"output_file_name   "<<output_file_name<<endl     ;
        cout<<"output_tree_name   "<<output_tree_name<<endl     ;
        cout<<"nVariables         "<<nVariables<<endl           ;
        for(int i1=0;i1<nVariables;i1++)
                cout<<variables[i1]<<" ";
        cout<<endl                      ;
        //Declare a new kMvaApplication obejct with database root file name
        kMvaApplication* kMva=new kMvaApplication(database_file_name);
        //Add the input variables' names
        for(int i1=0;i1<nVariables;i1++)
                kMva->addVariable(variables[i1])        ;                       //Add variable names
        vector<float> varVal(variables.size(),0)        ;       //store variable values;
        //Input root file name
        TFile* fin=TFile::Open(input_file_name.c_str()) ;
        if(!fin){
                cout<<"Can not open "<<input_file_name<<endl;
                return -1;
        }
        //Get input TTree
        TTree* treeIn=(TTree*) fin->Get(input_tree_name.c_str())        ;
        if(!treeIn){
                cout<<"Can not load "<<input_tree_name<<" from "<<database_file_name<<endl;
                return -1;
        }
        //Bind Input Variables' TBranch
        TBranch* br[1000];
        for(int i1=0;i1<nVariables;i1++){
                treeIn->SetBranchAddress(variables[i1].c_str(),&(varVal[i1]),&(br[i1]));
        }
        //Open Output file
        TFile* fout=TFile::Open(output_file_name.c_str(),"update");
        fout->cd();
        //Output TTree
        TTree* treeOut=new TTree(output_tree_name.c_str(),output_tree_name.c_str());
        for(int i1=0;i1<nVariables;i1++)
                        treeOut->Branch(variables[i1].c_str(),&(varVal[i1]),Form("%s/F",variables[i1].c_str()));
        double kLik;
        treeOut->Branch("kLik",&kLik,"kLik/D")  ;
        int entries=treeIn->GetEntries()                ;
        //Loop all the input data points and calculate the likelihood value
        for(int i1=0;i1<entries;i1++){
                for(int i2=0;i2<nVariables;i2++)        br[i2]->GetEntry(i1);
                kLik=kMva->eval(varVal); //calculate likelihood for varVal
                treeOut->Fill();
        }
        //Write the output
        fout->cd();
        treeOut->Write("",TObject::kOverwrite);
        fout->Close();
        return 0;
}
