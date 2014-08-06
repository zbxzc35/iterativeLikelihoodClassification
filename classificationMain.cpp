/* 
 * File:   classificationMain.cpp
 * Author: kaiwu
 *
 * Created on October 24, 2013, 10:15 AM
 */

#include <iostream>
#include <time.h>
#include "kMvaClassification.h"
using namespace std;
int main(int argc,char*argv[]){
        time_t start,end;
        time(&start);
        //Get the configuation file name for input
        string configurationFile=argv[1];
        //Declare a new kMvaClassification object
        kMvaClassification* kMva=new kMvaClassification(configurationFile);
        //Call start() procedure to start the training
        kMva->start();
        time(&end);
        cout<<"Time used "<<difftime(end,start)<<" seconds"<<endl;
        return 0;
}