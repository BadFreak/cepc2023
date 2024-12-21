#include <iostream>
#include<stdio.h>
#include<stdlib.h>
#include "Pedestal.h"
//#include "Calibration.h"
#include"TROOT.h"//https://root-forum.cern.ch/t/trees-with-vectors/8276/2

using namespace std;

int main(int argc,char *argv[])
{
    double start = clock();
    gROOT->ProcessLine("#include <vector>");//https://root-forum.cern.ch/t/trees-with-vectors/8276/2
    cout<<"start of main program"<<endl;
    treePedestal p;
    p.Pedestal(argv[1],argv[2],argv[3]);

    cout<<"end of main program"<<endl;
    double end = clock();
    cout<<"end of RawToRoot : Time : "<<(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
