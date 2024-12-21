#include <iostream>
#include <TROOT.h>
#include<stdio.h>
#include<stdlib.h>
#include "Extract.h"
#include"InterCalib.h"

using namespace std;

Int_t main(int argc,char *argv[])
{
    double start = clock();
  
	Extract * ex = new Extract();
	ex->PedestalExtract(argv[4]);
    interCalib calib(argv[1],argv[2],argv[3],ex);
    
    double end = clock();
    cout<<"end of HL_RATIO : Time : " << (end-start)/CLOCKS_PER_SEC << endl;
	delete ex;
  	return 0;
}
