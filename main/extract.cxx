#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include "SimpleCalibration.h"
#include "Extract.h"
#include "TROOT.h" //https://root-forum.cern.ch/t/trees-with-vectors/8276/2

using namespace std;

int main(int argc, char* argv[])
{
	double start = clock();
	gROOT->ProcessLine("#include <vector>"); //https://root-forum.cern.ch/t/trees-with-vectors/8276/2
	cout << "start of main program" << endl;

	Extract ex;
	ex.PedestalExtract(argv[1]);
	ex.HLGainRatioExtract(argv[2]);
	ex.MIPPeakValueExtract(argv[3]);
	ex.ReadAllDeadChannel();
	
	cout << "end of main program" << endl;
	double end = clock();
	cout << "end of RawToRoot : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
	return 0;
}
