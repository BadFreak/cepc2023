#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include "SimpleCalibration.h"
#include "MIPFit.h"
#include "Extract.h"
#include "TROOT.h" //https://root-forum.cern.ch/t/trees-with-vectors/8276/2

using namespace std;

int main(int argc, char* argv[])
{
	double start = clock();
	gROOT->ProcessLine("#include <vector>"); //https://root-forum.cern.ch/t/trees-with-vectors/8276/2
	cout << "start of main program" << endl;

	MIPCalibration t;
	Extract ex;

	//t.MIPFit(argv[5], argv[6], argv[7]);
	ex.PedestalExtract(argv[2]);
	ex.MIPPeakValueExtract(argv[3]);
	bool flag = (stoi(argv[6]) == 1);
	t.MIPFit(argv[1], ex, argv[4], argv[5], flag);

	cout << "end of main program" << endl;
	double end = clock();
	cout << "end of RawToRoot : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
	return 0;
}
