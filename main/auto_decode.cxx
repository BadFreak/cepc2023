#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "AutoGainDecode.h"

using namespace std;

Int_t main(int argc, char* argv[])
{
	double start = clock();

	raw2RootAuto tw;
	tw.AutoGainDecode(argv[1], argv[2]);

	double end = clock();
	cout << "end of RawToRoot under auto gain mode : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
}
