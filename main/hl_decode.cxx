#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "HLGainDecode.h"

using namespace std;

Int_t main(int argc, char* argv[])
{
	double start = clock();

	raw2RootHL tw;
	tw.HLGainDecode(argv[1], argv[2]);

	double end = clock();
	cout << "end of RawToRoot under hl gain mode : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
}
