double* EBUdecode(int LayerIDs, int ChipIDs, int ChannelIDs)
{
	// static const int layerNu = 32;
	static const int chipNu = 6;
	static const int chnNu = 36;
	int decodeID[chipNu][chnNu] = { 0, 42, 1, 43, 2, 44, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 54, 13, 55, 14, 56, 15, 57, 16, 58, 17, 59, 18, 60, 19, 61, 20, 62, 21, 22, 23,
									24, 66, 25, 67, 26, 68, 27, 69, 28, 70, 29, 71, 30, 72, 31, 73, 32, 74, 33, 75, 34, 76, 35, 77, 36, 78, 37, 79, 38, 80, 39, 81, 40, 82, 41, 83,
									149, 148, 147, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 63, 64, 65, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
									150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197, 156, 198, 157, 199, 158, 200, 159, 201, 160, 202, 161, 203, 162, 204, 163, 205, 164, 206, 165, 207, 166, 208, 167, 209,
									191, 190, 189, 188, 146, 187, 145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181, 139, 180, 138, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
									137, 136, 135, 134, 133, 132, 131, 130, 129, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 45, 46, 47, 48, 49, 50, 51, 52, 53, 210, 210, 210, 210, 210, 210 };

	int ScintillatorIDs = decodeID[ChipIDs][ChannelIDs];
	double layerZ;

	const double _xInterval = 5.3;	// 300 um gap in width direction
	const double _yInterval = 45.4; // 400 um gap in length direction
	const int rowNu = 42;
	const int columnNu = 5;
	int _xID = ScintillatorIDs % rowNu; // 0-41
	int _yID = ScintillatorIDs / rowNu; // 0-4
	static double _position[3] = { 0 };
	double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;//42 position
	double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;//5 position

	// for prototype test
	if (LayerIDs % 2 == 0) {
		_position[0] = -y0;
		_position[1] = -x0;
	}
	if (LayerIDs % 2 == 1) {
		_position[0] = -x0;
		_position[1] = -y0;
	}
	if (LayerIDs % 2 == 0)
		layerZ = 1 + LayerIDs / 2 * 19.9;
	//else layerZ = 12.95+(LayerIDs-1)/2*19.9;
	else
		layerZ = 12.2 + (LayerIDs - 1) / 2 * 19.9;
	_position[2] = layerZ;
	return _position;
}

void DrawBadPosition() {
    string inFileName = "/home/wangjx/HEP/CEPC2023/SiWECAL_analysis_2023_update/share/all_auto_muon_v4_Prime_update.root";
    TFile * mipFile = TFile::Open(TString(inFileName),"READ");
    if (!mipFile) cout << "=========== No MIP File ============" << endl;
    TTree * mipTree = (TTree*)mipFile->Get("MIP_Fit");
    if (!mipTree) cout << "=========== No MIP Tree ========" << endl;
    int _cellID;
    bool _isGoodChan;
    mipTree->SetBranchAddress("CellID", &_cellID);
    mipTree->SetBranchAddress("IsGoodChan", &_isGoodChan);
    TH2D* even = new TH2D("even","even",5,-111.5,111.5,42,-111.1,111.5);
    TH2D* odd= new TH2D("odd","odd",42,-111.1,111.5,5,-111.5,111.5);
    for (int i = 0; i < mipTree->GetEntries(); i++) {
        mipTree->GetEntry(i);
        int cellid = _cellID;
        int layerid = cellid/100000;
        int chipid = (cellid%100000)/10000;
        int channid = cellid%100;

        double *position = EBUdecode(layerid, chipid, channid);
        double xpos = *(position);
        double ypos = *(position+1);
        if (!_isGoodChan) {
            layerid % 2 == 0 ? even->Fill(xpos,ypos) : odd->Fill(xpos,ypos);
        }
    }
    
    TCanvas *c1 = new TCanvas("c1","c1",1300,600);
    c1->Divide(2,1);
    c1->cd(1);
    even->GetXaxis()->SetTitle("X Position (mm)");
    even->GetYaxis()->SetTitle("Y Position (mm)");
    even->GetXaxis()->SetTitleSize(0.045);
    even->GetYaxis()->SetTitleSize(0.045);
    even->GetXaxis()->SetLabelSize(0.045);
    even->GetYaxis()->SetLabelSize(0.045);;
    even->Draw("colz text");
    
    c1->cd(2);
    odd->GetXaxis()->SetTitle("X Position (mm)");
    odd->GetYaxis()->SetTitle("Y Position (mm)");
    odd->GetXaxis()->SetTitleSize(0.045);
    odd->GetYaxis()->SetTitleSize(0.045);
    odd->GetXaxis()->SetLabelSize(0.045);
    odd->GetYaxis()->SetLabelSize(0.045);
    odd->Draw("colz text");

}