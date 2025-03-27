#include <string>
#include <vector>
#include <map>

int TreeType(string inTreeName) { // Return inTree type Raw_Hit:0, Digi_Raw_Hit:1
	int treeType = -1;
	if (strcmp(inTreeName.c_str(), "Calib_Hit") == 0) {
		treeType = 0; // data 
	} 
	if (strcmp(inTreeName.c_str(), "Digi_Calib_Hit") == 0) {
		treeType = 1; // simu_digi
	}
	return treeType;
}

void ReadGoodChan(string inFileName, map<int,bool> &myMap) {
    TFile * mipFile = TFile::Open(TString(inFileName),"READ");
    if (!mipFile) cout << "=========== No MIP File ============" << endl;
    TTree * mipTree = (TTree*)mipFile->Get("MIP_Fit");
    if (!mipTree) cout << "=========== No MIP Tree ========" << endl;
    int _cellID;
    bool _isGoodChan;
    mipTree->SetBranchAddress("CellID", &_cellID);
    mipTree->SetBranchAddress("IsGoodChan", &_isGoodChan);
    for (int i = 0; i < mipTree->GetEntries(); i++) {
        mipTree->GetEntry(i);
        myMap[_cellID] = _isGoodChan;
    }
    mipFile->Close();
}

bool IsOverThreshold(double energy) {
    const double Thres = 0.8;
    const double MIPEnergy = 0.305;// MeV
    if (energy > Thres*MIPEnergy) 
        return true;
    else
        return false;
}

double TemperatureCorrection(double energy, double temp) {
    const double BaseTemp = 20.0;// degree
    const double TempCoef = -0.03;
    return energy * (1 - TempCoef * (temp - BaseTemp));
}

void DrawEnergy() {
    //string inFileName = "/home/wangjx/HEP/CEPC2023/Result_SPS/calib/e-/40GeV/40GeV_hl_electron.root";
    //string inMIPFileName = "/home/wangjx/HEP/CEPC2023/SiWECAL_analysis_2023_update/share/all_auto_muon_v4_Prime_update.root";
    string inFileName = "/mnt2/USTC/jxwang/CEPC_ScECAL/CEPC2023/Result/calib/e-/40GeV/40GeV_hl_electron.root";
    string inMIPFileName = "/mnt2/USTC/jxwang/CEPC_ScECAL/SiWECAL_analysis_2023_update/share/all_auto_muon_v4_Prime_update.root";
	TFile * inFile = TFile::Open(TString(inFileName),"READ");
    TTree * inTree = NULL;
    TIter next(inFile->GetListOfKeys());
	TKey *key;
	string inTreeName;
	while ((key = (TKey*)next())) {
		if (strcmp(key->GetClassName(), "TTree") == 0) { 
			inTree = (TTree*)inFile->Get(key->GetName());
			inTreeName = key->GetName();
		}
	}
	cout << "Tree name: " << inTree->GetName() << endl;
	int treeType = TreeType(inTreeName);
    cout << treeType << endl;	
	if (treeType != 0 && treeType != 1) {return;}
	if (!inTree) {cout << "!!!!! NO INPUT TREE !!!!!" << endl;}
    int Run_Num, Event_Time, Event_Num, DetectorID;
    std::vector<int> *CellID = nullptr;
    std::vector<double> *Hit_ADC = nullptr, *Hit_Energy = nullptr, *Hit_X = nullptr, *Hit_Y = nullptr, *Hit_Z = nullptr, *SiPM_Temp = nullptr;
    inTree->SetBranchAddress("Run_Num", &Run_Num);
    inTree->SetBranchAddress("Event_Time", &Event_Time);
    inTree->SetBranchAddress("Event_Num", &Event_Num);
    inTree->SetBranchAddress("DetectorID", &DetectorID);
    inTree->SetBranchAddress("CellID", &CellID);
    inTree->SetBranchAddress("Hit_ADC", &Hit_ADC);
    inTree->SetBranchAddress("Hit_Energy", &Hit_Energy);
    inTree->SetBranchAddress("Hit_X", &Hit_X);
    inTree->SetBranchAddress("Hit_Y", &Hit_Y);
    inTree->SetBranchAddress("Hit_Z", &Hit_Z);
    inTree->SetBranchAddress("SiPM_Temp", &SiPM_Temp);
    
    map<int,bool> isGoodChan;
    ReadGoodChan(inMIPFileName, isGoodChan);
   
    TH1D *EDep = new TH1D("EDep","Energy Deposition",500,400,2200);
    cout << "==================" << endl;
    for (int i=0; i<inTree->GetEntries(); i++) {
        inTree->GetEntry(i);
        double sumEnergy = 0;
        for (size_t j=0; j < Hit_Energy->size(); ++j) {
            int cellID = CellID->at(j);
            double hitEnergy = Hit_Energy->at(j);
            double SiPMTemp = SiPM_Temp->at(j);
            //if (!isGoodChan[cellID]) continue;
            if (!IsOverThreshold(hitEnergy)) continue;
            //hitEnergy = TemperatureCorrection(hitEnergy, SiPM_Temp->at(j));
            sumEnergy += hitEnergy;
        }
        EDep->Fill(sumEnergy);
        if (i%10000 == 0) {cout << "Entry No : " << i << " Sum Energy : "<< sumEnergy << endl;}
    }
    //Fit
    gStyle->SetOptFit(1111);
	double constant = EDep->GetBinContent(EDep->GetMaximumBin());
	double mean = EDep->GetMean();
	double sigma = EDep->GetRMS();
    double fitLeft, fitRight;
	TF1 * fCry = new TF1("fCry","crystalball");
	TF1 * fGaus = new TF1("fGaus","gaus(0)");
	EDep->Fit(fGaus,"0QR","",1000,2200);
	constant = fGaus->GetParameter(0);
	mean = fGaus->GetParameter(1);
	sigma = fGaus->GetParameter(2);
    fitLeft = mean - 1*sigma;
    fitRight = mean + 2*sigma;
	//cout << constant << "	" << mean << "	" << sigma << endl;
	//fCry->SetParameters(n,mu,sigma,alpha,N);
	fCry->SetParameters(constant,mean,30,1,10000);
	fCry->SetParLimits(0,0.7*constant,1.3*constant);
	fCry->SetParLimits(1,0.7*mean,1.3*mean);
	EDep->Fit(fCry,"R","",fitLeft,fitRight);
    //cout<<"Total Entries : "<<Entries<<" , Layer>=22 : "<<Entries_1<<" , Strips<=64 : "<<Entries_2<<endl;
    
    TCanvas * cEDep = new TCanvas("cEDep","Energy Deposition",800,600);
    cEDep->SetMargin(0.11,0.11,0.11,0.11);
    cEDep->SetGrid();
    cEDep->cd();
    EDep->SetLineWidth(2);
    EDep->GetXaxis()->SetTitle("Energy [MeV]");
    EDep->GetYaxis()->SetTitle("Counts");
    EDep->GetXaxis()->SetTitleSize(0.045);
    EDep->GetYaxis()->SetTitleSize(0.045);
    EDep->GetXaxis()->SetTitleOffset(0.9);
    EDep->GetYaxis()->SetTitleOffset(1.2);
    EDep->GetXaxis()->SetLabelSize(0.045);
    EDep->GetYaxis()->SetLabelSize(0.045);

    EDep->Draw();
    cEDep->SaveAs("EDep.png");

}
