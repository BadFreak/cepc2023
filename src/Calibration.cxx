#include <fstream>
#include <iostream>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TF1.h>
#include <TKey.h>
#include "Calibration.h"
#include "EBUdecode.h"
#include "Extract.h"

void processCali::Reset() {	
	if(_cellID_new.size()) 	_cellID_new.clear();
	if(_hitADC.size())		_hitADC.clear();
	if(_hitE.size())		_hitE.clear();
	if(_hitX.size())		_hitX.clear();
	if(_hitY.size())		_hitY.clear();
	if(_hitZ.size())		_hitZ.clear();
    if(_temp.size())		_temp.clear();
}

int processCali::TreeType(string inTreeName) { // Return tree type Raw_Hit:0, Digi_Raw_Hit:1
	int treeType = -1;
	if (strcmp(inTreeName.c_str(), "Raw_Hit") == 0) {
		treeType = 0; // data 
	} 
	if (strcmp(inTreeName.c_str(), "Digi_Raw_Hit") == 0) {
		treeType = 1; // simu_digi
	}
	return treeType;
}

void processCali::SetAllBranchAddress(int treeType) {
	inTree->SetBranchAddress("Run_Num",&_run_num);
	inTree->SetBranchAddress("Event_Time",&_event_time);
	inTree->SetBranchAddress("TriggerID",&_triggerID);
	inTree->SetBranchAddress("CellID",&_cellID);
	inTree->SetBranchAddress("BCID",&_bcID);
	inTree->SetBranchAddress("HitTag",&_hitTag);
	inTree->SetBranchAddress("GainTag",&_gainTag);
	inTree->SetBranchAddress("HG_Charge",&_hg_charge);
	inTree->SetBranchAddress("LG_Charge",&_lg_charge);
	inTree->SetBranchAddress("Hit_Time",&_hit_time);
	if(treeType==0) inTree->SetBranchAddress("Temperature",&_tempLayer);
};

void processCali::SetAllBranch() {
	outTree->Branch("Run_Num",&_run_num_new);
	outTree->Branch("Event_Time",&_event_time_new);
	outTree->Branch("Event_Num",&_event_num);
	outTree->Branch("DetectorID",&_detectorID);
	outTree->Branch("CellID",&_cellID_new);
	outTree->Branch("Hit_ADC",&_hitADC);
	outTree->Branch("Hit_Energy",&_hitE);
	outTree->Branch("Hit_X",&_hitX);
	outTree->Branch("Hit_Y",&_hitY);
	outTree->Branch("Hit_Z",&_hitZ);
	outTree->Branch("SiPM_Temp",&_temp);
}

double processCali::MIPCorrection(double MPV) {
	double MPVFitCT[60] = {
		71.365 , 70.7059 , 72.9961 , 72.2648 , 74.3081 , 76.8781 , 79.2888 , 81.4524 , 84.1601 , 85.3968 , 87.9643 , 
		88.5149 , 90.5417 , 92.2638 , 94.1959 , 94.6618 , 95.7191 , 96.0186 , 96.3926 , 97.7468 , 97.8199 , 97.8081 , 
		98.8043 , 99.2342 , 99.7387 , 101.198 , 102.17 , 103.785 , 105.827 , 106.424 , 108.777 , 110.714 , 111.86 , 113.437 ,
		115.041 , 116.128 , 117.451 , 118.914 , 119.554 , 120.591 , 120.931 , 121.436 , 122.4 , 122.773 , 123.095 , 124.73 , 
		125.317 , 126.593 , 127.866 , 129.154 , 130.597 , 132.147 , 132.925 , 134.907 , 136.073 , 137.594 , 138.826 , 140.341 , 141.392 , 142.424
	};
	int i= 0;
	while (MPV > MPVFitCT[i]) {
		i++;
		if (i==59) {
			break;
			return MPV/1.1;
		}
	}
	if (i==0) {return MPV/1.1;}
	if (MPVFitCT[i] < MPVFitCT[i-1]) {return MPV/1.1;}
	double ratio =  (MPVFitCT[i]- MPV) / (MPVFitCT[i] - MPVFitCT[i-1]);
	double MPVCor = (i+80) - ratio;
	return MPVCor;
}

bool processCali::HitCut(double xv, double yv, double slope, double intercept) {
	int interceptOffset = 8;
	double slopeCor = 1/ slope;
	double interceptCor = - intercept / slope;
	if (yv < xv * slopeCor  + interceptCor - interceptOffset) {
		return true;
	} else {
		return false;
	}
}

int processCali::Calibration(string inFileName, string outFileName, Extract & ex) {
	// Input settings
 	gStyle->SetOptFit(1111);
	inFile = TFile::Open(TString(inFileName),"READ");
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
	if (treeType != 0 && treeType != 1) {return 1;}
	if (!inTree) {cout << "!!!!! NO INPUT TREE !!!!!" << endl;}
	SetAllBranchAddress(treeType);
	// Output settings
    outFile = new TFile(TString(outFileName),"RECREATE");
	if (treeType == 0) {
    	outTree = new TTree("Calib_Hit","after calibration");
	} else if (treeType == 1) {
    	outTree = new TTree("Digi_Calib_Hit","after calibration");
	} else {
		return 1;
	}
	SetAllBranch();
	// ==============================================================================
	TH1D *EDep      = new TH1D("EDep","EDep(MeV)",1000,0,5000); // Event energy
	TH1I *hCellHit  = new TH1I("hCellHit","hCellHit",750,0,1500); // Cell hit number
	TH1D *EDepCut   = new TH1D("EDepCut","EDepCut(MeV)",1000,0,5000); // Event energy cut
	TH1I *hLayerHit = new TH1I("hLayerHit","hLayerHit",LAYERNO,0,LAYERNO); // Layer hit number
	
	TH2I *HitNoMap  = new TH2I("HitNoMap",  "HitNoMap",LAYERNO,0,LAYERNO,100,0,100); // Layer hitNo map
	TH2D *HitEMap   = new TH2D("HitEMap","HitEMap",LAYERNO,0,LAYERNO,500,0,500); // Layer Energy map
	TProfile *hHitPro = new TProfile("hHitPro","Profile Layer HitNo",LAYERNO,0,LAYERNO); // profile of HitNoMap
	TProfile *hEPro   = new TProfile("hEPro", "Profile Layer Energy",LAYERNO,0,LAYERNO); // profile of HitEMap
	TH2D *CellHit2EDep = new TH2D("CellHit2EDep","CellHit2EDep",750,0,3000,1000,0,1000); // Event CellHit vs EDep
	
	TH1I * habHit = new TH1I("abHit","abHit@ event level",100, 0, 100); // Event: hit drift number
	TH1I * habHitC = new TH1I("abHitC","abHitC@ event level",100, 0, 100); // Event: hit drift and in the center
	// ==============================================================================
    //                          Calibration signal 
    // ==============================================================================
	int Entries = inTree->GetEntries();
    cout<<" [ Prepare to processCali : ] Total Entries : "<<Entries<<" "<<endl;
	for (int entry=0; entry < Entries; entry++) {
	//for (int entry=0; entry < 10000; entry++) {
        inTree->GetEntry(entry);
        if (entry%10000==0)  
			cout<<" [ Process : Calibration : ]  event : "<<entry<<" "<<endl;
		int driftHitNu = 0;
		int centralChannelNu = 0;
		double EnergyDep = 0;
		double EnergyLayerLevel[LAYERNO] = {0};
		int hitChannelNuEventLevel = 0;
		int hitChannelNuLayerLevel[LAYERNO]={0};//hit channel Number: on each layer
		int hitFlagLayerLevel[LAYERNO]={0};//whether this layer has hit
		int hitLayerNuEventLevel = 0;

		// Reset branch and start to fill
		Reset();
		
		for (size_t ihit = 0; ihit != _hitTag->size(); ++ihit) {
            if (_hitTag->at(ihit)==0) continue;
			int layerid = _cellID->at(ihit)/1e5;
            int chipid	 = ( _cellID->at(ihit)%100000) / 1e4;
            int chanid   = _cellID->at(ihit)%100;	
			if (layerid>29) continue;// last two layers are double_side readout
            if (chipid==5 && chanid>29) {cout << "Wrong cellID decoding ..." << endl ;continue;}// only 210 channels per layer
			double caliCharge = _lg_charge->at(ihit);
            double caliTime   = _hg_charge->at(ihit);
			
			// reconstruction temperature distribution at each SiPM location
			double *_position = EBUdecode(layerid,chipid,chanid);
			double SiPMtemp;
			if (treeType == 0 && (_tempLayer->at(layerid).size()!=0)) {
				SiPMtemp = tempReconstruction(layerid, _position, _tempLayer->at(layerid));
			} else {
				SiPMtemp = 23;
			}
			// jiaxuan : temperature depedence correction
			// caliCharge = (1+(SiPMtemp-23)*0.03)*caliCharge;
            // caliTime = (1+(SiPMtemp-23)*0.03)*caliTime;
            caliCharge -= ex.pedMeanLow[layerid][chipid][chanid];
            caliTime   -= ex.pedMeanHigh[layerid][chipid][chanid];
            
			// jiaxuan :
			/*if (caliTime > 0 && caliTime <2000) {
				if (HitCut(caliTime, caliCharge, ex.HLRatio[layerid][chipid][chanid] ,ex.HLIntercept[layerid][chipid][chanid])) {
			//if(0){
					++driftHitNu;
					if ((fabs(*_position) < 67 && fabs(*_position+1) < 67)) {
						++centralChannelNu;
					}
				}
			}*/
			//if(caliTime<10 && _hit_time->at(ihit) < 0) continue; // exclude bug of SP2E of High Low Gain mode file

            // jiaxuan : high-low conversion without intercept
            caliCharge *= ex.HLRatio[layerid][chipid][chanid];
            //caliCharge += ex.HLIntercept[layerid][chipid][chanid];

			// MIP response calibration for Gain Mode
			int SiPMType = ((layerid < 4) || (layerid >27)) ? 1 : 0;
			// Exclude crosstalk: Not ever used here, change by using different MIP calibrated file
			// double MIPCT;
			// if (SiPMType == 0 ) {
			// 	MIPCT = ex.MIPPeakValue[layerid][chipid][chanid] / 1.1;
			// 	//MIPCT = (ex.MIPPeakValue[layerid][chipid][chanid] + 12.38) / 1.135;
			// 	//MIPCT = MIPCorrection(ex.MIPPeakValue[layerid][chipid][chanid]);
			// } else {
			// 	MIPCT = (ex.MIPPeakValue[layerid][chipid][chanid]) / 1.15;
			// 	//MIPCT = (ex.MIPPeakValue[layerid][chipid][chanid] + 16.58) / 1.162;
			// }

			// to clarify HL-gain mode or auto-gain mode
			double hitadc = 0;
			double hitenergy = 0;
			if (treeType == 0) {
				if (_hit_time->at(ihit)< 0) {// cout << "Processing HL-gain-mode file now ... ..." << endl;
					hitadc = _hg_charge->at(ihit)< ex.HighADCSwitch[layerid][chipid][chanid] ? caliTime: caliCharge; //ADC
				} else {// cout << "Processing auto-gain-mode file now ... ..." << endl;
					hitadc = _lg_charge->at(ihit)==-1 ? caliTime: caliCharge;
				}
			}	
			if (treeType == 1) {
				hitadc = _lg_charge->at(ihit)==-1 ? caliTime: caliCharge;
			}
			if (ex.MIPPeakValue[layerid][chipid][chanid] <= 50) continue;
			hitenergy = hitadc / ex.MIPPeakValue[layerid][chipid][chanid] * MIPEnergy;
			
			_hitADC.push_back(hitadc);
			_hitE.push_back(hitenergy);	
			_cellID_new.push_back(_cellID->at(ihit));
			_hitX.push_back(*_position);
			_hitY.push_back(*(_position+1));
			_hitZ.push_back(*(_position+2));
            _temp.push_back(SiPMtemp);
			
			EnergyDep += hitenergy ; // Event energy
			EnergyLayerLevel[layerid] += hitenergy; // Event Layer Energy
			
			hitChannelNuEventLevel++;
			hitChannelNuLayerLevel[layerid]++;
			hitFlagLayerLevel[layerid] = 1;
		}
		_run_num_new 	= _run_num;
		_event_time_new = _event_time;
		_event_num 		= _triggerID;
		_detectorID 	= 0;// 0 for Sci-ECAL
		outTree->Fill();
		
		habHit->Fill(driftHitNu);	
		habHitC->Fill(centralChannelNu);	
		for (int i_layer = 0; i_layer < LAYERNO; i_layer++) {	
			HitNoMap->Fill(i_layer,hitChannelNuLayerLevel[i_layer]);
			HitEMap->Fill(i_layer,EnergyLayerLevel[i_layer]);
			hitLayerNuEventLevel += hitFlagLayerLevel[i_layer];
		}
		hCellHit->Fill(hitChannelNuEventLevel);
		hLayerHit->Fill(hitLayerNuEventLevel);
		EDep->Fill(EnergyDep);
		if (driftHitNu < 2 && centralChannelNu < 2) EDepCut->Fill(EnergyDep);
		CellHit2EDep->Fill(EnergyDep,hitChannelNuEventLevel);
    }
	
	//CrystallBall Fit
	gStyle->SetOptFit(1111);
	double constant = EDep->GetBinContent(EDep->GetMaximumBin());
	double mean = EDep->GetMean();
	double sigma = EDep->GetRMS();
	TF1 * fCry = new TF1("fCry","crystalball");
	TF1 * fGaus = new TF1("fGaus","gaus(0)");
	EDep->Fit(fGaus);
	constant = fGaus->GetParameter(0);
	mean = fGaus->GetParameter(1);
	sigma = fGaus->GetParameter(2);
	cout << constant << "	" << mean << "	" << sigma << endl;
	//fCry->SetParameters(n,mu,sigma,alpha,N);
	fCry->SetParameters(constant,mean,30,1,10000);
	fCry->SetParLimits(0,0.7*constant,1.3*constant);
	fCry->SetParLimits(1,0.7*mean,1.3*mean);
	EDep->Fit(fCry);
    //cout<<"Total Entries : "<<Entries<<" , Layer>=22 : "<<Entries_1<<" , Strips<=64 : "<<Entries_2<<endl;
    
	outFile->cd();
	EDep->SetLineWidth(2);
	EDep->Write();
	EDepCut->SetLineWidth(2);
	EDepCut->Write();
	HitNoMap->Write();
	hCellHit->Write();
	hLayerHit->Write();
	CellHit2EDep->SetOption("colz");
	CellHit2EDep->GetXaxis()->SetTitle("Energy[MeV]");
	CellHit2EDep->GetYaxis()->SetTitle("HitNo");
	CellHit2EDep->Write();
	HitEMap->Write();
	
	habHit->Write();
	habHitC->Write();

	hHitPro = HitNoMap->ProfileX();
	hHitPro->Write();
	hEPro = HitEMap->ProfileX();
	hEPro->Write();
    
	outTree->Write("",TObject::kOverwrite);
    outFile->Close();
    inFile->Close();
    cout<<" [ End of Calibration Process ]  "<<endl;
	return 0;
}


