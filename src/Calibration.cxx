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

void processCali::Reset()
{	
	if(_cellID_new.size()) 	_cellID_new.clear();
	if(_hitADC.size())		_hitADC.clear();
	if(_hitE.size())		_hitE.clear();
	if(_hitX.size())		_hitX.clear();
	if(_hitY.size())		_hitY.clear();
	if(_hitZ.size())		_hitZ.clear();
    if(_temp.size())		_temp.clear();
}

int processCali::TreeType(string inFileName)
{
	int treeType = -1;
	inFile = TFile::Open(TString(inFileName),"READ");
    TIter next(inFile->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)next())) {
		if (strcmp(key->GetClassName(), "TTree") == 0) { 
			inTree = (TTree*)inFile->Get(key->GetName()); 
			cout << "Tree name: " << inTree->GetName() << endl;
			if (strcmp(inTree->GetName(), "Raw_Hit") == 0)  {
				treeType = 0; // data 
			} 
			if (strcmp(inTree->GetName(), "Digi_Raw_Hit") == 0)  {
				treeType = 1; // simu_digi
			}
		}
	}
	inFile->Close();
	return treeType;
}

void processCali::SetAllBranchAddress(int treeType)
{
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

void processCali::SetAllBranch()
{
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
	outTree->Branch("NewTemperature",&_temp);
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

int processCali::Calibration(string inFileName, string outFileName, Extract & ex)
{
 	gStyle->SetOptFit(1111);
	int treeType = -1;
	treeType = TreeType(inFileName);	
    inFile = TFile::Open(TString(inFileName),"READ");
	if (treeType == 0) {
		inTree = (TTree*)inFile->Get("Raw_Hit");
	} else if (treeType == 1) {
		inTree = (TTree*)inFile->Get("Digi_Raw_Hit");
		//inTree = (TTree*)inFile->Get("MC_Truth");
	} else {
		return 1;
	}
	if (!inTree) {cout << "!!!!! NO INPUT TREE !!!!!" << endl;}
	SetAllBranchAddress(treeType);

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
    //                          Calibration signal 
    // ==============================================================================
	TH1D *Edep = new TH1D("Edep","Edep(MeV)",1000,0,5000);
	TH1D *EdepCut = new TH1D("EdepCut","EdepCut(MeV)",1000,0,5000);
	TH1I *HitLayerNo = new TH1I("HitLayerNo","HitLayerNo",LAYERNO,0,LAYERNO);
	TH1I *HitNo = new TH1I("HitNo","HitNo",750,0,1500);
	TH2I *HitMap = new TH2I("HitMap","HitMap",LAYERNO,0,LAYERNO,100,0,100);
	TProfile *hitpro = new TProfile("hitpro","Profile Layer HitNo",30,0,30);
	TProfile *hpro = new TProfile("hpro","Profile Layer Energy",30,0,30);
	//TH2D *LayerEdep = new TH2D("LayerEdep","Layer Energy",30,0,30,100,0,400);
	TH2D *LayerEdep = new TH2D("LayerEdep","Layer Energy",30,0,30,500,0,500);
	TH2D *HitNo2EnergyDep = new TH2D("HitNo2EnergyDep","HitNo2EnergyDep",750,0,3000,1000,0,1000);

	TH1I * habHit = new TH1I("abHit","abHit@ event level",100, 0, 100);
	TH1I * habHitC = new TH1I("abHitC","abHitC@ event level",100, 0, 100);
	int Entries = inTree->GetEntries();
    cout<<" [ Prepare to processCali : ]oTotal Entries : "<<Entries<<" "<<endl;
	for(int entry=0; entry < Entries; entry++) {
	//for(int entry=0; entry < 10000; entry++) {
        inTree->GetEntry(entry);
        if (entry%10000==0)  
			cout<<" [ Process : Calibration : ]  event : "<<entry<<" "<<endl;
		int driftHitNu = 0;
		int positionHitNu = 0;
		double EnergyDep = 0;
		double EnergyLayerLevel[LAYERNO] = {0};
		int hitChannelNuEventLevel = 0;
		int hitChannelNuLayerLevel[LAYERNO]={0};//hit channel Number: on each layer
		int hitFlagLayerLevel[LAYERNO]={0};//whether this layer has hit
		int hitLayerNuEventLevel = 0;
		/////////////////////////////////
		//for	(size_t i_hit = 0; i_hit < _cellID->size(); i_hit++) {
		//	int layerID = _cellID->at(i_hit) / 1e5;
		//	double posX = _hitX->at(i_hit);
		//	double posY = _hitY->at(i_hit);
		//	double hitE = _hitE->at(i_hit);
		//	averX[layerID]  += posX * hitE;
		//	averY[layerID]  += posY * hitE;
		//	layerE[layerID] += hitE;
		//}
		//////////////////////////////////

		Reset();
		_run_num_new 	= _run_num;
		_event_time_new = _event_time;
		_event_num 		= _triggerID;
		_detectorID 	= 0;
		
		for (size_t ihit = 0; ihit != _hitTag->size(); ++ihit) {
            if(_hitTag->at(ihit)==0) continue;
			//if(IsDeadChannel(_cellID->at(ihit))) continue;

			int _layerID = _cellID->at(ihit)/1e5;
            int _chipID	 = ( _cellID->at(ihit)%100000) / 1e4;
            int _chnID   = _cellID->at(ihit)%100;
			if(_layerID>29) continue;// last two layers are double_side readout
            if(_chipID==5 && _chnID>29) continue;// only 210 channels per layer
            
			double caliCharge = _lg_charge->at(ihit);
            double caliTime   = _hg_charge->at(ihit);
			// reconstruction temperature distribution at each SiPM location
			double *_position = EBUdecode(_layerID,_chipID,_chnID);
			double SiPMtemp;
			if (treeType == 0) {	
				if(_tempLayer->at(_layerID).size()!=0) {
					SiPMtemp = tempReconstruction(_layerID, _position, _tempLayer->at(_layerID));
				} else {
					SiPMtemp = 23;// SiPMtemp=20
				}
			} else {
				SiPMtemp = 23;
			}
			// jiaxuan : temperature depedence correction
			// caliCharge = (1+(SiPMtemp-23)*0.03)*caliCharge;
            // caliTime = (1+(SiPMtemp-23)*0.03)*caliTime;
			
			// cout << "pedlow:	" << ex.pedMeanLow[_layerID][_chipID][_chnID]<< endl; 
			// cout << "pedhigh:	" << ex.pedMeanHigh[_layerID][_chipID][_chnID]<< endl; 
			// cout << "hlratio:	" << ex.HLRatio[_layerID][_chipID][_chnID]<< endl; 
			// cout << "MIPvalue:	" << ex.MIPPeakValue[_layerID][_chipID][_chnID]<< endl; 
			// cout << "ADCSwitch:	" << ex.HighADCSwitch[_layerID][_chipID][_chnID]<< endl; 
            caliCharge -= ex.pedMeanLow[_layerID][_chipID][_chnID];
            caliTime   -= ex.pedMeanHigh[_layerID][_chipID][_chnID];
            
			// jiaxuan :
			/*if (caliTime > 0 && caliTime <2000) {
				if (HitCut(caliTime, caliCharge, ex.HLRatio[_layerID][_chipID][_chnID] ,ex.HLIntercept[_layerID][_chipID][_chnID])) {
			//if(0){
					++driftHitNu;
					if ((fabs(*_position) < 67 && fabs(*_position+1) < 67)) {
						++positionHitNu;
					}
				}
			}*/
			//if(caliTime<10 && _hit_time->at(ihit) < 0) continue; // exclude bug of SP2E of High Low Gain mode file

            // jiaxuan : high-low conversion with intercept
            caliCharge *= ex.HLRatio[_layerID][_chipID][_chnID];
            //caliCharge += ex.HLIntercept[_layerID][_chipID][_chnID];

			// MIP response calibration for Gain Mode
			int SiPMType = ((_layerID < 4) || (_layerID >27)) ? 1 : 0;
			if (ex.MIPPeakValue[_layerID][_chipID][_chnID] <= 5) {
				ex.MIPPeakValue[_layerID][_chipID][_chnID] = (SiPMType == 0) ? 90 : 400;
			}
			// exclude crosstalk
			// double MIPCT;
			// if (SiPMType == 0 ) {
			// 	MIPCT = ex.MIPPeakValue[_layerID][_chipID][_chnID] / 1.1;
			// 	//MIPCT = (ex.MIPPeakValue[_layerID][_chipID][_chnID] + 12.38) / 1.135;
			// 	//MIPCT = MIPCorrection(ex.MIPPeakValue[_layerID][_chipID][_chnID]);
			// } else {
			// 	MIPCT = (ex.MIPPeakValue[_layerID][_chipID][_chnID]) / 1.15;
			// 	//MIPCT = (ex.MIPPeakValue[_layerID][_chipID][_chnID] + 16.58) / 1.162;
			// }
            
			//caliCharge /= ex.MIPPeakValue[_layerID][_chipID][_chnID];
            //caliTime   /= ex.MIPPeakValue[_layerID][_chipID][_chnID];
			
			// to clarify HL-gain mode or auto-gain mode
			double hitadc = 0;
			double hitenergy = 0;
			if ( treeType==0 ) {
				if (_hit_time->at(ihit)< 0) {
					// cout << "Processing HL-gain-mode file now ... ..." << endl;
					hitadc = _hg_charge->at(ihit)< ex.HighADCSwitch[_layerID][_chipID][_chnID] ? caliTime: caliCharge;//MeV
				} else {
					// cout << "Processing auto-gain-mode file now ... ..." << endl;
					hitadc = _lg_charge->at(ihit)==-1 ? caliTime: caliCharge;
				}
			}	
			if (treeType == 1) {
				hitadc = _lg_charge->at(ihit)==-1 ? caliTime: caliCharge;
			}
			hitenergy = hitadc / ex.MIPPeakValue[_layerID][_chipID][_chnID] * MIPEnergy;
			
			// jiaxuan mip cut
			double setThreshold = 0.8;
			//double threshold = setThreshold + gRandom->Gaus(setThreshold,0.05);
			if (hitenergy < setThreshold * 0.305) continue; 
			_hitADC.push_back(hitadc);
			_hitE.push_back(hitenergy);	
			_cellID_new.push_back(_cellID->at(ihit));
			_hitX.push_back(*_position);
			_hitY.push_back(*(_position+1));
			_hitZ.push_back(*(_position+2));
            _temp.push_back(SiPMtemp);
			//if(entry==0)cout << fixed << setprecision(2) << hitenergy << " " << *_position << " " <<  *(_position+1) <<  " " << *(_position+2)<< endl; // jiaxuan

			EnergyDep += hitenergy ;
			EnergyLayerLevel[_layerID] += hitenergy;
			
			hitChannelNuEventLevel++;
			hitChannelNuLayerLevel[_layerID]++;
			hitFlagLayerLevel[_layerID] = 1;
		}
		
		habHit->Fill(driftHitNu);	
		habHitC->Fill(positionHitNu);	
		
		for (int i_layer = 0; i_layer < LAYERNO; i_layer++) {	
			HitMap->Fill(i_layer,hitChannelNuLayerLevel[i_layer]);
			hitLayerNuEventLevel += hitFlagLayerLevel[i_layer];
			LayerEdep->Fill(i_layer,EnergyLayerLevel[i_layer]);
		}
		
		HitNo->Fill(hitChannelNuEventLevel);
		HitLayerNo->Fill(hitLayerNuEventLevel);
		outTree->Fill();
		//if(entry==0) cout<<"EnergyDep " << EnergyDep <<endl;// jiaxuan
		Edep->Fill(EnergyDep);
		if (driftHitNu < 2 && positionHitNu < 2) EdepCut->Fill(EnergyDep);
		HitNo2EnergyDep->Fill(EnergyDep,hitChannelNuEventLevel);
    }
	
	//CrystallBall Fit
	gStyle->SetOptFit(1111);
	double constant = Edep->GetBinContent(Edep->GetMaximumBin());
	double mean = Edep->GetMean();
	double sigma = Edep->GetRMS();
	TF1 * fCry = new TF1("fCry","crystalball");
	TF1 * fGaus = new TF1("fGaus","gaus(0)");
	Edep->Fit(fGaus);
	constant = fGaus->GetParameter(0);
	mean = fGaus->GetParameter(1);
	sigma = fGaus->GetParameter(2);
	cout << constant << "	" << mean << "	" << sigma << endl;
	//fCry->SetParameters(n,mu,sigma,alpha,N);
	fCry->SetParameters(constant,mean,30,1,10000);
	fCry->SetParLimits(0,0.7*constant,1.3*constant);
	fCry->SetParLimits(1,0.7*mean,1.3*mean);
	Edep->Fit(fCry);
    //cout<<"Total Entries : "<<Entries<<" , Layer>=22 : "<<Entries_1<<" , Strips<=64 : "<<Entries_2<<endl;
    
	outFile->cd();
	Edep->SetLineWidth(2);
	Edep->Write();
	EdepCut->SetLineWidth(2);
	EdepCut->Write();
	HitMap->Write();
	HitNo->Write();
	HitLayerNo->Write();
	HitNo2EnergyDep->SetOption("colz");
	HitNo2EnergyDep->GetXaxis()->SetTitle("Energy[MeV]");
	HitNo2EnergyDep->GetYaxis()->SetTitle("HitNo");
	HitNo2EnergyDep->Write();
	LayerEdep->Write();
	
	habHit->Write();
	habHitC->Write();

	hitpro = HitMap->ProfileX();
	hitpro->Write();
	hpro = LayerEdep->ProfileX();
	hpro->Write();
    
	outTree->Write("",TObject::kOverwrite);
    outFile->Close();
    inFile->Close();
    cout<<" [ End of Calibration Process ]  "<<endl;
	return 0;
}


