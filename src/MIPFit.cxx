/********************************************************
 * Author: Yazhou, Jiaxuan
 * Description : MIP File anlysis
 * Function_List ：
 * 		... ... 
 * Temperature correction
 * Last Edit : 2023/10/19
 *******************************************************/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TVirtualPad.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TError.h>
#include "MIPFit.h"
#include "EBUdecode.h"
#include <TKey.h>
#include <TMath.h>
#include <TLine.h>

// first layer (layer 0)
//		---------------------------------------------------------
// 		|		|		|		|		|		|		|		|	...
//		---------------------------------------------------------
// 		|		|		|		|		|		|		|		|	
//		---------------------------------------------------------
// 		|		|		|		|		|		|		|		|	
//		---------------------------------------------------------
// 		|		|		|		|		|		|		|		|	
//		---------------------------------------------------------
// 		|		|		|		|		|		|		|		|	
//		---------------------------------------------------------
// 		|		|		|		|		|		|		|		|	
//		---------------------------------------------------------
//								   ...   

using namespace std;

Double_t langaufun(Double_t* x, Double_t* par) {
	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation),
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.

	// Numeric constants
	Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
	Double_t mpshift = -0.22278298;			 // Landau maximum location

	// Control constants
	Double_t np = 100.0; // number of convolution steps
	Double_t sc = 5.0;	 // convolution extends to +-sc Gaussian sigmas

	// Variables
	Double_t xx;
	Double_t mpc;
	Double_t fland;
	Double_t sum = 0.0;
	Double_t xlow, xupp;
	Double_t step;
	Double_t i;

	// MP shift correction
	mpc = par[1] - mpshift * par[0];

	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];

	step = (xupp - xlow) / np;

	// Convolution integral of Landau and Gaussian by sum
	for (i = 1.0; i <= np / 2; i++) {
		xx = xlow + (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);

		xx = xupp - (i - .5) * step;
		fland = TMath::Landau(xx, mpc, par[0]) / par[0];
		sum += fland * TMath::Gaus(x[0], xx, par[3]);
	}
	return (par[2] * step * sum * invsq2pi / par[3]);

}

void MIPCalibration::ReadTree(int treeType) {
	inTree->SetBranchAddress("Run_Num",&_run_no);
	inTree->SetBranchAddress("Event_Time",&_event_time);
	inTree->SetBranchAddress("CycleID", &_cycleID);
	inTree->SetBranchAddress("TriggerID", &_triggerID);
	inTree->SetBranchAddress("CellID", &_cellID);
	inTree->SetBranchAddress("BCID", &_bcid);
	inTree->SetBranchAddress("HitTag", &_hitTag);
	inTree->SetBranchAddress("GainTag", &_gainTag);
	inTree->SetBranchAddress("HG_Charge", &_chargehi);//for high gain charge
	inTree->SetBranchAddress("LG_Charge", &_chargelo);//for low gain charge t
	inTree->SetBranchAddress("Hit_Time", &_hit_time);//for AutoGain mode time information
	if(treeType==0) inTree->SetBranchAddress("Temperature", &_temp);//ECAL temperature data is with scientific data
}

void MIPCalibration::SetTree() {
	mipTree->Branch("CellID", &_cellID_mip);
	mipTree->Branch("LandauMPV", &_landauMPV);
	mipTree->Branch("LandauWidth", &_landauWidth);
	mipTree->Branch("GausSigma", &_gausSigma);
	mipTree->Branch("ChiSquare", &_ChiSqr);
	mipTree->Branch("NDF", &_Ndf);
	mipTree->Branch("IsGoodChan", &_isGoodChan);
}

void MIPCalibration::DataClear() {
	if(_XPOS.size()) _XPOS.clear();
	if(_YPOS.size()) _YPOS.clear();
	if(_ZPOS.size()) _ZPOS.clear();
	if(_HADC.size()) _HADC.clear();
	if(_Layer.size()) _Layer.clear();
	if(_Chip.size()) _Chip.clear();
	if(_Channel.size()) _Channel.clear();
}

void MIPCalibration::ForceSetBranchValue(int layer, int chip, int channel){
	_cellID_mip = layer * 1e5 + chip * 1e4 + channel;
	_landauWidth = -10;
	_landauMPV = -10;  
	_gausSigma = -10;
	_ChiSqr = -10;
	_Ndf = -10;
	_isGoodChan = false;
}

int MIPCalibration::TreeType(string inFileName)
{
	int treeType = -1;
	TFile* testFile = TFile::Open(TString(inFileName),"READ");
    TIter next(testFile->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)next())) {
		if (strcmp(key->GetClassName(), "TTree") == 0) { 
			inTree = (TTree*)testFile->Get(key->GetName()); 
			cout << "Tree name: " << inTree->GetName() << endl;
			if (strcmp(inTree->GetName(), "Raw_Hit") == 0)  {
				treeType = 0; // data 
			} 
			if (strcmp(inTree->GetName(), "Digi_Raw_Hit") == 0)  {
				treeType = 1; // simu_digi
			}
		}
	}
	cout << "============ Function TreeType Finished ===============" << endl;
	return treeType;
}

/*void MIPCalibration::MuonTrackFit(int entry) { //entry just for checking
	gErrorIgnoreLevel = kFatal;
	gStyle->SetOptFit(1111);
	const double zerr = 1; // mm
	const double xerr = 2.5; // mm
	const double yerr = 2.5; // mm
	
	TF1 * lineXZ = new TF1("f1","pol1",0,300);
	TF1 * lineYZ = new TF1("f1","pol1",0,300);
	TGraphErrors* posXZ = new TGraphErrors();
	TGraphErrors* posYZ = new TGraphErrors();
	
	int i_XZ = 0;
	int i_YZ = 0;
	double averX = 0;
	double averY = 0;
	int allHitInOneEvent = _ZPOS.size();
	for(size_t i_pos = 0; i_pos < _ZPOS.size(); i_pos++)
	{
		if ( _Layer.at(i_pos)%2 == 1 ) { 
			posXZ->SetPoint(i_XZ, _ZPOS.at(i_pos), _XPOS.at(i_pos));
			posXZ->SetPointError(i_XZ, zerr, xerr);
			averX += _XPOS.at(i_pos); 
			i_XZ++;
		} else {
			posYZ->SetPoint(i_YZ, _ZPOS.at(i_pos), _YPOS.at(i_pos));
			posYZ->SetPointError(i_YZ, zerr, yerr);
			averY += _YPOS.at(i_pos);
			i_YZ++;
		}
	}
	averX /= i_XZ; averY /= i_YZ;
	//cout << "i_XZ : i_YZ " << i_XZ << " " << i_YZ << endl; 
	
	delete posYZ;
	delete posXZ;
	posXZ = new TGraphErrors();
	posYZ = new TGraphErrors();

	//////////////////////////////////////////////////////
	// 						fit first time 
	//////////////////////////////////////////////////////
	int i_Raw = 0;
	int i_XZ_Raw = 0;
	int i_YZ_Raw = 0;
	int cycle = 0 ;
	double firstThre = 40;
	while (cycle < allHitInOneEvent) 
	{
		cycle++;
		double residualXZ = _XPOS.at(i_Raw) - averX;
		double residualYZ = _YPOS.at(i_Raw) - averY;
		if ( _Layer.at(i_Raw)%2 == 1 ) {
			//cout << "layer : residualXZ " << _Layer.at(i_Raw) << residualXZ << endl;
			if (fabs(residualXZ) < firstThre) {
				posXZ->SetPoint(i_XZ_Raw, _ZPOS.at(i_Raw), _XPOS.at(i_Raw));
				posXZ->SetPointError(i_XZ_Raw, zerr, xerr);
				i_XZ_Raw++;
			} else {
				//cout << "====>Deleting ... " << fixed << setprecision(1) << " " <<_Layer.at(i_Raw) << " " <<_ZPOS.at(i_Raw) << " " <<_XPOS.at(i_Raw) <<  " " <<_YPOS.at(i_Raw) << endl; 
				_XPOS.erase(_XPOS.begin() + i_Raw);
				_YPOS.erase(_YPOS.begin() + i_Raw);
				_ZPOS.erase(_ZPOS.begin() + i_Raw);
				_HADC.erase(_HADC.begin() + i_Raw);
				_Layer.erase(_Layer.begin() + i_Raw);
				_Chip.erase(_Chip.begin() + i_Raw);
				_Channel.erase(_Channel.begin() + i_Raw);
				i_Raw--;
			}
		} else {
			//cout << "layer : residualYZ " << _Layer.at(i_Raw) << " " <<  residualYZ << endl;
			if (fabs(residualYZ) < firstThre) {
				posYZ->SetPoint(i_YZ_Raw, _ZPOS.at(i_Raw), _YPOS.at(i_Raw));
				posYZ->SetPointError(i_YZ_Raw, zerr, yerr);
				i_YZ_Raw++;
			} else {
				//cout << "====>Deleting ... " << fixed << setprecision(1) << " " <<_Layer.at(i_Raw) << " " <<_ZPOS.at(i_Raw) << " " <<_XPOS.at(i_Raw) <<  " " <<_YPOS.at(i_Raw) << endl; 
				_XPOS.erase(_XPOS.begin() + i_Raw);
				_YPOS.erase(_YPOS.begin() + i_Raw);
				_ZPOS.erase(_ZPOS.begin() + i_Raw);
				_HADC.erase(_HADC.begin() + i_Raw);
				_Layer.erase(_Layer.begin() + i_Raw);
				_Chip.erase(_Chip.begin() + i_Raw);
				_Channel.erase(_Channel.begin() + i_Raw);
				i_Raw--;	
			}
		}
		i_Raw++;
	}
	
	//cout << "i_Raw : i_XZ_Raw : i_YZ_Raw " << i_Raw << " " << i_XZ_Raw << " " << i_YZ_Raw << endl; 
	lineXZ->SetParameter(0,averX);
	lineXZ->SetParameter(1,0);
	lineXZ->SetParLimits(0,averX-15, averX+15);
	lineXZ->SetParLimits(1, -0.1, 0.1);
	posXZ->Fit(lineXZ,"Q");
	double interceptXZ = lineXZ->GetParameter(0);
	double slopeXZ = lineXZ->GetParameter(1);
	lineYZ->SetParameter(0,averY);
	lineYZ->SetParameter(1,0);
	lineYZ->SetParLimits(0,averY-15, averY+15);
	lineYZ->SetParLimits(1, -0.1, 0.1);
	posYZ->Fit(lineYZ,"Q");
	double interceptYZ = lineYZ->GetParameter(0);
	double slopeYZ = lineYZ->GetParameter(1);


	//////////////////////////////////////////////////////
	// 					fit second time
	//////////////////////////////////////////////////////
	int i_Track = 0;
	int i_XZ_Track = 0;
	int i_YZ_Track = 0;
	//cout << "===============================================" << endl;
	//cout << "_LAYER : _ZPOS : _XPOS : resXZ : _YPOS : resYZ		" << endl;
	cycle = 0 ;
	double secondThre = 7.5;
	while (cycle < i_Raw) {
		cycle++;
		double fitValueXZ = interceptXZ + slopeXZ * _ZPOS.at(i_Track);
		double residualXZ = fitValueXZ - _XPOS.at(i_Track);
		double fitValueYZ = interceptYZ + slopeYZ * _ZPOS.at(i_Track);
		double residualYZ = fitValueYZ - _YPOS.at(i_Track);
		//cout << _Layer.at(i_Track);
		//cout << fixed << setprecision(1) << "	" << _ZPOS.at(i_Track) ;
		//cout << fixed << setprecision(1) << "	" << _XPOS.at(i_Track) ;
		//cout << fixed << setprecision(1) << "	" << residualXZ ;
		//cout << fixed << setprecision(1) << "	" << _YPOS.at(i_Track) ;
		//cout << fixed << setprecision(1) << "	" << residualYZ << endl ;
		int layer = _Layer.at(i_Track);
		if ( layer%2 == 1 && fabs(residualXZ) > secondThre ) { 
			//cout << "==== abnormal ====" << endl;
			posXZ->RemovePoint(i_XZ_Track);
			i_XZ_Track--;
			// cout << "Deleting ... " << fixed << setprecision(1) << " " <<_Layer.at(i_Track) << " " <<_ZPOS.at(i_Track) << " " <<_XPOS.at(i_Track) <<  " " <<_YPOS.at(i_Track) << endl; 
			_XPOS.erase(_XPOS.begin() + i_Track);
			_YPOS.erase(_YPOS.begin() + i_Track);
			_ZPOS.erase(_ZPOS.begin() + i_Track);
			// for (int i = 0; i< _ZPOS.size(); i++){cout << _ZPOS.at(i) <<"	";if(i%5==4) cout << endl;}cout<<endl;
			_HADC.erase(_HADC.begin() + i_Track);
			_Layer.erase(_Layer.begin() + i_Track);
			_Chip.erase(_Chip.begin() + i_Track);
			_Channel.erase(_Channel.begin() + i_Track);
			i_Track--;
		} 
		if ( layer%2 == 0 && fabs(residualYZ) > secondThre ) {
			//cout << "==== abnormal ====" << endl;
			posYZ->RemovePoint(i_YZ_Track);
			i_YZ_Track--;
			//cout << "Deleting ... " << fixed << setprecision(1) << " " <<_Layer.at(i_Track) << " " <<_ZPOS.at(i_Track) << " " <<_XPOS.at(i_Track) <<  " " <<_YPOS.at(i_Track) << endl; 
			_XPOS.erase(_XPOS.begin() + i_Track);
			_YPOS.erase(_YPOS.begin() + i_Track);
			_ZPOS.erase(_ZPOS.begin() + i_Track);
			// for (int i = 0; i< _ZPOS.size(); i++){cout << _ZPOS.at(i) <<"	";if(i%5==4) cout << endl;}cout<<endl;
			_HADC.erase(_HADC.begin() + i_Track);
			_Layer.erase(_Layer.begin() + i_Track);
			_Chip.erase(_Chip.begin() + i_Track);
			_Channel.erase(_Channel.begin() + i_Track);
			i_Track--;	
		}
		//if(_Layer.at(i_Track)==0&&_Channel.at(i_Track)==0) cout << "======== Now Entry : " << entry << " =========" << endl;
		if ( layer%2 == 1 ) { 
			i_XZ_Track++;
		} else {
			i_YZ_Track++;
		}
		i_Track++;
	}
	//for (size_t i = 0; i< _ZPOS.size(); i++){cout << _ZPOS.at(i) <<"		";if(i%5==4) cout << endl;}cout << endl;cout<<"============"<<endl;
	//cout << "i_XZ_Track : i_YZ_Track : i_Track " << i_XZ_Track << " " << i_YZ_Track << " " << i_Track << endl; 
	
	lineXZ->SetParameter(0,interceptXZ);
	lineXZ->SetParameter(1,0);
	lineXZ->SetParLimits(0,interceptXZ-15, interceptXZ+15);
	lineXZ->SetParLimits(1, -0.1, 0.1);
	posXZ->Fit(lineXZ,"Q");
	interceptXZ = lineXZ->GetParameter(0);
	slopeXZ = lineXZ->GetParameter(1);
	lineYZ->SetParameter(0,interceptYZ);
	lineYZ->SetParameter(1,0);
	lineYZ->SetParLimits(0,interceptYZ-15, interceptYZ+15);
	lineYZ->SetParLimits(1, -0.1, 0.1);
	posYZ->Fit(lineYZ,"Q");
	interceptYZ = lineYZ->GetParameter(0);
	slopeYZ = lineYZ->GetParameter(1);
	
	//////////////////////////////////////////////////////
	// 				fit third time
	//////////////////////////////////////////////////////
	int i_Valid = 0;
	int i_XZ_Valid = 0;
	int i_YZ_Valid = 0;
	//cout << "===============================================" << endl;
	//cout << "_LAYER : _ZPOS : _XPOS : resXZ : _YPOS : resYZ		" << endl;
	int allHitInOneEvent_Track = _ZPOS.size();
	cycle = 0 ;
	double thirdThre = 3;
	while (cycle < allHitInOneEvent_Track) {
		cycle++;
		double fitValueXZ = interceptXZ + slopeXZ * _ZPOS.at(i_Valid);
		double residualXZ = fitValueXZ - _XPOS.at(i_Valid);
		double fitValueYZ = interceptYZ + slopeYZ * _ZPOS.at(i_Valid);
		double residualYZ = fitValueYZ - _YPOS.at(i_Valid);
		//cout << _Layer.at(i_Valid);
		//cout << fixed << setprecision(1) << "	" << _ZPOS.at(i_Valid) ;
		//cout << fixed << setprecision(1) << "	" << _XPOS.at(i_Valid) ;
		//cout << fixed << setprecision(1) << "	" << residualXZ ;
		//cout << fixed << setprecision(1) << "	" << _YPOS.at(i_Valid) ;
		//cout << fixed << setprecision(1) << "	" << residualYZ << endl ;
		int layer = _Layer.at(i_Valid);
		if ( layer%2 == 1 && fabs(residualXZ) > thirdThre ) { 
			//cout << "==== abnormal ====" << endl;
			posXZ->RemovePoint(i_XZ_Valid);
			i_XZ_Valid--;
			//cout << "Deleting ... " << fixed << setprecision(1) << " " <<_Layer.at(i_Valid) << " " <<_ZPOS.at(i_Valid) << " " <<_XPOS.at(i_Valid) <<  " " <<_YPOS.at(i_Valid) << endl; 
			_XPOS.erase(_XPOS.begin() + i_Valid);
			_YPOS.erase(_YPOS.begin() + i_Valid);
			_ZPOS.erase(_ZPOS.begin() + i_Valid);
			//for(int i=0; i<_ZPOS.size(); i++){cout <<_ZPOS.at(i)<<" ";if(i%5==4) cout<<endl;}cout<<endl;
			_HADC.erase(_HADC.begin() + i_Valid);
			_Layer.erase(_Layer.begin() + i_Valid);
			_Chip.erase(_Chip.begin() + i_Valid);
			_Channel.erase(_Channel.begin() + i_Valid);
			i_Valid--;
		} 
		if ( layer%2 == 0 && fabs(residualYZ) > thirdThre ) {
			//cout << "==== abnormal ====" << endl;
			posYZ->RemovePoint(i_YZ_Valid);
			i_YZ_Valid--;
			//cout << "Deleting ... " << fixed << setprecision(1) << " " <<_Layer.at(i_Valid) << " " <<_ZPOS.at(i_Valid) << " " <<_XPOS.at(i_Valid) <<  " " <<_YPOS.at(i_Valid) << endl; 
			_XPOS.erase(_XPOS.begin() + i_Valid);
			_YPOS.erase(_YPOS.begin() + i_Valid);
			_ZPOS.erase(_ZPOS.begin() + i_Valid);
			_HADC.erase(_HADC.begin() + i_Valid);
			_Layer.erase(_Layer.begin() + i_Valid);
			_Chip.erase(_Chip.begin() + i_Valid);
			_Channel.erase(_Channel.begin() + i_Valid);
			i_Valid--;
		}
		//if(_Layer.at(i_Valid)==0&&_Channel.at(i_Valid)==0) cout << "======== Now Entry : " << entry << " =========" << endl;
		
		if ( layer%2 == 1 ) { 
			i_XZ_Valid++;
		} else {
			i_YZ_Valid++;
		}
		i_Valid++;
	}

	allHit += allHitInOneEvent;
	validHit += i_Valid;
	validHit2All = i_Valid / allHitInOneEvent;
	//cout << "i_XZ_Valid: i_YZ_Valid: i_Valid " << i_XZ_Valid << " " << i_YZ_Valid << " " << i_Valid << endl; 
	//cout << "===============================================" << endl;
	//if(entry%100000==0)
	
	lineXZ->SetParameter(0,interceptXZ);
	lineXZ->SetParameter(1,0);
	lineXZ->SetParLimits(0,interceptXZ-10, interceptXZ+10);
	lineXZ->SetParLimits(1, -0.1, 0.1);
	posXZ->Fit(lineXZ,"Q");
	interceptXZ = lineXZ->GetParameter(0);
	slopeXZ = lineXZ->GetParameter(1);
	//cout << "===== XZ FIT FINISH ! ======" << endl;
	
	lineYZ->SetParameter(0,interceptYZ);
	lineYZ->SetParameter(1,0);
	lineYZ->SetParLimits(0,interceptYZ-10, interceptYZ+10);
	lineYZ->SetParLimits(1, -0.1, 0.1);
	posYZ->Fit(lineYZ,"Q");
	interceptYZ = lineYZ->GetParameter(0);
	slopeYZ = lineYZ->GetParameter(1);
	//cout << "===== YZ FIT FINISH ! ======" << endl;

	
	TCanvas *cXZ = new TCanvas(Form("cXZ_odd_%d",entry), Form("cXZ_odd_%d",entry) ,800,600);
	cXZ->cd();
	gPad->SetGridx();
	posXZ->SetTitle(Form("cXZ_odd_entry%d",entry));
	posXZ->GetXaxis()->SetRangeUser(0, 350);
	posXZ->GetYaxis()->SetRangeUser(-150, 150);
	posXZ->GetXaxis()->SetTitle("Z Direction (mm)");
	posXZ->GetYaxis()->SetTitle("X Direction (mm)");
	posXZ->GetYaxis()->SetTitleOffset(1);
	posXZ->Draw("AP");
	cXZ->Write();

	TCanvas *cYZ = new TCanvas(Form("cYZ_even_%d",entry), Form("cYZ_even_%d",entry),800,600);
	cYZ->cd();
	gPad->SetGridx();
	posYZ->SetTitle(Form("cYZ_even_entry%d",entry));
	posYZ->GetXaxis()->SetRangeUser(0, 350);
	posYZ->GetYaxis()->SetRangeUser(-150, 150);
	posYZ->GetXaxis()->SetTitle("Z Direction (mm)");
	posYZ->GetYaxis()->SetTitle("Y Direction (mm)");
	posYZ->GetYaxis()->SetTitleOffset(1);
	posYZ->Draw("AP");
	cYZ->Write();
	
	delete cXZ;
	delete cYZ;

	delete lineYZ;
	delete lineXZ;
	delete posYZ;
	delete posXZ;
}*/

// 1) 根据和 (averX, averY) 的距离来删除 (第一次)
void MIPCalibration::removeByAverage(std::vector<Hit>& hits, double threshold, double averX, double averY)
{
    hits.erase(std::remove_if(hits.begin(), hits.end(), 
                              [&](const Hit& h){
                                  if (h.layer % 2 == 1) {
                                      return (std::fabs(h.x - averX) > threshold);
                                  } else {
                                      return (std::fabs(h.y - averY) > threshold);
                                  }
                              }),
               hits.end());
}

// 2) 根据和线性拟合的残差来删除 (第二、第三次)
void MIPCalibration::removeByFitResidual(
    std::vector<Hit>& hits, double threshold,
    double interceptXZ, double slopeXZ,
    double interceptYZ, double slopeYZ )
{
    hits.erase(std::remove_if(hits.begin(), hits.end(), 
                              [&](const Hit& h){
                                  if (h.layer % 2 == 1) {
                                      // XZ平面: 拟合函数 f(z) = interceptXZ + slopeXZ * z
                                      double fitVal = interceptXZ + slopeXZ * h.z;
                                      return (std::fabs(fitVal - h.x) > threshold);
                                  } else {
                                      // YZ平面
                                      double fitVal = interceptYZ + slopeYZ * h.z;
                                      return (std::fabs(fitVal - h.y) > threshold);
                                  }
                              }),
               hits.end());
}

FitResult MIPCalibration::fitXZandYZ(const std::vector<Hit>& hits)
{
	gErrorIgnoreLevel = kFatal;
	gStyle->SetOptFit(1111);
    static const double zerr = 1;   // mm
    static const double xerr = 2.5; // mm
    static const double yerr = 2.5; // mm

    // 新建 TGraphErrors
    auto posXZ = std::make_unique<TGraphErrors>();
    auto posYZ = std::make_unique<TGraphErrors>();

    // 往里塞点
    int iXZ = 0, iYZ = 0;
    for (auto &h : hits) {
        if (h.layer % 2 == 1) {
            posXZ->SetPoint(iXZ, h.z, h.x);
            posXZ->SetPointError(iXZ, zerr, xerr);
            iXZ++;
        } else {
            posYZ->SetPoint(iYZ, h.z, h.y);
            posYZ->SetPointError(iYZ, zerr, yerr);
            iYZ++;
        }
    }

    // 准备 TF1
    TF1 lineXZ("lineXZ", "pol1", 0, 300);
    TF1 lineYZ("lineYZ", "pol1", 0, 300);

    // 简单先给一些初值 / 范围(也可以更精细)
    lineXZ.SetParameter(0, 0); // intercept
    lineXZ.SetParameter(1, 0); // slope
    lineXZ.SetParLimits(1, -0.5, 0.5); // example

    lineYZ.SetParameter(0, 0);
    lineYZ.SetParameter(1, 0);
    lineYZ.SetParLimits(1, -0.5, 0.5);

    // 拟合
    posXZ->Fit(&lineXZ, "Q"); 
    double interceptXZ = lineXZ.GetParameter(0);
    double slopeXZ     = lineXZ.GetParameter(1);

    posYZ->Fit(&lineYZ, "Q");
    double interceptYZ = lineYZ.GetParameter(0);
    double slopeYZ     = lineYZ.GetParameter(1);

    FitResult result { interceptXZ, slopeXZ, interceptYZ, slopeYZ };
    return result;
}

void MIPCalibration::DrawHits(const std::vector<Hit>& hits, int entry, int times) {
    TCanvas* c = new TCanvas(Form("c_entry%d_times%d",entry,times), Form("c_entry%d_times%d",entry,times), 1800, 900);
    c->Divide(2, 1);
	TGraphErrors* posXZ = new TGraphErrors();
	TGraphErrors* posYZ = new TGraphErrors();
	int i_XZ = 0;
	int i_YZ = 0;
	double xerr = 2.5; // mm
	double yerr = 2.5; // mm
	double zerr = 1;   // mm
	for (auto &h: hits) {
		if (h.layer%2 == 1) { 
			posXZ->SetPoint(i_XZ, h.z, h.x);
			posXZ->SetPointError(i_XZ, zerr, xerr); 
			i_XZ++;
		} else {
			posYZ->SetPoint(i_YZ, h.z, h.y);
			posYZ->SetPointError(i_YZ, zerr, yerr);
			i_YZ++;
		}
	}
	TLine * lineV[6];
	for (int i = 0; i < 6; i++) {
        lineV[i] = new TLine(45 + i*50, -120, 45 + i*50, 120);
        lineV[i]->SetLineColor(kBlack);
        lineV[i]->SetLineStyle(2);
    }
    c->cd(1);
	gPad->SetMargin(0.12,0.12,0.12,0.12);
    posXZ->SetTitle("XZ-plane;Z (mm);X (mm)");
    posXZ->SetMarkerStyle(20); // 圆点
    posXZ->SetMarkerColor(kBlue);
	posXZ->GetYaxis()->SetTitleOffset(1.2);
	posXZ->GetYaxis()->SetTitleSize(0.045);
	posXZ->GetXaxis()->SetTitleSize(0.045);
	posXZ->GetXaxis()->SetLabelSize(0.045);
	posXZ->GetYaxis()->SetLabelSize(0.045);
    posXZ->Draw("AP"); // A = Axis, P = Point
    // 根据需要可调范围
    // posXZ->GetXaxis()->SetRangeUser(0, 300); 
    posXZ->GetYaxis()->SetRangeUser(-120, 120);
    for (int i = 0; i < 6; i++) {lineV[i]->Draw();}

    // ===== 绘制 YZ 平面 =====
    c->cd(2);
	gPad->SetMargin(0.12,0.12,0.12,0.12);
    posYZ->SetTitle("YZ-plane;Z (mm);Y (mm)");
    posYZ->SetMarkerStyle(20);
    posYZ->SetMarkerColor(kRed);
	posYZ->GetYaxis()->SetTitleOffset(1.2);
	posYZ->GetYaxis()->SetTitleSize(0.045);
	posYZ->GetXaxis()->SetTitleSize(0.045);
	posYZ->GetXaxis()->SetLabelSize(0.045);
	posYZ->GetYaxis()->SetLabelSize(0.045);
    posYZ->Draw("AP");

    // 同理可调范围
    // posYZ->GetXaxis()->SetRangeUser(0, 300);
    posYZ->GetYaxis()->SetRangeUser(-120, 120);
    for (int i = 0; i < 6; i++) {lineV[i]->Draw();}

    // 刷新画布
    c->Update();
	c->SaveAs(Form("drawhits/c_entry%d_times%d.jpg",entry,times));
}

void MIPCalibration::MuonTrackFit(int entry) {
	// 假设本函数获取到所有 hit，放进 hits 这个 vector 中
	std::vector<Hit> hits;
	hits.reserve(_ZPOS.size());
	int allHitInOneEvent = _ZPOS.size();
	for (size_t i = 0; i < _ZPOS.size(); i++) {
		Hit h;
		h.x       = _XPOS[i];
		h.y       = _YPOS[i];
		h.z       = _ZPOS[i];
		h.hadc    = _HADC[i];
		h.layer   = _Layer[i];
		h.chip    = _Chip[i];
		h.channel = _Channel[i];
		hits.push_back(h);
	}
	

	// 计算 averX / averY（分别统计奇数层、偶数层）
	double sumX = 0.0, sumY = 0.0;
	int    countX = 0,   countY = 0;
	for (auto &h : hits) {
		if (h.layer % 2 == 1) { // odd
			sumX += h.x;
			countX++;
		} else {                // even
			sumY += h.y;
			countY++;
		}
	}
	double averX = (countX > 0) ? sumX / countX : 0;
	double averY = (countY > 0) ? sumY / countY : 0;
	// 2) 第一次筛选: 与平均值差大于 firstThre 就删除
    double firstThre = 40.0;
    removeByAverage(hits, firstThre, averX, averY);
    if (entry == 30) DrawHits(hits,entry,0);
    // 第一次拟合
    auto r1 = fitXZandYZ(hits);

    // 3) 第二次筛选: 用 r1 的拟合结果，残差阈值 secondThre
    double secondThre = 7.5;
    removeByFitResidual(hits, secondThre,
                        r1.interceptXZ, r1.slopeXZ,
                        r1.interceptYZ, r1.slopeYZ);

    // 第二次拟合
    auto r2 = fitXZandYZ(hits);

    // 4) 第三次筛选: 残差阈值 thirdThre
    double thirdThre = 3.0;
    removeByFitResidual(hits, thirdThre,
                        r2.interceptXZ, r2.slopeXZ,
                        r2.interceptYZ, r2.slopeYZ);

    // 第三次拟合
    //auto r3 = fitXZandYZ(hits);
	if (entry == 10) DrawHits(hits,entry,3);

    // 拟合完成后的最终参数
    //double interceptXZ = r3.interceptXZ;
    //double slopeXZ     = r3.slopeXZ;
    //double interceptYZ = r3.interceptYZ;
    //double slopeYZ     = r3.slopeYZ;

    // 统计保留下来的点数
    int validHitCount = hits.size();
    validHit2All = double(validHitCount) / double(allHitInOneEvent);
    allHit  += allHitInOneEvent;
    validHit += validHitCount;
	
	DataClear();
	for (auto &h : hits) {
		_XPOS.push_back(h.x);
		_YPOS.push_back(h.y);
		_ZPOS.push_back(h.z);
		_HADC.push_back(h.hadc);
		_Layer.push_back(h.layer);
		_Chip.push_back(h.chip);
		_Channel.push_back(h.channel);
	}
    // 如果需要，可以在这里把“保留下来的 hits”再写回 _XPOS, _YPOS, ... 
    // 或者把最终的 intercept/slope 放到某些成员变量里
    
    // Debug 输出
    // std::cout << "Entry " << entry
    //           << ": allHitInOneEvent=" << allHitInOneEvent 
    //           << ", validHitCount=" << validHitCount 
    //           << ", ratio=" << validHit2All << std::endl;
}

int MIPCalibration::MIPFit(string input_name, Extract ex, string hist_name, string mip_name, bool isTrackFit) 
{
	//PedestalExtract(pedfactor_name.c_str());
	// name operation
	std::cout << "********** MIPFit now **********" << std::endl;
	ofstream mipdeadtxt("DeadList.txt");
	string txtname(input_name);
	txtname.erase(0, txtname.find_last_of("/") + 1);
	string histname(hist_name);
	histname.insert(histname.length(), txtname);
	string mipname(mip_name);
	mipname.insert(mipname.length(), txtname);
	cout << "===========================" << endl;
	// *********************************************************************************
	cout << "MIP fit " << input_name << " ..." << endl;
	
	int treeType = -1;
	treeType = TreeType(input_name);	
    inFile = TFile::Open(TString(input_name),"READ");
	if (treeType == 0) {
		inTree = (TTree*)inFile->Get("Raw_Hit");
	} else if (treeType == 1) {
		inTree = (TTree*)inFile->Get("Digi_Raw_Hit");
		//inTree = (TTree*)inFile->Get("MC_Truth");
	} else {
		return 1;
	}
	ReadTree(treeType);
	mipFile = new TFile(mipname.c_str(), "RECREATE");
	mipTree = new TTree("MIP_Fit", "MIP_Calibration");
	SetTree();
	
	char char_tmp[100];
	TH1D* h_MIP[layerNu][chipNu][channelNu];
	
	// Set Format
	for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
				sprintf(char_tmp, "Layer%d_Chip%d_Chn%d", i_layer, i_chip, i_channel);
				if (i_layer > 3 && i_layer < 28) {
					h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 1000, -300, 3700);
				} else {
					h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 400, -300, 3700);
					//h_MIP[i_layer][i_chip][i_channel] = new TH1D(char_tmp, char_tmp, 800, -200, 3800);
				}
				h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetTitle("ADC");
				h_MIP[i_layer][i_chip][i_channel]->SetLineWidth(2);
				h_MIP[i_layer][i_chip][i_channel]->SetLineColor(4);
			}
		}
	}

	// Read Data
	int totalEntries = inTree->GetEntries();
	allHit = 0;
	validHit = 0;

	bool isLayerHit[layerNu];
	int layerHit[layerNu] = {0};
	int totalEffEvent = 0; 
	double hitEfficiency[layerNu] = {.0};

	if (isTrackFit) {
		cout << "========== TRACK FIT !!! ===========" << endl;
	} else {
		cout << "========== NO trackfit !!! ===========" << endl;
	}

	for (int entry = 0; entry < totalEntries; ++entry) {
	//for (int entry = 0; entry < 10; ++entry) {
		//cout << "======== Now Entry : " << entry << " =========" << endl;
		//if(entry%3==2)  continue; // jiaxuan : sample
		if (entry%100000==0) {	
			cout << " [ Process : MIP Fit : ]  event : " << entry << endl;
		}
		inTree->GetEntry(entry);
		// Hit level data in one event
		DataClear();
		int whichLayer[layerNu] = {0};
		hitLayerNo = 0;
		hitCellNo = 0;
		for (size_t i_hit = 0; i_hit != _cellID->size(); ++i_hit) {
			if(_hitTag->at(i_hit)==0)	continue;
			if(_gainTag->at(i_hit)==0)	continue; // for auto gain muon data
			int layerID = _cellID->at(i_hit) / 1e5;
			if( layerID>=layerNu || layerID<0) continue;
			int chipID = (_cellID->at(i_hit) % 100000) / 1e4;
			int channelID = _cellID->at(i_hit) % 100;
			whichLayer[layerID] = 1;
			if (entry%100000==0) cout << _cellID->at(i_hit) << " ";
			double *_position = EBUdecode(layerID,chipID,channelID);
			_XPOS.push_back(*_position);
			_YPOS.push_back(*(_position+1));
			_ZPOS.push_back(*(_position+2));
			_HADC.push_back(_chargehi->at(i_hit));
			_Layer.push_back(layerID);
			_Chip.push_back(chipID);
			_Channel.push_back(channelID);

			for(int i_layer = 0; i_layer < layerNu; i_layer++) { hitLayerNo += whichLayer[i_layer]; }
			hitCellNo++;
		}
		if (entry%100000==0) cout << endl;
		
		if (hitCellNo > 60 || hitLayerNo < 20) continue;
		//jiaxuan : if note , then all
		if (isTrackFit) {
			MuonTrackFit(entry);
		}
		for (size_t i_hit = 0; i_hit != _XPOS.size(); ++i_hit) {
			int layerID = _Layer.at(i_hit);
			int chipID = _Chip.at(i_hit);
			int channelID = _Channel.at(i_hit);

			//cout<<"i_hit "<<i_hit;
			//cout<<" layer : chip : channel "<<layerID<<" "<<chipID<<" "<<channelID<<endl; 
			//h_MIP[layerID][chipID][channelID]->Fill(_HADC.at(i_hit)-pedMeanTimesChnLevel[layerID][chipID][channelID]);
			double realADC = _HADC.at(i_hit)-ex.pedMeanHigh[layerID][chipID][channelID];
			if (treeType==0) {
				h_MIP[layerID][chipID][channelID]->Fill(realADC);
			} else if (treeType==1) {
				//if (realADC > 0.8 * ex.MIPPeakValue[layerID][chipID][channelID]) h_MIP[layerID][chipID][channelID]->Fill(realADC);
				h_MIP[layerID][chipID][channelID]->Fill(realADC);
			}
				
		}
		
		fill(isLayerHit,isLayerHit+layerNu,false);
		for (size_t i_hit = 0; i_hit != _XPOS.size(); ++i_hit) {
			int layerID = _Layer.at(i_hit);
			isLayerHit[layerID] = true;
		}

		for (int i_layer = 0; i_layer < layerNu; i_layer++) {
			if (entry%100000==0) cout << isLayerHit[i_layer] << " ";
			if (isLayerHit[i_layer])	layerHit[i_layer]++;
		}
		if (entry%100000==0) cout << endl;
		
		totalEffEvent++;

		//for (size_t i_hit = 0; i_hit != _cellID->size(); ++i_hit) {
			// double SiPMtemp;
            // if(_temp->at(layerID).size()!=0) 
			// 	SiPMtemp = tempReconstruction(layerID, _position, _temp->at(layerID));
			// else 
			// 	SiPMtemp = 23;
			// incident angle correction
			// time = 2. / _realLength * time;

			// jiaxuan: temperature dependence correction
			//time = (1+(SiPMtemp-23)*0.03)*time;
			
			//if (layerID < 3 || layerID > 28)
			//	time = (1 + (_temp_mip->at(i_hit) - 23) * 0.0122) * time;
			//else
			// 	time = (1 + (_temp_mip->at(i_hit) - 23) * 0.0073) * time;
			
			// MIP response decline correction
			//time = time/(1-0.17*(double)entry/(double)totalEntries);
			//int period1 = 20*13000;
			//int period2 = 75*13000;
			//if(entry<=period1) time = time/(1-0.1*(double)entry/(double)period1); // MIP response decline correction
			//else if(entry<=period2) time = time/(1-0.077*(double)(entry-period1)/(double)(period2-period1))*14./12.8; // MIP response decline correction
			//else  time = time/(1-0.017*(double)(entry-period2)/(double)(totalEntries-period2))*14/11.9; // MIP response decline correction
			//h_MIP[layerID][chipID][channelID]->Fill(time);
			//h_MIP_Chip[layerID][chipID]->Fill(time);
		//}
	}
	
	cout << "======== allHit : validHit  " << allHit << "  " << validHit << " ========== " << endl;	
	TH1D* hEfficiency = new TH1D("hEfficiency", "", 32, 0, 32);
	TH1D* hDeadCountLayer = new TH1D("hDeadCountLayer", "", 32, 0, 32);
	TH2D* hCountChn = new TH2D("hCountChn","",192,0,192,36,0,36);
	TH2D* hMIPChn = new TH2D("hMIPChn","",192,0,192,36,0,36);
	
	TFile* histfile = new TFile(histname.c_str(), "RECREATE");
	TCanvas* C_MIP[layerNu][chipNu];
	//TCanvas* C_MIP_chip[layerNu];
	
	string pdfname(txtname);
	pdfname.erase(pdfname.end()-5,pdfname.end());
	pdfname.append(".pdf");
	histname.erase(histname.begin()+histname.find("hist")-1,histname.end());
	pdfname = histname + "/pdf/" + pdfname; 
	string tmpp;
	cout <<  "================>>"<< pdfname << endl;
	
	//for Channel Level MIP
	TCanvas *Chn_MIP[2000];
	//TCanvas *Chn_MIP_Abnormal;
	int page_no = 0;
	int deadChannelNo = 0;
	cout << "Total Effective Event : " << totalEffEvent << endl;
	hEfficiency->SetOption("HIST");
	int MPVAbNu = 0;
	int SigmaAbNu = 0;
	int LackNu = 0;
	for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
	//for (int i_layer = 9; i_layer < 0; ++i_layer) {
		hitEfficiency[i_layer] = (layerHit[i_layer]+.0) / totalEffEvent;
		hEfficiency->Fill(i_layer,hitEfficiency[i_layer]);
		cout << " [ LandauGaus Fitting : ] Layer : " << i_layer << "	Hit : " << layerHit[i_layer] <<" Efficiency : " << hitEfficiency[i_layer] << endl;
		// jiaxuan : note for test
		for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
		//for (int i_chip = 2; i_chip < 3; ++i_chip) {
			histfile->cd();
			sprintf(char_tmp, "layer%d/chip%d", i_layer, i_chip);
			histfile->mkdir(char_tmp);
			histfile->cd(char_tmp);
			int chn_entries = 0;
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
				if (i_chip == 5 && i_channel > 29)	continue;
				if (i_channel%9 == 0) {
					page_no++;
					Chn_MIP[page_no] = new TCanvas(Form("page%d",page_no),Form("page%d",page_no),1);
					Chn_MIP[page_no]->Divide(3,3);
				}
				
				int chipID_mip = i_layer * 6 + i_chip ;
				chn_entries= h_MIP[i_layer][i_chip][i_channel]->GetEntries();
				if (chn_entries<10 && i_layer < layerNu) { 
					cout << i_layer*100000+i_chip*10000+i_channel <<", " ;
				}
				// jiaxaun : Smooth h_MIP
				//h_MIP[i_layer][i_chip][i_channel]->Smooth(3,"k5b");
				
				hCountChn->Fill(chipID_mip,i_channel,chn_entries);
				Chn_MIP[page_no]->cd(i_channel%9+1);
				gPad->SetGrid();
				
				TF1* fun = new TF1("langau", langaufun, -250, 1600, 4);
				fun->SetParName(0, "LandauWidth");
				fun->SetParName(1, "LandauMPV");
				fun->SetParName(2, "Area");
				fun->SetParName(3, "GausSigma");
			
				int SiPMType = -1;
				if (i_layer>=4 && i_layer <=27) {
					SiPMType = 10;
				} else if ((i_layer>=0 &&i_layer<=3) || (i_layer>=28 && i_layer<=31)) {
					SiPMType = 15;
				}
				double fitLeftThreshold = 0.6;
				if (SiPMType == 10) {
					fun->SetParameters(20, 120, 10000, 30);
					fun->SetParLimits(0, 10, 100);
					fun->SetParLimits(1, 20, 200);
					fun->SetParLimits(3, 5, 80);
					double initPeak = h_MIP[i_layer][i_chip][i_channel]->GetBinCenter(h_MIP[i_layer][i_chip][i_channel]->GetMaximumBin());
					double initLX = initPeak - 80;
					double initRX = 600;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", initLX, initRX);
					double peakX = fun->GetMaximumX(initLX, initRX);
					// jiaxuan : range change
					double rangeLX = fitLeftThreshold * peakX ;
					double rangeRX = initRX;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", rangeLX, rangeRX);
					peakX = fun->GetMaximumX(rangeLX, rangeRX);
					rangeLX = fitLeftThreshold * peakX ;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", rangeLX, rangeRX);
					
					h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(-10,600);
					h_MIP[i_layer][i_chip][i_channel]->Draw();
					fun->Draw("same");
					gStyle->SetOptFit(1111);

				} else if (SiPMType == 15) {
					fun->SetParameters(64, 500, 25000, 100);
					fun->SetParLimits(0, 20, 200);
					fun->SetParLimits(1, 20, 800);
					fun->SetParLimits(3, 5, 250);
					double initPeak = h_MIP[i_layer][i_chip][i_channel]->GetBinCenter(h_MIP[i_layer][i_chip][i_channel]->GetMaximumBin());
					double initLX = initPeak - 300;
					double initRX = 1200; 
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", initLX, initRX);
					double peakX = fun->GetMaximumX(initLX,initRX);
					// jiaxuan : range change
					double rangeLX = fitLeftThreshold * peakX;
					double rangeRX = initRX;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", rangeLX, rangeRX);
					peakX = fun->GetMaximumX(rangeLX, rangeRX);
					rangeLX = fitLeftThreshold * peakX ;
					h_MIP[i_layer][i_chip][i_channel]->Fit(fun, "QR", "", rangeLX, rangeRX);
					
					h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(-10,2000);
					h_MIP[i_layer][i_chip][i_channel]->Draw();
					fun->Draw("same");
					gStyle->SetOptFit(1111);
				} else {
					delete fun;
					continue;
				}	
				
				// Set Branch Value
				double landauMPV = fun->GetParameter(1);
				double gausSigma = fun->GetParameter(3);
				//if (landauMPV < 50 || chn_entries < 500 || (gausSigma < 25 && SiPMType == 10) || (gausSigma < 50 && SiPMType == 15) ) {
				if (landauMPV < 70 || chn_entries < 500 || gausSigma < 10) {
					deadChannelNo++;
					if (landauMPV < 70 && SiPMType==10 ) MPVAbNu++;
					if (gausSigma < 10 && SiPMType==10  ) SigmaAbNu++;
					if (chn_entries < 500&& SiPMType==10 ) LackNu++;
					hDeadCountLayer->Fill(i_layer);
					hMIPChn->Fill(chipID_mip,i_channel,0);
					//cout << i_layer << " " << i_chip << " " << i_channel << endl;
					mipdeadtxt << fixed << setprecision(0) << i_layer*1e5 + i_chip*1e4 + i_channel << ", ";
					ForceSetBranchValue(i_layer,i_chip,i_channel);
					_landauMPV = fun->GetParameter(1);
					sprintf(char_tmp, "Layer%d_Chip%d_Chn%d_Abandoned", i_layer, i_chip, i_channel);
					h_MIP[i_layer][i_chip][i_channel]->SetTitle(char_tmp);
					h_MIP[i_layer][i_chip][i_channel]->SetName(char_tmp);
				} else {
					hMIPChn->Fill(chipID_mip,i_channel,fun->GetParameter(1));
					_cellID_mip = i_layer * 1e5 + i_chip * 1e4 + i_channel;
					_landauWidth = fun->GetParameter(0);
					_landauMPV = fun->GetParameter(1);
					_gausSigma = fun->GetParameter(3);
					_ChiSqr = fun->GetChisquare();
					_Ndf = fun->GetNDF();
					_isGoodChan = true;
				}
					
				gPad->Update();
				mipTree->Fill();
				h_MIP[i_layer][i_chip][i_channel]->Write();
				
				if (i_layer < layerNu) {
					if ((i_channel%9 == 8) || (i_channel==29 && i_chip==5)) {
						if(tmpp.size())	tmpp.clear();
						if (i_layer==0 && i_chip==0 && i_channel==8) {
							tmpp = pdfname + "(";
							Chn_MIP[page_no]->Print(tmpp.c_str());
						} else if (i_layer==(layerNu) && i_chip==5 && i_channel==29) {
							tmpp = pdfname + ")";
							Chn_MIP[page_no]->Print(tmpp.c_str());
						} else {
							tmpp = pdfname + "";
							Chn_MIP[page_no]->Print(tmpp.c_str());
						}
					}
				}
			}
			
			sprintf(char_tmp, "layer%d", i_layer);
			histfile->cd(char_tmp);
			sprintf(char_tmp, "layer%d_chip%d", i_layer, i_chip);
			C_MIP[i_layer][i_chip] = new TCanvas(char_tmp, char_tmp, 1);
			C_MIP[i_layer][i_chip]->Divide(6, 6);
			for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
				C_MIP[i_layer][i_chip]->cd(i_channel + 1);
				if (i_layer>=4 && i_layer <=27) {
					h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(-50,250);
				} else {
					h_MIP[i_layer][i_chip][i_channel]->GetXaxis()->SetRangeUser(-50,1200);
				}
				h_MIP[i_layer][i_chip][i_channel]->Draw();
			}
			C_MIP[i_layer][i_chip]->Write();
		}
	
	}
	cout << "MPV ABNORMAL Nu: " << MPVAbNu <<endl;
	cout << "SIGMA ABNORMAL Nu: " << SigmaAbNu <<endl;
	cout << "LACK STATISTICS Nu: " << LackNu <<endl;

	
	cout << "Total Dead Chn : " << deadChannelNo << " **************" << endl;
	histfile->cd();
	hEfficiency->Write();
	hDeadCountLayer->Write();
	hMIPChn->Draw("colz");
	hCountChn->Draw("colz");
	hMIPChn->Write();
	hCountChn->Write();
	mipFile->cd();
	mipTree->Write();
	
	inFile->Close();
	histfile->Close();
	mipFile->Close();

	//delete hEfficiency;
	//delete hDeadCountLayer;
	//delete hMIPChn;
	//delete hCountChn;
	//delete mipFile;
	//delete mipTree;
	
	cout << "========== MIP Fit Done ================="<< endl;
	return 0;
}


