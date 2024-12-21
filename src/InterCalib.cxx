/////////////////////////////////////////////////////
////////////////				/////////////////////
////			binFit with TH2D				/////
////////////////				/////////////////////
/////////////////////////////////////////////////////
#include<fstream>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include<TGraph.h>
#include<TRandom.h>
#include<TSpectrum.h>
#include<TVirtualFitter.h>
#include<TH2F.h>
#include"InterCalib.h"
#include"TPaveText.h"
#include"TText.h"
#include"TProfile.h"
#include"TPaveStats.h"
#include"EBUdecode.h"
#include"Extract.h"
#include <cmath>

using namespace std;

void interCalib::Init() {
	cout<<"=========== Initiation ============"<<endl;
	char fname[50];
	for (int i_layer = 0; i_layer < LayerNo; ++i_layer) {
        for (int i_chip = 0; i_chip < chipNu; ++i_chip)   {
            for (int i_chn = 0; i_chn < channelNu; ++i_chn) {
				sprintf(fname,"Layer%d_Chip%d_Chn%d",i_layer,i_chip,i_chn);
				//hisg[i_layer][i_chip][i_chn] = new TH2F(Form("%d_%d_%d",i_layer,i_chip,i_chn),Form("%d_%d_%d",i_layer,i_chip,i_chn),2000,-100.5,3899.5,400,-100.5,299.5);
				hisg[i_layer][i_chip][i_chn] = new TH2F("hl_ratio","hl_ratio",1000,-300,3700,400,-100,300);
				hisgRo[i_layer][i_chip][i_chn] = new TH2F("hl_rotate","hl_rotate",100,-10,190,1000,0,4000);
				hisg[i_layer][i_chip][i_chn]->SetName(fname);
				hisg[i_layer][i_chip][i_chn]->SetTitle(fname);
				sprintf(fname,"Layer%d_Chip%d_Chn%d_Rotate",i_layer,i_chip,i_chn);
				hisgRo[i_layer][i_chip][i_chn]->SetName(fname);
				hisgRo[i_layer][i_chip][i_chn]->SetTitle(fname);
				hisg[i_layer][i_chip][i_chn]->GetXaxis()->SetRangeUser(-300,3700);
				hisg[i_layer][i_chip][i_chn]->GetXaxis()->SetTitle("HighGain [ADC]");
				hisg[i_layer][i_chip][i_chn]->GetYaxis()->SetRangeUser(-100,300);
				hisg[i_layer][i_chip][i_chn]->GetYaxis()->SetTitle("LowGain [ADC]");
				hisg[i_layer][i_chip][i_chn]->SetOption("colz");
				/*
				grap[i_layer][i_chip][i_chn] = new TGraph();
				grap[i_layer][i_chip][i_chn]->SetName(fname);
				grap[i_layer][i_chip][i_chn]->SetTitle(fname);
				grap[i_layer][i_chip][i_chn]->GetXaxis()->SetLimits(-100,3900);
				grap[i_layer][i_chip][i_chn]->GetXaxis()->SetTitle("HighGain [ADC]");
				grap[i_layer][i_chip][i_chn]->GetYaxis()->SetRangeUser(-100,300);
				grap[i_layer][i_chip][i_chn]->GetYaxis()->SetTitle("LowGain [ADC]");
				*/	
				eventNo[i_layer][i_chip][i_chn]=0;
				Xmax[i_layer][i_chip][i_chn]=0;
			}
		}
	}

}

void interCalib::HistCut(TH2F* &hist2d, double slope, double intercept) {
	int interceptOffset = 8;
	for (int ix = 1; ix <= hist2d->GetNbinsX(); ++ix) {
		for (int iy = 1; iy <= hist2d->GetNbinsY(); ++iy) {
			double xv = hist2d->GetXaxis()->GetBinCenter(ix);
			double yv = hist2d->GetYaxis()->GetBinCenter(iy);
			if (yv < xv*slope + intercept - interceptOffset) {
				hist2d->SetBinContent(ix, iy, 0);
			}
		}
	}	
}

interCalib::interCalib(string root_name, string spshist_name, string spsfit_name, Extract *& ex) {
    Init();
	ofstream hlratiodeadtxt("temp/hlratiodead.txt");
    // ================ Read DATA =======================    
    TFile* infile = TFile::Open(root_name.c_str(),"READ");
    if(!infile)cout<<"read file failed"<<endl;
    TTree* mipTreeRead = (TTree*)infile->Get("Raw_Hit");
    //TTree* mipTreeRead = (TTree*)infile->Get("Digi_Raw_Hit");
    if(!mipTreeRead)cout<<"Get tree failed"<<endl;
    mipTreeRead->SetBranchAddress("CellID",&_cellID);
    mipTreeRead->SetBranchAddress("HitTag",&_hitTag);
    mipTreeRead->SetBranchAddress("LG_Charge",&_lg_charge);
    mipTreeRead->SetBranchAddress("HG_Charge",&_hg_charge);
    mipTreeRead->SetBranchAddress("Temperature",&_tempLayer);

    for(int entry=0; entry<mipTreeRead->GetEntries(); ++entry) {
    //for (int entry=0; entry<100000; ++entry) {
		int badHitNo = 0;
        mipTreeRead->GetEntry(entry);
        if (entry%10000==0)	cout<<" [ Process ]  event : "<< entry << endl;
        for (unsigned i_hit=0; i_hit!=_cellID->size(); ++i_hit) {
			if(_hitTag->at(i_hit)==0) continue;
            int layerID = _cellID->at(i_hit)/1e5;
            int chipID = _cellID->at(i_hit)%100000/1e4;
            int chnID = _cellID->at(i_hit)%100;
			double hg_charge = _hg_charge->at(i_hit);
            double lg_charge = _lg_charge->at(i_hit);
			double *_position = EBUdecode(layerID,chipID,chnID);
			double SiPMtemp;
			if(_tempLayer->at(layerID).size()!=0) {
				SiPMtemp = tempReconstruction(layerID, _position, _tempLayer->at(layerID));
			} else {
				SiPMtemp = 23;// SiPMtemp=20
			}
			//lg_charge *= (1+ (SiPMtemp-23)*0.03);
			//hg_charge *= (1+ (SiPMtemp-23)*0.03);
			lg_charge -= ex->pedMeanLow[layerID][chipID][chnID];
			hg_charge -= ex->pedMeanHigh[layerID][chipID][chnID];
			hisg[layerID][chipID][chnID]->Fill(hg_charge, lg_charge);
			//grap[layerID][chipID][chnID]->AddPoint(hg_charge-pedMeanTimesChnLevel[layerID][chipID][chnID],lg_charge-pedMeanChargesChnLevel[layerID][chipID][chnID]);
			if (Xmax[layerID][chipID][chnID] < hg_charge ) Xmax[layerID][chipID][chnID] = hg_charge;
			++eventNo[layerID][chipID][chnID];
        }
    }
    infile->Close();

    // ================ LINEAR FITTING =======================    
    TFile* spsfile = new TFile(spshist_name.c_str(),"RECREATE");
	TFile* spsfilefit = new TFile(spsfit_name.c_str(),"RECREATE");
    TTree* fNtuple = new TTree("InterCalib_Cut","fit SPS with multigaus");
    fNtuple->Branch("CellID",&_hlCellID); //i_chn + i_chip*36 + i_layer*210
    fNtuple->Branch("Slope",&_slope);
    fNtuple->Branch("Intercept",&_intercept);
    fNtuple->Branch("SlopeError",&_slopeError);
    fNtuple->Branch("InterceptError",&_interceptError);
    fNtuple->Branch("ChiSquare",&_chiSquare);
    fNtuple->Branch("NDF",&_NDF);
	fNtuple->Branch("XMax",&_xmax);
	fNtuple->Branch("Statistics",&_statistics);
	fNtuple->Branch("BadHitNo",&_badHitNo);

    char fname[50];
    int badChanNo=0;
	TH2D *hlratiocount = new TH2D("Statistics","Statistics",192,0,192,36,0,36);
	TH2D *hlratiovalue = new TH2D("High-low gain ratio","High-low gain ratio",192,0,192,36,0,36);
    //for (int i_layer = 0; i_layer < LayerNo; ++i_layer) {
    for (int i_layer = 0; i_layer < 30; ++i_layer) {
        cout<<" -------- "<< i_layer <<" -------- "<<endl;
        for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
    		TCanvas* c1 = new TCanvas(Form("c_%d_%d",i_layer,i_chip),Form("c_%d_%d",i_layer,i_chip),10,10,1500,1200);
			c1->Divide(6,6);
			int pad_no = 0;
            sprintf(fname,"Layer%d/Chip%d",i_layer,i_chip);
            spsfile->mkdir(fname);
            spsfile->cd(fname);
            for (int i_chn = 0; i_chn < channelNu; ++i_chn) {
				bool deadchannelyet = 0;
    			TCanvas* c2 = new TCanvas(Form("c_%d_%d_%d",i_layer,i_chip,i_chn),Form("c_%d_%d_%d",i_layer,i_chip,i_chn),800,600);
				//if(Xmax[i_layer][i_chip][i_chn] < 900 && Xmax[i_layer][i_chip][i_chn] != 0) cout << "Xmax : " << Xmax[i_layer][i_chip][i_chn] <<" layer:" << i_layer <<" chip:" << i_chip <<" chn:"<< i_chn << endl;
				pad_no++;
				// jiaxuan : position
				//double *posi = EBUdecode(i_layer, i_chip, i_chn);
				if((i_chip==5&&i_chn>29) && pad_no != 36) continue;
				if((i_chip==5&&i_chn>29) && pad_no == 36) {c1->Write(); continue;}
				if(hisg[i_layer][i_chip][i_chn]->GetEntries()<50 || Xmax[i_layer][i_chip][i_chn] < 900 ) {
				//if(grap[i_layer][i_chip][i_chn]->GetN()<50 || Xmax[i_layer][i_chip][i_chn] < 900 ) {
					cout << "===>"<< i_layer << " "<< i_chip << " "<< i_chn << endl;
					hlratiodeadtxt << fixed << setprecision(0) << i_layer*1e5 + i_chip*1e4 + i_chn << ", ";
					deadchannelyet = 1;
					badChanNo++;
					if(pad_no!=36) continue;
					if(pad_no==36) {c1->Write(); continue;}
				}
				
				int highGainOffset = 600; 
				int Xmax_Fit = Xmax[i_layer][i_chip][i_chn] - highGainOffset;
				Xmax_Fit = Xmax_Fit > 2000 ? Xmax_Fit : 2000; 
				int Xmin_Fit = 0;
				if(hisg[i_layer][i_chip][i_chn]->GetEntries()>100000) {
					Xmin_Fit = 800;
				} else if(hisg[i_layer][i_chip][i_chn]->GetEntries()<10000) {
					Xmin_Fit = 200;
				} else {
					Xmin_Fit = 500;
				}
				TF1* fun = new TF1("fun","[0]*x+[1]",Xmin_Fit,Xmax_Fit);
				fun->SetParameter(0,0.03);
				fun->SetParameter(1,0.0);
                fun->SetParLimits(0,0.020,0.037);//27~43.4
                fun->SetParLimits(1,-100,100);
                fun->SetLineColor(2);
                fun->SetLineWidth(2);
                
				hisg[i_layer][i_chip][i_chn]->Fit(fun,"qn","",Xmin_Fit,Xmax_Fit);
				slope[i_layer][i_chip][i_chn] = fun->GetParameter(0);
				intercept[i_layer][i_chip][i_chn] = fun->GetParameter(1);

				//HistCut(hisg[i_layer][i_chip][i_chn], slope[i_layer][i_chip][i_chn], intercept[i_layer][i_chip][i_chn]);
				
				fun->SetParameter(0,fun->GetParameter(0));
				fun->SetParameter(1,fun->GetParameter(1));
                fun->SetParLimits(0,0.020,0.037);//27~43.4
                fun->SetParLimits(1,-100,100);
				hisg[i_layer][i_chip][i_chn]->Fit(fun,"qn","",Xmin_Fit,Xmax_Fit);
				fitLeftLimits[i_layer][i_chip][i_chn] = Xmin_Fit;
				fitRightLimits[i_layer][i_chip][i_chn] = Xmax_Fit;
				slope[i_layer][i_chip][i_chn] = fun->GetParameter(0);
				intercept[i_layer][i_chip][i_chn] = fun->GetParameter(1);

				/*	
				if(grap[i_layer][i_chip][i_chn]->GetN()>40000) {
					grap[i_layer][i_chip][i_chn]->Fit(fun,"qn","",300,Xmax_Fit);
				} else {
					grap[i_layer][i_chip][i_chn]->Fit(fun,"qn","",200,Xmax_Fit);
				}
				*/	
				TPaveText *pt1 = new TPaveText(.16,.5,.42,.65,"blNDC");
				TText *t1 = pt1->AddText(Form("L = %f * H + %f", fun->GetParameter(0), fun->GetParameter(1)));
				TText *t2 = pt1->AddText(Form("H = %f * L + %f", 1/fun->GetParameter(0), -fun->GetParameter(1)/fun->GetParameter(0)));
                // ==============================================================
                c1->cd(pad_no);
				gStyle->SetOptStat(0000);
                hisg[i_layer][i_chip][i_chn]->Draw();
                //grap[i_layer][i_chip][i_chn]->Draw("AP");
                fun->Draw("same");
				pt1->Draw("same");
				gPad->Modified();
				gPad->Update();
				if(pad_no==36) {c1->Write();}
					
				c2->cd();
				gStyle->SetOptStat(0000);
				hisg[i_layer][i_chip][i_chn]->Draw();
                //grap[i_layer][i_chip][i_chn]->Draw("AP");
                fun->Draw("same");
				pt1->Draw("same");
				c2->Modified();
				c2->Update();
				c2->Write();
                // ==============================================================
				
				hlratiovalue->Fill(i_layer*6+i_chip, i_chn, 1/(fun->GetParameter(0)));
				hlratiocount->Fill(i_layer*6+i_chip, i_chn, hisg[i_layer][i_chip][i_chn]->GetEntries());
				//hlratiocount->Fill(i_layer*6+i_chip, i_chn, grap[i_layer][i_chip][i_chn]->GetN());
                _hlCellID = i_chn + i_chip*1e4 + i_layer*1e5;
				_slope = fun->GetParameter(0);
				if (_hlCellID==1020014) cout << "1020014 slope : " << _slope << endl;
                _intercept = fun->GetParameter(1);
                _slopeError = fun->GetParError(0);
                _interceptError = fun->GetParError(1);
                _chiSquare = fun->GetChisquare();
                _NDF = fun->GetNDF();
				_xmax = Xmax[i_layer][i_chip][i_chn];
				_statistics = hisg[i_layer][i_chip][i_chn]->GetEntries();
				if(_slope < 0.022 || _slope > 0.035) {	
					cout << "slope beyond limits (0.022, 0.035) " << i_layer << " " << i_chip << " " << i_chn << endl;
					if(!deadchannelyet) hlratiodeadtxt << fixed << setprecision(0) << i_layer*1e5 + i_chip*1e4 + i_chn << ", "; 
				}
                fNtuple->Fill();
            }
        }
    }
	
	FillRotateHist(root_name, ex);
	TH1D * hDriftRatio = new TH1D("driftRatio","driftRatio",100,0,0.05);
	TH1D * hDriftRatioCenter = new TH1D("driftRatioC","driftRatioC",100,0,0.05);
	TH1D * hDriftAmplitude = new TH1D("driftMaxAmp","driftMaxAmp",300,0,300);
	TH1D * hDriftAmplitudeCenter = new TH1D("driftMaxAmpC","driftMaxAmpC",300,0,300);
	for (int i_layer = 0; i_layer < 30; ++i_layer) {
        cout<<" -------- "<< i_layer <<" -------- "<<endl;
        for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
            sprintf(fname,"Layer%d/Chip%d",i_layer,i_chip);
            spsfile->cd(fname);
            for (int i_chn = 0; i_chn < channelNu; ++i_chn) {
				TH1D * hProj = hisgRo[i_layer][i_chip][i_chn]->ProjectionX();
				TF1 * fProj = new TF1("fProj","gaus(0)");
				hProj->Fit(fProj,"Q");
				double mean = fProj->GetParameter(1);
				double sigma = fProj->GetParameter(2);
				int binIndex1 = hProj->GetXaxis()->FindBin(mean-3*sigma);
				int binIndex2 = hProj->GetXaxis()->FindBin(mean+3*sigma);
				int normalEntry = hProj->Integral(binIndex1,binIndex2);
				double ratio = 1- normalEntry / hProj->GetEntries();
				hDriftRatio->Fill(ratio);
				if (hProj->GetEntries() > 20000) hDriftRatioCenter->Fill(ratio);
    			TCanvas* c3 = new TCanvas(Form("crot_%d_%d_%d",i_layer,i_chip,i_chn),Form("crot_%d_%d_%d",i_layer,i_chip,i_chn),1300,600);
				c3->Divide(2,1);
				gStyle->SetOptStat(1111);
				c3->cd(1);
				hisgRo[i_layer][i_chip][i_chn]->Draw("colz");
				c3->cd(2);
				gPad->SetLogy();
				hProj->Draw();
				c3->Modified();
				c3->Update();
				c3->Write();
				hDriftAmplitude->Fill(maxLowDrift[i_layer][i_chip][i_chn]);
				if (hProj->GetEntries() > 20000) hDriftAmplitudeCenter->Fill(maxLowDrift[i_layer][i_chip][i_chn]);
			}
		}
	}
	TCanvas * cDrift = new TCanvas("CDrift","CDrift",1000,1000);
	gStyle->SetOptStat(1111);
	cDrift->Divide(2,2);
	spsfile->cd();
	cDrift->cd(1);
	hDriftRatio->Draw();
	cDrift->cd(2);
	hDriftRatioCenter->Draw();
	cDrift->cd(3);
	hDriftAmplitude->Draw();
	cDrift->cd(4);
	hDriftAmplitudeCenter->Draw();
	cDrift->Write();
	hlratiodeadtxt.close();
	spsfile->Close();

	TCanvas *coverview = new TCanvas("","",1400,1400);
	coverview->Divide(1,2);
	coverview->cd(1); hlratiovalue->Draw("colz");
	coverview->cd(2); hlratiocount->Draw("colz");
    // *********************************************************************************
    spsfilefit->cd();
	coverview->Write();
    fNtuple->Write();
    spsfilefit->Close();
	
}

void interCalib::RotatePoint(double &xv, double &yv, double slopek) {
	double slopeRad = atan(slopek);
	double thetaRad = M_PI/2 - slopeRad;  
    double xPrime = xv * cos(thetaRad) - yv * sin(thetaRad);
    double yPrime = xv * sin(thetaRad) + yv * cos(thetaRad);
	xv = xPrime;
	yv = yPrime;
}

void interCalib::FillRotateHist(string root_name, Extract *&ex) {
	TFile* infile = TFile::Open(root_name.c_str(),"READ");
    if(!infile)cout<<"read file failed"<<endl;
    TTree* mipTreeRead = (TTree*)infile->Get("Raw_Hit");
    //TTree* mipTreeRead = (TTree*)infile->Get("Digi_Raw_Hit");
    if(!mipTreeRead)cout<<"Get tree failed"<<endl;
    mipTreeRead->SetBranchAddress("CellID",&_cellID);
    mipTreeRead->SetBranchAddress("HitTag",&_hitTag);
    mipTreeRead->SetBranchAddress("LG_Charge",&_lg_charge);
    mipTreeRead->SetBranchAddress("HG_Charge",&_hg_charge);
    mipTreeRead->SetBranchAddress("Temperature",&_tempLayer);

    for(int entry=0; entry<mipTreeRead->GetEntries(); ++entry) {
    //for (int entry=0; entry<100000; ++entry) {
		int badHitNo = 0;
        mipTreeRead->GetEntry(entry);
        if (entry%10000==0)	cout<<" [ Process ]  event : "<< entry << endl;
        for (unsigned i_hit=0; i_hit!=_cellID->size(); ++i_hit) {
			if(_hitTag->at(i_hit)==0) continue;
            int layerID = _cellID->at(i_hit)/1e5;
            int chipID = _cellID->at(i_hit)%100000/1e4;
            int chnID = _cellID->at(i_hit)%100;
			double hg_charge = _hg_charge->at(i_hit);
            double lg_charge = _lg_charge->at(i_hit);
			double *_position = EBUdecode(layerID,chipID,chnID);
			double SiPMtemp;
			//if(_tempLayer->at(layerID).size()!=0) {
			//	SiPMtemp = tempReconstruction(layerID, _position, _tempLayer->at(layerID));
			//} else {
			//	SiPMtemp = 23;// SiPMtemp=20
			//}
			//lg_charge *= (1+ (SiPMtemp-23)*0.03);
			//hg_charge *= (1+ (SiPMtemp-23)*0.03);
			lg_charge -= ex->pedMeanLow[layerID][chipID][chnID];
			hg_charge -= ex->pedMeanHigh[layerID][chipID][chnID];
			if (hg_charge < fitRightLimits[layerID][chipID][chnID] && hg_charge > fitLeftLimits[layerID][chipID][chnID]) {
				double deltaDrift = hg_charge * slope[layerID][chipID][chnID] - lg_charge;
				maxLowDrift[layerID][chipID][chnID] = maxLowDrift[layerID][chipID][chnID] > deltaDrift ? maxLowDrift[layerID][chipID][chnID] : deltaDrift;
				RotatePoint(hg_charge, lg_charge, slope[layerID][chipID][chnID]);
				hisgRo[layerID][chipID][chnID]->Fill(hg_charge, lg_charge);
        	}
		}
    }
    infile->Close();
}

 
