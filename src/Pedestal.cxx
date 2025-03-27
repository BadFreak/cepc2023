#include<fstream>
#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<vector>
#include<TROOT.h>
#include<TStyle.h>
#include<TMath.h>
#include<TMath.h>
#include<TH2.h>
#include"TCanvas.h"
#include"../include/Pedestal.h"

using namespace std;

void MyFormatChan(TH2D* &h2d)
{
	h2d->GetXaxis()->SetTitle("chipNu");
	h2d->GetXaxis()->SetTitleOffset(2);
	h2d->GetYaxis()->SetTitle("channelNu");
	h2d->GetXaxis()->SetRangeUser(0,192.);
	h2d->GetYaxis()->SetRangeUser(0,36.);
}

void treePedestal::Reset()
{
    _pedCellID.clear();
    _pedMeanCharges.clear();
    _pedSigmaCharges.clear();
    _pedMeanTimes.clear();
    _pedSigmaTimes.clear();
}


void treePedestal::pedInit() 
{
    char hChargeName[200];
    //Hlist1 = new TObjArray(10);
    //Hlist3 = new TObjArray(10);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
    		Hlist2l[l][i] = new TObjArray(10);
			Hlist2h[l][i] = new TObjArray(10);
            for(int j=0; j!=channelNu; ++j)
            {
				/*
                for(int k=0; k!=memoryNu; ++k)
                {
                    sprintf(hChargeName,"Charges_layer%d_chip%d_chan%d_cell%d",l,i,j,k);
                    hCharge1[l][i][j][k] = new TH1D(hChargeName,hChargeName,4096,0,4095);
                    Hlist1->Add(hCharge1[l][i][j][k]);
                    sprintf(hChargeName,"Times_layer%d_chip%d_chan%d_cell%d",l,i,j,k);

                    hTime1[l][i][j][k] = new TH1D(hChargeName,hChargeName,4096,0,4095);
                    Hlist1->Add(hTime1[l][i][j][k]);
                }
				*/
                sprintf(hChargeName,"lg_layer%d_chip%d_chan%d",l,i,j);
                hCharge2[l][i][j] = new TH1D(hChargeName,hChargeName,4096,0,4095);
                Hlist2l[l][i]->Add(hCharge2[l][i][j]);                
                sprintf(hChargeName,"hg_layer%d_chip%d_chan%d",l,i,j);
                hTime2[l][i][j] = new TH1D(hChargeName,hChargeName,4096,0,4095);
                Hlist2h[l][i]->Add(hTime2[l][i][j]);
            }
            //sprintf(hChargeName,"Charges_layer%d_chip%d",l,i);
            //hCharge3[l][i] = new TH1D(hChargeName,hChargeName,4096,0,4095);
            //Hlist3->Add(hCharge3[l][i]);
            //sprintf(hChargeName,"Times_layer%d_chip%d",l,i);
            //hTime3[l][i] = new TH1D(hChargeName,hChargeName,4096,0,4095);
            //Hlist3->Add(hTime3[l][i]);
        }
    }
}


treePedestal::~treePedestal(){}


int treePedestal::Pedestal(string root_name,string pedestal_name,string fitdata_name)
{
    pedInit();
    string txtname(root_name);
	cout << "txtname  : " << txtname << endl;
    txtname.erase(0,txtname.find_last_of("/")+1);
	cout << "txtname  : " << txtname << endl;

    // *********************************************************************************
    //TFile *f =  TFile::Open("test.root","RECREATE");
    rawfile =  TFile::Open(root_name.c_str(),"READ");
    if(!rawfile){ cout<<" Get open eventAlign file failed !!!"<<endl;return 0;}
    rawfile->GetObject("Raw_Hit",tree);
    if(!tree){ cout<<" Get tree object failed !!!"<<endl;return 0;}
    int Entries = tree->GetEntries();

    cout<<" [ Prepare to PedestalExtraction : ] Total Entries : "<<Entries<<" "<<endl;
	//tree ->SetBranchAddress("Run_Num",&_run_num);
	tree ->SetBranchAddress("Event_Time",&_event_time);
    //tree ->SetBranchAddress("CycleID",&_cycleID);
	//tree ->SetBranchAddress("TriggerID",&_triggerID);
    tree ->SetBranchAddress("CellID",&_cellID);
    tree ->SetBranchAddress("HitTag",&_hitTag);
    tree ->SetBranchAddress("GainTag",&_gainTag);
    tree ->SetBranchAddress("LG_Charge",&_lg_charge);
    tree ->SetBranchAddress("HG_Charge",&_hg_charge);

    // *********************************************************************************
    pedestal_name.insert(pedestal_name.length(),txtname);
    fitdata_name.insert(fitdata_name.length(), txtname);
    // *********************************************************************************
	// drop some pedestal
	// jiaxuan : note
	/*tree->GetEntry(tree->GetEntries()-1);
	const int event_time_max = _event_time + 1 ;
	cout << "event_time_max : " << event_time_max << endl;
	vector<vector<int> > count(192,vector<int>(event_time_max));
	for(int entry=0; entry!=Entries; ++entry)
    {
        tree->GetEntry(entry);
		int event_time = _event_time;
		for(unsigned i_hit=0; i_hit!=_hitTag->size(); ++i_hit)
        {
			int channelID = _cellID->at(i_hit)%100;
			if(channelID==0)//channel 0 on each chip is used
			{
				int layerID = _cellID->at(i_hit)/1e5;
				int chipID = _cellID->at(i_hit)%100000/1e4;
				count[layerID*6+chipID][event_time]++;
			}
		}
	}*/
    // ==============================================================================
    //                          Extracting pedestal 
    // ==============================================================================
	//int skip_flag[192] = {0};
	//int last_skip_event_time[192] = {0};
    for(int entry=0; entry!=Entries; ++entry)
    {
        tree->GetEntry(entry);
        if(entry%10000==0)  cout<<" [ Process : PedestalExtraction : ]  event : "<< entry <<" "<<endl;
		int event_time = _event_time;
		if(event_time==0) continue;
		
		for(unsigned i_hit=0; i_hit!=_hitTag->size(); ++i_hit)
        {
            int layerID = _cellID->at(i_hit)/1e5;
            int chipID = _cellID->at(i_hit)%100000/1e4;
            //int memoryCells = _cellID->at(i_hit)%10000/1e2;
            int channelID = _cellID->at(i_hit)%100;
            double  lg_charge = _lg_charge->at(i_hit);
            double  hg_charge  = _hg_charge->at(i_hit);
            //if(_hitTag->at(i_hit)!=0) continue;// only save hit=0 events
            if(_hitTag->at(i_hit)==0) continue;// only save hit=1 events
			
			//if(_cellID->at(i_hit)==600000)	cout << "Event time : " << event_time << endl;	
			
			// drop some pedestal
			// jiaxuan : note
			/*if(count[layerID*6+chipID][event_time-1]==0 && count[layerID*6+chipID][event_time]!=0 && last_skip_event_time[layerID*6+chipID] < event_time)	
			{
				skip_flag[layerID*6+chipID] = 1; 
				//if(_cellID->at(i_hit)==600000) cout << " Drop Event time : " << event_time << endl;	
				continue;
			}*/

		    //hCharge1[layerID][chipID][channelID][memoryCells] ->Fill(lg_charge);
            //hTime1[layerID][chipID][channelID][memoryCells] ->Fill(hg_charge);
            hCharge2[layerID][chipID][channelID] ->Fill(lg_charge);
            hTime2[layerID][chipID][channelID] ->Fill(hg_charge);
            //hCharge3[layerID][chipID] ->Fill(lg_charge);
            //hTime3[layerID][chipID] ->Fill(hg_charge);
        }
		// jiaxuan : forget the meaning, are intend to do the drop
		//for(int j=0; j<192; j++)
		//{
		//	if(skip_flag[j])
		//	{
		//		skip_flag[j] = 0;
		//		last_skip_event_time[j] = event_time;
		//	}
		//}
    }
	
	// ==============================================================================
    //                          Fitting pedestal 
    // ==============================================================================
    gausfun = new TF1("gausfun","gaus");
//    cout<<fitdata_name<<endl;
    fitfile = new TFile(fitdata_name.c_str(),"RECREATE");
    
	cout<<" ============== Begin Channel Level pedestal fitting  Process==============="<<endl;
    Reset();
    fNtuple2 = new TTree("ChnLevel","fit pedestal at channel Level");
    fNtuple2->Branch("CellID",&_pedCellID);
    fNtuple2->Branch("PedLowMean",&_pedMeanCharges);
    fNtuple2->Branch("PedLowSig",&_pedSigmaCharges);
    fNtuple2->Branch("PedHighMean",&_pedMeanTimes);
    fNtuple2->Branch("PedHighSig",&_pedSigmaTimes);
	low_mean = new TH2D("low_gain_mean","low_gain_mean",192,0,192,36,0,36);
	low_rms = new TH2D("low_gain_rms","low_gain_rms",192,0,192,36,0,36);
	high_mean = new TH2D("high_gain_mean","high_gain_mean",192,0,192,36,0,36);
	high_rms = new TH2D("high_gain_rms","high_gain_rms",192,0,192,36,0,36);
    for(int l=0; l!=layerNu; ++l)
    {
        for(int i=0; i!=chipNu; ++i)
        {
            for(int j=0; j!=channelNu; ++j)
            {
                // Channel Level
                if(hCharge2[l][i][j]->GetEntries()==0) {cout <<l<<" "<<i<<" "<<j<<" abnormal"<< endl; continue;}
                if(hTime2[l][i][j]->GetEntries()==0) {cout <<l<<" "<<i<<" "<<j<<" abnormal"<< endl; continue;}
                _pedCellID.push_back(l*pow(10,5)+i*pow(10,4)+0*pow(10,2)+j);
				if(j == (channelNu-1) ) cout << endl;
				//  The ADC distribution of hits with hitTag == 0 shows several small peaks in larger region.
				//  So the fit based on the gaussian mean and RMS works not well for the beam test data.
				
                //double meanCharges = hCharge2[l][i][j]->GetMean();
                //double rmsCharges = hCharge2[l][i][j]->GetRMS();
                double meanCharges = hCharge2[l][i][j]->GetBinCenter(hCharge2[l][i][j]->GetMaximumBin());
                double rmsCharges = 5;
                hCharge2[l][i][j] ->Fit(gausfun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = gausfun->GetParameter(1);
                rmsCharges = gausfun->GetParameter(2);
                hCharge2[l][i][j] ->Fit(gausfun,"Q","",meanCharges-2*rmsCharges,meanCharges+2*rmsCharges);
                meanCharges = gausfun->GetParameter(1);
                rmsCharges = gausfun->GetParameter(2);
                _leftC = (meanCharges-10*rmsCharges)>0?(meanCharges-10*rmsCharges):0;
                _righC = (meanCharges+10*rmsCharges)<4096?(meanCharges+10*rmsCharges):4096;
				if(rmsCharges>10) cout << "Low layer:chip:chan " << l << " " << i << " " << j << endl; 
                //hCharge2[l][i][j] ->GetXaxis()->SetRangeUser(_leftC,_righC);
                hCharge2[l][i][j]->GetXaxis()->SetRangeUser(300,700);
				
				//cout << "layerNu, chipNu, channelNu:" << l <<" "<< i<<" " << j << endl;
				low_mean->Fill((i+l*6), j, meanCharges ? meanCharges : 1000);
				low_rms->Fill((i+l*6), j, rmsCharges ? rmsCharges : 1000);
				_pedMeanCharges.push_back(meanCharges);
                _pedSigmaCharges.push_back(rmsCharges);

                //double meanTimes = hTime2[l][i][j]->GetMean();
                //double rmsTimes = hTime2[l][i][j]->GetRMS();
                double meanTimes = hTime2[l][i][j]->GetBinCenter(hTime2[l][i][j]->GetMaximumBin());
                double rmsTimes = 3;
                hTime2[l][i][j] ->Fit(gausfun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = gausfun->GetParameter(1);
                rmsTimes = gausfun->GetParameter(2);
                hTime2[l][i][j] ->Fit(gausfun,"Q","",meanTimes-2*rmsTimes,meanTimes+2*rmsTimes);
                meanTimes = gausfun->GetParameter(1);
                rmsTimes = gausfun->GetParameter(2);
                _leftC = (meanTimes-10*rmsTimes)>0 ? (meanTimes-10*rmsTimes) : 0;
                _righC = (meanTimes+10*rmsTimes)<4096 ? (meanTimes+10*rmsTimes) : 4096;
				if(rmsTimes>10) cout << "High layer:chip:chan " << l << " " << i << " " << j << endl; 
                //hTime2[l][i][j]->GetXaxis()->SetRangeUser(_leftC,_righC);
                hTime2[l][i][j]->GetXaxis()->SetRangeUser(300,700);
				high_mean->Fill((i+l*6), j, meanTimes ? meanTimes : 0);
				high_rms->Fill((i+l*6), j, rmsTimes ? rmsTimes : 0);
				_pedMeanTimes.push_back(meanTimes);
                _pedSigmaTimes.push_back(rmsTimes);

/*				char tmpchar[100];
				sprintf(tmpchar,"Charges_layer%d_chip%d_chan%d.png",l,i,j);
				if(meanCharges>500 || meanCharges < 400 || rmsCharges < 0.1)
					hCharge2[l][i][j]->Print(tmpchar);
					
				sprintf(tmpchar,"Times_layer%d_chip%d_chan%d.png",l,i,j);
				if(meanTimes > 500 || meanTimes < 400 || rmsTimes < 0.1)
					hTime2[l][i][j]->Print(tmpchar);
*/					
            }
        }
    }
    fNtuple2->Fill();
    cout<<" ============== Finish Channel Level pedestal fitting Process==============="<<endl;
    

// **************************** Chip Cells Level *****************************************************
// *********************************************************************************
	MyFormatChan(low_mean);
	MyFormatChan(low_rms);
	MyFormatChan(high_mean);
	MyFormatChan(high_rms);
	low_mean->GetZaxis()->SetRangeUser(0,550.);
	high_mean->GetZaxis()->SetRangeUser(0,550.);
	low_rms->GetZaxis()->SetRangeUser(0,30.);
	high_rms->GetZaxis()->SetRangeUser(0,30.);

	low_mean->Write();	
 	low_rms->Write();	
 	high_mean->Write();	
 	high_rms->Write();	
    
//	fNtuple1->Write();
    fNtuple2->Write();
//  fNtuple3->Write();
    fitfile->Close();

// **********************************************************************************

    pedestalfile = new TFile(pedestal_name.c_str(),"RECREATE");
    //pedestalfile ->mkdir("MemCellLevel");
    //pedestalfile ->cd("MemCellLevel");
    //Hlist1 ->Write();

	for(int i_layer=0; i_layer<layerNu; i_layer++)
	{
		for(int i_chip=0; i_chip<chipNu; i_chip++)
		{
			pedestalfile->mkdir(Form("ChnLevel/Layer%d/Chip%d",i_layer,i_chip));
			pedestalfile->cd(Form("ChnLevel/Layer%d/Chip%d",i_layer,i_chip));
			Hlist2h[i_layer][i_chip]->Write();
			Hlist2l[i_layer][i_chip]->Write();
		}
	}
    //pedestalfile ->mkdir("ChipLevel");
    //pedestalfile ->cd("ChipLevel");
    //Hlist3 ->Write();
    pedestalfile->Close();
    
    cout<<"Ending of the process!!!!"<<endl;
	exit(0);
}
