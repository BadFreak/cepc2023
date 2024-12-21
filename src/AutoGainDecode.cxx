/********************************************************
 * Filename : AutoGainDecode.cxx
 * Author : Jiaxuan
 * Description : Decode  mode file
 * Function_List ：
 * 		... ... 
 * Last_Edit : 2023/3/21 
 *******************************************************/

#include "AutoGainDecode.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <TFile.h>
#include "TTree.h"
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TStyle.h>
#include "TMath.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
using namespace std;

int small_abnormal_bag_auto = 0;
int temp_abnormal_bag_auto = 0;
int data_abnormal_bag_auto = 0;
int empty_bag_auto = 0;
int good_bag_auto = 0;

//AutoGainTag
TH1F *hgaintag = new TH1F("","",2,0,2);
TH1F *hhittag = new TH1F("","",2,0,2);


/*int raw2RootAuto::IsDataAvailable(int layerID, int chipID)
{
	if(_chip_v[layerID][chipID].size())
		return 1;
	else
		return 0;
}
*/

// judge if element in vector are same
bool raw2RootAuto::IsSame(vector<int> vec)
{
	int stdvalue = vec[0];
	for(unsigned int i = 0; i < vec.size(); i++)
	{
		if(vec[i]!=stdvalue) return false;
	}
	return true;
}

// branch of vector clear
void raw2RootAuto::BranchClear()
{
	_cellID.clear();
	_bcid.clear();
	_hitTag.clear();
	_gainTag.clear();
	_chargehi.clear();
	_chargelo.clear();
	_hit_time.clear();
	for (int i_layer = 0; i_layer < layer_No; i_layer++)
		_temp[i_layer].clear();
}

void raw2RootAuto::SetTreeBranch(TTree* &tree)
{
	tree->Branch("Run_Num",&_run_no);
	tree->Branch("Event_Time",&_event_time);
	tree->Branch("CycleID", &_cycleID);
	tree->Branch("TriggerID", &_triggerID);
	tree->Branch("CellID", &_cellID);
	tree->Branch("BCID", &_bcid);
	tree->Branch("HitTag", &_hitTag);
	tree->Branch("GainTag", &_gainTag);
	tree->Branch("HG_Charge", &_chargehi);//for high gain charge
	tree->Branch("LG_Charge", &_chargelo);//for low gain charge 
	tree->Branch("Hit_Time", &_hit_time);//for AutoGain mode time information
	tree->Branch("Temperature", &_temp);//ECAL temperature data is with scientific data
}


int raw2RootAuto::CatchSPIROCBag(ifstream& f_in, vector<int>& buffer_v, int& cycleID, int& triggerID)
{
	bool event_begin_flag = 0;
	bool event_end_flag = 0;
	bool layer_begin_flag = 0;
	bool layer_end_flag = 0;
	int buffer = 0;
	int layerID;
	vector<int> testbuffer;
	if(buffer_v.size()) buffer_v.clear();
	// find event bag head "fbee fbee"
	while(!event_begin_flag && f_in.read((char*)(&buffer),1) && !f_in.eof())//each time 1 byte, that is 2 0x number
	{
		buffer_v.push_back(buffer);
		if(buffer_v.size()>4)	{buffer_v.erase(buffer_v.begin(),buffer_v.end()-4);}
		if(buffer_v[0]==0xfb &&buffer_v[1]==0xee && buffer_v[2]==0xfb && buffer_v[3]==0xee && buffer_v.size()==4)
		{
			//cout << "******* Find the event bag head *******" << endl;
			event_begin_flag = 1;
			buffer_v.clear();
			LayerTriggerID.clear();
		}
	}
		
	// start collect each layer 
	//************************************************************************************************
	//*** fa5a fa5a cycleid(32bit) sci_dat(36+36+bcid+chip) triggerid f001 temp(16) f002 feee feee ***
	//************************************************************************************************
	while(event_begin_flag && !event_end_flag && !f_in.eof())
	{
		while(!layer_begin_flag && !f_in.eof())
		{
			f_in.read((char*)(&buffer),1);
			buffer_v.push_back(buffer);
			if(buffer_v.size()>4)	buffer_v.erase(buffer_v.begin(),buffer_v.end()-4);
			if(buffer_v[0]==0xfa && buffer_v[1]==0x5a && buffer_v[2]==0xfa && buffer_v[3] == 0x5a && buffer_v.size()==4)	layer_begin_flag = 1;
		}
		// -------------	 find bag at layer level "fa5a fa5a"	-----------------
		while(layer_begin_flag && !layer_end_flag && !f_in.eof())
		{  
			f_in.read((char*)(&buffer),1);
			buffer_v.push_back(buffer);
			if(buffer_v.size()>=12 && *(buffer_v.end()-1)==0xee && *(buffer_v.end()-2)==0xfe && *(buffer_v.end()-3)==0xee && *(buffer_v.end()-4)==0xfe) // "fa5a fa5a cycleID feee feee" has 12 bytes 
				layer_end_flag = 1;
		}
		
		//dddddddddddddddddddddddddddddddddddddddddddddddddddd
		//for(auto it=buffer_v.begin(); it!=buffer_v.end(); ++it)	{cout<<hex << *it <<" " ;}
		//cout << endl;

		f_in.read((char*)(&buffer), 1);// read "ff"
		//cout << "buffer(ff) : " << buffer << endl;
		f_in.read((char*)(&buffer), 1);// read layer number
		//cout << "buffer(layerID) : " << buffer << endl;
		if(buffer<0 || buffer>31)	
		{
			cout << "******* Layer No is not in [0,31] *******" << endl;
			for(auto it=buffer_v.begin(); it!=buffer_v.end(); ++it)	{cout<<hex << *it <<" " ;} cout << endl;
			buffer_v.clear();
			layer_begin_flag = 0;
			layer_end_flag = 0;
			continue;
		}
		
		layerID = buffer;		
		
		if(layerID==0x1f)// event_time and event bag tail "fedd fedd "
		{
			event_time = 0;
			for(auto i_time : {0x1000000, 0x10000, 0x100, 0x1})
			{ 
				f_in.read((char*)(&buffer),1);
				if(buffer != 0)	event_time += buffer * i_time;
			}
			
			for(int read_time = 0 ; read_time < 4 ; read_time++)
			{
				f_in.read((char*)(&buffer),1);
				testbuffer.push_back(buffer);
			}
			
			/*while(testbuffer[0]!=0xfe || testbuffer[1]!=0xdd || testbuffer[2]!=0xfe || testbuffer[3]!= 0xdd)	
			{
				testbuffer.erase(testbuffer.begin());
				f_in.read((char*)(&buffer),1);
				testbuffer.push_back(buffer);
				//cout << "testbuffer_v : " << hex <<  testbuffer[0] << " " << testbuffer[1] << " " << testbuffer[2] << " " << testbuffer[3] << endl;
			}*/
			while(testbuffer[0]!=0xfe || testbuffer[1]!=0xdd || testbuffer[2]!=0xfe || testbuffer[3]!= 0xdd)	
			{
				cout << "Wrong event bag tail !!" << endl;
				cout << "testbuffer_v : " << hex <<  testbuffer[0] << " " << testbuffer[1] << " " << testbuffer[2] << " " << testbuffer[3] << endl;
				f_in.read((char*)(&buffer),1);
				testbuffer.push_back(buffer);
				testbuffer.erase(testbuffer.begin());
			}
			event_begin_flag = 1;
			event_end_flag = 1;
			//cout << "******* Find the event bag tail *******" << endl;
		}


		// empty bag && abnormal bag 
		if(buffer_v.size()<198) // 99 * 2 (2fa5a + 2cycl + 72+1+1 + 1trig + 1+16+1 + 2feee) 
		{	
			layer_begin_flag = 0;
			layer_end_flag = 0;
			if(buffer_v.size()!=12) 
			{
				//cout << "SMALL abnormal Bag first selected ! Skipping ..." << endl;
				//for(auto it=buffer_v.begin(); it!=buffer_v.end(); ++it)	{cout<<hex << *it <<" " ;}cout << endl;
				small_abnormal_bag_auto++;
			}
			else {empty_bag_auto++;}
			buffer_v.clear();
			layer_begin_flag = 0;
			layer_end_flag = 0;
			continue;
		}

		// -------------	merge 2 0x number of data bag	-----------------
		if(buffer_v.size()%2)
		{
			cout << "******* When merging buffer_v.size()%2 equal 1 *******" << endl;
			cout << "******* And layerID : " << hex << layerID << " *******" << endl;
			buffer_v.clear(); 
			layer_begin_flag = 0;
			layer_end_flag = 0;
			continue;
		}
		
		for(unsigned int i_size=0; i_size < buffer_v.size()/2; i_size++)
		{
			buffer_v[i_size] = buffer_v[i_size*2] * 0x100 + buffer_v[i_size*2+1];
		}
		buffer_v.erase(buffer_v.end()-buffer_v.size()/2,buffer_v.end());
	
		// extract cycleID form buffer_v and remove it from buffer_v
		cycleID = buffer_v[3] + buffer_v[2] * 0x10000;
		buffer_v.erase(buffer_v.begin() + 2, buffer_v.begin() + 4);
		// triggerID
		triggerID = *(buffer_v.end()-(2+16+3));
		buffer_v.erase(buffer_v.end()-(2+16+3));
		LayerTriggerID.push_back(triggerID);
		
		if(buffer_v.size()>=96) // "fa5a fa5a 36+36 1bcid 1chipid f001 16 f002 feee feee"	
		{
			if(*(buffer_v.end()-3)!=0xf002 && *(buffer_v.end()-20)!=0xf001)
			{
				cout << "******** Temperature Abormal********" << endl;
				cout << "******** triggerID : " << hex << triggerID << " *******" << endl;
				cout << "******** buffer size : " << dec << buffer_v.size() << " ********" << endl;
				std::for_each(buffer_v.begin(), buffer_v.end(), [](int elem) {printf("%x, ", elem); });cout << endl;
				temp_abnormal_bag_auto++;
				buffer_v.clear();
				layer_begin_flag = 0;
				layer_end_flag = 0;
				continue;
			}
		}
		/*else 
		{
			cout << "******** Small Abnormal Bag ********" << endl;
			cout << "******** And triggerID : " << hex << triggerID << " *******" << endl;
			cout << "******** And buffer size : " << dec << buffer_v.size() << " ********" << endl;
			std::for_each(buffer_v.begin(), buffer_v.end(), [](int elem) {printf("%x, ", elem); });cout << endl;
			buffer_v.clear();
			layer_begin_flag = 0;
			layer_end_flag = 0;
			continue;
		}*/
		///////////////////////////////////////////////////////////////////////////////////////////
		
		buffer_v.erase(buffer_v.begin(),buffer_v.begin()+2);	// erase "fa5a fa5a"
		buffer_v.erase(buffer_v.end()-2,buffer_v.end());	// erase "feee feee"
		
		buffer_v.erase(buffer_v.end()-(1+16+1)); // erase "f001"
		buffer_v.erase(buffer_v.end()-1);// erase "f002"
	
		//for(auto it=buffer_v.begin(); it!=buffer_v.end(); ++it)	{cout<<hex << *it <<" " ;}
		//cout << endl;
		//cout << "buffer_v.size " << dec << buffer_v.size() << endl;
		
		const int sci_data_len = 74;
		if(buffer_v.size()%sci_data_len!=16)
		{
			cout << "******* SCI DATA OR TEMP DATA Abnormal Bag  *******" << endl;
			cout << "******** triggerID : " << hex << triggerID << " *******" << endl;
			//dddddddddddddddddddddddddddddddddddddddddddddddddddd
			for(auto it=buffer_v.begin(); it!=buffer_v.end(); ++it)	{cout<<hex << *it <<" " ;}
			cout << endl;
			cout <<	"buffer_v.size " << dec << buffer_v.size() << endl;
			data_abnormal_bag_auto++;
			buffer_v.clear();
			layer_begin_flag = 0;
			layer_end_flag = 0;
			continue;
		}
		
		bool chipID_bad_flag=0;
		while(buffer_v.size()!=16)
		{
			int chipID = buffer_v[sci_data_len-1] - 1;
			if(chipID>5||chipID<0) 
			{
				chipID_bad_flag = 1;
				cout << "******* chipID is " << dec << chipID<< " not in [0,5], continue ********" << endl;
				cout << "buffer_v size : " << buffer_v.size() << endl;
				buffer_v.erase(buffer_v.begin(),buffer_v.begin()+sci_data_len);
				continue;
			}
			_chip_v[layerID][chipID].assign(buffer_v.begin(),buffer_v.begin()+sci_data_len);

			buffer_v.erase(buffer_v.begin(),buffer_v.begin()+sci_data_len);
		}
		_tmpr_v[layerID].assign(buffer_v.begin(), buffer_v.end()); // only temperature dat left
		
		if(chipID_bad_flag==1) 
			data_abnormal_bag_auto++;		
		else
			good_bag_auto++;
		buffer_v.clear();
		layer_end_flag = 0;
		layer_begin_flag = 0;
		

	} // end collect each layer
	//cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ " << endl;
	return 0;
}

int raw2RootAuto::DecodeAEvent(vector<int>& chip_v, vector<int>& tmpr_v, int layerID)
{
	int size = chip_v.size();
	//int Memo_No = size / channel_FEE;
	//int offset = 72 * (Memo_No - 1);// No Memo_ID , thus, offset = 0 
	//cout << "offset : " << offset << endl; 
	
	//for(auto it=chip_v.begin(); it!=chip_v.end(); ++it)	{cout<<hex << *it <<" " ;}
	//cout << endl;
	//cout <<	"chip_v.size " << dec << chip_v.size() << endl;
	
	for (int i_ch = 0; i_ch < channel_No; ++i_ch)  // 72 + BCID + Chip
	{
		////int tdc = chip_v.at(i_ch + offset) & 0x0fff;	
		////int adc = chip_v.at(i_ch + offset + channel_No) & 0x0fff;
		//int chanID = channel_No - 1 - i_ch;// chanID from 35 to 0
		//int gain = chip_v.at(i_ch + offset) & 0x2000; // bit operation
		//int hit = chip_v.at(i_ch + offset) & 0x1000;
		//int time = chip_v.at(i_ch + offset) & 0x0fff;
		//int high,low;
		int chanID = channel_No - 1 - i_ch;// chanID from 35 to 0
		int gain = chip_v.at(i_ch + channel_No) & 0x2000; // bit operation
		int hit = chip_v.at(i_ch + channel_No) & 0x1000;
		int gaintest = chip_v.at(i_ch) & 0x2000; // bit operation
		int hittest = chip_v.at(i_ch) & 0x1000;
		int time = chip_v.at(i_ch) & 0x0fff;
		int high,low;
		if(gain>0) 	
		{
			high = chip_v.at(i_ch + channel_No) & 0x0fff;
			low = -1;
		}
		else
		{
			low = chip_v.at(i_ch + channel_No) & 0x0fff;
			high = -1;
		}

		int chip = chip_v[size-1] - 1;
		int bcid = chip_v[size-2];

		_cellID.push_back(layerID * 1E5 + chip * 1E4 + 0 * 1E2 + chanID);// 0 is for memoryID
		_bcid.push_back(bcid);
		_hitTag.push_back(hit>0 ? 1 : 0 );
		_gainTag.push_back(gain>0 ? 1 : 0);
		_chargehi.push_back(high);
		_chargelo.push_back(low);
		_hit_time.push_back(time);
	}

	if(_temp_checkfirst[layerID]) 
	{
		for (int i_temp = 0; i_temp < temp_No; i_temp++) 
		{
			double tmpr = tmpr_v.at(i_temp) / 128.;
			_temp[layerID].push_back(tmpr);
		}
		_temp_checkfirst[layerID] = 0;
	}

	//chip_v.erase(chip_v.end() - 2, chip_v.end() - 1);
	//chip_v.erase(chip_v.begin() + offset, chip_v.begin() + offset + 72);
	chip_v.erase(chip_v.begin(), chip_v.begin() + 74);
	if(chip_v.size()) 
	{
		cout << "chip_v.size() != 0" << endl;
		chip_v.clear();
	}
	return 0;
}


int raw2RootAuto::AutoGainDecode(string input_file, string output_file)
{
	ifstream f_in;
	int cycleID, triggerID;
	//int BCID[layer_No][chip_No];
	//int Memo_ID[layer_No][chip_No];
	int preID = -1;
	//int layerID;
	int decode_begin_flag;
	
	trigger_period = 0; // When triggerID of DAQ turns back to one, trigger_period ++
	Event_No = 0;
	
	f_in.open(input_file, ios::in);
	if(!f_in)	{cout << "cant open " << input_file << endl; return 0;}
	// Init 
	/*for (int i_layer = 0; i_layer < layer_No; ++i_layer) 
	{
		if(_buffer_v.size()) _buffer_v.clear();
		for (int i_chip = 0; i_chip < chip_No; ++i_chip) 
		{
			//BCID[i_layer][i_chip] = -1;
			//Memo_ID[i_layer][i_chip] = 0;
			if(_chip_v[i_layer][i_chip].size())	_chip_v[i_layer][i_chip].clear();
		}
	}*/

	// TFile and TTree
	string str_out = output_file;
	TFile* fout =  TFile::Open(str_out.c_str(), "RECREATE");
	if (!fout) {cout << "cant create " << str_out << endl; return 0;}
	//str_out = str_out.substr(str_out.rfind("Run")+3,3);
	str_out.erase(str_out.begin(),str_out.begin()+str_out.find("Run")+3);
	str_out.erase(str_out.begin()+str_out.find("_"),str_out.end());
	run_no = atoi(str_out.c_str());
	//run_no = -1;
	cout << "******** Run Number is " << str_out << " ********"<< endl; 
	TTree* tree = new TTree("Raw_Hit", "data from binary file");
	tree->SetAutoSave(32 * 1024 * 102e4);
	SetTreeBranch(tree);
	
	for(int i_layer = 0; i_layer < layer_No; i_layer++)
	{
		_temp.push_back(vector<double>());
		_temp_checkfirst[i_layer] = 1;
	} 
	
	while (!(f_in.eof())) 
	{
		decode_begin_flag = 0;
		if(_buffer_v.size())	_buffer_v.clear();
		if(LayerTriggerID.size())	LayerTriggerID.clear();
		for(int i_layer = 0 ; i_layer < layer_No; i_layer++)
		{
			for(int i_chip = 0 ; i_chip < chip_No; i_chip++)
			{
				if(_chip_v[i_layer][i_chip].size())	{_chip_v[i_layer][i_chip].clear();}
			}
		}

		CatchSPIROCBag(f_in, _buffer_v, cycleID, triggerID);
		
		/*for(int i_layer = 0 ; i_layer < layer_No; i_layer++)
		{
			for(int i_chip = 0 ; i_chip < chip_No; i_chip++)
			{
				if(_chip_v[i_layer][i_chip].size())	{decode_begin_flag = 1;	break;}
			}
		}*/
		
		if(!LayerTriggerID.size())
			continue;
		else
		{
			decode_begin_flag = 1;
		}
		
		if(!IsSame(LayerTriggerID))
		{
			cout << "!!!!!!!!!!!!!!!!!! TriggerID NOT SAME !!!!!!!!!!!!!!!!!!!" << endl; 
			for(unsigned int i = 0 ; i< LayerTriggerID.size() ; i++)
			{
				cout << hex << LayerTriggerID.at(i) << " ";
			}
			cout << endl;
			cout << "DELETING ... ..." << endl; 
			cout << "NEXT EVENT ... ..." << endl;
			continue;
		}
				
		if(decode_begin_flag && IsSame(LayerTriggerID))
		{
			//cout << "==================== TriggerID SAME ===================" << endl; 
			if(preID == triggerID) 
			{
				cout << "triggerID of last event and this event is same !" << endl;
				//cout <<	"cycleID : " << hex << cycleID << endl; 
				cout << "preID : triggerID " << hex << preID  << " " << triggerID << endl;
				continue;
			}
			if(preID > triggerID)
			{
				if(abs(preID-triggerID)>30000 && preID > 0xf000)	
					trigger_period += 65536;
				else 
				{
					cout << "triggerID of last event large than this events'" << endl;
					//cout <<	"cycleID : " << cycleID << endl; 
					cout << "preID : triggerID " << hex << preID  << " " << triggerID << endl;
					continue;// some disorder of trigger ID 
				}
			}
			preID = triggerID;
			BranchClear();
			bool event_test_flag = false ; 
			for (int i_layer = 0; i_layer < layer_No; ++i_layer) 
			{
				_temp_checkfirst[i_layer] = 1;
				for (int i_chip = 0; i_chip < chip_No; ++i_chip) 
				{
					if(_chip_v[i_layer][i_chip].size() == 0)	continue;
					event_test_flag = true;
					_run_no = run_no;
					_event_time = event_time;
					_cycleID = cycleID;
					_triggerID = triggerID+ trigger_period;
					DecodeAEvent(_chip_v[i_layer][i_chip], _tmpr_v[i_layer], i_layer);
					//Memo_ID[i_layer][i_chip]--;
				}
			}
			if(event_test_flag) 
			{
				Event_No++;
				if(Event_No%5000==0)
				//if(1)
				{
					cout << "triggerID : " << hex << triggerID << endl;
					cout << "******* Event_No is " << dec <<  Event_No << " *******" << endl;
				}
			}
			tree->Fill();
			//if(Event_No==1) 
			if(0) 
				{tree->Write();fout->Write();fout->Close();return 1;}

	//		if (abs(tmp_ID - preID) > 30000) {trigger_period += 65536;}
		}
	}
	cout << "good		empty		small_abnormal		temp_abnormal		data_abnormal" << endl;
	cout << dec << good_bag_auto << "		" << empty_bag_auto << "		" << small_abnormal_bag_auto << "		" << temp_abnormal_bag_auto << "		" << data_abnormal_bag_auto << endl;

	f_in.close();
	tree->Write();
	fout->Write();
	fout->Close();
	return 0;
}

