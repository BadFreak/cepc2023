#ifndef TreeCalibration_h_
#define TreeCalibration_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>
#include"Extract.h"

using namespace std;

class MIPCalibration
{
    // for draw some histrogram
    private:
		// from auto gain mode decode 
    	Int_t _run_no, _event_time, _cycleID, _triggerID;
    	vector< int >* _cellID=nullptr;
    	vector< int >* _bcid=nullptr;
    	vector< int >* _hitTag=nullptr;
    	vector< int >* _gainTag=nullptr;
    	vector< double >* _chargehi=nullptr;
    	vector< double >* _chargelo=nullptr;
    	vector< double >* _hit_time=nullptr;
    	std::vector<vector<double> >* _temp=nullptr;

		double _landauMPV, _landauWidth, _gausSigma, _ChiSqr;
		int _Ndf, _cellID_mip, _chn_entries;
		TFile * inFile;
		TTree * inTree;
		TFile * mipFile;
		TTree * mipTree;

    private:
		static const int layerNu = 32;
		static const int chipNu = 6;
		static const int channelNu = 36;
		static const int memoryNu = 16;
		
		// for fitting parameters
		TFile* pedFile;
		TTree* fNtuple1;
		std::vector<int>* _pedCellID=nullptr;
		std::vector<float>* _pedMeanTimes=nullptr;
		//double pedMeanTimesChnLevel[layerNu][chipNu][channelNu]={0};
		vector<double> _XPOS;
		vector<double> _YPOS;
		vector<double> _ZPOS;
		vector<double> _HADC;
		vector<int>	   _Layer;
		vector<int>	   _Chip;
		vector<int>	   _Channel;

	public:
		int hitLayerNo;
		int hitCellNo;
		double validHit2All;
		int allHit, validHit;
    	
		MIPCalibration(){};
		int TreeType(string inFileName);
		void MuonTrackFit(int entry);
    	int MIPFit(string input_name,Extract ex,string hist_name,string mip_name, bool isTrackFit);
    	//int MIPFit(string input_name,string pedfactor_name,string hist_name,string mip_name);
		void ReadTree(int treeType);
		void SetTree();
		void DataClear();
		void ForceSetBranchValue(int layer, int chip, int channel);
		//void Init();
		//void Reset();
		void PedestalExtract(string pedestal_name);
};

#endif
