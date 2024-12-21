#ifndef Calibration_h_
#define Calibration_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>
#include "Extract.h"

using namespace std;

class processCali
{
public:
	static const int LAYERNO = 30;
    static const int CHIPNO = 6;
    static const int CHANNO = 36;
	const double MIPEnergy = 0.305; // MeV
    
	int allHit;
	int validHit;
	double validHit2All;
	int hitCellNo;
	int hitLayerNo;
	vector<double> _XPOS;
	vector<double> _YPOS;
	vector<double> _ZPOS;
	vector<double> _HADC;
	vector<int>	   _Layer;
	vector<int>	   _Chip;
	vector<int>	   _Channel;
	void Init();
	void DataClear();
	void MuonTrackFit(int entry);
    
	void Reset();
	int TreeType(string inFileName);
	void SetAllBranchAddress(int treeType);
	void SetAllBranch();
	double MIPCorrection(double MPV);
    int Calibration(string inFileName, string outFileName, Extract& ex);
	bool HitCut(double xv, double yv, double slope, double intercept);

private:   
    int _hlCellID;
    double _slope,_intercept,_xmax;
       
	// for Event Energy Reconstruction
    TFile* rawfile;
    TTree* fNtuple2;
    int _run_num, _event_time, _triggerID;
    std::vector<int>*    _cellID=nullptr;
    std::vector<int>*    _bcID=nullptr;
    std::vector<int>*    _hitTag=nullptr;
    std::vector<int>*    _gainTag=nullptr;
    std::vector<double>* _hg_charge=nullptr;
    std::vector<double>* _lg_charge=nullptr;
    std::vector<double>* _hit_time=nullptr;
    std::vector<vector<double> >* _tempLayer=nullptr;

    // for Calibration
    TFile* inFile;
	TFile* outFile;
    TTree* inTree;
    TTree* outTree;
	int _run_num_new, _event_time_new, _event_num, _detectorID;
    std::vector<int> _cellID_new;
    std::vector<double> _hitADC;
    std::vector<double> _hitE;
    std::vector<double> _hitX;
    std::vector<double> _hitY;
    std::vector<double> _hitZ;
    std::vector<double> _temp;
    //std::vector<double> _hitTag_new;
    //std::vector<double> _hit_time_new;
    //std::vector<double> _hg_cali;
    //std::vector<double> _lg_cali;
};
#endif
