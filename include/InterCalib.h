#ifndef TreeSelectRoot_h_
#define TreeSelectRoot_h_

#include<TROOT.h>
#include<TChain.h>
#include<TFile.h>
#include<TH1.h>
#include<TH2.h>
#include<TF1.h>
#include<TCanvas.h>
#include<vector>
#include<TGraph.h>
#include"Extract.h"

using namespace std;

class interCalib {
    public:
    void Init();
	void PedestalExtract(string pedestal_name);
	void RotatePoint(double &xv, double &yv, double slopek); 
	void FillRotateHist(string root_name, Extract* &ex); 
    interCalib(string root_name,string spshist_name,string spsfit_name,Extract *&ex);

    static const int chipNu = 6;
    static const int channelNu = 36;
    static const int memoryNu = 2;
    static const int LayerNo=32;
	void HistCut(TH2F*&hist2d, double slope, double intercept);

    protected:
    // Channels Level
    TObjArray *Hlist;
    TObjArray *Glist = new TObjArray(210);
    TGraph *gp[LayerNo][chipNu][channelNu];
	TH2F* hisg[LayerNo][chipNu][channelNu];
	TH2F* hisgRo[LayerNo][chipNu][channelNu];
	TGraph* grap[LayerNo][chipNu][channelNu];
    int eventNo[LayerNo][chipNu][channelNu];
	int Xmax[LayerNo][chipNu][channelNu];

    int _hlCellID, _Npeaks, _DAC, _statistics, _badHitNo;
    double _Interval;
    vector<double> _spsPar;
    double _slope, _slopeError;
    double _intercept, _interceptError;
    double _chiSquare,_NDF,_xmax;

	int _event_time , _triggerID;
   	std::vector<int>* _cellID=nullptr;
   	std::vector<int>* _hitTag=nullptr;
    std::vector<double>* _lg_charge=nullptr;
    std::vector<double>* _hg_charge=nullptr;
    std::vector<vector<double> >* _tempLayer=nullptr;

    private:
    TFile* rawfile;
	double slope[LayerNo][chipNu][channelNu] = {0};
	double intercept[LayerNo][chipNu][channelNu] = {0};
	double fitLeftLimits[LayerNo][chipNu][channelNu] = {0};
	double fitRightLimits[LayerNo][chipNu][channelNu] = {0};
	double maxLowDrift[LayerNo][chipNu][channelNu] = {0};


};

#endif
