#ifndef _EXTRACT_H
#define _EXTRACT_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
// #include "MCDigi.h"

using namespace std;

class Extract
{
	public:
		static constexpr const double MIPEnergy = 0.305;
		static const int LAYERNO = 30;
		static const int CHIPNO = 6;
		static const int CHANNO = 36;
		static const int CELLNO = 6300;

		int pedDeadChannel;
		int ratioDeadChannel; 
		int mipDeadChannel;
		int totalDeadChannel;
		
		double pedMeanHigh[LAYERNO][CHIPNO][CHANNO]; 
		double pedSigmaHigh[LAYERNO][CHIPNO][CHANNO]; 
		double pedMeanLow[LAYERNO][CHIPNO][CHANNO]; 
		double pedSigmaLow[LAYERNO][CHIPNO][CHANNO]; 

		double HLRatio[LAYERNO][CHIPNO][CHANNO];
		double HLIntercept[LAYERNO][CHIPNO][CHANNO];
		double HighADCSwitch[LAYERNO][CHIPNO][CHANNO];
		
		double MIPPeakValue[LAYERNO][CHIPNO][CHANNO]; 
		double SiPMGain[LAYERNO][CHIPNO][CHANNO]; 
		double SiPMGainError[LAYERNO][CHIPNO][CHANNO];

		bool pedOn[LAYERNO][CHIPNO][CHANNO];
		bool ratioOn[LAYERNO][CHIPNO][CHANNO];
		bool mipOn[LAYERNO][CHIPNO][CHANNO];
		bool channelOn[LAYERNO][CHIPNO][CHANNO];
	
	public:
		Extract(){
			for (int il=0; il<LAYERNO; il++) {
				for (int ip=0; ip<CHIPNO; ip++) {
					for (int ic=0; ic<CHANNO; ic++) {
						pedOn[il][ip][ic]=1;
						ratioOn[il][ip][ic]=1;
						mipOn[il][ip][ic]=1;
						channelOn[il][ip][ic]=1;
					}
				}
			}
		}
		void PedestalExtract(string inPedFileName);
		void HLGainRatioExtract(string inHLGainRatioFileName);
		void MIPPeakValueExtract(string inMIPFileName);
		void SiPMGainExtract(string SiPMGainFileName);
		void ReadPedDeadChannel();
		void ReadRatioDeadChannel();
		void ReadMIPDeadChannel();
		void ReadAllDeadChannel();
	
};

#endif
