#ifndef _MUONTRACKFIT_H
#define	_MUONTRACKFIT_H
#include <iostream>
#include <vector>

using namespace std;

class MuonTrackFit
{
	public:
		MuonTrackFit(vector<double> xpos, vector<double> ypos, vector<double> zpos)
		{
			_XPOS.assign(xpos.begin(), xpos.end());
			_YPOS.assign(ypos.begin(), ypos.end());
			_ZPOS.assign(zpos.begin(), zpos.end());
		};
		void MuonTrackFitXZ();
		void MuonTrackFitYZ();

	private:
		vector<double> 	_XPOS;
		vector<double> 	_YPOS;
		vector<double> 	_ZPOS;
	//	vector<int>		_HADC;

};

#endif
