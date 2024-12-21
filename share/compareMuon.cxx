void compareMuon()
{
	TFile *f1 = TFile::Open("all_auto_muon_v4.root","READ");
	TTree *t1 = (TTree*)f1->Get("MIP_Fit");
	
	TFile *f2 = TFile::Open("all_auto_muon_v3_2023.root","READ");
	TTree *t2 = (TTree*)f2->Get("MIP_Fit");

	t1->AddFriend(t2);
	t1->Draw("t2.LandauMPV");

}
