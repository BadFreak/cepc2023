#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

void plot_LandauMPV() {
    // 打开两个 ROOT 文件
    TFile *file1 = new TFile("../share/20201216_20210308_MIPdeclineCor_Update.root");
    TFile *file2 = new TFile("../share/all_auto_muon_v4_Prime_update.root");

    // 读取 TTree
    TTree *tree1 = (TTree*)file1->Get("MIP_Fit");
    TTree *tree2 = (TTree*)file2->Get("MIP_Fit");

    // 创建 TH1D 直方图
    TH1D *h1 = new TH1D("h1", "LandauMPV Comparison;LandauMPV;Counts", 200, 0, 1000);
    TH1D *h2 = new TH1D("h2", "LandauMPV Comparison;LandauMPV;Counts", 200, 0, 1000);
	TH1D *h3 = new TH1D("h3", "Residual: MIP_{CosmicRay}-MIP_{TestBeam}", 100, -100, 500);
	TH1D *h4 = new TH1D("h4", "15 um chns", 100, -100, 500);
	TH1D *h5 = new TH1D("h5", "10 um chns", 100, -100, 500);

    // 填充直方图
    tree1->Draw("LandauMPV>>h1");
    tree2->Draw("LandauMPV>>h2");
	tree1->AddFriend(tree2,"tree2");
	
    int CellID1, CellID2;
    double LandauMPV1, LandauMPV2;

    tree1->SetBranchAddress("CellID", &CellID1);
    tree1->SetBranchAddress("LandauMPV", &LandauMPV1);
    tree2->SetBranchAddress("CellID", &CellID2);
    tree2->SetBranchAddress("LandauMPV", &LandauMPV2);

    // 使用 map 存储 tree2 的 CellID -> LandauMPV
    std::map<int, double> cellMap;
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; i++) {
        tree2->GetEntry(i);
        if (LandauMPV2 >0)
			cellMap[CellID2] = LandauMPV2;
	}

    // 遍历 tree1 并匹配
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; i++) {
        tree1->GetEntry(i);
        if (cellMap.find(CellID1) != cellMap.end()) {
            double diff = LandauMPV1 - cellMap[CellID1];
            h3->Fill(diff);
			if (CellID1 >= 400000 && CellID1 < 2800000){
				h4->Fill(diff);
			} else {
				h5->Fill(diff);
			}
        }
    }

    // 设置不同颜色
	h1->GetXaxis()->SetTitle("LandauMPV[ADC]");
    h1->SetLineColor(kBlack);
	h1->SetLineWidth(3);
	h1->SetLineStyle(2);
    h2->SetLineColor(kBlue);
	h2->SetLineWidth(2);
	h2->SetFillStyle(3002);
	h2->SetFillColor(kBlue);

    // 创建画布并绘制
    TCanvas *c1 = new TCanvas("c1", "LandauMPV Comparison", 800, 600);
	gPad->SetLogy();
    h1->Draw();
    h2->Draw("SAME");

    // 添加图例
    auto legend = new TLegend(0.62, 0.68, 0.85, 0.85);
    legend->AddEntry(h1, "Cosmic ray", "l");
    legend->AddEntry(h2, "Test beam", "l");
    legend->Draw();

	h3->GetXaxis()->SetRangeUser(-100, 500);
	h3->SetLineColor(kBlue);
	h3->SetLineWidth(4);
	h3->SetFillStyle(3002);
	h3->SetFillColor(kBlue);
	h3->GetXaxis()->SetTitle("MPV_Diff[ADC]");
	h3->GetYaxis()->SetTitle("Counts");
	h4->SetLineColor(kRed);
	h4->SetLineWidth(2);
	h5->SetLineColor(kGreen);
	h5->SetLineWidth(2);
	
	TCanvas *c2 = new TCanvas("c2", "Residual: MIP_{CosmicRay}-MIP_{TestBeam} ", 800, 600);
	c2->cd();
	gPad->SetLogy();
	h3->SetStats(0);
	h3->Draw();
	h4->Draw("SAME");
	h5->Draw("SAME");

	auto Legend2 = new TLegend(0.62, 0.68, 0.85, 0.85);
    Legend2->AddEntry(h3, "All  chns", "l");
    Legend2->AddEntry(h4, "10um chns", "l");
    Legend2->AddEntry(h5, "15um chns", "l");
    Legend2->Draw();

    // 保存图片
    c1->SaveAs("LandauMPV_Comparison.png");
	c2->SaveAs("LandauMPV_Comparison_diff.png");
    // 关闭文件
//     file1->Close();
//     file2->Close();
 }

