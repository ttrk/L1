#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void DrawRateAtFullEfficiency(int centrality=0)
{
	const Int_t INPUTFILES = 2;

	int sampleIndex = 1;
	TString sample[2] = {"Hydjet502Dijet30", "Hydjet502Dijet80"};
	TString outFileTag = Form("RateAtFullEfficiency_%s", sample[sampleIndex].Data());
	TString inFileName[INPUTFILES] = {"twoByTwoANDzeroWalls", "twoByTwoANDzeroWallsANDsigmaSubtraction"};
	TString dirName[INPUTFILES] = {"/export/d00/scratch/tatar/output/out_L1EmulatorMacros_v4", "/export/d00/scratch/tatar/output/out_L1EmulatorMacros_vFinerBin"};

	TFile *inFile[INPUTFILES];
	const Int_t COLORS[INPUTFILES] = {kRed, kBlue};		// {kBlack, kRed, kBlue};
	TGraphAsymmErrors *asymm[INPUTFILES];

	for(int i = 0; i < INPUTFILES; i++){
//		inFile[i] = TFile::Open(Form("results/filerate_Hydjet502_%s_cent%d.root",inFileName[i].Data(),centrality));
		inFile[i] = TFile::Open(Form("%s/filerate_%s_%s_cent%d.root",dirName[i].Data(),sample[sampleIndex].Data(),inFileName[i].Data(),centrality));
		asymm[i] = (TGraphAsymmErrors*)inFile[i]->Get("Graph");
		asymm[i]->SetMarkerColor(COLORS[i]);
		asymm[i]->SetLineColor(COLORS[i]);
		asymm[i]->SetLineWidth(2);
	}

	TH1D *hEmpty = new TH1D("hEmpty",Form(";Offline Jet p_{T} (GeV);rate at full efficiency"),100,0,150);

	TCanvas *c1 = new TCanvas();
	hEmpty->SetMinimum(0.1);
	hEmpty->SetMaximum(400000);
	hEmpty->Draw();
	c1->SetLogy();

	for(int i = 0; i < INPUTFILES; i++)
	{
		asymm[i]->Draw("lp");
	}

	TLegend *leg = new TLegend(0.55,0.2,0.9,0.5,"|#eta| < 2");
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetTextSizePixels(18);

	for(int i = 0; i < INPUTFILES; i++)
	{ 
		leg->AddEntry(asymm[i], Form("%s", inFileName[i].Data()), "lp");
	}

	leg->Draw();

	c1->SaveAs(Form("%s/%s_cent%d.pdf",dirName[1].Data(),outFileTag.Data(),centrality));	// save results to the "*vFinerBin" directory
	c1->SaveAs(Form("%s/%s_cent%d.root",dirName[1].Data(),outFileTag.Data(),centrality));

}
