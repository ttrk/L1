#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void DrawRateAtFullEfficiency(int centrality=0)
{
	const Int_t INPUTFILES = 3;


	TString outFileTag = "RateAtFullEfficiency";
	TString inFileName[INPUTFILES] = {"zeroWalls","twoByTwoANDzeroWalls", "twoByTwoANDzeroWallsANDsigmaSubtraction"};

	TFile *inFile[INPUTFILES];
	const Int_t COLORS[INPUTFILES] = {kBlack, kRed, kBlue};
	TGraphAsymmErrors *asymm[INPUTFILES];

	for(int i = 0; i < INPUTFILES; i++){
		inFile[i] = TFile::Open(Form("results/filerate_Hydjet502_%s_cent%d.root",inFileName[i].Data(),centrality));
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

	c1->SaveAs(Form("%s_cent%d.pdf",outFileTag.Data(),centrality));
	c1->SaveAs(Form("%s_cent%d.root",outFileTag.Data(),centrality));

}
