#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void plotTurnOn(TString inFileName, TString outFileTag)
{
  TFile *inFile = TFile::Open(inFileName);

  const Int_t THRESHOLDS = 4;
  const Double_t L1_THRESHOLD[THRESHOLDS] = {4, 12, 20, 29};
  const Int_t COLORS[THRESHOLDS] = {kBlack, kRed, kBlue, kGreen+3};//, kMagenta+3};
  TGraphAsymmErrors *asymm[THRESHOLDS];//[2];

  for(int i = 0; i < THRESHOLDS; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      asymm[i] = (TGraphAsymmErrors*)inFile->Get(Form("asymm_pt_%d_0",(int)L1_THRESHOLD[i]));
      asymm[i]->SetMarkerColor(COLORS[i]);
      asymm[i]->SetLineColor(COLORS[i]);
    }
    //asymm[i][1]->SetMarkerStyle(25);
  }

  // these values MUST MATCH those used in makeTurnOn.C
  const int nBins = 200;
  const double maxPt = 200;

  TH1D *hEmpty = new TH1D("hEmpty",Form(";Photon p_{T} (GeV);Efficiency"),nBins,0,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  TLine *line = new TLine(0,1,maxPt,1);
  line->Draw();

  for(int i = 0; i < THRESHOLDS; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      asymm[i]->Draw("p");
    }
  }

  TLegend *leg = new TLegend(0.55,0.2,0.9,0.5,"|#eta| < 2");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(18);

  for(int i = 0; i < THRESHOLDS; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      leg->AddEntry(asymm[i], Form("L1_BarrelRegion%d", (int)L1_THRESHOLD[i]), "lp");
    }
  }

  leg->Draw();

  c1->SaveAs(Form("%s_turnon.pdf",outFileTag.Data()));

  //////// Kaya's modificiation ////////
  //  plot L1-Offline correlation
  TH2D* corr = (TH2D*)inFile->Get("corr");
  TCanvas* c2 = new TCanvas();
  corr->Draw("COLZ");
  gPad->SetLogz();
  c2->SaveAs(Form("%s_corr.pdf",outFileTag.Data()));
  //////// Kaya's modificiation - END ////////
}

int main(int argc, char **argv)
{
  if(argc == 3)
  {
    plotTurnOn(argv[1], argv[2]);
    return 0;
  }
  else
  {
    std::cout << "Usage:\nplotTurnOn.exe <input_filename> <output_tag>" << std::cout;
    std::cout << "An output pdf will be named <output_tag>_turnon.pdf" << std::cout;
    return 1;
  }
}
