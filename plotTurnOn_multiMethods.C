#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void plotTurnOn()
{
  const Int_t INPUTFILES = 4;


  // TString outFileTag = "MBData_200Hz";
  // TString inFileName[INPUTFILES] = {"hist_MBData_zeroWalls.root",
  // 				    "hist_MBData_sigmaSubtraction.root",
  // 				    "hist_MBData_oneByOneAndzeroWalls.root",
  // 				    "hist_MBData_twoByTwoANDzeroWalls.root"
  // };
  // TString labels[INPUTFILES] = {"zeroWalls, L1_40, 206Hz",
  // 				"sigma subtracted, L1_24, 187Hz",
  // 				"1x1 jets +ZW, L1_16, 191Hz",
  // 				"2x2 jets +ZW, L1_52, 228Hz"};
  // const Double_t L1_THRESHOLD[INPUTFILES] = {40, 24, 16, 32};
  // // TString labels[INPUTFILES] = {"zeroWalls, L1_44, 550Hz",
  // // 				"sigma subtracted, L1_20, 505Hz",
  // // 				"1x1 jets +ZW, L1_20, 441Hz",
  // // 				"2x2 jets +ZW, L1_36, 444Hz"};
  // // const Double_t L1_THRESHOLD[INPUTFILES] = {44, 20, 20, 36};

  // TString outFileTag = "Hydjet276_270Hz";
  // TString inFileName[INPUTFILES] = {"hist_Hydjet276_zeroWalls.root",
  // 				    "hist_Hydjet276_sigmaSubtraction.root",
  // 				    "hist_Hydjet276_oneByOneAndzeroWalls.root",
  // 				    "hist_Hydjet276_twoByTwoANDzeroWalls.root"};
  // // const TString labels[INPUTFILES] = {"zeroWalls, L1_60, 104Hz",
  // // 				      "sigma subtracted, L1_28, 102Hz",
  // // 				      "1x1 jets, L1_28, 112Hz",
  // // 				      "2x2 jets, L1_48, 107Hz"};
  // // const Double_t L1_THRESHOLD[INPUTFILES] = {60, 28, 28, 48};
  // const TString labels[INPUTFILES] = {"zeroWalls, L1_40, 276Hz",
  // 				      "sigma subtracted, L1_24, 280Hz",
  // 				      "1x1 jets +ZW, L1_20, 270Hz",
  // 				      "2x2 jets +ZW, L1_32, 272Hz"};
  // const Double_t L1_THRESHOLD[INPUTFILES] = {40, 24, 20, 32};


  TString outFileTag = "Hydjet502_200Hz";
  TString inFileName[INPUTFILES] = {"hist_Hydjet502_zeroWalls.root",
  				    "hist_Hydjet502_sigmaSubtraction.root",
  				    "hist_Hydjet502_oneByOneAndzeroWalls.root",
  				    "hist_Hydjet502_twoByTwoANDzeroWalls.root"
  };
  // // TString labels[INPUTFILES] = {"zeroWalls, L1_80, 150Hz",
  // // 				"sigma subtracted, L1_36, 163Hz",
  // // 				"1x1 jets, L1_36, 167Hz",
  // // 				"2x2 jets, L1_88, 155Hz"};
  // // const Double_t L1_THRESHOLD[INPUTFILES] = {80, 36, 36, 64};
  TString labels[INPUTFILES] = {"zeroWalls, L1_56, 164Hz",
  				"sigma subtracted, L1_36, 162Hz",
  				"1x1 jets +ZW, L1_28, 154Hz",
  				"2x2 jets +ZW, L1_44, 195Hz"};
  const Double_t L1_THRESHOLD[INPUTFILES] = {56, 36, 28, 44};

  TFile *inFile[INPUTFILES];
  const Int_t COLORS[INPUTFILES] = {kBlack, kRed, kBlue, kGreen+3};//, kMagenta+3};
  TGraphAsymmErrors *asymm[INPUTFILES];//[2];

  for(int i = 0; i < INPUTFILES; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      inFile[i] = TFile::Open(inFileName[i]);
      asymm[i] = (TGraphAsymmErrors*)inFile[i]->Get(Form("asymm_pt_%d_0",(int)L1_THRESHOLD[i]));
      asymm[i]->SetMarkerColor(COLORS[i]);
      asymm[i]->SetLineColor(COLORS[i]);
    }
    //asymm[i][1]->SetMarkerStyle(25);
  }

  // these values MUST MATCH those used in makeTurnOn.C
  const int nBins = 75;
  const double maxPt = 300;

  TH1D *hEmpty = new TH1D("hEmpty",Form(";Jet p_{T} (GeV);Efficiency"),nBins,0,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  TLine *line = new TLine(0,1,maxPt,1);
  line->Draw();

  for(int i = 0; i < INPUTFILES; i++)
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

  for(int i = 0; i < INPUTFILES; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      leg->AddEntry(asymm[i], Form("%s", labels[i].Data()), "lp");
    }
  }

  leg->Draw();

  c1->SaveAs(Form("%s_turnon.pdf",outFileTag.Data()));

}

int main(int argc, char **argv)
{
  if(argc == 1)
  {
    plotTurnOn();
    return 0;
  }
  else
  {
    std::cout << "Usage:\nplotTurnOn_multiMethods.exe <output_tag>" << std::cout;
    std::cout << "An output pdf will be named <output_tag>_turnon.pdf" << std::cout;
    return 1;
  }
}
