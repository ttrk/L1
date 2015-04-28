#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void plotTurnOn()
{
  //const Int_t INPUTFILES = 11;
  const Int_t INPUTFILES = 2;

/*
  TString outFileTag = "Hydjet502_500Hz";
  TString inFileName[INPUTFILES] = {"hist_Hydjet502_twoByTwoANDzeroWalls.root",
  				    "hist_Hydjet502_twoByTwoANDzeroWallsANDsigmaSubtraction.root"};
  TString labels[INPUTFILES] = {"twoByTwoANDzeroWalls, L1_38.5, 499.487Hz",
				"twoByTwoANDzeroWallsANDsigmaSubtraction, L1_21.0, 486.777Hz"
  };
  const Float_t L1_THRESHOLD[INPUTFILES] = {38.5, 21.0};

  

  TString outFileTag = "Hydjet502_200Hz";
  TString inFileName[INPUTFILES] = {"hist_Hydjet502_twoByTwoANDzeroWalls.root",
  				    "hist_Hydjet502_twoByTwoANDzeroWallsANDsigmaSubtraction.root"};
  TString labels[INPUTFILES] = {"twoByTwoANDzeroWalls, L1_44, Hz 196.471",
				"twoByTwoANDzeroWallsANDsigmaSubtraction, L1_24.5, 208.2 Hz"
  };
  const Float_t L1_THRESHOLD[INPUTFILES] = {44,24.5};
*/
  TString outFileTag = "Hydjet502_100Hz";
  TString inFileName[INPUTFILES] = {"hist_Hydjet502_twoByTwoANDzeroWalls.root",
  				    "hist_Hydjet502_twoByTwoANDzeroWallsANDsigmaSubtraction.root"};
  TString labels[INPUTFILES] = {"twoByTwoANDzeroWalls, L1_49, 99.701 Hz",
				"twoByTwoANDzeroWallsANDsigmaSubtraction, L1_28, 104.589 Hz"
  };
  const Float_t L1_THRESHOLD[INPUTFILES] = {49,28};


  
  
  
  
  TFile* inFile[INPUTFILES];
  const Int_t COLORS[INPUTFILES] = {kViolet+1, kBlue};
  TGraphAsymmErrors *asymm[INPUTFILES];//[2];
  
  for(int i=0; i<INPUTFILES; i++)
  {
    inFile[i] = TFile::Open(inFileName[i]);
    //asymm[i] = (TGraphAsymmErrors*)inFile[i]->Get(Form("asymm_pt_%d_0",(int)L1_THRESHOLD[i]));
    asymm[i] = (TGraphAsymmErrors*)inFile[i]->Get(Form("asymm_pt_%.1f_0",L1_THRESHOLD[i]));
    asymm[i]->SetMarkerColor(COLORS[i]);
    asymm[i]->SetLineColor(COLORS[i]);
    asymm[i]->SetLineWidth(2);
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

  for(int i=0; i<INPUTFILES; i++)
  {
    asymm[i]->Draw("p");
  }

  TLegend *leg = new TLegend(0.55,0.2,0.9,0.5,"|#eta| < 2");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(18);

  for(int i=0; i<INPUTFILES; i++)
    {
      leg->AddEntry(asymm[i], Form("%s", labels[i].Data()), "lp");
    }
  leg->Draw();
  c1->SaveAs(Form("%s_turnon.pdf",outFileTag.Data()));
}

/*
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
*/
