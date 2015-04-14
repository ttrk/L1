/*
 * detect L1_THRESHOLD pt for which turn-on curve becomes 100% at pT = 60 GeV
 *
 * Alex : For the single region, 2x2, and 3x3 seeds, we should show the turn-on curve that reaches 100% efficiency at photon pt = 60 GeV all on the same plot and state the rates. I can make the combo plot and I h
 * */
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>

#include <iostream>
#include <stdlib.h>     /* atof */

void plotTurnOn_detect100Efficency(TString inFileName, double pt_photon = 60);

/*
 * detect L1_THRESHOLD pt for which
 * */
void plotTurnOn_detect100Efficency(TString inFileName, double pt_photon /* = 60 */)
{
	std::cout << "inFileName = " << inFileName << std::endl;
	std::cout << "pt_photon  = " << pt_photon  << std::endl;

  TFile *inFile = TFile::Open(inFileName);

  Int_t THRESHOLDS = 80;
  if(inFileName.Contains("HydjetMB") 				||
		  inFileName.Contains("AllQCDPhoton30")  	||
		  inFileName.Contains("PyquenUnquenched")	)
  {
	  THRESHOLDS=60;
  }
  TGraphAsymmErrors *asymm[THRESHOLDS];//[2];

  int fNpoints;
  int first_threshold = 0;
  bool first_threshold_FOUND = false;
  for(int i = THRESHOLDS-1; i >=0 ; i--)	// we want the highest threshold value possible.
  {
	  if(first_threshold_FOUND) break;

      asymm[i]  = (TGraphAsymmErrors*)inFile->Get(Form("asymm_pt_%d_0",i));

      fNpoints   = asymm[i]->GetN();
      double* fX = asymm[i]->GetX();
      double* fY = asymm[i]->GetY();
      for(int j=fNpoints -1 ; j >= 0; j--)
      {
    	  if(fX[j] <= pt_photon)	// we passed the bin where fX[bin]=pt_photon
    	  {
    		  if (fY[j] >= 1)	// 100% efficieny at pt_photon
    		  {
    			  // this is the highest possible L1_THRESHOLD where efficiency is still 100% at or below pt_photon (60 GeV)
    			  // keep index "i"
    			  first_threshold = i;
    			  first_threshold_FOUND = true;
    			  std::cout << "j =" << j <<std::endl;
    			  std::cout << "fX[j] =" << fX[j] <<std::endl;

    		  }
    		  break;
    	  }
      }
  }

  std::cout << "first_threshold = " << first_threshold << std::endl;

 /*
  //      const Int_t THRESHOLDS_plot = 4;
  //  const Double_t L1_THRESHOLD[THRESHOLDS_plot] = {4, 12, 20, 29};
  //  const Int_t COLORS[THRESHOLDS_plot] = {kBlack, kRed, kBlue, kGreen+3};//, kMagenta+3};

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

  c1->SaveAs(Form("%s_turnon.pdf",outFileTag.Data()));
*/
}

int main(int argc, char **argv)
{
  if(argc == 4)
  {
	plotTurnOn_detect100Efficency(argv[1], atof(argv[2]));
    return 0;
  }
  else
  {
	  std::cout << "wrong usage" << std::endl;
//    std::cout << "Usage:\plotTurnOn_detect100Efficency.exe <input_filename> <output_tag>" << std::cout;
//    std::cout << "An output pdf will be named <output_tag>_turnon.pdf" << std::cout;
    return 1;
  }
}
