/*
 * find L1_THRESHOLD pt for which turn-on curve becomes 100% at pT = 60 GeV
 *
 * Alex : For the single region, 2x2, and 3x3 seeds, we should show the turn-on curve that reaches 100% efficiency at photon pt = 60 GeV all on the same plot and state the rates. I can make the combo plot and I h
 * */
#include <TFile.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <stdlib.h>     /* atof */

void plotTurnOn_detect100Efficency(TString inFileName, double pt_photon = 60);

/*
 * find the first (highest) L1_THRESHOLD pt for which the efficiency becomes at or before "pt_photon"
 * */
void plotTurnOn_detect100Efficency(TString inFileName, double pt_photon /* = 60 */)
{
	std::cout << "inFileName = " << inFileName << std::endl;
	std::cout << "pt_photon  = " << pt_photon  << std::endl;

  TFile *inFile = TFile::Open(inFileName);

  Int_t THRESHOLDS = 80;

  // an ugly way to determine the number of "THRESHOLDS".
  // usually histogram files for MC samples have 60 thresholds.  see "makeTurnOn_fromSameFile_photons.C"
  // usually histogram files for data samples have 80 thresholds.  see "makeTurnOn_fromSameFile_photons.C"
  if(inFileName.Contains("HydjetMB") 			||
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
}

int main(int argc, char **argv)
{
	if(argc == 2)
	{
		plotTurnOn_detect100Efficency(argv[1]);
		return 0;
	}
	else if(argc == 3)
	{
		plotTurnOn_detect100Efficency(argv[1], atof(argv[2]));
		return 0;
	}
	else
	{
		std::cout << "wrong usage" << std::endl;
		return 1;
	}
}
