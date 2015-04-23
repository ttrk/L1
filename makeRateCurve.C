#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>

void makeRateCurve(TString inFileName = "Hydjet502_JetResults.root", bool secondaryFile = false)
{
  std::cout << "Analyzing " << inFileName << std::endl;
  TH1::SetDefaultSumw2();

  const int MAXJETS = 8;
  const int nBins = 300;
  const int maxPt = 300; // make sure that maxPt/nBins = 4.

  TFile *inFile = TFile::Open(inFileName);
  TTree *inTree;
  if(secondaryFile)
    inTree = (TTree*)inFile->Get("L1UpgradeTree");
  else
    inTree = (TTree*)inFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_pt[MAXJETS];
  inTree->SetBranchAddress("jet_pt",l1_pt);

  TH1D *counts = new TH1D("counts","counts;Leading L1 Jet p_{T};Count",nBins,0,maxPt);

  long long entries = inTree->GetEntries();
  for(long long i = 0; i < entries; ++i)
  {
    inTree->GetEntry(i);

    double maxl1pt = 0;
    double maxCenJet = l1_pt[0];
    double maxForJet = l1_pt[4];
    maxl1pt = std::max(maxCenJet, maxForJet);

    counts->Fill(maxl1pt);
  }

  //TCanvas *c0 = new TCanvas();
  //counts->Draw("e");
  //c0->SetLogy();

  TH1D *rate;
  rate = new TH1D("rate",";L1 p_{T};Rate",nBins,0,maxPt);
  double total_integral = counts->Integral();

  std::cout << "Trigger Value \t Rate @ 30kHz" << std::endl;
  for(int i = 0; i < nBins; i++)
  {
    double j = (double)i*(double)maxPt/(double)nBins;
    double integral = counts->Integral(i+1, nBins);
    rate->Fill(j, (double)integral/total_integral);
    std::cout << "L1_SingleJet" << j << "\t" << integral/total_integral*30000 << std::endl;
    //std::cout << integral/total_integral*30000 << std::endl;
  }

  //TCanvas *c1 = new TCanvas();
  //rate->Draw("hist");
  //c1->SetLogy();
}

int main(int argc, char **argv)
{
  if(argc == 1){
    makeRateCurve();
    return 0;
  }
  else if ( argc == 2 )
  {
    makeRateCurve(argv[1]);
    return 0;
  }
  else if ( argc == 3 )
  {
    makeRateCurve(argv[1], atoi(argv[2]));
  }
  else
  {
    std::cout << "Usage: \nmakeRateCurve.exe <input_file_name>" << std::endl;
    return 1;
  }
}
