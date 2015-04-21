#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <iostream>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <cmath>
#include <iostream>
using namespace std;


#include <TH1D.h>
#include <TCanvas.h>
#include <TFrame.h>

#define BIN_NUM 40;
const int MAXJETS = 8;
const int nBins = 64;
const int maxPt = 256; // make sure that maxPt/nBins = 4.


int find(TString infname, Double_t pTthes, Double_t effthes, int cent)
{
  TFile* inf = new TFile(infname);
  int i=0,j=0;
  Bool_t flag=false;
  TString ingname;
  for(i=116;i>=0;i-=4)
    {
      ingname = Form("asymm_pt_%i_%d",i,cent);
      TGraphAsymmErrors* ga = (TGraphAsymmErrors*)inf->Get(ingname);
      if(!ga) break;
      Double_t vx,vy,interx,intermin=1000000.;
      Bool_t flagx=false;
      for(j=0;j<ga->GetN();j++)
	{
	  ga->GetPoint(j,vx,vy);
	  if(vx==pTthes)
	    {
	      //cout<<vy<<endl;
	      flagx=true;
	      if(vy>=effthes)
		{
		  flag=true;
		}
	      break;
	    }
	  if(TMath::Abs(vx-pTthes)<=intermin)
	    {
	      intermin = TMath::Abs(vx-pTthes);
	      interx = vx;
	    }
	}
      if(!flagx)
	{
	  cout<<endl;
	  cout<<">>>> ERROR"<<endl;
	  cout<<">>>> Graph <"<<ingname<<"> has no point at "<<pTthes<<"GeV"<<endl;
	  cout<<">>>> The closest point is"<<interx<<endl;;
	  cout<<">>>> ERROR ENDS"<<endl;
	  cout<<endl;
	  return -1;
	}
      if(flag) break;
    }
  if(!flag)
    {
      cout<<endl;
      cout<<">>>> ERROR"<<endl;
      cout<<"ERROR: File<"<<infname<<"> has no thredshold for "<<effthes*100<<"% at "<<pTthes<<" GeV/c"<<endl;
      cout<<">>>> ERROR ENDS"<<endl;
      cout<<endl;
	  return -1;

    }
  else
    {
      //cout<<endl;
      //cout<<">>>> RESULT"<<endl;
      //cout<<">>>> File <"<<infname<<"> has the thredshold <"<<ingname<<"> for "<<effthes*100<<"% at "<<pTthes<<" GeV/c"<<endl;
      //cout<<">>>> RESULT ENDS"<<endl;
      //cout<<endl;
      return i;
    }
  
}


void findthes(TString inFileName = "Hydjet502_JetResults_zeroWalls.root",TString infn = "hist_Hydjet502_zeroWalls.root",TString outfile = "rate_Hydjet502_zeroWalls", int centrality=0)
{
  TH1::SetDefaultSumw2();

  TFile *inFile = TFile::Open(inFileName);
  TTree *inTree;
  inTree = (TTree*)inFile->Get("L1UpgradeTree");

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
  
  TH1D *rate;
  rate = new TH1D("rate",";L1 p_{T};Rate",nBins,0,maxPt);
  double total_integral = counts->Integral();

  std::cout << "Trigger Value \t Rate @ 30kHz" << std::endl;
  for(int i = 0; i < nBins; i++)
  {
    double j = (double)i*(double)maxPt/(double)nBins;
    double integral = counts->Integral(i+1, nBins);
    rate->Fill(j, (double)integral/total_integral);
    //std::cout << "L1_SingleJet" << j << "\t" << integral/total_integral*30000 << std::endl;
  }

  const int Nthresholds=11;
  double offlinethresholds[Nthresholds]={26.,34.,42.,50.,62.,74.,86.,98.,110.,122.,130.};
  int L1thresholds[Nthresholds]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  double rates[Nthresholds]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  
  for(int m=0; m<Nthresholds;m++){
    L1thresholds[m]=find(infn, offlinethresholds[m], 1.,centrality);
    std::cout<<"threshold"<<L1thresholds[m]<<std::endl;
    rates[m]=rate->GetBinContent(int(L1thresholds[m]/4)+1)*30000;
  }  
   TCanvas* c1 = new TCanvas("c1","A Simple Graph with assymetric error bars",200,10,700,500);
   c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   
   Double_t exl[Nthresholds] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t eyl[Nthresholds] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t exh[Nthresholds] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   Double_t eyh[Nthresholds] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
   TGraphAsymmErrors *gr = new TGraphAsymmErrors(Nthresholds,offlinethresholds,rates,exl,exh,eyl,eyh);
   gr->SetTitle("TGraphAsymmErrors Example");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");
   
   TFile*foutput=new TFile(Form("%s_cent%d.root",outfile.Data(),centrality),"recreate");
   foutput->cd();
   gr->Write();
   foutput->Close();
   
  
}

int main(int argc, char **argv)
{
  if(argc == 5)
  {
    findthes(argv[1], argv[2], argv[3], atoi(argv[4]));
    return 0;
  }else  {
    std::cout << "Usage: \nmakeTurnOn_fromSameFile.exe <input_HiForest_file> <output_file>" << std::endl;
    return 1;
  }
}
