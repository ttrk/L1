#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include <vector>
#include <iostream>

//#include "EventMatchingCMS.h"

const int MAXL1JETS = 8;
const int MAXJETS = 500;
const Int_t THRESHOLDS = 200; // This will correspond to 0 to 199.5 GeV in 0.5 GeV increments


void makeTurnOn(TString inL1Name, TString inHiForestFileName, TString outFileName, bool montecarlo = false, bool genJets = false)
{

  Double_t L1_THRESHOLD[THRESHOLDS];
  for (int i=0;i<60;i++) L1_THRESHOLD[i]=i*2;

  TFile *outFile = new TFile(outFileName,"RECREATE");

  //TFile *lFile = TFile::Open(inL1FileName);
  //TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeTree");

  // TChain *f1Tree = new TChain("akPu3CaloJetAnalyzer/t","f1Tree");
  // TChain *fEvtTree = new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
  // TChain *fSkimTree = new TChain("skimanalysis/HltTree","fSkimTree");
  // TChain *l1Tree = new TChain("L1UpgradeAnalyzer/L1UpgradeTree","l1Tree");

  // fEvtTree->Add(inHiForestFileName);
  // fSkimTree->Add(inHiForestFileName);
  // f1Tree->Add(inHiForestFileName);
  // l1Tree->Add(inHiForestFileName);

  TFile *inFile = TFile::Open(inHiForestFileName);
  TFile *inL1File = TFile::Open(inL1Name);
  TTree *f1Tree = (TTree*)inFile->Get("akPu3CaloJetAnalyzer/t");
  TTree *fEvtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *fSkimTree = (TTree*)inFile->Get("skimanalysis/HltTree");
  //TTree *l1Tree = (TTree*)inFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");
  TTree *l1Tree = (TTree*)inL1File->Get("L1UpgradeTree");

  //f1Tree->AddFriend(l1Tree);

  Int_t l1_event, l1_run, l1_lumi;
  Int_t l1_hwPt[MAXL1JETS], l1_hwEta[MAXL1JETS], l1_hwPhi[MAXL1JETS];
  Float_t l1_pt[MAXL1JETS];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  l1Tree->SetBranchAddress("jet_pt",l1_pt);

  // f1Tree->SetBranchAddress("event",&l1_event);
  // f1Tree->SetBranchAddress("run",&l1_run);
  // f1Tree->SetBranchAddress("lumi",&l1_lumi);
  // f1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  // f1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  // f1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  // f1Tree->SetBranchAddress("jet_pt",l1_pt);

  Int_t f_evt, f_run, f_lumi;
  Float_t vz;
  Int_t hiBin;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  fEvtTree->SetBranchAddress("vz",&vz);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t f_num;
  Float_t f_pt[MAXJETS];
  Float_t f_eta[MAXJETS];
  Float_t f_phi[MAXJETS];
  Float_t f_rawpt[MAXJETS];

  Int_t num_gen;
  Float_t genpt[MAXJETS], geneta[MAXJETS], genphi[MAXJETS];

  f1Tree->SetBranchAddress("nref",&f_num);
  f1Tree->SetBranchAddress("jtpt",f_pt);
  f1Tree->SetBranchAddress("jteta",f_eta);
  f1Tree->SetBranchAddress("jtphi",f_phi);
  f1Tree->SetBranchAddress("rawpt",f_rawpt);

  if(montecarlo)
  {
    f1Tree->SetBranchAddress("ngen",&num_gen);
    f1Tree->SetBranchAddress("genpt",genpt);
    f1Tree->SetBranchAddress("geneta",geneta);
    f1Tree->SetBranchAddress("genphi",genphi);
  }

  const int nBins = 600;
  const double maxPt = 300;

  TH1D *l1Pt = new TH1D("l1Pt",";L1 p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fPt[3];
  fPt[0] = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_1");
  fPt[2] = (TH1D*)fPt[0]->Clone("fPt_2");
  TH1D *accepted[THRESHOLDS][3];

  for(int i = 0; i < THRESHOLDS; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      accepted[i][j] = new TH1D(Form("accepted_pt%.1f_%d",(i*0.5),j),";offline p_{T}",nBins,0,maxPt);
    }
  }

  TH2D *corr = new TH2D("corr",";offline p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);

  Long64_t entries = f1Tree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    // Only use good collision events ********
    fEvtTree->GetEntry(j);
    fSkimTree->GetEntry(j);

    //std::cout << accepted[0][0]->GetName() << std::endl;
    l1Tree->GetEntry(j);
    f1Tree->GetEntry(j);

    //std::cout << accepted[0][0]->GetName() << std::endl;

    double maxl1pt = 0;
    for(int i = 0; i < MAXL1JETS; ++i)
    {
      if(l1_pt[i] > maxl1pt)
	maxl1pt = l1_pt[i];
    }

    double maxfpt = 0;
    if(!genJets)
    {
      for(int i = 0; i < f_num; ++i)
      {
	if(TMath::Abs(f_eta[i]) > 2.0) continue;
	if(f_pt[i] > maxfpt) {
	  maxfpt = f_pt[i];
	}
      }
    }
    else
    {
      for(int i = 0; i < num_gen; ++i)
      {
	if(TMath::Abs(geneta[i]) > 2.0) continue;
	if(genpt[i] > maxfpt)
	  maxfpt = genpt[i];
      }
    }
    l1Pt->Fill(maxl1pt);

    bool goodEvent = false;
    if((pcollisionEventSelection == 1) && (montecarlo || (pHBHENoiseFilter == 1)) && (TMath::Abs(vz) < 15))
    {
      goodEvent = true;
    }
    if(!goodEvent) continue;

    fPt[0]->Fill(maxfpt);
    if(hiBin < 60)
      fPt[1]->Fill(maxfpt);
    else if (hiBin >= 100)
      fPt[2]->Fill(maxfpt);

    corr->Fill(maxfpt,maxl1pt);

    if((maxfpt > 65) && (maxl1pt < 25))
      std::cout << "Event: " << l1_event << " Lumi: " << l1_lumi << " Run: " << l1_run << std::endl;

    for(int i = 0; i < THRESHOLDS; ++i)
    {
      if(maxl1pt>(i*0.5))
      {
	accepted[i][0]->Fill(maxfpt);
	if(hiBin < 60)
	  accepted[i][1]->Fill(maxfpt);
	else if (hiBin >= 100)
	  accepted[i][2]->Fill(maxfpt);
      }
    }
  }

  TGraphAsymmErrors *a[THRESHOLDS][3];
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
      a[k][l]->SetName(Form("asymm_pt_%.1f_%d",(k*0.5),l));
    }
  }

  outFile->cd();

  l1Pt->Write();
  fPt[0]->Write();
  fPt[1]->Write();
  fPt[2]->Write();
  corr->Write();
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      accepted[k][l]->Write();
      a[k][l]->Write();
    }
  }

  outFile->Close();
}

int main(int argc, char **argv)
{
  if(argc == 4)
  {
    makeTurnOn(argv[1], argv[2], argv[3]);
    return 0;
  }
  else if(argc == 6)
  {
    makeTurnOn(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
  }
  else
  {
    std::cout << "Usage: \nmakeTurnOn_fromSameFile.exe <input_HiForest_file> <output_file>" << std::endl;
    return 1;
  }
}
