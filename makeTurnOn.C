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

#include "EventMatchingCMS.h"

const int MAXL1JETS = 8;
const int MAXJETS = 500;
const Int_t THRESHOLDS = 11;
const Double_t L1_THRESHOLD[THRESHOLDS] = {16,20,32,36,40,44,52,68,80,92,128};

void makeTurnOn(TString inL1FileName, TString inHiForestFileName, TString outFileName, bool montecarlo = false)
{
  TFile *lFile = TFile::Open(inL1FileName);
  TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeTree");

  Int_t l1_event, l1_run, l1_lumi;
  Int_t l1_hwPt[MAXL1JETS], l1_hwEta[MAXL1JETS], l1_hwPhi[MAXL1JETS];
  Int_t l1_pt[MAXL1JETS];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  l1Tree->SetBranchAddress("jet_pt",l1_pt);

  TChain *f1Tree = new TChain("akPu3CaloJetAnalyzer/t","f1Tree");
  TChain *fEvtTree = new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
  TChain *fSkimTree = new TChain("skimanalysis/HltTree","fSkimTree");

  fEvtTree->Add(inHiForestFileName);
  fSkimTree->Add(inHiForestFileName);
  f1Tree->Add(inHiForestFileName);

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

  TFile *outFile = new TFile(outFileName,"RECREATE");

  const int nBins = 75;
  const double maxPt = 300;

  TH1D *l1Pt = new TH1D("l1Pt",";L1 p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fPt[3];
  fPt[0] = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_1");
  fPt[2] = (TH1D*)fPt[0]->Clone("fPt_2");
  TH1D *accepted[THRESHOLDS][3];

  for(int i = 0; i < THRESHOLDS; ++i)
    for(int j = 0; j < 3; ++j)
    {
      accepted[i][j] = new TH1D(Form("accepted_pt%d_%d",(int)L1_THRESHOLD[i],j),";offline p_{T}",nBins,0,maxPt);
    }

  TH2D *corr = new TH2D("corr",";offline p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);

  // Make the event-matching map ************
  EventMatchingCMS *matcher = new EventMatchingCMS();
  std::cout << "Begin making map." << std::endl;
  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);
    matcher->addEvent(l1_event, l1_lumi, l1_run, j);
    //matcher->addEvent(l1_event, 0, l1_run, j);
  }
  std::cout << "Finished making map." << std::endl;
  // **********************

  outFile->cd();

  int count = 0;

  Long64_t entries = f1Tree->GetEntries();
  std::cout << "L1 entries: " << l_entries << std::endl;
  std::cout << "HiForest entries: " << entries << std::endl;
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    // Only use good collision events ********
    fEvtTree->GetEntry(j);
    fSkimTree->GetEntry(j);
    bool goodEvent = false;
    // 5.02 TeV Hydjet missing pcollisionEventSelection for now
    //if((pcollisionEventSelection == 1) && (montecarlo || (pHBHENoiseFilter == 1)) && (TMath::Abs(vz) < 15))
    if((montecarlo || (pHBHENoiseFilter == 1)) && (TMath::Abs(vz) < 15))
    {
      goodEvent = true;
    }
    if(!goodEvent) continue;
    //**************

    // retrieve the matching entry from the l1 input ***********
    long long l1Entry = matcher->retrieveEvent(f_evt, f_lumi, f_run);
    //long long l1Entry = matcher->retrieveEvent(f_evt, 0, f_run);
    if( l1Entry == -1 )  continue;
    //**************

    l1Tree->GetEntry(l1Entry);
    f1Tree->GetEntry(j);
    count++;

    double maxl1pt = 0;
    for(int i = 0; i < MAXL1JETS; ++i)
    {
      if(l1_pt[i] > maxl1pt)
	maxl1pt = l1_pt[i];
    }

    double maxfpt = 0;
    if(!montecarlo)
    {
      for(int i = 0; i < f_num; ++i)
      {
	//if(TMath::Abs(PuJet_eta[i]) > 2) continue;
	if(f_pt[i] > maxfpt) {
	  maxfpt = f_pt[i];
	}
      }
    }
    else
    {
      for(int i = 0; i < num_gen; ++i)
      {
	//if(TMath::Abs(genJet_eta[i]) > 2.0) continue;
	if(genpt[i] > maxfpt)
	  maxfpt = genpt[i];
      }
    }
    l1Pt->Fill(maxl1pt);


    fPt[0]->Fill(maxfpt);
    if(hiBin < 60)
      fPt[1]->Fill(maxfpt);
    else if (hiBin >= 100)
      fPt[2]->Fill(maxfpt);

    corr->Fill(maxfpt,maxl1pt);

    for(int k = 0; k < THRESHOLDS; ++k)
    {
      if(maxl1pt>L1_THRESHOLD[k])
      {
	accepted[k][0]->Fill(maxfpt);
	if(hiBin < 60)
	  accepted[k][1]->Fill(maxfpt);
	else if (hiBin >= 100)
	  accepted[k][2]->Fill(maxfpt);
      }
    }
  }

  std::cout << "Matching entries: " << count << std::endl;

  TGraphAsymmErrors *a[THRESHOLDS][3];
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
      a[k][l]->SetName(Form("asymm_pt_%d_%d",(int)L1_THRESHOLD[k],l));
    }
  }

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

  lFile->Close();
  outFile->Close();
}

int main(int argc, char **argv)
{
  if(argc == 4)
  {
    makeTurnOn(argv[1], argv[2], argv[3]);
    return 0;
  }
  else if(argc == 5)
  {
    makeTurnOn(argv[1], argv[2], argv[3], atoi(argv[4]));
  }
  else
  {
    std::cout << "Usage: \nmakeTurnOn.exe <input_l1_file> <input_HiForest_file> <output_file>" << std::endl;
    return 1;
  }
}
