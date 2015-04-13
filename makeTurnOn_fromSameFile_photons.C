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
#include "L1EmulatorSimulator.h"

const int MAXL1EMCANDS = 144;
const int MAXL1REGIONS = 396;
const int MAXL1JETS = 8;
const int MAXPHOTONS = 500;
const Int_t THRESHOLDS = 60;
// const Double_t L1_THRESHOLD[THRESHOLDS] = {0, 4, 8, 12, 16, 20, 24,
// 					   28, 32, 36, 40, 44, 48,
// 					   52, 56, 60, 64, 68, 72,
// 					   76, 80, 84, 88, 92, 96,
// 					   100, 104, 108, 112, 116};


void makeTurnOn(/*TString inL1Name,*/ TString inHiForestFileName, TString outFileName)
{
	//////// Kaya's modificiation ////////
	// zero out regions with 0 <= region.eta <= 5 or 16 <= region.eta <= 21
	// decide zeroOut from file name
	bool zeroOut = ( outFileName.Contains("2x2") ||  outFileName.Contains("3x3") );
	std::cout << "zeroOut = " << zeroOut << std::endl;
	//////// Kaya's modificiation - END ////////

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
  //TFile *inL1File = TFile::Open(inL1Name);
  //TTree *f1Tree = (TTree*)inFile->Get("akPu3CaloJetAnalyzer/t");
  TTree *f1Tree = (TTree*)inFile->Get("multiPhotonAnalyzer/photon");
  TTree *fEvtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *fSkimTree = (TTree*)inFile->Get("skimanalysis/HltTree");
  TTree *l1Tree = (TTree*)inFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");
  //TTree *l1Tree = (TTree*)inL1File->Get("L1UpgradeTree");

  //f1Tree->AddFriend(l1Tree);

  Int_t l1_event, l1_run, l1_lumi;
  Int_t l1_hwPt[MAXL1JETS], l1_hwEta[MAXL1JETS], l1_hwPhi[MAXL1JETS];
  Double_t l1_pt[MAXL1JETS];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  l1Tree->SetBranchAddress("jet_pt",l1_pt);

  Int_t emcand_hwPt[MAXL1EMCANDS], emcand_hwEta[MAXL1EMCANDS], emcand_hwPhi[MAXL1EMCANDS];
  Int_t region_hwPt[MAXL1REGIONS], region_hwEta[MAXL1REGIONS], region_hwPhi[MAXL1REGIONS];

  l1Tree->SetBranchAddress("emcand_hwPt",emcand_hwPt);
  l1Tree->SetBranchAddress("emcand_hwEta",emcand_hwEta);
  l1Tree->SetBranchAddress("emcand_hwPhi",emcand_hwPhi);
  l1Tree->SetBranchAddress("region_hwPt",region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta",region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi",region_hwPhi);

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

  Int_t nPhoton;
  Float_t photon_pt[MAXPHOTONS];
  Float_t photon_eta[MAXPHOTONS];
  Float_t photon_phi[MAXPHOTONS];
  Float_t cc4[MAXPHOTONS];
  Float_t cr4[MAXPHOTONS];
  Float_t ct4PtCut20[MAXPHOTONS];
  Float_t trkSumPtHollowConeDR04[MAXPHOTONS];
  Float_t hcalTowerSumEtConeDR04[MAXPHOTONS];
  Float_t ecalRecHitSumEtConeDR04[MAXPHOTONS];
  Float_t hadronicOverEm[MAXPHOTONS];
  Float_t sigmaIetaIeta[MAXPHOTONS];
  Int_t isEle[MAXPHOTONS];
  Float_t sigmaIphiIphi[MAXPHOTONS];
  Float_t swissCrx[MAXPHOTONS];
  Float_t seedTime[MAXPHOTONS];

  f1Tree->SetBranchAddress("nPhotons",&nPhoton);
  f1Tree->SetBranchAddress("pt",photon_pt);
  f1Tree->SetBranchAddress("eta",photon_eta);
  f1Tree->SetBranchAddress("phi",photon_phi);

  f1Tree->SetBranchAddress("cc4",cc4);
  f1Tree->SetBranchAddress("cr4",cr4);
  f1Tree->SetBranchAddress("ct4PtCut20",ct4PtCut20);
  f1Tree->SetBranchAddress("trkSumPtHollowConeDR04",trkSumPtHollowConeDR04);
  f1Tree->SetBranchAddress("hcalTowerSumEtConeDR04",hcalTowerSumEtConeDR04);
  f1Tree->SetBranchAddress("ecalRecHitSumEtConeDR04",ecalRecHitSumEtConeDR04);
  f1Tree->SetBranchAddress("hadronicOverEm",hadronicOverEm);
  f1Tree->SetBranchAddress("sigmaIetaIeta",sigmaIetaIeta);
  f1Tree->SetBranchAddress("isEle",isEle);
  f1Tree->SetBranchAddress("sigmaIphiIphi",sigmaIphiIphi);
  f1Tree->SetBranchAddress("swissCrx",swissCrx);
  f1Tree->SetBranchAddress("seedTime",seedTime);

  const int nBins = 200;
  const double maxPt = 200;

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
      accepted[i][j] = new TH1D(Form("accepted_pt%d_%d",i,j),";offline p_{T}",nBins,0,maxPt);
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
    bool goodEvent = false;
    // 5.02 TeV Hydjet missing pcollisionEventSelection for now
    //if((pcollisionEventSelection == 1) && (montecarlo || (pHBHENoiseFilter == 1)) && (TMath::Abs(vz) < 15))
    if((TMath::Abs(vz) < 15))
    {
      goodEvent = true;
    }
    if(!goodEvent) continue;

    //std::cout << accepted[0][0]->GetName() << std::endl;
    l1Tree->GetEntry(j);
    f1Tree->GetEntry(j);

    double maxl1pt = 0;
    // for(int i = 0; i < MAXL1EMCANDS; ++i)
    // {
    //   if(emcand_hwPt[i] > maxl1pt)
    // 	maxl1pt = emcand_hwPt[i];
    // }
    L1EmulatorSimulator::cand regions[396];
    L1EmulatorSimulator::cand subregions_tmp[396];
    for(int i = 0; i < MAXL1REGIONS; ++i)
    {
      regions[i].pt = region_hwPt[i];
      regions[i].eta = region_hwEta[i];
      regions[i].phi = region_hwPhi[i];
    }
    //    L1EmulatorSimulator::CaloRingBackgroundSubtraction(regions, subregions);
    L1EmulatorSimulator::CaloRingBackgroundSubtraction(regions, subregions_tmp);

    //////// Kaya's modificiation ////////
    /*
     *  run 2x2 and 3x3 jet finder algorithms with forward regions (eta<=5 or eta >=16) zeroed out.
     *  zero out regions with 0 <= region.eta <= 5 or 16 <= region.eta <= 21
    */
    int newMAXL1REGIONS = MAXL1REGIONS;
    if(zeroOut)
    {
    	newMAXL1REGIONS = 8;

		// zero out forward regions
    	for(int i = 0; i < MAXL1REGIONS; i++) {
    		if(subregions_tmp[i].eta <= 5 || subregions_tmp[i].eta >= 16)
    		{
    			subregions_tmp[i].pt = 0;
    		}
    	}
    }
    L1EmulatorSimulator::cand subregions[newMAXL1REGIONS];
    if(zeroOut)
    {
    	if(outFileName.Contains("2x2"))
    	{
    		L1EmulatorSimulator::SlidingWindowJetFinder(subregions_tmp,subregions, L1EmulatorSimulator::twoByTwoANDzeroWalls);
    	}
    	else if(outFileName.Contains("3x3"))
    	{
    		L1EmulatorSimulator::SlidingWindowJetFinder(subregions_tmp,subregions, L1EmulatorSimulator::nominal);
    	}
    	else
    	{
    		std::cout << "Use either 2x2 or 3x3" << std::endl;
    		exit(1);
    	}
    }
    else
    {
    	for(int i = 0; i < newMAXL1REGIONS; i++)
    	{
        	subregions[i].pt = subregions_tmp[i].pt;
        	subregions[i].eta = subregions_tmp[i].eta;
        	subregions[i].phi = subregions_tmp[i].phi;
    	}
    }
    //////// Kaya's modificiation - END ////////

    //    for(int i = 0; i < MAXL1REGIONS; ++i)
    for(int i = 0; i < newMAXL1REGIONS; ++i)
    {
    	if(subregions[i].eta < 6 || subregions[i].eta > 15) continue;
    	//////// Kaya's modificiation ////////
    	//      if(subregions[i].pt > maxl1pt)
    	//	maxl1pt = subregions[i].pt;
    	// take hardware scale change into account
    	if(zeroOut)
    	{
    		if(subregions[i].pt * 4 > maxl1pt)
    			maxl1pt = subregions[i].pt * 4 ;
    	}
    	else
    	{
    		if(subregions[i].pt * 0.5 > maxl1pt)
    			maxl1pt = subregions[i].pt * 0.5;
    	}
    	//////// Kaya's modificiation - END ////////
    }

    double maxfpt = 0;
    for(int i = 0; i < nPhoton; ++i)
    {
      if((cc4[i] + cr4[i] + ct4PtCut20[i]) < 0.9)
      if(TMath::Abs(photon_eta[i]) < 1.4791)
      if(!isEle[i])
      if(TMath::Abs(seedTime[i])<3)
      if(swissCrx[i] < 0.9)
      if(sigmaIetaIeta[i] > 0.002)
      if(sigmaIphiIphi[i] > 0.002)
      if(hadronicOverEm[i] < 0.1)
      if(photon_pt[i] > maxfpt) {
	maxfpt = photon_pt[i];
      }
    }
    l1Pt->Fill(maxl1pt);

    fPt[0]->Fill(maxfpt);
    if(hiBin < 60)
      fPt[1]->Fill(maxfpt);
    else if (hiBin >= 100)
      fPt[2]->Fill(maxfpt);

    corr->Fill(maxfpt,maxl1pt);

    for(int i = 0; i < THRESHOLDS; ++i)
    {
      if(maxl1pt> i)
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
      a[k][l]->SetName(Form("asymm_pt_%d_%d",k,l));
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
  if(argc == 3)
  {
    makeTurnOn(argv[1], argv[2]);
    return 0;
  }
  else
  {
    std::cout << "Usage: \nmakeTurnOn_fromSameFile.exe <input_HiForest_file> <output_file>" << std::endl;
    return 1;
  }
}
