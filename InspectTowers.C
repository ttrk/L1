#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TH2F.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include "L1EmulatorSimulator.h"

void L1JetEmulator(TString forest_input = "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v5.root", TString outFileName = "Hydjet502_JetResults.root")
{
  std::cout << "Processing file: " << forest_input << std::endl;
  std::cout << "Saving to: " << outFileName << std::endl;

  L1EmulatorSimulator::cand regions[396];

  // ************ this block makes regions out of towers

  const int nEta = 22;
  const int nPhi = 18;

  TChain *fTowerTree = new TChain("rechitanalyzer/tower","fTowerTree");
  TChain *fEvtTree =  new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
  fTowerTree->Add(forest_input);
  fEvtTree->Add(forest_input);
  TChain *l1Tree = new TChain("L1UpgradeAnalyzer/L1UpgradeTree","l1Tree");
  l1Tree->Add(forest_input);
  TChain *fSkimTree = new TChain("skimanalysis/HltTree");
  fSkimTree->Add(forest_input);

  Int_t l1_event, l1_run, l1_lumi;
  Int_t region_hwPtOriginal[396], region_hwEtaOriginal[396], region_hwPhiOriginal[396];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("region_hwPt", region_hwPtOriginal);
  l1Tree->SetBranchAddress("region_hwEta", region_hwEtaOriginal);
  l1Tree->SetBranchAddress("region_hwPhi", region_hwPhiOriginal);

  Int_t f_evt, f_run, f_lumi,f_pcollisionEventSelection,f_pHBHENoiseFilter;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&f_pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&f_pHBHENoiseFilter);


  Int_t f_nTow = -1;
  Float_t f_eta[10000];
  Float_t f_phi[10000];
  Float_t f_et[10000];

  fTowerTree->SetBranchAddress("eta"                     , &f_eta);
  fTowerTree->SetBranchAddress("phi"                     , &f_phi);
  fTowerTree->SetBranchAddress("et"                      , &f_et);
  fTowerTree->SetBranchAddress("n"                       , &f_nTow);

  Double_t rctEtaMap[23] = {-5.100, -4.500, -4.000, -3.500, -3.000, -2.172, -1.740, -1.392, -1.044, -0.696,
			    -0.348,  0.000,  0.348,  0.696,  1.044,  1.392,  1.740,  2.172,  3.000,  3.500,
			    4.000,  4.500,  5.10};  
  Double_t rctPhiMap[20] = {-1*TMath::Pi(), -2.967060, -2.617990, -2.268930, -1.919860, -1.570800,
			    -1.221730, -0.872665, -0.523599, -0.174533,  0.174533,
			    0.523599,  0.872665,  1.221730,  1.570800,  1.919860,
			    2.268930,  2.617990,  2.96706, TMath::Pi()};

  TH2F* hFtmp    = new TH2F("Ftmp",    "F_tmp"   , nEta , rctEtaMap, nPhi+1, rctPhiMap);
  TH2F* hFtmpReb = new TH2F("FtmpReb", "F_tmpReb", nEta , 0, nEta  , nPhi  , 0, nPhi);
  // **********************************************

  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  TTree *outTree = new TTree("L1UpgradeTree","L1UpgradeTree");

  Int_t run, lumi, evt,pcollisionEventSelection,pHBHENoiseFilter;

  Int_t region_hwPt_[396], region_hwPhi_[396], region_hwEta_[396];
  Int_t region_hwPtOriginal_[396], region_hwPhiOriginal_[396], region_hwEtaOriginal_[396];

  outTree->Branch("run",&run,"run/I");
  outTree->Branch("lumi",&lumi,"lumi/I");
  outTree->Branch("event",&evt,"event/I");
  outTree->Branch("pcollisionEventSelection",&pcollisionEventSelection,"pcollisionEventSelection/I");
  outTree->Branch("pHBHENoiseFilter",&pHBHENoiseFilter,"pHBHENoiseFilter/I");

  outTree->Branch("region_hwPt",region_hwPt_,"region_hwPt[396]/I");
  outTree->Branch("region_hwPhi",region_hwPhi_,"region_hwPhi[396]/I");
  outTree->Branch("region_hwEta",region_hwEta_,"region_hwEta[396]/I");
  outTree->Branch("region_hwPtOriginal",region_hwPtOriginal_,"region_hwPtOriginal[396]/I");
  outTree->Branch("region_hwPhiOriginal",region_hwPhiOriginal_,"region_hwPhiOriginal[396]/I");
  outTree->Branch("region_hwEtaOriginal",region_hwEtaOriginal_,"region_hwEtaOriginal[396]/I");

  Long64_t l_entries = fTowerTree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    fTowerTree->GetEntry(j);
    fEvtTree->GetEntry(j);
    l1Tree->GetEntry(j);
        
    for (int i=0;i<396;i++){
      region_hwPtOriginal_[i]=region_hwPtOriginal[i];
      region_hwPhiOriginal_[i]=region_hwPhiOriginal[i];
      region_hwEtaOriginal_[i]=region_hwEtaOriginal[i];
    }
    
    run = f_run;
    evt = f_evt;
    lumi = f_lumi;
    pcollisionEventSelection=(int)(f_pcollisionEventSelection);
    pHBHENoiseFilter=(int)(f_pHBHENoiseFilter);

    for(int itow = 0; itow < f_nTow; ++itow)
    {
      hFtmp->Fill(f_eta[itow], f_phi[itow], f_et[itow]);
    }

    for(int ieta = 0; ieta < nEta; ++ieta)
    {
      double content = 0.;
      for(int iphi = 0; iphi < nPhi; ++iphi)
      {
	if(iphi == 0 ) content = hFtmp->GetBinContent(ieta+1,1) + hFtmp->GetBinContent(ieta+1,19);
	else           content = hFtmp->GetBinContent(ieta+1,iphi+1);

	if(iphi < 9)         	       hFtmpReb->SetBinContent(ieta+1, iphi+10, content);
	else if(iphi >= 9 && iphi < 18) hFtmpReb->SetBinContent(ieta+1, iphi-8, content);
      }
    }

    int nReg = 0;
    for(int ieta = 0; ieta < nEta; ++ieta)
    {
      for(int iphi = 0; iphi < nPhi; ++iphi)
      {
	regions[nReg].pt = hFtmpReb->GetBinContent(ieta+1, iphi+1) * 2.0;
	regions[nReg].eta = ieta;
	regions[nReg].phi = iphi;

	region_hwPt_[nReg] = regions[nReg].pt;
	region_hwEta_[nReg] = regions[nReg].eta;
	region_hwPhi_[nReg] = regions[nReg].phi;
	nReg++;
      }
    }

    hFtmpReb->Reset();
    hFtmp->Reset();

    outTree->Fill();
  }

  outTree->Write();

  //lFile->Close();
  outFile->Close();
}

int main(int argc, char **argv)
{
  if(argc == 1)
  {
    L1JetEmulator();
    return 0;
  }
  else if (argc == 3)
  {
    L1JetEmulator(argv[1], argv[2]);
    return 0;
  }
  else
  {
    std::cout << "Usage: \nL1JetEmulator.exe <input_file_name> <output_file_name>" << std::endl;
    return 1;
  }
}
