#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "L1EmulatorSimulator.h"

void L1JetEmulator(TString l1_input = "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root", TString outFileName = "Hydjet502_JetResults.root")
{
  std::cout << "Processing file: " << l1_input << std::endl;
  std::cout << "Saving to: " << outFileName << std::endl;

  TFile *lFile = TFile::Open(l1_input);
  TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_event, l1_run, l1_lumi;
  Int_t region_hwPt[396], region_hwEta[396], region_hwPhi[396];
  L1EmulatorSimulator::cand regions[396];
  L1EmulatorSimulator::cand subRegions[396];
  L1EmulatorSimulator::cand outJets[8];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("region_hwPt", region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta", region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi", region_hwPhi);

  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  TTree *outTree = new TTree("L1UpgradeTree","L1UpgradeTree");

  Int_t run, lumi, evt;

  Int_t region_hwPt_[396], region_hwPhi_[396], region_hwEta_[396];
  Int_t subregion_hwPt[396], subregion_hwPhi[396], subregion_hwEta[396];
  Int_t jet_hwPt[8], jet_hwEta[8], jet_hwPhi[8], jet_pt[8];

  outTree->Branch("run",&run,"run/I");
  outTree->Branch("lumi",&lumi,"lumi/I");
  outTree->Branch("event",&evt,"evt/I");

  outTree->Branch("region_hwPt",region_hwPt_,"region_hwPt[396]/I");
  outTree->Branch("region_hwPhi",region_hwPhi_,"region_hwPhi[396]/I");
  outTree->Branch("region_hwEta",region_hwEta_,"region_hwEta[396]/I");
  outTree->Branch("subregion_hwPt",subregion_hwPt,"subregion_hwPt[396]/I");
  outTree->Branch("subregion_hwPhi",subregion_hwPhi,"subregion_hwPhi[396]/I");
  outTree->Branch("subregion_hwEta",subregion_hwEta,"subregion_hwEta[396]/I");

  outTree->Branch("jet_hwPt",jet_hwPt,"jet_hwPt[8]/I");
  outTree->Branch("jet_pt",jet_pt,"jet_pt[8]/I");
  outTree->Branch("jet_hwPhi",jet_hwPhi,"jet_hwPhi[8]/I");
  outTree->Branch("jet_hwEta",jet_hwEta,"jet_hwEta[8]/I");

  Long64_t l_entries = l1Tree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    l1Tree->GetEntry(j);

    run = l1_run;
    lumi = l1_lumi;
    evt = l1_event;

    // copy input to output tree, also
    // pack regions into array for subtraction method.
    for(int i = 0; i < 396; ++i)
    {
      region_hwPt_[i] = region_hwPt[i];
      region_hwEta_[i] = region_hwEta[i];
      region_hwPhi_[i] = region_hwPhi[i];

      regions[i].pt = region_hwPt_[i];
      regions[i].eta = region_hwEta_[i];
      regions[i].phi = region_hwPhi_[i];
    }

    // perform phi-ring avg background subtraction
    L1EmulatorSimulator::CaloRingBackgroundSubtraction(regions, subRegions);

    // copy sub regions to output tree
    for(int i = 0; i < 396; ++i)
    {
      subregion_hwPt[i] =  subRegions[i].pt;
      subregion_hwEta[i] = subRegions[i].eta;
      subregion_hwPhi[i] = subRegions[i].phi;
    }

    // run the jet finder on the input regions
    L1EmulatorSimulator::SlidingWindowJetFinder(subRegions, outJets);
    //RegionJetFinder(subRegions, outJets);

    // copy the jets to the output tree
    // The factor of 4 is the conversion factor from hardware value to
    // physical GeV.
    for(int i = 0; i < 8; i++)
    {
      jet_hwPt[i] = outJets[i].pt;
      jet_pt[i] = outJets[i].pt * 4;
      jet_hwEta[i] = outJets[i].eta;
      jet_hwPhi[i] = outJets[i].phi;
    }

    outTree->Fill();
  }

  outTree->Write();

  lFile->Close();
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
