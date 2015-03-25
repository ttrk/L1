#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <algorithm>

struct cand{
  int pt;
  int eta;
  int phi;
};

void CaloRingBackgroundSubtraction(cand region[396], cand subregion[396]);

void SlidingWindowJetFinder(cand input[396], cand output[8]);
void RegionJetFinder(cand input[396], cand output[8]);

int deltaGctPhi(int phi1, int phi2);

void L1JetEmulator(TString l1_input = "/export/d00/scratch/luck/HydjetMB_740pre8_MCHI2_74_V3_53XBS_L1UpgradeAnalyzer_GT_MCHI2_74_V3.root", TString outFileName = "Hydjet502_JetResults.root")
{
  std::cout << "Processing file: " << l1_input << std::endl;
  std::cout << "Saving to: " << outFileName << std::endl;
  
  TFile *lFile = TFile::Open(l1_input);
  //TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");
  
  // Int_t l1_event, l1_run, l1_lumi;
  // Int_t region_hwPt[396], region_hwEta[396], region_hwPhi[396];
  cand regions[396];

  // ************ this block makes regions out of towers

  const int nEta = 22;
  const int nPhi = 18;

  TTree *fTowerTree;
  fTowerTree = (TTree*)lFile->Get("rechitanalyzer/tower");
  
  Int_t f_nTow = -1;
  Float_t f_eta[10000];
  Float_t f_phi[10000];
  Float_t f_et[10000];

  fTree->SetBranchAddress("eta"                     , &f_eta);
  fTree->SetBranchAddress("phi"                     , &f_phi);
  fTree->SetBranchAddress("et"                      , &f_et);
  fTree->SetBranchAddress("n"                       , &f_nTow);

  Double_t etaMap[65] = {-5.100, -4.500, -4.000, -3.500, -3.000, -2.725, -2.500, -2.336, -2.172, -2.051,
		       	 -1.930, -1.835, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, 
		       	 -1.044, -0.957, -0.870, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, 
		       	 -0.174, -0.087,  0.000,  0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609, 
			 0.696,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,  1.305,  1.392,  1.479, 
			 1.566,  1.653,  1.740,  1.835,  1.930,  2.051,  2.172,  2.336,  2.500,  2.725, 
			 3.000,  3.500,  4.000,  4.500,  5.100};
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
  cand subRegions[396];
  cand outJets[8];

  // l1Tree->SetBranchAddress("event",&l1_event);
  // l1Tree->SetBranchAddress("lumi",&l1_lumi);
  // l1Tree->SetBranchAddress("run",&l1_run);
  // l1Tree->SetBranchAddress("region_hwPt", region_hwPt);
  // l1Tree->SetBranchAddress("region_hwEta", region_hwEta);
  // l1Tree->SetBranchAddress("region_hwPhi", region_hwPhi);

  TFile *outFile = TFile::Open(outFileName,"RECREATE");
  TTree *outTree = new TTree("L1UpgradeTree","L1UpgradeTree");

  Int_t run, lumi, evt;

  Int_t region_hwPt_[396], region_hwPhi_[396], region_hwEta_[396];
  Int_t subregion_hwPt[396], subregion_hwPhi[396], subregion_hwEta[396];
  Int_t jet_hwPt[8], jet_hwEta[8], jet_hwPhi[8], jet_pt[8];

  outTree->Branch("run",&run,"run/I");
  outTree->Branch("lumi",&lumi,"lumi/I");
  outTree->Branch("evt",&evt,"evt/I");

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

  Long64_t l_entries = fTowerTree->GetEntries();
  for(Long64_t j = 0; j < l_entries; ++j)
  {
    fTowerTree->GetEntry(j);

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

	region_hwPt_[i] = regions[nReg].pt;
	region_hwEta_[i] = regions[nReg].eta;
	region_hwPhi_[i] = regions[nReg].phi = iphi;
	nReg++;
      }
    }
    
    hFtmpReb->Reset();
    hL1tmp->Reset();

    // perform phi-ring avg background subtraction
    CaloRingBackgroundSubtraction(regions, subRegions);

    // copy sub regions to output tree
    for(int i = 0; i < 396; ++i)
    {
      subregion_hwPt[i] =  subRegions[i].pt;
      subregion_hwEta[i] = subRegions[i].eta; 
      subregion_hwPhi[i] = subRegions[i].phi;
    }

    // run the jet finder on the input regions
    SlidingWindowJetFinder(subRegions, outJets);
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

void CaloRingBackgroundSubtraction(cand region[396], cand subregion[396])
{
  int etaCount[22];
  int puLevelHI[22];
  float r_puLevelHI[22];

  // 22 values of eta
  for(unsigned i = 0; i < 22; ++i)
  {
    puLevelHI[i] = 0;
    r_puLevelHI[i] = 0.0;
    etaCount[i] = 0;
  }

  for(int i = 0; i < 396; ++i){
    r_puLevelHI[region[i].eta] += region[i].pt;
    etaCount[region[i].eta]++;
  }

  for(unsigned i = 0; i < 22; ++i)
  {
    if(etaCount[i] != 18)
      std::cout << "ERROR: wrong number of regions in phi ring." << std::endl;
    puLevelHI[i] = floor(r_puLevelHI[i]/18. + 0.5); // this floating point operation should probably be replaced
  }

  for(int i = 0; i < 396; ++i){
    subregion[i].pt = std::max(0, region[i].pt - puLevelHI[region[i].eta]);
    subregion[i].eta = region[i].eta;
    subregion[i].phi = region[i].phi;
  }
}

void SlidingWindowJetFinder(cand region[396], cand output[8])
{
  std::vector<cand> forjets;
  std::vector<cand> cenjets;

  for(int i = 0; i < 396; i++) {
    int regionET = region[i].pt;
    int regionEta = region[i].eta;
    int regionPhi = region[i].phi;
    //if(regionEta == 4 || regionEta == 17) regionET =0;
    int neighborN_et = 0;
    int neighborS_et = 0;
    int neighborE_et = 0;
    int neighborW_et = 0;
    int neighborNE_et = 0;
    int neighborSW_et = 0;
    int neighborNW_et = 0;
    int neighborSE_et = 0;
    unsigned int nNeighbors = 0;
    for(int j = 0; j < 396; j++) {
      int neighborET = region[j].pt;
      int neighborEta = region[j].eta;
      //if(neighborEta == 4 || neighborEta == 17) neighborET =0;
      int neighborPhi = region[j].phi;
      if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	 (regionEta ) == neighborEta) {
	neighborN_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	      (regionEta    ) == neighborEta) {
	neighborS_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
	      (regionEta + 1) == neighborEta) {
	neighborE_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
	      (regionEta - 1) == neighborEta) {
	neighborW_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	      (regionEta + 1) == neighborEta) {
	neighborNE_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	      (regionEta - 1) == neighborEta) {
	neighborSW_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	      (regionEta - 1) == neighborEta) {
	neighborNW_et = neighborET;
	nNeighbors++;
	continue;
      }
      else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	      (regionEta + 1) == neighborEta) {
	neighborSE_et = neighborET;
	nNeighbors++;
	continue;
      }
    }
    if(regionET > neighborN_et &&
       regionET > neighborNW_et &&
       regionET > neighborW_et &&
       regionET > neighborSW_et &&
       regionET >= neighborNE_et &&
       regionET >= neighborE_et &&
       regionET >= neighborSE_et &&
       regionET >= neighborS_et) {
      unsigned int jetET = regionET +
	neighborN_et + neighborS_et + neighborE_et + neighborW_et +
	neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;

      int jetPhi = regionPhi;
      int jetEta = regionEta;

      bool neighborCheck = (nNeighbors == 8);
      // On the eta edge we only expect 5 neighbors
      if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
	neighborCheck = true;

      if (!neighborCheck) {
	std::cout << "ERROR (jet finder): Wrong number of neighbor regions." << std::endl;
	std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
      }

      cand theJet;
      theJet.pt = jetET / 8; // factor of 8 comes from hardware scale change
      theJet.eta = jetEta;
      theJet.phi = jetPhi;

      const bool forward = (jetEta < 4 || jetEta > 17);

      if(forward)
	forjets.push_back(theJet);
      else
	cenjets.push_back(theJet);
    }
  }

  auto comp = [&](cand i, cand j)-> bool {
    return (i.pt > j.pt );
  };


  // sort the jet collections and only take the largest 4
  // emulator only outputs top 4 in each category.
  std::sort(forjets.begin(), forjets.end(), comp);
  std::sort(cenjets.begin(), cenjets.end(), comp);
  forjets.resize(4);
  cenjets.resize(4);

  output[0]=cenjets.at(0);
  output[1]=cenjets.at(1);
  output[2]=cenjets.at(2);
  output[3]=cenjets.at(3);
  output[4]=forjets.at(0);
  output[5]=forjets.at(1);
  output[6]=forjets.at(2);
  output[7]=forjets.at(3);
}

void RegionJetFinder(cand region[396], cand output[8])
{
  std::vector<cand> forjets;
  std::vector<cand> cenjets;

  for(int i = 0; i < 396; i++) {
    region[i].pt = region[i].pt /8;
    if(region[i].eta < 4 || region[i].eta > 18) {
      forjets.push_back(region[i]);
    } else {
      cenjets.push_back(region[i]);
    }
  }

  auto comp = [&](cand i, cand j)-> bool {
    return (i.pt > j.pt );
  };

  // sort the jet collections and only take the largest 4
  // emulator only outputs top 4 in each category.
  std::sort(forjets.begin(), forjets.end(), comp);
  std::sort(cenjets.begin(), cenjets.end(), comp);
  forjets.resize(4);
  cenjets.resize(4);

  output[0]=cenjets.at(0);
  output[1]=cenjets.at(1);
  output[2]=cenjets.at(2);
  output[3]=cenjets.at(3);
  output[4]=forjets.at(0);
  output[5]=forjets.at(1);
  output[6]=forjets.at(2);
  output[7]=forjets.at(3);

}


int deltaGctPhi(int phi1, int phi2)
{
  int diff = phi1 - phi2;
  if (std::abs(phi1 - phi2) == 17) { //18 regions in phi
    diff = -diff/std::abs(diff);
  }
  return diff;
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
