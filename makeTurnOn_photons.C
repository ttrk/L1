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

const int MAXPHOTONS = 500;

void makeTurnOn_photons(TString inHiForestFileName, TString outFileName, double offlineEtaCut=1.44, bool offlineIso=0)
{
  TFile *inFile = TFile::Open(inHiForestFileName);
  TTree *photonTree = (TTree*)inFile->Get("multiPhotonAnalyzer/photon");
  TTree *hltTree = (TTree*)inFile->Get("hltanalysis/HltTree");
  TTree *skimTree = (TTree*)inFile->Get("skimanalysis/HltTree");
  TTree *evtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
//  photonTree->AddFriend(hltTree);

  Int_t triggerFlag;
  const char* triggerName = "HLT_PAPhoton30_NoCaloIdVL_v1";
  hltTree->SetBranchAddress(triggerName,&triggerFlag);

  Float_t vz;
  Int_t hiBin;
  evtTree->SetBranchAddress("vz",&vz);
  evtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  skimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  skimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

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

  photonTree->SetBranchAddress("nPhotons",&nPhoton);
  photonTree->SetBranchAddress("pt",photon_pt);
  photonTree->SetBranchAddress("eta",photon_eta);
  photonTree->SetBranchAddress("phi",photon_phi);

  photonTree->SetBranchAddress("cc4",cc4);
  photonTree->SetBranchAddress("cr4",cr4);
  photonTree->SetBranchAddress("ct4PtCut20",ct4PtCut20);
  photonTree->SetBranchAddress("trkSumPtHollowConeDR04",trkSumPtHollowConeDR04);
  photonTree->SetBranchAddress("hcalTowerSumEtConeDR04",hcalTowerSumEtConeDR04);
  photonTree->SetBranchAddress("ecalRecHitSumEtConeDR04",ecalRecHitSumEtConeDR04);
  photonTree->SetBranchAddress("hadronicOverEm",hadronicOverEm);
  photonTree->SetBranchAddress("sigmaIetaIeta",sigmaIetaIeta);
  photonTree->SetBranchAddress("isEle",isEle);
  photonTree->SetBranchAddress("sigmaIphiIphi",sigmaIphiIphi);
  photonTree->SetBranchAddress("swissCrx",swissCrx);
  photonTree->SetBranchAddress("seedTime",seedTime);

  TFile* outFile = new TFile(outFileName, "RECREATE");
  outFile->cd();

  const int nBins = 200;
  const double maxPt = 100;

  TH1D* allPt      = new TH1D("allPt",";reco p_{T} (GeV)",nBins,0,maxPt);
  TH1D* allPt_1    = new TH1D("allPt",";reco p_{T} (GeV)",nBins,0,maxPt);
  TH1D* allPt_2    = new TH1D("allPt",";reco p_{T} (GeV)",nBins,0,maxPt);
  TH1D* accepted = new TH1D("accepted",";accepted reco p_{T}",nBins,0,maxPt);

//  photonTree->Draw("pt[0]>>allPt");
//  photonTree->Draw("pt[0]>>accepted", Form("%s > 0", triggerName));

  Long64_t entries = photonTree->GetEntries();
  std::cout << "HiForest entries: " << entries << std::endl;

  for(Long64_t j = 0; j < entries; ++j)
  {
      if(j % 10000 == 0)
        printf("%lld / %lld\n",j,entries);

      photonTree->GetEntry(j);
      hltTree->GetEntry(j);

      double maxfpt = 0;
      double maxfeta = -10;
      double maxfphi = -10;
      for(int i = 0; i < nPhoton; ++i)
      {
        if(TMath::Abs(photon_eta[i]) < offlineEtaCut)
      if(!isEle[i])
        if(TMath::Abs(seedTime[i])<3)
          if(swissCrx[i] < 0.9)
            if(sigmaIetaIeta[i] > 0.002)
          if(sigmaIphiIphi[i] > 0.002)
            if(photon_pt[i] > maxfpt) {
              if(offlineIso){
                if((cc4[i] + cr4[i] + ct4PtCut20[i]) < 0.9)
              if(hadronicOverEm[i] < 0.1)
              {
                maxfpt = photon_pt[i];
                maxfeta = photon_eta[i];
                maxfphi = photon_phi[i];
              }
              } else {
                maxfpt = photon_pt[i];
                maxfeta = photon_eta[i];
                maxfphi = photon_phi[i];
              }
            }
      }

      bool goodEvent = false;
//      if((pcollisionEventSelection == 1) && (TMath::Abs(vz) < 15))
      if((TMath::Abs(vz) < 15))
      {
        goodEvent = true;
      }
      if(!goodEvent) continue;

      allPt->Fill(maxfpt);
      if(triggerFlag>0)
      {
          accepted->Fill(maxfpt);
      }
/*
      if(hiBin < 60)
          allPt_1->Fill(maxfpt);
      else if (hiBin >= 100)
          allPt_2->Fill(maxfpt);
      */
  }

  TGraphAsymmErrors* a = new TGraphAsymmErrors();
  a->BayesDivide(accepted,allPt);
  a->SetName("asymm_pt");



  outFile->Write();
  outFile->Close();
}

int main(int argc, char **argv)
{
  if(argc == 3)
  {
    makeTurnOn_photons(argv[1], argv[2]);
    return 0;
  }
  else if(argc == 4)
  {
    makeTurnOn_photons(argv[1], argv[2], atoi(argv[3]));
    return 0;
  }
  else if(argc == 5)
  {
    makeTurnOn_photons(argv[1], argv[2], atof(argv[3]), atoi(argv[4]));
    return 0;
  }
  else
  {
    return 1;
  }
}
