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

const int MAXL1EMCANDS = 144;
const int MAXL1REGIONS = 396;
const int MAXL1JETS = 8;
const int MAXPHOTONS = 500;
const Int_t THRESHOLDS = 1;

bool compareHistograms(TH1* h1, TH1* h2);

void makeTurnOn(TString inHiForestFileName, TString outFileName, double offlineEtaCut=1.44, bool offlineIso=0)
{
  TFile *outFile = new TFile(outFileName,"RECREATE");

  TFile *inFile = TFile::Open(inHiForestFileName);
  TTree *f1Tree = (TTree*)inFile->Get("multiPhotonAnalyzer/photon");
  TTree *fEvtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *fSkimTree = (TTree*)inFile->Get("skimanalysis/HltTree");
  TTree *hltTree = (TTree*)inFile->Get("hltanalysis/HltTree");
  f1Tree->AddFriend(fEvtTree);

  Int_t HLT_PAPhoton30_NoCaloIdVL_v1;
  hltTree->SetBranchAddress("HLT_PAPhoton30_NoCaloIdVL_v1",&HLT_PAPhoton30_NoCaloIdVL_v1);

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

  const int nBins = 100;
  const int maxPt = 100;

  TH1D *fPt[3];
  fPt[0] = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_1");
  fPt[2] = (TH1D*)fPt[0]->Clone("fPt_2");

  TH1D *fPt_vTest[3];
  fPt_vTest[0] = new TH1D("fPt_vTest_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt_vTest[1] = (TH1D*)fPt_vTest[0]->Clone("fPt_vTest_1");
  fPt_vTest[2] = (TH1D*)fPt_vTest[0]->Clone("fPt_vTest_2");

  TString formula = "pt";
  TString condition = Form("abs(eta)<%f", offlineEtaCut);
  condition += " && isEle==0";
  condition += " && abs(seedTime)<3";
  condition += " && swissCrx<0.9";
  condition += " && sigmaIetaIeta>0.002";
  condition += " && sigmaIphiIphi>0.002";
  if (offlineIso)
  {
      condition += " && (cc4 + cr4 + ct4PtCut20) < 0.9";
      condition += " && hadronicOverEm<0.01";
  }
  const char* cut_hiBin60 ="hiBin<60";
  const char* cut_hiBin100="hiBin>=100";

//  f1Tree->Draw(Form("MaxIf$(%s,%s) >> %s",formula.Data(),condition.Data(),fPt_vTest[0]->GetName()));
//  f1Tree->Draw(Form("MaxIf$(%s,nPhotons>0) >> %s",formula.Data(),fPt_vTest[0]->GetName()));
//  f1Tree->Draw(Form("Max$(%s) >> %s",formula.Data(),fPt_vTest[0]->GetName()));
  f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula.Data(), condition.Data(), fPt_vTest[0]->GetName()));
  f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula.Data(), condition.Data(), fPt_vTest[1]->GetName()), cut_hiBin60);
  f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula.Data(), condition.Data(), fPt_vTest[2]->GetName()), cut_hiBin100);
  std::cout<<"condition = "<<condition.Data()<<std::endl;

  TH1D *fPtSigma[3];
  fPtSigma[0] = new TH1D("fPtSigma_0",";#sigma_{#eta #eta}", 100, 0, 0.1);
  fPtSigma[1] = (TH1D*)fPtSigma[0]->Clone("fPtSigma_1");
  fPtSigma[2] = (TH1D*)fPtSigma[0]->Clone("fPtSigma_2");

  TH1D *accepted[THRESHOLDS][3];
  TH1D *accepted_sieie[THRESHOLDS][3];

  for(int i = 0; i < THRESHOLDS; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      accepted[i][j] = new TH1D(Form("accepted_pt%d_%d",i,j),";offline p_{T}",nBins,0,maxPt);
      accepted_sieie[i][j] = new TH1D(Form("accepted_sieie%d_%d",i,j), ";#sigma_{#eta #eta}",100,0,0.1);
    }
  }

  Long64_t entries = f1Tree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    fEvtTree->GetEntry(j);
    fSkimTree->GetEntry(j);

    // bool goodEvent = false;
    // if((pcollisionEventSelection == 1) && (TMath::Abs(vz) < 15))
    // {
    //   goodEvent = true;
    // }
    // if(!goodEvent) continue;

    hltTree->GetEntry(j);
    f1Tree->GetEntry(j);


    double maxfpt = 0;
    //double maxfeta = -10;
    //double maxfphi = -10;
    double sigmaietaieta = -1;
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
			  //maxfeta = photon_eta[i];
			  //maxfphi = photon_phi[i];
			  sigmaietaieta = sigmaIetaIeta[i];
			}
		    } else {
		      maxfpt = photon_pt[i];
		      //maxfeta = photon_eta[i];
		      //maxfphi = photon_phi[i];
		      sigmaietaieta = sigmaIetaIeta[i];
		    }
		  }
    }


    fPt[0]->Fill(maxfpt);
    if(hiBin < 60){
      fPt[1]->Fill(maxfpt);
    }
    else if (hiBin >= 100) {
      fPt[2]->Fill(maxfpt);
    }

    if(maxfpt > 40)
    {
      fPtSigma[0]->Fill(sigmaietaieta);
      if(hiBin < 60){
	fPtSigma[1]->Fill(sigmaietaieta);
      }
      else if (hiBin >= 100) {
	fPtSigma[2]->Fill(sigmaietaieta);
      }
    }

    for(int i = 0; i < THRESHOLDS; ++i)
    {
      if(HLT_PAPhoton30_NoCaloIdVL_v1)
      {
	accepted[i][0]->Fill(maxfpt);
	if(hiBin < 60){
	  accepted[i][1]->Fill(maxfpt);
	}
	else if (hiBin >= 100){
	  accepted[i][2]->Fill(maxfpt);
	}

	if(maxfpt > 40)
	{
	  accepted_sieie[i][0]->Fill(sigmaietaieta);
	  if(hiBin < 60){
	    accepted_sieie[i][1]->Fill(sigmaietaieta);
	  }
	  else if (hiBin >= 100){
	    accepted_sieie[i][2]->Fill(sigmaietaieta);
	  }
	}
      }
    }
  }

  TGraphAsymmErrors *a[THRESHOLDS][3];
  TGraphAsymmErrors *a_sieie[THRESHOLDS][3];
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
      a[k][l]->SetName(Form("asymm_pt_%d_%d",k,l));

      a_sieie[k][l] = new TGraphAsymmErrors();
      a_sieie[k][l]->BayesDivide(accepted_sieie[k][l], fPtSigma[l]);
      a_sieie[k][l]->SetName(Form("asymm_sieie_%d_%d",k,l));
    }
  }

  outFile->cd();

  fPt[0]->Write();
  fPt[1]->Write();
  fPt[2]->Write();
  fPt_vTest[0]->Write();
  fPt_vTest[1]->Write();
  fPt_vTest[2]->Write();
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      accepted[k][l]->Write();
      accepted_sieie[k][l]->Write();
      a[k][l]->Write();
      a_sieie[k][l]->Write();
    }
  }

  outFile->Close();

//  bool same=compareHistograms(fPt[0],fPt_vTest[0]);
  std::cout << compareHistograms(fPt[0],fPt_vTest[0]) << std::endl;
  std::cout << compareHistograms(fPt[1],fPt_vTest[1]) << std::endl;
  std::cout << compareHistograms(fPt[2],fPt_vTest[2]) << std::endl;
}

/*
 * compare two histograms bin by bin.
 * */
bool compareHistograms(TH1* h1, TH1* h2)
{
    int numBins=h1->GetNbinsX();
    if(numBins != h2->GetNbinsX())
        return false;

    for(int i=0; i<numBins; i++)
    {
        if(h1->GetBinContent(i)!=h2->GetBinContent(i))
            return false;
    }
    return true;
}

int main(int argc, char **argv)
{
    if(argc == 3)
    {
        makeTurnOn(argv[1], argv[2]);
        return 0;
    }
    else if(argc == 4)
    {
        makeTurnOn(argv[1], argv[2], atoi(argv[3]));
        return 0;
    }
    else if(argc == 5)
    {
        makeTurnOn(argv[1], argv[2], atof(argv[3]), atoi(argv[4]));
        return 0;
    }
    else
    {
        return 1;
    }
}
