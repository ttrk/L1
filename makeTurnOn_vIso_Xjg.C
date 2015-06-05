/*
 * turn-on curves for the pPb MC samples.
 *      1) turn-on curves for photon pT
 *      2) turn-on curves for sigmaIetaIeta

 * Xjg distributions for different isolation/sideband selections and signal/background regions
 *      1) isolation selection - signal region
 *      2) isolation selection - background region
 *      3) sideband selection - signal region
 *      4) sideband selection - background region
 *
 * Samples are listed here :
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/PhotonAnalyses2014#pAWinter13_pPb_pythia_HIJING
 *
 * New Samples are listed here :
 * https://twiki.cern.ch/twiki/bin/view/CMS/PhotonAnalyses2015#pPb_centrally_produced
 *
 * The selections are based on : 1) isolation
 *                               2) sideband
 *
 * */
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

// cut for signal/background
const float cut_sigmaIetaIeta = 0.01;

// cuts for jets
const float cut_jet_pt = 30;
const float cut_jet_eta = 1.6;
const float cut_jet_photon_deltaPhi = (7 * TMath::Pi()) / 8 ;   // 7/8 * TMath::Pi() --> evaluates to 0;

//const int MAXL1EMCANDS = 144;
//const int MAXL1REGIONS = 396;
//const int MAXL1JETS = 8;
const int MAXPHOTONS = 500;
const int MAXJETS = 500;
const Int_t THRESHOLDS = 1;
enum process {
    AllQCDPhotons,  //forward
    EmEnrichedDijets,    // forward
    AllQCDPhotons_NEW,  //forward
    EmEnrichedDijets_NEW    // forward
};
const char* processNames[] = {"AllQCDPhotons", "EmEnrichedDijets", "AllQCDPhotons v3", "EmEnrichedDijets v3"};

using namespace std;

Double_t getDPHI( Double_t phi1, Double_t phi2);

void makeTurnOn(process hiForestProcess, TString outFileName, double offlineEtaCut=1.44)
{
    TFile *outFile = new TFile(outFileName,"RECREATE");

    TChain *f1Tree = new TChain("multiPhotonAnalyzer/photon","f1Tree");
    TChain *fEvtTree = new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
    TChain *fSkimTree = new TChain("skimanalysis/HltTree","fSkimTree");
    TChain *hltTree = new TChain("hltanalysis/HltTree","hltTree");
    TChain *f1JetTree = new TChain("akPu3PFJetAnalyzer/t","f1JetTree");   // must decide on the tree accroding to collision type

    const char* processName = processNames[hiForestProcess];

    // prepare file names and weights according to process type
    int MAXFILES=5;
    const char* fileNames[MAXFILES];
    double weights[MAXFILES];
    if(hiForestProcess == AllQCDPhotons)
    {
        fileNames[0]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton30.root";
        fileNames[1]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton50.root";
        fileNames[2]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton80.root";
        fileNames[3]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton120.root";
        fileNames[4]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_AllQCDPhoton170.root";

        weights[0]=62744./62744. ;
        weights[1]=29499./107309.;
        weights[2]=7640. /106817.;
        weights[3]=1868. /104443.;
        weights[4]=649.  /139647.;
    }
    else if(hiForestProcess == EmEnrichedDijets)
    {
        fileNames[0]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_EmEnriched30.root";
        fileNames[1]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_EmEnriched50.root";
        fileNames[2]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_EmEnriched80.root";
        fileNames[3]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_EmEnriched120.root";
        fileNames[4]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v1/HiForest_pPb_MIX_EmEnriched170.root";

        weights[0]=62744./62744. ;
        weights[1]=29499./107309.;
        weights[2]=7640. /106817.;
        weights[3]=1868. /104443.;
        weights[4]=649.  /139647.;
    }
    else if(hiForestProcess == AllQCDPhotons_NEW)
    {
        fileNames[0]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_AllQCDPhoton30_localJEC_v3.root";
        fileNames[1]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_AllQCDPhoton50_localJEC_v3.root";
        fileNames[2]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_AllQCDPhoton80_localJEC_v3.root";
        fileNames[3]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_AllQCDPhoton120_localJEC_v3.root";
        fileNames[4]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_AllQCDPhoton170_localJEC_v3.root";

        weights[0]=62744./62744. ;
        weights[1]=29499./107309.;
        weights[2]=7640. /106817.;
        weights[3]=1868. /104443.;
        weights[4]=649.  /139647.;
    }
    else if(hiForestProcess == EmEnrichedDijets_NEW)
    {
        fileNames[0]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_EmEnriched30_localJEC_v3.root";
        fileNames[1]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_EmEnriched50_localJEC_v3.root";
        fileNames[2]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_EmEnriched80_localJEC_v3.root";
        fileNames[3]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_EmEnriched120_localJEC_v3.root";
        fileNames[4]="/mnt/hadoop/cms/store/user/luck/2014-photon-forests/pPb_MIX_localJEC_v3/pPb_MIX_EmEnriched170_localJEC_v3.root";

        weights[0]=62744./62744. ;
        weights[1]=29499./107309.;
        weights[2]=7640. /106817.;
        weights[3]=1868. /104443.;
        weights[4]=649.  /139647.;
    }
    else
    {
        std::cout << "Please enter a valid process type"<<std::endl;
        exit(1);
    }

    for(int i=0; i<MAXFILES;i++)
    {
        f1Tree->Add(fileNames[i]);
        fEvtTree->Add(fileNames[i]);
        fSkimTree->Add(fileNames[i]);
        hltTree->Add(fileNames[i]);
        f1JetTree->Add(fileNames[i]);
    }
    // files are added to the chain.

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

    Int_t nJet;
    Float_t jet_pt[MAXJETS];
    Float_t jet_eta[MAXJETS];
    Float_t jet_phi[MAXJETS];

    f1JetTree->SetBranchAddress("nref",&nJet);
    f1JetTree->SetBranchAddress("jtpt",jet_pt);
    f1JetTree->SetBranchAddress("jteta",jet_eta);
    f1JetTree->SetBranchAddress("jtphi",jet_phi);

    const int nBins = 100;
    const int maxPt = 100;

    TH1::SetDefaultSumw2();
    // tree names that go like fPt* : events that pass selection

    //  TH1D *fPt[3];
    //  fPt[0] = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
    //  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_1");
    //  fPt[2] = (TH1D*)fPt[0]->Clone("fPt_2");

    //  TH1D *fPtSigma[3];
    //  fPtSigma[0] = new TH1D("fPtSigma_0",";#sigma_{#eta #eta}", 100, 0, 0.1);
    //  fPtSigma[1] = (TH1D*)fPtSigma[0]->Clone("fPtSigma_1");
    //  fPtSigma[2] = (TH1D*)fPtSigma[0]->Clone("fPtSigma_2");

    //  TH1D *accepted[THRESHOLDS][3];
    //  TH1D *accepted_sieie[THRESHOLDS][3];

    TH1D *fPt_iso[3];
    fPt_iso[0] = new TH1D("fPt_iso_0","isolated;offline p_{T} (GeV)",nBins,0,maxPt);
    fPt_iso[1] = (TH1D*)fPt_iso[0]->Clone("fPt_iso_1");
    fPt_iso[2] = (TH1D*)fPt_iso[0]->Clone("fPt_iso_2");

    TH1D *fPtSigma_iso[3];
    fPtSigma_iso[0] = new TH1D("fPtSigma_iso_0","isolated;#sigma_{#eta #eta}", 100, 0, 0.1);
    fPtSigma_iso[1] = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_1");
    fPtSigma_iso[2] = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_2");

    TH1D *fPt_sideBand[3];
    fPt_sideBand[0] = new TH1D("fPt_sideBand_0","sideband;offline p_{T} (GeV)",nBins,0,maxPt);
    fPt_sideBand[1] = (TH1D*)fPt_sideBand[0]->Clone("fPt_sideBand_1");
    fPt_sideBand[2] = (TH1D*)fPt_sideBand[0]->Clone("fPt_sideBand_2");

    TH1D *fPtSigma_sideBand[3];
    fPtSigma_sideBand[0] = new TH1D("fPtSigma_sideBand_0","sideband;#sigma_{#eta #eta}", 100, 0, 0.1);
    fPtSigma_sideBand[1] = (TH1D*)fPtSigma_sideBand[0]->Clone("fPtSigma_sideBand_1");
    fPtSigma_sideBand[2] = (TH1D*)fPtSigma_sideBand[0]->Clone("fPtSigma_sideBand_2");

    // Xjg distributions
    // isolation selection - signal region
    TH1D *fXjg_iso_signal[3];
    fXjg_iso_signal[0] = new TH1D("fXjg_iso_signal_0",";X_{J #gamma}", 100, 0, 2);
    fXjg_iso_signal[0]->SetTitle(Form("%s - Xjg - isolation selection in signal region",processName));
    fXjg_iso_signal[1] = (TH1D*)fXjg_iso_signal[0]->Clone("fXjg_iso_signal_1");
    fXjg_iso_signal[2] = (TH1D*)fXjg_iso_signal[0]->Clone("fXjg_iso_signal_2");

    // isolation selection - background region
    TH1D *fXjg_iso_bkg[3];
    fXjg_iso_bkg[0] = new TH1D("fXjg_iso_bkg_0",";X_{J #gamma}", 100, 0, 2);
    fXjg_iso_bkg[0]->SetTitle(Form("%s - Xjg - isolation selection in background region",processName));
    fXjg_iso_bkg[1] = (TH1D*)fXjg_iso_bkg[0]->Clone("fXjg_iso_bkg_1");
    fXjg_iso_bkg[2] = (TH1D*)fXjg_iso_bkg[0]->Clone("fXjg_iso_bkg_2");

    // sideband selection - signal region
    TH1D *fXjg_sideBand_signal[3];
    fXjg_sideBand_signal[0] = new TH1D("fXjg_sideBand_signal_0",";X_{J #gamma}", 100, 0, 2);
    fXjg_sideBand_signal[0]->SetTitle(Form("%s - Xjg - sideband selection in signal region",processName));
    fXjg_sideBand_signal[1] = (TH1D*)fXjg_sideBand_signal[0]->Clone("fXjg_sideBand_signal_1");
    fXjg_sideBand_signal[2] = (TH1D*)fXjg_sideBand_signal[0]->Clone("fXjg_sideBand_signal_2");

    // sideband selection - background region
    TH1D *fXjg_sideBand_bkg[3];
    fXjg_sideBand_bkg[0] = new TH1D("fXjg_sideBand_bkg_0",";X_{J #gamma}", 100, 0, 2);
    fXjg_sideBand_bkg[0]->SetTitle(Form("%s - Xjg - sideband selection in background region",processName));
    fXjg_sideBand_bkg[1] = (TH1D*)fXjg_sideBand_bkg[0]->Clone("fXjg_sideBand_bkg_1");
    fXjg_sideBand_bkg[2] = (TH1D*)fXjg_sideBand_bkg[0]->Clone("fXjg_sideBand_bkg_2");

    // isolation selection - signal region - HLT fired
    TH1D *fXjg_iso_signal_HLT1[3];
    fXjg_iso_signal_HLT1[0] = new TH1D("fXjg_iso_signal_0_HLT1",";X_{J #gamma}", 100, 0, 2);
    fXjg_iso_signal_HLT1[0]->SetTitle(Form("%s - Xjg - isolation selection in signal region - HLT triggered",processName));

    // isolation selection - background region - HLT fired
    TH1D *fXjg_iso_bkg_HLT1[1];
    fXjg_iso_bkg_HLT1[0] = new TH1D("fXjg_iso_bkg_0_HLT1",";X_{J #gamma}", 100, 0, 2);
    fXjg_iso_bkg_HLT1[0]->SetTitle(Form("%s - Xjg - isolation selection in background region - HLT triggered",processName));

    // sideband selection - signal region - HLT fired
    TH1D *fXjg_sideBand_signal_HLT1[1];
    fXjg_sideBand_signal_HLT1[0] = new TH1D("fXjg_sideBand_signal_0_HLT1",";X_{J #gamma}", 100, 0, 2);
    fXjg_sideBand_signal_HLT1[0]->SetTitle(Form("%s - Xjg - sideband selection in signal region - HLT triggered",processName));

    // sideband selection - background region - HLT fired
    TH1D *fXjg_sideBand_bkg_HLT1[1];
    fXjg_sideBand_bkg_HLT1[0] = new TH1D("fXjg_sideBand_bkg_0_HLT1",";X_{J #gamma}", 100, 0, 2);
    fXjg_sideBand_bkg_HLT1[0]->SetTitle(Form("%s - Xjg - sideband selection in background region - HLT triggered",processName));

    // correlation plots
    TH2D* fPt_Sigma_iso_corr = new TH2D("fPt_Sigma_iso_corr", "isolated;p_{T}^{#gamma};#sigma_{#eta #eta}", 60,40,100, 40,0,0.04);
    TH2D* fPt_Sigma_iso_corr_HLT1 = new TH2D("fPt_Sigma_iso_corr_HLT1", "isolated - HLTtriggered;p_{T}^{#gamma};#sigma_{#eta #eta}", 60,40,100, 40,0,0.04);

    TH2D* fPt_Sigma_iso_corr_jet = new TH2D("fPt_Sigma_iso_corr_jet", "isolated;p_{T}^{Jet};#sigma_{#eta #eta}", 170,30,200, 40,0,0.04);
    TH2D* fPt_Sigma_iso_corr_jet_HLT1 = new TH2D("fPt_Sigma_iso_corr_jet_HLT1", "isolated - HLTtriggered;p_{T}^{Jet};#sigma_{#eta #eta}", 170,30,200, 40,0,0.04);

    TH2D* fPt_Sigma_sideBand_corr = new TH2D("fPt_Sigma_sideBand_corr", "sideband;p_{T}^{#gamma};#sigma_{#eta #eta}", 60,40,100, 40,0,0.04);
    TH2D* fPt_Sigma_sideBand_corr_HLT1 = new TH2D("fPt_Sigma_sideBand_corr_HLT1", "sideband - HLTtriggered;p_{T}^{#gamma};#sigma_{#eta #eta}", 60,40,100, 40,0,0.04);

    TH2D* fPt_Sigma_sideBand_corr_jet = new TH2D("fPt_Sigma_sideBand_corr_jet", "sideband;p_{T}^{Jet};#sigma_{#eta #eta}", 170,30,200, 40,0,0.04);
    TH2D* fPt_Sigma_sideBand_corr_jet_HLT1 = new TH2D("fPt_Sigma_sideBand_corr_jet_HLT1", "sideband - HLTtriggered;p_{T}^{Jet};#sigma_{#eta #eta}", 170,30,200, 40,0,0.04);

    // tree names that go like accepted* : events that pass selection and that are triggered
    TH1D *accepted_iso[THRESHOLDS][3];
    TH1D *accepted_sieie_iso[THRESHOLDS][3];

    TH1D *accepted_sideBand[THRESHOLDS][3];
    TH1D *accepted_sieie_sideBand[THRESHOLDS][3];

    for(int i = 0; i < THRESHOLDS; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            //      accepted[i][j] = new TH1D(Form("accepted_pt%d_%d",i,j),";offline p_{T}",nBins,0,maxPt);
            //      accepted_sieie[i][j] = new TH1D(Form("accepted_sieie%d_%d",i,j), ";#sigma_{#eta #eta}",100,0,0.1);

            accepted_iso[i][j] = new TH1D(Form("accepted_pt_iso%d_%d",i,j),"isolated;offline p_{T}",nBins,0,maxPt);
            accepted_sieie_iso[i][j] = new TH1D(Form("accepted_sieie_iso%d_%d",i,j), "isolated;#sigma_{#eta #eta}",100,0,0.1);

            accepted_sideBand[i][j] = new TH1D(Form("accepted_pt_sideBand%d_%d",i,j),"sideband;offline p_{T}",nBins,0,maxPt);
            accepted_sieie_sideBand[i][j] = new TH1D(Form("accepted_sieie_sideBand%d_%d",i,j), "sideband;#sigma_{#eta #eta}",100,0,0.1);
        }
    }

    // loop over all events in all files
    int fileEntries[MAXFILES];
    int treeIndex;
    double weight;
    Long64_t entries = f1Tree->GetEntries();
    for(Long64_t j = 0; j < entries; ++j)
    {
        if(j % 10000 == 0)
            printf("%lld / %lld\n",j,entries);

        fEvtTree->GetEntry(j);
        fSkimTree->GetEntry(j);

        treeIndex=fEvtTree->GetTreeNumber();
        weight = weights[treeIndex];
        fileEntries[treeIndex]=fEvtTree->GetTree()->GetEntries();

        // bool goodEvent = false;
        // if((pcollisionEventSelection == 1) && (TMath::Abs(vz) < 15))
        // {
        //   goodEvent = true;
        // }
        // if(!goodEvent) continue;

        hltTree->GetEntry(j);
        f1Tree->GetEntry(j);
        f1JetTree->GetEntry(j);

        //    double maxfpt = 0;
        //    //double maxfeta = -10;
        //    //double maxfphi = -10;
        //    double sigmaietaieta = -1;

        double maxfpt_iso = 0;
        double maxfpt_sideBand = 0;

        double sigmaietaieta_iso = -1;
        double sigmaietaieta_sideBand = -1;

        int index_maxfpt_iso = -1;
        int index_maxfpt_sideBand = -1;

        for(int i = 0; i < nPhoton; ++i)
        {
            if(TMath::Abs(photon_eta[i]) < offlineEtaCut)
            if(!isEle[i])
            if(TMath::Abs(seedTime[i])<3)
            if(swissCrx[i] < 0.9)
            if(sigmaIetaIeta[i] > 0.002)
            if(sigmaIphiIphi[i] > 0.002)
            if(hadronicOverEm[i] < 0.1)
            {
                // isolation selection
                if(photon_pt[i] > maxfpt_iso) {
                    if(ecalRecHitSumEtConeDR04[i] < 4.2 && hcalTowerSumEtConeDR04[i] < 2.2 && trkSumPtHollowConeDR04[i] < 2)
                    {
                        maxfpt_iso=photon_pt[i];
                        sigmaietaieta_iso = sigmaIetaIeta[i];
                        index_maxfpt_iso=i;
                    }
                }
                // sideband selection
                if (photon_pt[i] > maxfpt_sideBand)  {
                    if(((cc4[i] + cr4[i] + ct4PtCut20[i])/0.9) > 10 && ((cc4[i] + cr4[i] + ct4PtCut20[i])/0.9) < 20)
                    {
                        maxfpt_sideBand=photon_pt[i];
                        sigmaietaieta_sideBand = sigmaIetaIeta[i];
                        index_maxfpt_sideBand=i;
                    }
                }
            }
        }

        double maxfpt_jet_iso_signal = 0;   // sigmaIetaIeta < 0.01
        double maxfpt_jet_iso_bkg = 0;      // sigmaIetaIeta > 0.01
        double maxfpt_jet_iso = 0;          // will be used in jet-sigmaIetaIeta TH2 plots

        double Xjg_iso_signal = -1;   // sigmaIetaIeta < 0.01
        double Xjg_iso_bkg = -1;      // sigmaIetaIeta > 0.01

        // match jets to photons from isolation selection
        for(int i = 0; i < nJet; ++i)
        {
            if(index_maxfpt_iso < 0)
                break;      // do not loop over jets, no photon was selected

            if(TMath::Abs(jet_eta[i]) < cut_jet_eta)
            if(TMath::Abs(getDPHI(jet_phi[i],photon_phi[index_maxfpt_iso])) > cut_jet_photon_deltaPhi)
            if(jet_pt[i] > cut_jet_pt)
            {
                // signal region
                if(jet_pt[i] > maxfpt_jet_iso_signal) {
                    if(sigmaIetaIeta[index_maxfpt_iso] < cut_sigmaIetaIeta){
                        maxfpt_jet_iso_signal=jet_pt[i];
                        Xjg_iso_signal = maxfpt_jet_iso_signal / maxfpt_iso;
                    }
                }
                // background region
                if (jet_pt[i] > maxfpt_jet_iso_bkg)  {
                    if(sigmaIetaIeta[index_maxfpt_iso] > cut_sigmaIetaIeta){
                        maxfpt_jet_iso_bkg=jet_pt[i];
                        Xjg_iso_bkg = maxfpt_jet_iso_bkg / maxfpt_iso;
                    }
                }
           }
        }
        if(maxfpt_jet_iso_signal >= maxfpt_jet_iso_bkg) {
            maxfpt_jet_iso = maxfpt_jet_iso_signal;
        }
        else {
            maxfpt_jet_iso = maxfpt_jet_iso_bkg;
        }

        // Xjg
        double maxfpt_jet_sideBand_signal = 0;   // sigmaIetaIeta < 0.01
        double maxfpt_jet_sideBand_bkg = 0;      // sigmaIetaIeta > 0.01
        double maxfpt_jet_sideBand = 0;          // will be used in jet-sigmaIetaIeta TH2 plots

        double Xjg_sideBand_signal = -1;   // sigmaIetaIeta < 0.01
        double Xjg_sideBand_bkg = -1;      // sigmaIetaIeta > 0.01

        // match jets to photons from sideBand selection
        for(int i = 0; i < nJet; ++i)
        {
            if(index_maxfpt_sideBand < 0)
                break;      // do not loop over jets, no photon was selected

            if(TMath::Abs(jet_eta[i]) < cut_jet_eta)
            if(TMath::Abs(getDPHI(jet_phi[i],photon_phi[index_maxfpt_sideBand])) > cut_jet_photon_deltaPhi)
            if(jet_pt[i] > cut_jet_pt)
            {
                // signal region
                if(jet_pt[i] > maxfpt_jet_sideBand_signal) {
                    if(sigmaIetaIeta[index_maxfpt_sideBand] < cut_sigmaIetaIeta){
                        maxfpt_jet_sideBand_signal=jet_pt[i];
                        Xjg_sideBand_signal = maxfpt_jet_sideBand_signal / maxfpt_sideBand;
                    }
                }
                // background region
                if (jet_pt[i] > maxfpt_jet_sideBand_bkg)  {
                    if(sigmaIetaIeta[index_maxfpt_sideBand] > cut_sigmaIetaIeta){
                        maxfpt_jet_sideBand_bkg=jet_pt[i];
                        Xjg_sideBand_bkg = maxfpt_jet_sideBand_bkg / maxfpt_sideBand;
                    }
                }
           }
        }
        if(maxfpt_jet_sideBand_signal >= maxfpt_jet_sideBand_bkg) {
            maxfpt_jet_sideBand = maxfpt_jet_sideBand_signal;
        }
        else {
            maxfpt_jet_sideBand = maxfpt_jet_sideBand_bkg;
        }

        // fill Xjg distributions
        if(!(index_maxfpt_iso < 0) && maxfpt_iso > 40)  // there is an isolated photon in the event and it passes pt cut
        {
            if(Xjg_iso_signal!=-1) fXjg_iso_signal[0]->Fill(Xjg_iso_signal,weight);
            if(Xjg_iso_bkg!=-1)    fXjg_iso_bkg[0]->Fill(Xjg_iso_bkg,weight);
            if(maxfpt_jet_iso!=0)  fPt_Sigma_iso_corr_jet->Fill(maxfpt_jet_iso,sigmaIetaIeta[index_maxfpt_iso],weight);
            if(HLT_PAPhoton30_NoCaloIdVL_v1)
            {
                if(Xjg_iso_signal!=-1) fXjg_iso_signal_HLT1[0]->Fill(Xjg_iso_signal,weight);
                if(Xjg_iso_bkg!=-1)    fXjg_iso_bkg_HLT1[0]->Fill(Xjg_iso_bkg,weight);
                if(maxfpt_jet_iso!=0)  fPt_Sigma_iso_corr_jet_HLT1->Fill(maxfpt_jet_iso,sigmaIetaIeta[index_maxfpt_iso],weight);
            }
            if(hiBin < 60){
                if(Xjg_iso_signal!=-1) fXjg_iso_signal[1]->Fill(Xjg_iso_signal,weight);
                if(Xjg_iso_bkg!=-1)    fXjg_iso_bkg[1]->Fill(Xjg_iso_bkg,weight);
            }
            else if (hiBin >= 100) {
                if(Xjg_iso_signal!=-1) fXjg_iso_signal[2]->Fill(Xjg_iso_signal,weight);
                if(Xjg_iso_bkg!=-1)    fXjg_iso_bkg[2]->Fill(Xjg_iso_bkg,weight);
            }
        }
        if(!(index_maxfpt_sideBand < 0) && maxfpt_sideBand > 40)  // there is a sideband photon in the event and it passes pt cut
        {
            if(Xjg_sideBand_signal!=-1) fXjg_sideBand_signal[0]->Fill(Xjg_sideBand_signal,weight);
            if(Xjg_sideBand_bkg!=-1)    fXjg_sideBand_bkg[0]->Fill(Xjg_sideBand_bkg,weight);
            if(maxfpt_jet_sideBand!=0)  fPt_Sigma_sideBand_corr_jet->Fill(maxfpt_jet_sideBand,sigmaIetaIeta[index_maxfpt_sideBand],weight);
            if(HLT_PAPhoton30_NoCaloIdVL_v1)
            {
                if(Xjg_sideBand_signal!=-1) fXjg_sideBand_signal_HLT1[0]->Fill(Xjg_sideBand_signal,weight);
                if(Xjg_sideBand_bkg!=-1)    fXjg_sideBand_bkg_HLT1[0]->Fill(Xjg_sideBand_bkg,weight);
                if(maxfpt_jet_sideBand!=0)  fPt_Sigma_sideBand_corr_jet_HLT1->Fill(maxfpt_jet_sideBand,sigmaIetaIeta[index_maxfpt_sideBand],weight);
            }
            if(hiBin < 60){
                if(Xjg_sideBand_signal!=-1) fXjg_sideBand_signal[1]->Fill(Xjg_sideBand_signal,weight);
                if(Xjg_sideBand_bkg!=-1)    fXjg_sideBand_bkg[1]->Fill(Xjg_sideBand_bkg,weight);
            }
            else if (hiBin >= 100) {
                if(Xjg_sideBand_signal!=-1) fXjg_sideBand_signal[2]->Fill(Xjg_sideBand_signal,weight);
                if(Xjg_sideBand_bkg!=-1)    fXjg_sideBand_bkg[2]->Fill(Xjg_sideBand_bkg,weight);
            }
        }

        if(!(index_maxfpt_iso<0))   //  there is an isolated photon in the event
        {
            fPt_iso[0]->Fill(maxfpt_iso,weight);
            if(hiBin < 60){
                fPt_iso[1]->Fill(maxfpt_iso,weight);
            }
            else if (hiBin >= 100) {
                fPt_iso[2]->Fill(maxfpt_iso,weight);
            }
        }

        if(!(index_maxfpt_sideBand<0))   //  there is a sideband photon in the event
        {
            fPt_sideBand[0]->Fill(maxfpt_sideBand,weight);
            if(hiBin < 60){
                fPt_sideBand[1]->Fill(maxfpt_sideBand,weight);
            }
            else if (hiBin >= 100) {
                fPt_sideBand[2]->Fill(maxfpt_sideBand,weight);
            }
        }

        if(maxfpt_iso > 40)
        {
            fPtSigma_iso[0]->Fill(sigmaietaieta_iso,weight);
            fPt_Sigma_iso_corr->Fill(maxfpt_iso, sigmaietaieta_iso, weight);    // corr : photon pt - sigmaIetaIeta
            if(hiBin < 60){
                fPtSigma_iso[1]->Fill(sigmaietaieta_iso,weight);
            }
            else if (hiBin >= 100) {
                fPtSigma_iso[2]->Fill(sigmaietaieta_iso,weight);
            }
        }

        if(maxfpt_sideBand > 40)
        {
            fPtSigma_sideBand[0]->Fill(sigmaietaieta_sideBand,weight);
            fPt_Sigma_sideBand_corr->Fill(maxfpt_sideBand, sigmaietaieta_sideBand, weight);    // corr : photon pt - sigmaIetaIeta
            if(hiBin < 60){
                fPtSigma_sideBand[1]->Fill(sigmaietaieta_sideBand,weight);
            }
            else if (hiBin >= 100) {
                fPtSigma_sideBand[2]->Fill(sigmaietaieta_sideBand,weight);
            }
        }

        for(int i = 0; i < THRESHOLDS; ++i)
        {
            if(HLT_PAPhoton30_NoCaloIdVL_v1)
            {
                if(!(index_maxfpt_iso<0)) accepted_iso[i][0]->Fill(maxfpt_iso,weight);
                if(!(index_maxfpt_sideBand<0)) accepted_sideBand[i][0]->Fill(maxfpt_sideBand,weight);
                if(hiBin < 60){
                    if(!(index_maxfpt_iso<0)) accepted_iso[i][1]->Fill(maxfpt_iso,weight);
                    if(!(index_maxfpt_sideBand<0)) accepted_sideBand[i][1]->Fill(maxfpt_sideBand,weight);
                }
                else if (hiBin >= 100){
                    if(!(index_maxfpt_iso<0)) accepted_iso[i][2]->Fill(maxfpt_iso,weight);
                    if(!(index_maxfpt_sideBand<0)) accepted_sideBand[i][2]->Fill(maxfpt_sideBand,weight);
                }

                if(maxfpt_iso > 40)
                {
                    accepted_sieie_iso[i][0]->Fill(sigmaietaieta_iso,weight);
                    fPt_Sigma_iso_corr_HLT1->Fill(maxfpt_iso, sigmaietaieta_iso, weight);    // corr : photon pt - sigmaIetaIeta
                    if(hiBin < 60){
                        accepted_sieie_iso[i][1]->Fill(sigmaietaieta_iso,weight);
                    }
                    else if (hiBin >= 100){
                        accepted_sieie_iso[i][2]->Fill(sigmaietaieta_iso,weight);
                    }
                }

                if(maxfpt_sideBand > 40)
                {
                    accepted_sieie_sideBand[i][0]->Fill(sigmaietaieta_sideBand,weight);
                    fPt_Sigma_sideBand_corr_HLT1->Fill(maxfpt_sideBand, sigmaietaieta_sideBand, weight);    // corr : photon pt - sigmaIetaIeta
                    if(hiBin < 60){
                        accepted_sieie_sideBand[i][1]->Fill(sigmaietaieta_sideBand,weight);
                    }
                    else if (hiBin >= 100){
                        accepted_sieie_sideBand[i][2]->Fill(sigmaietaieta_sideBand,weight);
                    }
                }

            }
        }
    }
    // print weights and number of events in files
    double totWeightedEntries=0;
    for(int i=0; i<MAXFILES; i++)
    {
        totWeightedEntries +=fileEntries[i]*weights[i];
        std::cout<<"i = "<<i<<" , entries = "<<fileEntries[i]<<" , weight = "<<weights[i]<<std::endl;
    }
    std::cout<<"total weighted entries = "<<totWeightedEntries<<std::endl;

    std::cout<<"writing files"<<std::endl;

    // add full photon tree chain
    f1Tree->SetNameTitle("photon","full Photon TChain");
    f1Tree->Write();

    //  TGraphAsymmErrors *a[THRESHOLDS][3];
    //  TGraphAsymmErrors *a_sieie[THRESHOLDS][3];

    TGraphAsymmErrors *a_iso[THRESHOLDS][3];
    TGraphAsymmErrors *a_sieie_iso[THRESHOLDS][3];

    TGraphAsymmErrors *a_sideBand[THRESHOLDS][3];
    TGraphAsymmErrors *a_sieie_sideBand[THRESHOLDS][3];
    for(int k = 0; k < THRESHOLDS; ++k){
        for(int l = 0; l < 3; ++l)
        {
            //      a[k][l] = new TGraphAsymmErrors();
            //      a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
            //      a[k][l]->SetName(Form("asymm_pt_%d_%d",k,l));
            //
            //      a_sieie[k][l] = new TGraphAsymmErrors();
            //      a_sieie[k][l]->BayesDivide(accepted_sieie[k][l], fPtSigma[l]);
            //      a_sieie[k][l]->SetName(Form("asymm_sieie_%d_%d",k,l));

            a_iso[k][l] = new TGraphAsymmErrors();
            a_iso[k][l]->BayesDivide(accepted_iso[k][l],fPt_iso[l]);
            a_iso[k][l]->SetName(Form("asymm_pt_iso_%d_%d",k,l));
            a_iso[k][l]->SetTitle(Form("%s - turn on for pT - isolated",processName));

            a_sieie_iso[k][l] = new TGraphAsymmErrors();
            a_sieie_iso[k][l]->BayesDivide(accepted_sieie_iso[k][l], fPtSigma_iso[l]);
            a_sieie_iso[k][l]->SetName(Form("asymm_sieie_iso_%d_%d",k,l));
            a_sieie_iso[k][l]->SetTitle(Form("%s - turn on for sigmaIEtaEta - isolated", processName));

            a_sideBand[k][l] = new TGraphAsymmErrors();
            a_sideBand[k][l]->BayesDivide(accepted_sideBand[k][l],fPt_sideBand[l]);
            a_sideBand[k][l]->SetName(Form("asymm_pt_sideBand_%d_%d",k,l));
            a_sideBand[k][l]->SetTitle(Form("%s - turn on for pT - sideband", processName));

            a_sieie_sideBand[k][l] = new TGraphAsymmErrors();
            a_sieie_sideBand[k][l]->BayesDivide(accepted_sieie_sideBand[k][l], fPtSigma_sideBand[l]);
            a_sieie_sideBand[k][l]->SetName(Form("asymm_sieie_sideBand_%d_%d",k,l));
            a_sieie_sideBand[k][l]->SetTitle(Form("%s - turn on for sigmaIEtaEta - sideband",processName));
        }
    }
    TGraphAsymmErrors* a_Xjg_iso_bkg = new TGraphAsymmErrors();
    a_Xjg_iso_bkg->SetName("asymm_Xjg_iso_bkg");
    a_Xjg_iso_bkg->SetTitle(Form("%s - turn on for Xjg - isolated - background region",processName));
    TGraphAsymmErrors* a_Xjg_sideBand_bkg = new TGraphAsymmErrors();
    a_Xjg_sideBand_bkg->SetName("asymm_Xjg_sideBand_bkg");
    a_Xjg_sideBand_bkg->SetTitle(Form("%s - turn on for Xjg - sideband - background region",processName));

    TH1D* clone1_fXjg_iso_bkg=(TH1D*)fXjg_iso_bkg[0]->Clone();
    TH1D* clone1_fXjg_iso_bkg_HLT1=(TH1D*)fXjg_iso_bkg_HLT1[0]->Clone();
    TH1D* clone1_fXjg_sideBand_bkg=(TH1D*)fXjg_sideBand_bkg[0]->Clone();
    TH1D* clone1_fXjg_sideBand_bkg_HLT1=(TH1D*)fXjg_sideBand_bkg_HLT1[0]->Clone();

    clone1_fXjg_iso_bkg->Rebin(2);
    clone1_fXjg_iso_bkg_HLT1->Rebin(2);
    clone1_fXjg_sideBand_bkg->Rebin(2);
    clone1_fXjg_sideBand_bkg_HLT1->Rebin(2);

    a_Xjg_iso_bkg->BayesDivide(clone1_fXjg_iso_bkg_HLT1, clone1_fXjg_iso_bkg);
    a_Xjg_sideBand_bkg->BayesDivide(clone1_fXjg_sideBand_bkg_HLT1, clone1_fXjg_sideBand_bkg);

    TCanvas* c1 = new TCanvas();
    c1->SetTitle(outFileName.Data());
    c1->Divide(2,2);
    c1->cd(1);
    a_sieie_iso[0][0]->Draw("a p");
    c1->cd(3);
    a_sieie_sideBand[0][0]->Draw("a p");
    c1->cd(2);
    fPtSigma_iso[0]->Draw();
    c1->cd(4);
    fPtSigma_sideBand[0]->Draw();

    outFile->cd();

    //  fPt[0]->Write();
    //  fPt[1]->Write();
    //  fPt[2]->Write();
    c1->Write();

    fPt_iso[0]->Write();
    fPt_iso[1]->Write();
    fPt_iso[2]->Write();

    fPt_sideBand[0]->Write();
    fPt_sideBand[1]->Write();
    fPt_sideBand[2]->Write();

    fPtSigma_iso[0]->Write();
    fPtSigma_iso[1]->Write();
    fPtSigma_iso[2]->Write();

    fPtSigma_sideBand[0]->Write();
    fPtSigma_sideBand[1]->Write();
    fPtSigma_sideBand[2]->Write();

    for(int k = 0; k < THRESHOLDS; ++k){
        for(int l = 0; l < 3; ++l)
        {
            //      accepted[k][l]->Write();
            //      accepted_sieie[k][l]->Write();

            accepted_iso[k][l]->Write();
            accepted_sieie_iso[k][l]->Write();

            accepted_sideBand[k][l]->Write();
            accepted_sieie_sideBand[k][l]->Write();
        }
    }

    for(int k = 0; k < THRESHOLDS; ++k){
        for(int l = 0; l < 3; ++l)
        {
            //      a[k][l]->Write();
            //      a_sieie[k][l]->Write();

            a_iso[k][l]->Write();
            a_sieie_iso[k][l]->Write();

            a_sideBand[k][l]->Write();
            a_sieie_sideBand[k][l]->Write();
        }
    }
    a_Xjg_iso_bkg->Write();
    a_Xjg_sideBand_bkg->Write();

    // write Xjg distributions
    fXjg_iso_signal[0]->Write();
    fXjg_iso_signal[1]->Write();
    fXjg_iso_signal[2]->Write();
    fXjg_iso_bkg[0]->Write();
    fXjg_iso_bkg[1]->Write();
    fXjg_iso_bkg[2]->Write();
    fXjg_sideBand_signal[0]->Write();
    fXjg_sideBand_signal[1]->Write();
    fXjg_sideBand_signal[2]->Write();
    fXjg_sideBand_bkg[0]->Write();
    fXjg_sideBand_bkg[1]->Write();
    fXjg_sideBand_bkg[2]->Write();

    fXjg_iso_signal_HLT1[0]->Write();
    fXjg_iso_bkg_HLT1[0]->Write();
    fXjg_sideBand_signal_HLT1[0]->Write();
    fXjg_sideBand_bkg_HLT1[0]->Write();

    fPt_Sigma_iso_corr->Write();
    fPt_Sigma_iso_corr_HLT1->Write();
    fPt_Sigma_iso_corr_jet->Write();
    fPt_Sigma_iso_corr_jet_HLT1->Write();

    fPt_Sigma_sideBand_corr->Write();
    fPt_Sigma_sideBand_corr_HLT1->Write();
    fPt_Sigma_sideBand_corr_jet->Write();
    fPt_Sigma_sideBand_corr_jet_HLT1->Write();

    outFile->Close();
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;

  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }

  return dphi;
}

int main(int argc, char **argv)
{
    if(argc == 3)
    {
        makeTurnOn((process)atoi(argv[1]), argv[2]);
        return 0;
    }
    else if(argc == 4)
    {
        makeTurnOn((process)atoi(argv[1]), argv[2], atof(argv[3]));
        return 0;
    }
    else
    {
        return 1;
    }
}
