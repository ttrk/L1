/*
 * turn-on curves for the pPb MC samples. Samples are listed here :
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/PhotonAnalyses2014#pAWinter13_pPb_pythia_HIJING
 *
 * The selections are based on : 1) isolation
 *                               2) sideband
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

const int MAXL1EMCANDS = 144;
const int MAXL1REGIONS = 396;
const int MAXL1JETS = 8;
const int MAXPHOTONS = 500;
const Int_t THRESHOLDS = 1;
enum process {
    AllQCDPhotons,  //forward
    EmEnrichedDijets    // forward
};

void makeTurnOn(process hiForestProcess, TString outFileName, double offlineEtaCut=1.44)
{
    TFile *outFile = new TFile(outFileName,"RECREATE");

    TChain *f1Tree = new TChain("multiPhotonAnalyzer/photon","f1Tree");
    TChain *fEvtTree = new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
    TChain *fSkimTree = new TChain("skimanalysis/HltTree","fSkimTree");
    TChain *hltTree = new TChain("hltanalysis/HltTree","hltTree");

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

        //    double maxfpt = 0;
        //    //double maxfeta = -10;
        //    //double maxfphi = -10;
        //    double sigmaietaieta = -1;

        double maxfpt_iso = 0;
        double maxfpt_sideBand = 0;

        double sigmaietaieta_iso = -1;
        double sigmaietaieta_sideBand = -1;

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
                    }
                }
                // sideband selection
                if (photon_pt[i] > maxfpt_sideBand)  {
                    if(((cc4[i] + cr4[i] + ct4PtCut20[i])/0.9) > 10 && ((cc4[i] + cr4[i] + ct4PtCut20[i])/0.9) < 20)
                    {
                        maxfpt_sideBand=photon_pt[i];
                        sigmaietaieta_sideBand = sigmaIetaIeta[i];
                    }
                }
            }
        }

        fPt_iso[0]->Fill(maxfpt_iso,weight);
        fPt_sideBand[0]->Fill(maxfpt_sideBand,weight);
        if(hiBin < 60){
            fPt_iso[1]->Fill(maxfpt_iso,weight);
            fPt_sideBand[1]->Fill(maxfpt_sideBand,weight);
        }
        else if (hiBin >= 100) {
            fPt_iso[2]->Fill(maxfpt_iso,weight);
            fPt_sideBand[2]->Fill(maxfpt_sideBand,weight);
        }

        if(maxfpt_iso > 40)
        {
            fPtSigma_iso[0]->Fill(sigmaietaieta_iso,weight);
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
                accepted_iso[i][0]->Fill(maxfpt_iso,weight);
                accepted_sideBand[i][0]->Fill(maxfpt_sideBand,weight);
                if(hiBin < 60){
                    accepted_iso[i][1]->Fill(maxfpt_iso,weight);
                    accepted_sideBand[i][1]->Fill(maxfpt_sideBand,weight);
                }
                else if (hiBin >= 100){
                    accepted_iso[i][2]->Fill(maxfpt_iso,weight);
                    accepted_sideBand[i][2]->Fill(maxfpt_sideBand,weight);
                }

                if(maxfpt_iso > 40)
                {
                    accepted_sieie_iso[i][0]->Fill(sigmaietaieta_iso,weight);
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
            a_iso[k][l]->SetTitle("turn on for pT - isolated");

            a_sieie_iso[k][l] = new TGraphAsymmErrors();
            a_sieie_iso[k][l]->BayesDivide(accepted_sieie_iso[k][l], fPtSigma_iso[l]);
            a_sieie_iso[k][l]->SetName(Form("asymm_sieie_iso_%d_%d",k,l));
            a_sieie_iso[k][l]->SetTitle("turn on for sigmaIEtaEta - isolated");

            a_sideBand[k][l] = new TGraphAsymmErrors();
            a_sideBand[k][l]->BayesDivide(accepted_sideBand[k][l],fPt_sideBand[l]);
            a_sideBand[k][l]->SetName(Form("asymm_pt_sideBand_%d_%d",k,l));
            a_sideBand[k][l]->SetTitle("turn on for pT - sideband");

            a_sieie_sideBand[k][l] = new TGraphAsymmErrors();
            a_sieie_sideBand[k][l]->BayesDivide(accepted_sieie_sideBand[k][l], fPtSigma_sideBand[l]);
            a_sieie_sideBand[k][l]->SetName(Form("asymm_sieie_sideBand_%d_%d",k,l));
            a_sieie_sideBand[k][l]->SetTitle("turn on for sigmaIEtaEta - sideband");
        }
    }

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

    outFile->Close();
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
