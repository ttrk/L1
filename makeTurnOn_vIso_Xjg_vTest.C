/*
 * The whole purpose of this file is to test whether it is possible to plot maximum element properties
 *  without looping over a tree. That is whether it is possible to implement the selections in the loop
 *  inside TTree::Draw() function.
 *
 * Samples are listed here :
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
#include "/net/hisrv0001/home/tatar/code/HIUtils/treeUtil.h"  // tested and works

#include <vector>
#include <iostream>

// cut for signal/background
const float cut_sigmaIetaIeta = 0.01;

// cuts for jets
const float cut_jet_pt = 30;
const float cut_jet_eta = 1.6;
const float cut_jet_photon_deltaPhi = (7 * TMath::Pi()) / 8 ;   // 7/8 * TMath::Pi() --> evaluates to 0;

const int MAXPHOTONS = 500;
const int MAXJETS = 500;
const Int_t THRESHOLDS = 1;
enum process {
    AllQCDPhotons,  //forward
    EmEnrichedDijets    // forward
};
const char* processNames[] = {"AllQCDPhotons", "EmEnrichedDijets"};

using namespace std;

Double_t getDPHI( Double_t phi1, Double_t phi2);
bool compareHistograms(TH1* h1, TH1* h2);

void makeTurnOn(process hiForestProcess, TString outFileName, double offlineEtaCut=1.44)
{
    TFile *outFile = new TFile(outFileName,"RECREATE");

    TChain *f1Tree = new TChain("multiPhotonAnalyzer/photon","f1Tree");
    TChain *fEvtTree = new TChain("hiEvtAnalyzer/HiTree","fEvtTree");
    TChain *fSkimTree = new TChain("skimanalysis/HltTree","fSkimTree");
    TChain *hltTree = new TChain("hltanalysis/HltTree","hltTree");
    TChain *f1JetTree = new TChain("akPu3PFJetAnalyzer/t","f1JetTree");   // must decide on the tree accroding to collision type

//    const char* processName = processNames[hiForestProcess];

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

    ///////////// TEST
    TH1D *fPt_iso_test_Cond = (TH1D*)fPt_iso[0]->Clone("fPt_iso_test_Cond");
    TH1D *fPt_iso_test_Loop = (TH1D*)fPt_iso[0]->Clone("fPt_iso_test_Loop");
    TH1D *fPt_iso_test_Cond_subLeading = (TH1D*)fPt_iso[0]->Clone("fPt_iso_test_Cond_subLeading");
    TH1D *fPt_iso_test_Loop_subLeading = (TH1D*)fPt_iso[0]->Clone("fPt_iso_test_Loop_subLeading");

    TH1D *fPtSigma_iso[3];
    fPtSigma_iso[0] = new TH1D("fPtSigma_iso_0","isolated;#sigma_{#eta #eta}", 100, 0, 0.1);
    fPtSigma_iso[1] = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_1");
    fPtSigma_iso[2] = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_2");

    ///////////// TEST
    TH1D *fPtSigma_iso_test_Cond = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_test_Cond");
    TH1D *fPtSigma_iso_test_Loop = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_test_Loop");
    TH1D *fPtSigma_iso_test_Cond_subLeading = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_test_Cond_subLeading");
    TH1D *fPtSigma_iso_test_Loop_subLeading = (TH1D*)fPtSigma_iso[0]->Clone("fPtSigma_iso_test_Loop_subLeading");

//    // tree names that go like accepted* : events that pass selection and that are triggered
//    TH1D *accepted_iso[THRESHOLDS][3];
//    TH1D *accepted_sieie_iso[THRESHOLDS][3];
//
//    TH1D *accepted_sideBand[THRESHOLDS][3];
//    TH1D *accepted_sieie_sideBand[THRESHOLDS][3];

    TString formula_pt = "pt";
    TString formula_sigma = "sigmaIetaIeta";
    TString condition = Form("abs(eta)<%f", offlineEtaCut);
    condition += " && isEle==0";
    condition += " && abs(seedTime)<3";
    condition += " && swissCrx<0.9";
    condition += " && sigmaIetaIeta>0.002";
    condition += " && sigmaIphiIphi>0.002";
    condition += " && hadronicOverEm<0.1";
    condition += " && ecalRecHitSumEtConeDR04<4.2 && hcalTowerSumEtConeDR04<2.2 && trkSumPtHollowConeDR04<2";
//    const char* cut_hiBin60 ="hiBin<60";
//    const char* cut_hiBin100="hiBin>=100";

  //  f1Tree->Draw(Form("MaxIf$(%s,%s) >> %s",formula.Data(),condition.Data(),fPt_vTest[0]->GetName()));
  //  f1Tree->Draw(Form("MaxIf$(%s,nPhotons>0) >> %s",formula_pt.Data(),fPt_vTest[0]->GetName()));
  //  f1Tree->Draw(Form("Max$(%s) >> %s",formula_pt.Data(),fPt_vTest[0]->GetName()));


    // leading photon pT
    // if no photon satisfies condition, then 0 will be plotted.
//    f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula_pt.Data(), condition.Data(), fPt_iso_test_Cond->GetName()));      // tested and works
//    drawMaximum(f1Tree, formula_pt, condition, fPt_iso_test_Cond, true);  // tested and works

    // if no photon satisfies condition, then nothing will be plotted.
//    f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula_pt.Data(), condition.Data(), fPt_iso_test_Cond->GetName()),condition.Data());   // does not work
//    f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula_pt.Data(), condition.Data(), fPt_iso_test_Cond->GetName()), Form("pt == Max$(pt*(%s))", condition.Data()) );   // tested and works
    drawMaximum(f1Tree, formula_pt, condition, fPt_iso_test_Cond);     // tested and works
//    drawMaximumGeneral(f1Tree, formula_pt, formula_pt, condition, fPt_iso_test_Cond);   // tested and works

     // subleading photon pT
    // if no subleading photon found, then 0 will be plotted.
//    f1Tree->Draw(Form("Max$(%s*(pt < Max$(pt*(%s)))*(%s)) >> %s",formula_pt.Data(), condition.Data(), condition.Data(), fPt_iso_test_Cond_subLeading->GetName()));  // tested and works
//    drawMaximum2nd(f1Tree, formula_pt, condition, fPt_iso_test_Cond_subLeading, true);  // tested and works

    // if no subleading photon found, then nothing will be plotted.
    drawMaximum2nd(f1Tree, formula_pt, condition, fPt_iso_test_Cond_subLeading);   // tested and works
//    drawMaximum2nd(f1Tree, formula_pt, condition, "1==1", fPt_iso_test_Cond_subLeading);   // tested and works
//    drawMaximum2ndGeneral(f1Tree, formula_pt, formula_pt, condition, fPt_iso_test_Cond_subLeading);    // tested and works

     // leading photon sigmaIetaIeta with leading photon pT > 40
//    f1Tree->Draw(Form("%s >> %s",formula_sigma.Data(), fPtSigma_iso_test_Cond->GetName()), Form("pt > 40 && pt == Max$(%s*(%s))",formula_pt.Data(), condition.Data()));
     drawMaximumGeneral(f1Tree, formula_sigma, formula_pt, condition, Form("%s > 40", formula_pt.Data()) ,fPtSigma_iso_test_Cond);     // tested and works

     // subleading photon sigmaIetaIeta with leading photon pT > 40 : following two lines give the same output
//     f1Tree->Draw(Form("%s >> %s",formula_sigma.Data(), fPtSigma_iso_test_Cond_subLeading->GetName()), Form("Max$(pt*(%s)) > 40 && pt == Max$(pt*(pt<Max$(pt*(%s)))*(%s))", condition.Data(), condition.Data(), condition.Data())); // works
     drawMaximum2ndGeneral(f1Tree, formula_sigma, formula_pt, condition, Form("Max$(pt*(%s)) > 40", condition.Data()) ,fPtSigma_iso_test_Cond_subLeading);   // tested and works


//    f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula_pt.Data(), condition.Data(), fPt_vTest[1]->GetName()), cut_hiBin60);
//    f1Tree->Draw(Form("Max$(%s*(%s)) >> %s",formula_pt.Data(), condition.Data(), fPt_vTest[2]->GetName()), cut_hiBin100);
    std::cout<<"condition = "<<condition.Data()<<std::endl;

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

        hltTree->GetEntry(j);
        f1Tree->GetEntry(j);
        f1JetTree->GetEntry(j);

        //    double maxfpt = 0;
        //    //double maxfeta = -10;
        //    //double maxfphi = -10;
        //    double sigmaietaieta = -1;

        double maxfpt_iso = 0;
        double maxfpt_iso_subLeading = 0;        //TEST
        int index_maxfpt_iso = -1;
        int index_maxfpt_iso_subLeading = -1;

        double sigmaietaieta_iso = -1;
        double sigmaietaieta_iso_subLeading = -1;        //TEST

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
                        maxfpt_iso_subLeading = maxfpt_iso;   //TEST
                        sigmaietaieta_iso_subLeading = sigmaietaieta_iso;   //TEST

                        maxfpt_iso=photon_pt[i];
                        sigmaietaieta_iso = sigmaIetaIeta[i];
                        index_maxfpt_iso=i;
                    }
                }
                // TEST
                else if(photon_pt[i] > maxfpt_iso_subLeading)
                {
                    if(ecalRecHitSumEtConeDR04[i] < 4.2 && hcalTowerSumEtConeDR04[i] < 2.2 && trkSumPtHollowConeDR04[i] < 2)
                    {
                        maxfpt_iso_subLeading = photon_pt[i];
                        sigmaietaieta_iso_subLeading = sigmaIetaIeta[i];
                        index_maxfpt_iso_subLeading = i;
                    }
                }
                // TEST - END
            }
        }

        if(!(index_maxfpt_iso<0))
        {
            fPt_iso[0]->Fill(maxfpt_iso,weight);
            fPt_iso_test_Loop->Fill(maxfpt_iso);     //TEST
            if(!(index_maxfpt_iso_subLeading<0)) fPt_iso_test_Loop_subLeading->Fill(maxfpt_iso_subLeading);     //TEST
//            fPt_iso_test_Loop_subLeading->Fill(maxfpt_iso_subLeading);     //TEST
            if(hiBin < 60){
                fPt_iso[1]->Fill(maxfpt_iso,weight);
            }
            else if (hiBin >= 100) {
                fPt_iso[2]->Fill(maxfpt_iso,weight);
            }
        }

        if(maxfpt_iso > 40)
        {
            if(sigmaietaieta_iso<=0)
                cout << sigmaietaieta_iso << endl;
            fPtSigma_iso[0]->Fill(sigmaietaieta_iso,weight);
            fPtSigma_iso_test_Loop->Fill(sigmaietaieta_iso);    //TEST
            if(sigmaietaieta_iso_subLeading>-1) fPtSigma_iso_test_Loop_subLeading->Fill(sigmaietaieta_iso_subLeading);    //TEST
//            fPtSigma_iso_test_Loop_subLeading->Fill(sigmaietaieta_iso_subLeading);    //TEST
            if(hiBin < 60){
                fPtSigma_iso[1]->Fill(sigmaietaieta_iso,weight);
            }
            else if (hiBin >= 100) {
                fPtSigma_iso[2]->Fill(sigmaietaieta_iso,weight);
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

    std::cout<<compareHistograms(fPt_iso_test_Cond, fPt_iso_test_Loop)<<std::endl;
    std::cout<<compareHistograms(fPtSigma_iso_test_Cond, fPtSigma_iso_test_Loop)<<std::endl;
    std::cout<<compareHistograms(fPt_iso_test_Cond_subLeading, fPt_iso_test_Loop_subLeading)<<std::endl;
    std::cout<<compareHistograms(fPtSigma_iso_test_Cond_subLeading, fPtSigma_iso_test_Loop_subLeading)<<std::endl;

    outFile->cd();
    fPt_iso_test_Cond->Write();
    fPt_iso_test_Loop->Write();
    fPtSigma_iso_test_Cond->Write();
    fPtSigma_iso_test_Loop->Write();
    //
    fPt_iso_test_Cond_subLeading->Write();
    fPt_iso_test_Loop_subLeading->Write();
    fPtSigma_iso_test_Cond_subLeading->Write();
    fPtSigma_iso_test_Loop_subLeading->Write();

    // add full photon tree chain
    f1Tree->SetNameTitle("photon","full Photon TChain");
    f1Tree->Write();

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
