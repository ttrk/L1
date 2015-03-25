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
  Double_t rctEtaMap[23] = {-5.100, -4.500, -4.000, -3.500, -3.000, -2.172, -1.740, -1.392, -1.044, -0.696,-0.348,  0.000,  0.348,  0.696,  1.044,  1.392,  1.740,  2.172,  3.000,  3.500, 4.000,  4.500,  5.10}; 


void kinematicCorrelator(TString inL1FileName, TString inHiForestFileName, TString outFileName)
{
	TFile *lFile = TFile::Open(inL1FileName);
	TTree *l1Tree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

	Int_t l1_event, l1_run, l1_lumi;
	Int_t l1_hwPt[MAXL1JETS], l1_hwEta[MAXL1JETS], l1_hwPhi[MAXL1JETS];
	Double_t l1_pt[MAXL1JETS];

	Int_t region_hwPt[396], region_hwEta[396], region_hwPhi[396];

	l1Tree->SetBranchAddress("event",&l1_event);
	l1Tree->SetBranchAddress("run",&l1_run);
	l1Tree->SetBranchAddress("lumi",&l1_lumi);
	l1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
	l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
	l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
	l1Tree->SetBranchAddress("jet_pt",l1_pt);
	l1Tree->SetBranchAddress("region_hwPt", region_hwPt);
	l1Tree->SetBranchAddress("region_hwEta", region_hwEta);
	l1Tree->SetBranchAddress("region_hwPhi", region_hwPhi);

	TFile *fFile = TFile::Open(inHiForestFileName);
	TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
	TTree *fGenParticle = (TTree*)fFile->Get("HiGenParticleAna/hi");
	fEvtTree->AddFriend(fGenParticle);

	Int_t f_evt, f_run, f_lumi;
	Float_t vz;
	Int_t hiBin;
	Int_t hiNpix;
	Float_t pt[200000];
	Float_t eta[200000];
	Int_t n;
	Float_t npart;


	fEvtTree->SetBranchAddress("evt",&f_evt);
	fEvtTree->SetBranchAddress("run",&f_run);
	fEvtTree->SetBranchAddress("lumi",&f_lumi);
	fEvtTree->SetBranchAddress("vz",&vz);
	fEvtTree->SetBranchAddress("hiBin",&hiBin);
	fEvtTree->SetBranchAddress("hiNpix",&hiNpix);
	fEvtTree->SetBranchAddress("pt",pt);
	fEvtTree->SetBranchAddress("eta",eta);
	fEvtTree->SetBranchAddress("n",&n);
	fEvtTree->SetBranchAddress("npart",&npart);


	TFile *outFile = new TFile(outFileName,"RECREATE");

	const int NHISTS = 22;
	TH2D *aveptregion_aveptgen[NHISTS];
	TH2D *aveptregion_hiNpix[NHISTS];
	
	TH2D *aveptregion_multgen_lowpt[NHISTS];
	TH2D *aveptregion_multgen_middlept[NHISTS];
	TH2D *aveptregion_multgen_highpt[NHISTS];

	aveptregion_aveptgen[0]=new TH2D("aveptregion_aveptgen_0",";<pt_{gen}>;<pt_{region}(#eta=0)>",100,0,2,50,0,100);
	aveptregion_hiNpix[0]=new TH2D("aveptregion_hiNpix_0",";N_{pixels};<pt_{region}(#eta=0)>",500,0,10000,50,0,100);
	aveptregion_multgen_lowpt[0]=new TH2D("aveptregion_multgen_lowpt_0",";N_{part}(p_{t}<0.6 GeV);<pt_{region}(#eta=0)>",500,0,500,50,0,200);
	aveptregion_multgen_middlept[0]=new TH2D("aveptregion_multgen_middlept_0",";N_{part}(0.6<p_{t}<1.0 GeV);<pt_{region}(#eta=0)>",500,0,500,50,0,200);
	aveptregion_multgen_highpt[0]=new TH2D("aveptregion_multgen_highpt_0",";N_{part}(p_{t}>1.0 GeV);<pt_{region}(#eta=0)>",500,0,500,50,0,200);

	for(int i = 0; i < NHISTS; i++)
	{   
		if(i ==0 ) continue;
		aveptregion_aveptgen[i]=(TH2D*)aveptregion_aveptgen[0]->Clone(Form("aveptregion_aveptgen_%i",i));
		aveptregion_hiNpix[i]=(TH2D*)aveptregion_hiNpix[0]->Clone(Form("aveptregion_hiNpix_%i",i));
		aveptregion_multgen_lowpt[i]=(TH2D*)aveptregion_multgen_lowpt[0]->Clone(Form("aveptregion_multgen_lowpt_%i",i));
		aveptregion_multgen_middlept[i]=(TH2D*)aveptregion_multgen_middlept[0]->Clone(Form("aveptregion_multgen_middlept_%i",i));
		aveptregion_multgen_highpt[i]=(TH2D*)aveptregion_multgen_highpt[0]->Clone(Form("aveptregion_multgen_highpt_%i",i));

	}

	// Make the event-matching map ************
	EventMatchingCMS *matcher = new EventMatchingCMS();
	std::cout << "Begin making map." << std::endl;
	Long64_t l_entries = l1Tree->GetEntries();
	for(Long64_t j = 0; j < l_entries; ++j)
	{
		l1Tree->GetEntry(j);
		matcher->addEvent(l1_event, l1_lumi, l1_run, j);
	}
	std::cout << "Finished making map." << std::endl;
	// **********************

	outFile->cd();

	int count = 0;

	Long64_t entries = fEvtTree->GetEntries();
	std::cout << "L1 entries: " << l_entries << std::endl;
	std::cout << "HiForest entries: " << entries << std::endl;
	for(Long64_t j = 0; j < entries; ++j)
	{
		if(j % 10000 == 0)
			printf("%lld / %lld\n",j,entries);

		fEvtTree->GetEntry(j);
		// retrieve the matching entry from the l1 input ***********
		long long l1Entry = matcher->retrieveEvent(f_evt, f_lumi, f_run);
		//long long l1Entry = matcher->retrieveEvent(f_evt, 0, f_run);
		if( l1Entry == -1 )  continue;
		//**************

		l1Tree->GetEntry(l1Entry);
		count++;

		double sums[NHISTS];
		double sumsptgen=0;
		int multgenlowpt[NHISTS];
		int multgenmiddlept[NHISTS];
		int multgenhighpt[NHISTS];

		for(int i = 0; i < NHISTS; i++)
		{
			sums[i] = 0.;
			multgenlowpt[i] = 0.;
			multgenmiddlept[i] = 0.;
			multgenhighpt[i] = 0.;
		}
		for(int i = 0; i < 396; i++)
		{
			sums[region_hwEta[i]] += region_hwPt[i];
		}

		for(int i = 0; i < n; i++){
			sumsptgen+=pt[i];
			for(int m = 0; m < 22; m++){
			  if(eta[i]>rctEtaMap[m] && eta[i]<rctEtaMap[m+1]){
			    if(pt[i]<0.6) multgenlowpt[m]++;
			    else if(pt[i]>0.6 && pt[i]<1.0) multgenmiddlept[m]++;
			    else if (pt[i]>1.0) multgenhighpt[m]++;
			  }
			}
		}

		for(int i = 0; i < NHISTS; i++)
		{
			aveptregion_aveptgen[i]->Fill(sumsptgen/n,sums[i]/18.);
			aveptregion_hiNpix[i]->Fill(hiNpix,sums[i]/18.);
			aveptregion_multgen_lowpt[i]->Fill(multgenlowpt[i],sums[i]/18.);
			aveptregion_multgen_middlept[i]->Fill(multgenmiddlept[i],sums[i]/18.);
			aveptregion_multgen_highpt[i]->Fill(multgenhighpt[i],sums[i]/18.);
		}
	}

	outFile->cd();

	for(int i = 0; i < NHISTS; i++)
	{
		aveptregion_aveptgen[i]->Write();
		aveptregion_hiNpix[i]->Write();
		aveptregion_multgen_lowpt[i]->Write();
		aveptregion_multgen_middlept[i]->Write();
		aveptregion_multgen_highpt[i]->Write();
	}

	std::cout << "Matching entries: " << count << std::endl;

	lFile->Close();
	outFile->Close();
}

int main(int argc, char **argv)
{
	if(argc == 4)
	{
		kinematicCorrelator(argv[1], argv[2], argv[3]);
		return 0;
	}
	else
	{
		std::cout << "Usage: \nkinematicPlotterCentrality.exe <input_l1_file> <input_HiForest_file> <output_file>" << std::endl;
		return 1;
	}
}
