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

void kinematicPlotterCentrality(TString inL1FileName, TString inHiForestFileName, TString outFileName)
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
	const int CENTBINS = 4;
	const int CENTCUT[CENTBINS] = {0, 10000, 20000, 30000};

	TH1D *ebeSigma[CENTBINS][NHISTS];
	TH1D *range[CENTBINS][NHISTS];
	TH1D *phiDelta[CENTBINS][NHISTS];
	TH1D *ptDists[CENTBINS][NHISTS];

	range[0][0] = new TH1D("range_0_0","range;max pT - min pT",1024,0,1024);
	phiDelta[0][0] = new TH1D("phiDelta_0_0","phiDelta;PhiMax - PhiMin",10,0,10);
	ebeSigma[0][0] = new TH1D("ebeSigma_0_0","ebeSigma;#sigma;count",300,0,50);
	ptDists[0][0] = new TH1D("ptDists_0_0",";region_hwPt",1024,0,1024);

	TH1D *distSigma[CENTBINS];
	TH1D *distPt[CENTBINS];
	TH2D *aveptregion_aveptgen[NHISTS];
	TH2D *aveptregion_hiNpix[NHISTS];

	distSigma[0]= new TH1D("distSigma_0",";#phi index;#sigma", NHISTS,0,NHISTS);
	distPt[0]= new TH1D("distPt_0",";#phi index;<p_{T}>",NHISTS,0,NHISTS);
	aveptregion_aveptgen[0]=new TH2D("aveptregion_aveptgen_0",";<pt_{gen}>;<pt_{region}(#eta=0)>",100,0,2,50,0,100);
	aveptregion_hiNpix[0]=new TH2D("aveptregion_hiNpix_0",";N_{pixels};<pt_{region}(#eta=0)>",500,0,10000,50,0,100);


	for(int cent = 0; cent < CENTBINS; cent++)
	{
		for(int i = 0; i < NHISTS; i++)
		{
			if(i == 0 && cent ==0) continue;
			ebeSigma[cent][i] = (TH1D*)ebeSigma[0][0]->Clone(Form("ebeSigma_%i_%i",cent,i));
			ptDists[cent][i] = (TH1D*)ptDists[0][0]->Clone(Form("ptDists_%i_%i",cent,i));
			range[cent][i] = (TH1D*)range[0][0]->Clone(Form("range_%i_%i",cent,i));
			phiDelta[cent][i] = (TH1D*)phiDelta[0][0]->Clone(Form("phiDelta_%i_%i",cent,i));
		}
		if(cent ==0 ) continue;
		distSigma[cent] = (TH1D*)distSigma[0]->Clone(Form("distSigma_%i", cent));
		distPt[cent] = (TH1D*)distPt[0]->Clone(Form("distPt_%i", cent));
	}

	for(int i = 0; i < NHISTS; i++)
	{   
		if(i ==0 ) continue;
		aveptregion_aveptgen[i]=(TH2D*)aveptregion_aveptgen[0]->Clone(Form("aveptregion_aveptgen_%i",i));
		aveptregion_hiNpix[i]=(TH2D*)aveptregion_hiNpix[0]->Clone(Form("aveptregion_hiNpix_%i",i));

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

		int centBin = 0;
		for(int i = 0; i < CENTBINS; i++)
		{
			if(hiNpix > CENTCUT[i])
				centBin = i;
		}

		double sums[NHISTS];
		double sums2[NHISTS];
		double maxLocation[NHISTS];
		double maxValue[NHISTS];
		double minLocation[NHISTS];
		double minValue[NHISTS];
		double sumsptgen=0;

		for(int i = 0; i < NHISTS; i++)
		{
			sums[i] = 0.;
			sums2[i] = 0.;
			maxLocation[i] = -1;
			maxValue[i] = -1;
			minLocation[i] = -1;
			minValue[i] = -1;
		}
		for(int i = 0; i < 396; i++)
		{
			ptDists[centBin][region_hwEta[i]]->Fill(region_hwPt[i]);
			sums[region_hwEta[i]] += region_hwPt[i];
			sums2[region_hwEta[i]] += (region_hwPt[i] * region_hwPt[i]);

			if(maxValue[region_hwEta[i]] < region_hwPt[i])
			{
				maxValue[region_hwEta[i]] = region_hwPt[i];
				maxLocation[region_hwEta[i]] = region_hwPhi[i];
			}
			if(minValue[region_hwEta[i]] > region_hwPt[i])
			{
				minValue[region_hwEta[i]] = region_hwPt[i];
				minLocation[region_hwEta[i]] = region_hwPhi[i];
			}
		}

		for(int i = 0; i < n; i++){
			sumsptgen+=pt[i];
		}

		for(int i = 0; i < NHISTS; i++)
		{
			double ebesigma = TMath::Sqrt( (sums2[i]/18.) - ((sums[i]/18.0)*(sums[i]/18.0)) );
			ebeSigma[centBin][i]->Fill(ebesigma);
			range[centBin][i]->Fill(maxValue[i] - minValue[i]);
			double phidelta = TMath::Abs(maxLocation[i] - minLocation[i]);
			if(phidelta > 9)
				phidelta = 18 - phidelta;
			phiDelta[centBin][i]->Fill(phidelta);

			aveptregion_aveptgen[i]->Fill(sumsptgen/n,sums[i]/18.);
			aveptregion_hiNpix[i]->Fill(hiNpix,sums[i]/18.);
		}
	}

	for(int cent = 0; cent < CENTBINS; cent++)
	{
		for(int i = 0; i < NHISTS; i++)
		{
			distSigma[cent]->SetBinContent(i+1, ptDists[cent][i]->GetStdDev());
			distSigma[cent]->SetBinError(i+1, ptDists[cent][i]->GetStdDevError());

			distPt[cent]->SetBinContent(i+1, ptDists[cent][i]->GetMean());
			distPt[cent]->SetBinError(i+1, ptDists[cent][i]->GetMeanError());
		}
	}

	outFile->cd();
	for(int cent = 0; cent < CENTBINS; cent++)
	{
		for(int i = 0; i < NHISTS; i++)
		{
			ptDists[cent][i]->Write();
			ebeSigma[cent][i]->Write();
			range[cent][i]->Write();
			phiDelta[cent][i]->Write();
		}
		distSigma[cent]->Write();
		distPt[cent]->Write();
	}

	for(int i = 0; i < NHISTS; i++)
	{
		aveptregion_aveptgen[i]->Write();
		aveptregion_hiNpix[i]->Write();
	}

	std::cout << "Matching entries: " << count << std::endl;

	lFile->Close();
	outFile->Close();
}

int main(int argc, char **argv)
{
	if(argc == 4)
	{
		kinematicPlotterCentrality(argv[1], argv[2], argv[3]);
		return 0;
	}
	else
	{
		std::cout << "Usage: \nkinematicPlotterCentrality.exe <input_l1_file> <input_HiForest_file> <output_file>" << std::endl;
		return 1;
	}
}
