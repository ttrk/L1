void plotCorrelationL1(){


  TFile *filPP=new TFile("testingoutput.root");
 
  const int NHISTS = 22;
  TH2D *aveptregion_aveptgen[NHISTS];
  TH2D *aveptregion_hiNpix[NHISTS];
  TH2D *aveptregion_multgen_lowpt[NHISTS];
  TH2D *aveptregion_multgen_middlept[NHISTS];
  TH2D *aveptregion_multgen_highpt[NHISTS];


  for(int i = 0; i < NHISTS; i++){
    aveptregion_aveptgen[i] = (TH2D*)filPP->Get(Form("aveptregion_aveptgen_%d",i));
    aveptregion_hiNpix[i] = (TH2D*)filPP->Get(Form("aveptregion_hiNpix_%d",i));  
	aveptregion_aveptgen[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));
	aveptregion_hiNpix[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));
	aveptregion_hiNpix[i]->GetXaxis()->SetLabelSize(0.03);
	
    aveptregion_multgen_lowpt[i] = (TH2D*)filPP->Get(Form("aveptregion_multgen_lowpt_%d",i));
    aveptregion_multgen_middlept[i] = (TH2D*)filPP->Get(Form("aveptregion_multgen_middlept_%d",i));
    aveptregion_multgen_highpt[i] = (TH2D*)filPP->Get(Form("aveptregion_multgen_highpt_%d",i));
    
	aveptregion_multgen_lowpt[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));
	aveptregion_multgen_middlept[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));
	aveptregion_multgen_highpt[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));

	aveptregion_multgen_middlept[i]->GetXaxis()->SetRangeUser(0,100);
	aveptregion_multgen_highpt[i]->GetXaxis()->SetRangeUser(0,100);


  }

  TCanvas*canvasptgen=new TCanvas("canvasptgen","canvasptgen",1200,800);
  canvasptgen->Divide(3,2);
  canvasptgen->cd(1);
  aveptregion_aveptgen[0]->Draw("colz");
  canvasptgen->cd(2);
  aveptregion_aveptgen[2]->Draw("colz");
  canvasptgen->cd(3);
  aveptregion_aveptgen[4]->Draw("colz");
  canvasptgen->cd(4);
  aveptregion_aveptgen[8]->Draw("colz");
  canvasptgen->cd(5);
  aveptregion_aveptgen[15]->Draw("colz");
  canvasptgen->cd(6);
  aveptregion_aveptgen[20]->Draw("colz");
  canvasptgen->SaveAs("canvasptgen.pdf");

  TCanvas*canvashiNpix=new TCanvas("canvashiNpix","canvashiNpix",1200,800);
  canvashiNpix->Divide(3,2);
  canvashiNpix->cd(1);
  aveptregion_hiNpix[0]->Draw("colz");
  canvashiNpix->cd(2);
  aveptregion_hiNpix[2]->Draw("colz");
  canvashiNpix->cd(3);
  aveptregion_hiNpix[4]->Draw("colz");
  canvashiNpix->cd(4);
  aveptregion_hiNpix[8]->Draw("colz");
  canvashiNpix->cd(5);
  aveptregion_hiNpix[15]->Draw("colz");
  canvashiNpix->cd(6);
  aveptregion_hiNpix[20]->Draw("colz");
  canvashiNpix->SaveAs("canvashiNpix.pdf");

  TCanvas*canvashimultgen_lowpt=new TCanvas("canvashimultgen_lowpt","canvashimultgen_lowpt",1200,800);
  canvashimultgen_lowpt->Divide(3,2);
  canvashimultgen_lowpt->cd(1);
  aveptregion_multgen_lowpt[0]->Draw("colz");
  canvashimultgen_lowpt->cd(2);
  aveptregion_multgen_lowpt[2]->Draw("colz");
  canvashimultgen_lowpt->cd(3);
  aveptregion_multgen_lowpt[4]->Draw("colz");
  canvashimultgen_lowpt->cd(4);
  aveptregion_multgen_lowpt[8]->Draw("colz");
  canvashimultgen_lowpt->cd(5);
  aveptregion_multgen_lowpt[15]->Draw("colz");
  canvashimultgen_lowpt->cd(6);
  aveptregion_multgen_lowpt[20]->Draw("colz");
  canvashimultgen_lowpt->SaveAs("canvashimultgen_lowpt.pdf");

  TCanvas*canvashimultgen_middlept=new TCanvas("canvashimultgen_middlept","canvashimultgen_middlept",1200,800);
  canvashimultgen_middlept->Divide(3,2);
  canvashimultgen_middlept->cd(1);
  aveptregion_multgen_middlept[0]->Draw("colz");
  canvashimultgen_middlept->cd(2);
  aveptregion_multgen_middlept[2]->Draw("colz");
  canvashimultgen_middlept->cd(3);
  aveptregion_multgen_middlept[4]->Draw("colz");
  canvashimultgen_middlept->cd(4);
  aveptregion_multgen_middlept[8]->Draw("colz");
  canvashimultgen_middlept->cd(5);
  aveptregion_multgen_middlept[15]->Draw("colz");
  canvashimultgen_middlept->cd(6);
  aveptregion_multgen_middlept[20]->Draw("colz");
  canvashimultgen_middlept->SaveAs("canvashimultgen_middlept.pdf");



  TCanvas*canvashimultgen_highpt=new TCanvas("canvashimultgen_highpt","canvashimultgen_highpt",1200,800);
  canvashimultgen_highpt->Divide(3,2);
  canvashimultgen_highpt->cd(1);
  aveptregion_multgen_highpt[0]->Draw("colz");
  canvashimultgen_highpt->cd(2);
  aveptregion_multgen_highpt[2]->Draw("colz");
  canvashimultgen_highpt->cd(3);
  aveptregion_multgen_highpt[4]->Draw("colz");
  canvashimultgen_highpt->cd(4);
  aveptregion_multgen_highpt[8]->Draw("colz");
  canvashimultgen_highpt->cd(5);
  aveptregion_multgen_highpt[15]->Draw("colz");
  canvashimultgen_highpt->cd(6);
  aveptregion_multgen_highpt[20]->Draw("colz");
  canvashimultgen_highpt->SaveAs("canvashimultgen_highpt.pdf");

  
}