void plotCorrelationL1(){


  TFile *filPP=new TFile("testingoutput.root");
 
  const int NHISTS = 22;
  TH2D *aveptregion_aveptgen[NHISTS];
  TH2D *aveptregion_hiNpix[NHISTS];

  for(int i = 0; i < NHISTS; i++){
    aveptregion_aveptgen[i] = (TH2D*)filPP->Get(Form("aveptregion_aveptgen_%d",i));
    aveptregion_hiNpix[i] = (TH2D*)filPP->Get(Form("aveptregion_hiNpix_%d",i));  
	aveptregion_aveptgen[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));
	aveptregion_hiNpix[i]->GetYaxis()->SetTitle(Form("<pt_{region} (#eta=%d)>",i));
	aveptregion_hiNpix[i]->GetXaxis()->SetLabelSize(0.03);

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

  
}