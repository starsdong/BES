#include "style.C+"
#include "draw.C+"

Double_t xgx(Double_t *x, Double_t *par)
{
  double A1 = par[0];
  double A2 = par[1];
  double mu = par[2];
  double s1 = par[3];
  double s2 = par[4];
  //  double s3 = par[5];
  double xl = TMath::Log10(x[0]) - TMath::Log10(mu);
  if(x[0]<mu) {
    return -A1*(xl+A2*xl*xl)*TMath::Exp(-0.5*TMath::Power(xl/s1,2.0)); // - A2*xl*TMath::Exp(-0.25*TMath::Power(xl/s2,4.0));
  } else {
    return -A1*xl*TMath::Exp(-0.5*TMath::Power(xl/s2,2.0));
  }
}

void kappasigma2(){
  style();

  // parametrization from Andronic: 1710.09425
  TF1 *fTCF = new TF1("TCF","158.4/(1+exp(2.60-log(x))/0.45)",2.,1e5);
  TF1 *fMuB = new TF1("MuB","1307.5/(1+0.288*x)",2.,1e5);
  
  const Int_t NE = 10;
  const Double_t E[NE] = {3.0, 7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 54.4, 62.4, 200.};  // 3.0 - -0.5<y<0,  others - |y|<0.5
  Double_t MuB[NE];
  Double_t C4C2[NE] = {-0.846, 1.767, 0.696, 1.469, 0.141, 0.196, 0.740, 0.633, 0.793, 0.901};
  Double_t C4C2_e[NE] = {0.086, 1.151, 0.595, 0.396, 0.329, 0.219, 0.147, 0.055, 0.251, 0.209};
  Double_t C4C2_se[NE] = {0.822, 0.414, .267, 0.254, 0.156, 0.143, 0.136, 0.139, 0.122, 0.139};
  Double_t C4C2_e_tot[NE];

  const Int_t NE_ref = 5;
  const Double_t E_ref[NE_ref] = {3.0, 8.8, 17.27, 27.0, 62.4};  // 3.0 - UrQMD, others - HRG - CE
  const Double_t C4C2_ref[NE_ref] = {-0.774, 0.43147, 0.65726, 0.73934, 0.81176};
  TGraph *gr_ref = new TGraph(NE_ref, E_ref, C4C2_ref);

  for(int i=0;i<NE;i++) {
    C4C2[i] -= gr_ref->Eval(E[i]);
    C4C2_e_tot[i] = sqrt(C4C2_e[i]*C4C2_e[i]+C4C2_se[i]*C4C2_se[i]);
  }
  TGraphErrors *gr_e = new TGraphErrors(NE, E, C4C2, 0, C4C2_e);
  TGraphErrors *gr_se = new TGraphErrors(NE, E, C4C2, 0, C4C2_se);
  TGraphErrors *gr_e_tot = new TGraphErrors(NE, E, C4C2, 0, C4C2_e_tot);
  
  
  const Int_t NEP = 4;
  const Double_t EP[NEP] = {7.7, 14.5, 19.6, 27};
  Double_t MuBP[NEP];
  const Double_t y_0[NEP] = {0.0, 0.0, 0.0, 0.0};
  const Double_t y_e[NEP] = {0.3, 0.13, 0.061, 0.065};
  TGraphErrors *gr_p = new TGraphErrors(NEP, EP, y_0, 0, y_e);
  
  const Int_t NEL = 4;
  const Double_t EL[NEL] = {3.0, 7.7, 19.6, 200.};
  const Char_t* TextL[NEL] = {"3.0", "7.7", "19.6", "200"};
  Double_t MuBL[NEL];
  for(int i=0;i<NEL;i++) {
      MuBL[i] = fMuB->Eval(EL[i]);
  }

  for(int i=0;i<NEP;i++) {
    MuBP[i] = fMuB->Eval(EP[i]);
  }
  for(int i=0;i<NE;i++) {
    MuB[i] = fMuB->Eval(E[i]);
  }
  TGraphErrors *gr_e_muB = new TGraphErrors(NE, MuB, C4C2, 0, C4C2_e);
  TGraphErrors *gr_se_muB = new TGraphErrors(NE, MuB, C4C2, 0, C4C2_se);
  TGraphErrors *gr_p_muB = new TGraphErrors(NEP, MuBP, y_0, 0, y_e);
  
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->SetLeftMargin(0.14);
  c1->cd()->SetLogx();

  double x1 = 1.9;
  double x2 = 4e2;
  double y1 = -1.3;
  double y2 = 2.9;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->SetYTitle("Net-proton #kappa#sigma^{2}_{mea} - #kappa#sigma^{2}_{ref}");
  h0->GetYaxis()->SetTitleSize(0.07);
  h0->GetYaxis()->SetTitleOffset(0.8);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw("c");

  //  drawColorBox(2.0, y1, 5.0, y2, 5, 0.5);
  
  gr_p->SetFillColor(kGreen-3);
  gr_p->SetLineColor(kGreen-3);
  gr_p->SetLineWidth(15);
  //  gr_p->Draw("e3");
  drawLine(x1,0,x2,0,2,9,1);

  TF1 *fun = new TF1("fun",xgx,x1,x2,5);
  fun->SetParameters(2.915,-18.51,13.,0.177,0.243);
  fun->FixParameter(2,13.0);
  fun->SetLineWidth(2);
  fun->SetLineStyle(2);
  fun->SetLineColor(4);
  //  gr_e_tot->Fit("fun","R");
  fun->Draw("same");
  
	   
  drawSysBox(gr_se, 0.05, 16, 1);
  setGraphMarker(gr_e, 20, 1, 2.2);
  setGraphLine(gr_e, 1, 1, 3);
  gr_e->Draw("p");
  TGraphErrors *gr_e_2 = (TGraphErrors *)gr_e->Clone();
  setGraphMarker(gr_e_2, 20, 2, 1.5);
  gr_e_2->Draw("p");
  
  
  drawHistBox(x1,x2,y1,y2);
  drawText(2.85, y1-0.25, "3");
  drawText(8.5, y1-0.25, "10");
  drawText(26.5, y1-0.25, "30");
  drawText(85, y1-0.25, "100");
  drawText(2.6, -0.8, "-0.5<y<0",42,0.035,90);

  drawText(23, y2*0.85, "Central Au+Au Collisions",22,0.055);
  drawText(42, y2*0.72, "0.4<p_{T}<2.0 GeV/c, |y|<0.5",42,0.04);
  drawText(60, y2*0.6, " ref: HRG CE/UrQMD",12,0.04);
  TLegend *leg = new TLegend(0.67, 0.66, 0.9, 0.74);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  leg->SetTextFont(12);
  leg->AddEntry(gr_p,"  BES-II proj.","l");
  //  leg->Draw();
  
  c1->Update();
  c1->SaveAs("kappasigma2VsE.pdf");
  c1->SaveAs("kappasigma2VsE.png");


  TCanvas *c2 = new TCanvas("c2","c2",800,0,800,600);
  c2->SetLeftMargin(0.14);
  c2->SetTickx(0);
  c2->cd();

  x1 = -50;
  x2 = 850;
  h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Baryon Chemical Potential #mu_{B} (MeV)");
  h0->GetXaxis()->SetLabelOffset(0.008);
  h0->GetXaxis()->SetTitleOffset(1.1);
  h0->SetYTitle("Net-proton #kappa#sigma^{2}_{mea} - #kappa#sigma^{2}_{ref}");
  h0->GetYaxis()->SetTitleSize(0.07);
  h0->GetYaxis()->SetTitleOffset(0.8);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw("c");

  //  drawColorBox(fMuB->Eval(19.6), y1, fMuB->Eval(3.0), y2, kGreen-4, 0.3);
  drawColorBox(fMuB->Eval(5.0), y1, fMuB->Eval(2.0), y2, 5, 0.5);

  gr_p_muB->SetFillColor(kGreen-3);
  gr_p_muB->SetLineColor(kGreen-3);
  gr_p_muB->SetLineWidth(15);
  gr_p_muB->Draw("e3");
  drawLine(x1,0,x2,0,2,9,1);
	   
  drawSysBox(gr_se_muB, 10, 16);
  setGraphMarker(gr_e_muB, 20, 1, 2.2);
  setGraphLine(gr_e_muB, 1, 1, 3);
  gr_e_muB->Draw("p");
  TGraphErrors *gr_e_muB_2 = (TGraphErrors *)gr_e_muB->Clone();
  setGraphMarker(gr_e_muB_2, 20, 2, 1.5);
  gr_e_muB_2->Draw("p");
  
  
  drawHistBox(x1,x2,y1,y2);
  drawText(740, -0.8, "-0.5<y<0",42,0.035,90);

  drawText(0, y2*0.85, "Central Au+Au Collisions",22,0.055);
  drawText(0, y2*0.72, "0.4<p_{T}<2.0 GeV/c, |y|<0.5",42,0.04);
  drawText(0, y2*0.6, " ref: HRG CE/UrQMD",12,0.04);
  leg = new TLegend(0.18, 0.66, 0.4, 0.74);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  leg->SetTextFont(12);
  leg->AddEntry(gr_p_muB,"  BES-II proj.","l");
  leg->Draw();

  c2->Update();
  c2->SaveAs("kappasigma2VsMuB.pdf");
  c2->SaveAs("kappasigma2VsMuB.png");
}


