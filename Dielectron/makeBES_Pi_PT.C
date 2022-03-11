#include "draw.C+"
#include "style.C+"

void makeBES_Pi_PT()
{
  style();

  // parametrization from Andronic: 1710.09425
  TF1 *fTCF = new TF1("TCF","158.4/(1+exp(2.60-log(x))/0.45)",2.,1e5);
  TF1 *fMuB = new TF1("MuB","1307.5/(1+0.288*x)",2.,1e5);

  const Int_t NE = 7;
  const Double_t E[NE] = {2.42, 17.3, 19.6, 27.0, 39.0, 62.4, 200};
  // const Double_t Y[NE] = {7.44, 7.00, 5.65, 9.51, 7.85, 8.29};
  // const Double_t Ye[NE] = {0.26, 2.83, 2.07, 1.28, 1.78, 1.43};
  // const Double_t Yes[NE] = {0.00, 1.19, 0.72, 1.02, 1.22, 1.67};
  // const Double_t Yec[NE] = {0.00, 0.00, 0.64, 0.75, 0.85, 0.00};
  // Excess yield / pi+- 0.3-0.7
  const Double_t Y[NE] = {9.12, 12.95, 12.49, 8.85, 14.92, 12.22, 14.22};
  const Double_t Ye[NE] = {0.48, 1.00, 2.80, 1.97, 1.17, 1.54, 1.40};
  const Double_t Yes[NE] = {1.74, 0.00, 1.94, 1.54, 1.74, 1.30, 1.67};
  Double_t muB[NE];
  Double_t Yet[NE], Yest[NE];
  for(int i=0;i<NE;i++) {
    muB[i] = fMuB->Eval(E[i]);
    Yet[i] = sqrt(Ye[i]*Ye[i]+Yes[i]*Yes[i]);
  }

  TGraph *gr_S = new TGraph("ExcessPi_RR.txt","%lg %lg");
  const Int_t NM = 2;  // Number of models
  const Int_t NPMAX = 20;
  const Int_t NPM[NM] = {11, 11};
  const Char_t *NameM[NM] = {"Rapp-noPT","Rapp-PT"};
  const Double_t EM[NM][NPMAX] = {{4.75, 6.32, 7.69, 9.18, 11.51, 14.5, 19.57, 27.0, 39.0, 62.4, 200},
				  {4.75, 6.32, 7.69, 9.18, 11.51, 14.5, 19.57, 27.0, 39.0, 62.4, 200}};
  const Double_t YMo[NM][NPMAX] = {{11.29, 11.35, 11.88, 12.29, 12.76, 13.29, 13.88, 1., 1., 1., 1.},
				   {16.35, 13.41, 13.41, 13.41, 13.52, 13.70, 13.94, 1., 1., 1., 1.}}; // 0.3-0.7 Rapp-2021
  Double_t muBM[NM][NPMAX];
  Double_t YM[NM][NPMAX];
  for(int i=0;i<NM;i++) {
    for(int j=0;j<NPM[i];j++) {
      YM[i][j] = YMo[i][j] * gr_S->Eval(EM[i][j]) / YMo[0][j];
      muBM[i][j] = fMuB->Eval(EM[i][j]);
    }
  }

  const Int_t NEL = 4;
  const Double_t EL[NEL] = {2.4, 7.7, 19.6, 200.};
  const Char_t* TextL[NEL] = {"2.4", "7.7", "19.6", "200"};
  Double_t MuBL[NEL];
  for(int i=0;i<NEL;i++) {
      MuBL[i] = fMuB->Eval(EL[i]);
  }

  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->SetLogx();
  c1->SetLeftMargin(0.15);
  c1->Draw();

  double x1 = 1.9;
  double x2 = 300.;
  double y1 = 0;
  double y2 = 19.9;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->GetYaxis()->SetNdivisions(104);
  h0->SetYTitle("N_{ee}/N_{#pi^{+}#pi^{-}} #times 10^{6}");
  h0->GetYaxis()->SetTitleOffset(0.9);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();


  TGraph *gr_M[NM];
  const Int_t lineColor[NM] = { 1, 1};
  const Int_t lineStyle[NM] = { 1, 2};
  
  for(int i=0;i<NM;i++) {
    gr_M[i] = new TGraph(NPM[i], EM[i], YM[i]);
    gr_M[i]->SetName(NameM[i]);
    gr_M[i]->SetLineWidth(2);
    gr_M[i]->SetLineColor(lineColor[i]);
    gr_M[i]->SetLineStyle(lineStyle[i]);
    gr_M[i]->Draw("L");
  }

  TGraphErrors *gr_Ds = new TGraphErrors(NE, E, Y, 0, Yes);
  drawSysError(gr_Ds, 0.05, 0.3, 1, 1, 0);

  TGraphErrors *gr_D = new TGraphErrors(NE, E, Y, 0, Ye);
  gr_D->SetMarkerStyle(20);
  gr_D->SetMarkerSize(1.8);
  gr_D->SetLineWidth(2);
  gr_D->Draw("p");

  drawHistBox(x1,x2,y1,y2);
  drawText(2.8, -1.2, "3");
  drawText(9., -1.2, "10");
  drawText(27, -1.2, "30");
  drawText(85, -1.2, "100");

  TLegend *leg = new TLegend(0.2, 0.2, 0.52, 0.36);
  leg->SetLineColor(10);
  leg->AddEntry(gr_M[1], " 1st-order PT", "l");
  leg->AddEntry(gr_M[0], " no 1st-order PT", "l");
  leg->Draw();

  c1->Update();
  c1->SaveAs("DielectronPiVsE_Rapp.pdf");
  c1->SaveAs("DielectronPiVsE_Rapp.png");

  
  TCanvas *c2 = new TCanvas("c2","c2",800,0,800,600);
  c2->SetTopMargin(0.08);
  c2->SetTickx(0);
  c2->cd();

  x1 = -50;
  x2 = 850;
  y1 = 0;
  y2 = 19.9;
  h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Baryon Chemical Potential #mu_{B} (MeV)");
  h0->GetXaxis()->SetLabelOffset(0.008);
  h0->GetXaxis()->SetTitleOffset(1.1);
  h0->GetYaxis()->SetNdivisions(104);
  h0->SetYTitle("N_{ee}/N_{#pi^{+}#pi^{-}} #times 10^{6}");
  h0->GetYaxis()->SetTitleOffset(0.9);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();

  drawColorBox(fMuB->Eval(5.0), y1, fMuB->Eval(2.0), y2, 5, 0.5);

  TLine *ll[NEL];
  for(int i=0;i<NEL;i++) {
    ll[i] = new TLine(MuBL[i], y2*0.98, MuBL[i], y2);
    ll[i]->SetLineWidth(2);
    ll[i]->Draw("same");
    drawText(MuBL[i]-20, y2*1.02, TextL[i], 42, 0.04);
  }
  drawText(x1-50, y2*1.03, "#sqrt{s_{NN}}", 42, 0.045);
  

  TGraph *gr_muB_M[NM];
  
  for(int i=0;i<NM;i++) {
    gr_muB_M[i] = new TGraph(NPM[i], muBM[i], YM[i]);
    gr_muB_M[i]->SetName(NameM[i]);
    gr_muB_M[i]->SetLineWidth(2);
    gr_muB_M[i]->SetLineColor(lineColor[i]);
    gr_muB_M[i]->SetLineStyle(lineStyle[i]);
    gr_muB_M[i]->Draw("L");
  }

  TGraphErrors *gr_Ds_muB = new TGraphErrors(NE, muB, Y, 0, Yes);
  drawSysError(gr_Ds_muB, 12, 0.3, 1, 0, 0);

  TGraphErrors *gr_D_muB = new TGraphErrors(NE, muB, Y, 0, Ye);
  gr_D_muB->SetMarkerStyle(20);
  gr_D_muB->SetMarkerSize(1.8);
  gr_D_muB->SetLineWidth(2);
  gr_D_muB->Draw("p");

  leg = new TLegend(0.22, 0.22, 0.52, 0.36);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_M[1], " 1st-order PT", "l");
  leg->AddEntry(gr_M[0], " no 1st-order PT", "l");
  leg->Draw();

  drawHistBox(x1,x2,y1,y2);

  c2->Update();
  c2->SaveAs("DielectronPiVsMuB_Rapp.pdf");
  c2->SaveAs("DielectronPiVsMuB_Rapp.png");

}
