#include "draw.C+"
#include "style.C+"

void makeBES_PT()
{
  style();

  // parametrization from Andronic: 1710.09425
  TF1 *fTCF = new TF1("TCF","158.4/(1+exp(2.60-log(x))/0.45)",2.,1e5);
  TF1 *fMuB = new TF1("MuB","1307.5/(1+0.288*x)",2.,1e5);

  const Int_t NE = 6;
  const Double_t E[NE] = {17.3, 19.6, 27.0, 39.0, 62.4, 200};
  const Double_t Y[NE] = {7.44, 7.00, 5.65, 9.51, 7.85, 8.29};
  const Double_t Ye[NE] = {0.26, 2.83, 2.07, 1.28, 1.78, 1.43};
  const Double_t Yes[NE] = {0.00, 1.19, 0.72, 1.02, 1.22, 1.67};
  const Double_t Yec[NE] = {0.00, 0.00, 0.64, 0.75, 0.85, 0.00};
  Double_t muB[NE];
  Double_t Yet[NE], Yest[NE];
  for(int i=0;i<NE;i++) {
    muB[i] = fMuB->Eval(E[i]);
    Yet[i] = sqrt(Ye[i]*Ye[i]+Yes[i]*Yes[i]+Yec[i]*Yec[i]);
    Yest[i] = sqrt(Yes[i]*Yes[i]+Yec[i]*Yec[i]);
  }

  const Int_t NM = 5;  // Number of models
  const Int_t NPM[NM] = {5, 5, 5, 7, 7};
  const Char_t *NameM[NM] = {"Rapp","Endres","PHSD","Rapp-21-PT","Rapp-21-noPT"};
  const Double_t EM[NM][10] = {{19.6, 27.0, 39.0, 62.4, 200., 0, 0, 0, 0},
			       {19.6, 27.0, 39.0, 62.4, 200., 0, 0, 0, 0},
			       //  const Double_t EM[NM][10] = {{17.3, 19.6, 27.0, 39.0, 62.4, 200., 0, 0, 0, 0},
			       //			       {17.3, 19.6, 27.0, 39.0, 62.4, 200., 0, 0, 0, 0},
			       {19.6, 27.0, 39.0, 62.4, 200., 0, 0, 0, 0, 0},
			       {4.75, 6.32, 7.69, 9.18, 11.51, 14.5, 19.57, 0, 0, 0},
			       {4.75, 6.32, 7.69, 9.18, 11.51, 14.5, 19.57, 0, 0, 0}};
  const Double_t YMo[NM][10] = {{7.81, 8.02, 8.15, 8.86, 10.0, 0, 0, 0, 0, 0},      // 0.4 - 0.75
				{7.09, 7.13, 7.47, 7.78, 11.4, 0, 0, 0, 0, 0},
				{6.23, 7.15, 7.69, 7.39, 9.64, 0, 0, 0, 0, 0},
				{16.35, 13.41, 13.41, 13.41, 13.52, 13.70, 13.94, 0, 0, 0},  // 0.3-0.7 Rapp-2021
				{11.29, 11.35, 11.88, 12.29, 12.76, 13.29, 13.88, 0, 0, 0}}; // 0.3-0.7 Rapp-2021
  const Double_t sc[NM] = {1., 1., 1., 0.56, 0.56};   // 0.56 arbitrary scaling factor to match 19.6 calculations
  Double_t muBM[NM][10];
  Double_t YM[NM][10];
  for(int i=0;i<NM;i++) {
    for(int j=0;j<NPM[i];j++) {
      YM[i][j] = YMo[i][j] * sc[i];
      muBM[i][j] = fMuB->Eval(EM[i][j]);
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->SetLogx();
  c1->SetLeftMargin(0.15);
  c1->Draw();

  double x1 = 1.9;
  double x2 = 300.;
  double y1 = 0;
  double y2 = 15;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->SetYTitle("N_{ee}/N_{ch} #times 10^{6}");
  h0->GetYaxis()->SetTitleOffset(0.9);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();

  TGraph *gr_M[NM];
  const Int_t lineColor[NM] = {2, 3, 4, 1, 1};
  const Int_t lineStyle[NM] = {2, 2, 3, 1, 2};
  
  for(int i=0;i<NM;i++) {
    if(i==1 || i==2) continue;
    gr_M[i] = new TGraph(NPM[i], EM[i], YM[i]);
    gr_M[i]->SetName(NameM[i]);
    gr_M[i]->SetLineWidth(2);
    gr_M[i]->SetLineColor(lineColor[i]);
    gr_M[i]->SetLineStyle(lineStyle[i]);
    gr_M[i]->Draw("L");
  }

  TGraphErrors *gr_Ds = new TGraphErrors(NE, E, Y, 0, Yest);
  gr_Ds->RemovePoint(0);
  drawSysError(gr_Ds, 0.05, 0.3, 1, 1, 0);

  TGraphErrors *gr_D = new TGraphErrors(NE, E, Y, 0, Ye);
  gr_D->RemovePoint(0); // remove NA60 data points, STAR only
  gr_D->SetMarkerStyle(20);
  gr_D->SetMarkerSize(1.8);
  gr_D->SetLineWidth(2);
  gr_D->Draw("p");

  TGraphErrors *gr_NA60 = new TGraphErrors(1, E, Y, 0, Ye);
  gr_NA60->SetMarkerStyle(24);
  gr_NA60->SetMarkerSize(1.8);
  gr_NA60->SetLineWidth(2);
  gr_NA60->Draw("p");

  drawHistBox(x1,x2,y1,y2);
  drawText(2.8, -0.9, "3");
  drawText(9., -0.9, "10");
  drawText(27, -0.9, "30");
  drawText(85, -0.9, "100");

  TLegend *leg = new TLegend(0.2, 0.2, 0.52, 0.36);
  leg->SetLineColor(10);
  leg->AddEntry(gr_M[3], " 1st-order PT", "l");
  leg->AddEntry(gr_M[4], " no 1st-order PT", "l");
  leg->Draw();

  c1->Update();
  c1->SaveAs("DielectronVsE_Rapp.pdf");
  c1->SaveAs("DielectronVsE_Rapp.png");
  
  TCanvas *c2 = new TCanvas("c2","c2",800,0,800,600);
  c2->cd();

  x1 = -50;
  x2 = 990;
  y1 = 0;
  y2 = 15;
  h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("#mu_{B} (MeV)");
  h0->GetXaxis()->SetLabelOffset(0.005);
  h0->SetYTitle("N_{ee}/N_{ch} #times 10^{6}");
  h0->GetYaxis()->SetTitleOffset(0.9);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();

  TGraph *gr_muB_M[NM];
  
  for(int i=0;i<NM;i++) {
    if(i==1 || i==2) continue;
    gr_muB_M[i] = new TGraph(NPM[i], muBM[i], YM[i]);
    gr_muB_M[i]->SetName(NameM[i]);
    gr_muB_M[i]->SetLineWidth(2);
    gr_muB_M[i]->SetLineColor(lineColor[i]);
    gr_muB_M[i]->SetLineStyle(lineStyle[i]);
    gr_muB_M[i]->Draw("L");
  }

  TGraphErrors *gr_Ds_muB = new TGraphErrors(NE, muB, Y, 0, Yest);
  gr_Ds_muB->RemovePoint(0);
  drawSysError(gr_Ds_muB, 12, 0.3, 1, 0, 0);

  TGraphErrors *gr_D_muB = new TGraphErrors(NE, muB, Y, 0, Ye);
  gr_D_muB->RemovePoint(0); // remove NA60 data points, STAR only
  gr_D_muB->SetMarkerStyle(20);
  gr_D_muB->SetMarkerSize(1.8);
  gr_D_muB->SetLineWidth(2);
  gr_D_muB->Draw("p");

  TGraphErrors *gr_NA60_muB = new TGraphErrors(1, muB, Y, 0, Ye);
  gr_NA60_muB->SetMarkerStyle(24);
  gr_NA60_muB->SetMarkerSize(1.8);
  gr_NA60_muB->SetLineWidth(2);
  gr_NA60_muB->Draw("p");

  leg = new TLegend(0.6, 0.2, 0.92, 0.36);
  leg->SetLineColor(10);
  leg->AddEntry(gr_M[3], " 1st-order PT", "l");
  leg->AddEntry(gr_M[4], " no 1st-order PT", "l");
  leg->Draw();

  drawHistBox(x1,x2,y1,y2);

  c2->Update();
  c2->SaveAs("DielectronVsMuB_Rapp.pdf");
  c2->SaveAs("DielectronVsMuB_Rapp.png");

}
