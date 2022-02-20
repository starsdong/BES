#include "style.C+"
#include "draw.C+"

void s3(){
  style();

  // parametrization from Andronic: 1710.09425
  TF1 *fTCF = new TF1("TCF","158.4/(1+exp(2.60-log(x))/0.45)",2.,1e5);
  TF1 *fMuB = new TF1("MuB","1307.5/(1+0.288*x)",2.,1e5);
  

  // const Int_t NM = 5;
  // const Char_t *FileName[NM] = {"hybridurqmd","stringmeltingampt","defaultampt","coalescencedcm","thermalmodelgsi"};
  // const Char_t *LegendName[NM] = {"Hybrid UrQMD","AMPT (string melt.)","AMPT (default)","Coal. (DCM)","Thermal - GSI"};
  const Int_t NM = 4;
  const Char_t *FileName[NM] = {"hybridurqmd","stringmeltingampt","defaultampt","thermalmodelgsi"};
  const Char_t *LegendName[NM] = {"Hybrid UrQMD","AMPT (string melt.)","AMPT (default)","Thermal - GSI"};
  TGraph *gr_M[NM];
  TGraph *gr_M_muB[NM];
  for(int i=0;i<NM;i++) {
    gr_M[i] = new TGraph(Form("xin_s3/s3%s.txt",FileName[i]),"%lg %lg");
    double muB[100];
    for(int j=0;j<gr_M[i]->GetN();j++) muB[j] = fMuB->Eval(gr_M[i]->GetX()[j]);
    gr_M_muB[i] = new TGraph(gr_M[i]->GetN(), muB, gr_M[i]->GetY());
  }

  const Int_t ND = 6;
  const Char_t *DataName[ND] = {"HypHI","E864_AuPt_0_10","STAR_AuAu_0_10","STAR_AuAu_0_80","ALICE_PbPb_0_10","ALICE_pPb_0_40"};
  const Bool_t sys[ND] = {0, 0, 1, 1, 1, 1};
  TGraphErrors *gr_D[ND];
  TGraphErrors *gr_D_sys[ND];  
  TGraphErrors *gr_D_muB[ND];
  TGraphErrors *gr_D_muB_sys[ND];  
  TFile *fin = new TFile("s3_data.root");
  for(int i=0;i<ND;i++) {
    gr_D[i] = (TGraphErrors *)fin->Get(DataName[i]);
    gr_D[i]->Print();
    double muB[10];
    for(int j=0;j<gr_D[i]->GetN();j++) muB[j] = fMuB->Eval(gr_D[i]->GetX()[j]);
    gr_D_muB[i] = new TGraphErrors(gr_D[i]->GetN(), muB, gr_D[i]->GetY(), 0, gr_D[i]->GetEY());
    
    if(!sys[i]) continue;
    
    gr_D_sys[i] = (TGraphErrors *)fin->Get(Form("%s_sys",DataName[i]));
    for(int j=0;j<gr_D_sys[i]->GetN();j++) muB[j] = fMuB->Eval(gr_D_sys[i]->GetX()[j]);
    gr_D_muB_sys[i] = new TGraphErrors(gr_D_sys[i]->GetN(), muB, gr_D_sys[i]->GetY(), 0, gr_D_sys[i]->GetEY());
  }

  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd()->SetLogx();

  double x1 = 1.5;
  double x2 = 6e3;
  double y1 = 0;
  double y2 = 1.2;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->SetYTitle("S_{3} #equiv ({}^{3}_{#Lambda}H/^{3}He)/(#Lambda/p)");
  h0->GetYaxis()->SetTitleOffset(1.1);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();

  //  drawColorBox(3.0, y1, 19.6, y2, kGreen-4, 0.3);
  drawColorBox(2.0, y1, 5.0, y2, 5, 0.5);
  
  // const Int_t lineColor[NM] = {kBlack,kRed,kGreen+2,kBlue,kBlack};
  // const Int_t lineStyle[NM] = {1,2,3,4,2};
  const Int_t lineColor[NM] = {kBlack,kRed,kGreen+2,kBlack};
  const Int_t lineStyle[NM] = {1,2,3,2};
  for(int i=0;i<NM;i++) {
    gr_M[i]->SetLineColor(lineColor[i]);
    gr_M[i]->SetLineStyle(lineStyle[i]);
    gr_M[i]->SetLineWidth(2);
    if(i==0) gr_M[i]->SetLineWidth(1);
    gr_M[i]->Draw("c");
  }

  const Int_t markerColor[ND] = {1, 1, 1, 1, 1, 1};
  const Int_t markerStyle[ND] = {24, 26, 20, 20, 21, 21};
  const Double_t markerSize[ND] = {1.8,1.8, 2.2, 2.2, 1.8, 1.8};
  for(int i=0;i<ND;i++) {
    if(sys[i]) {
      setGraphMarker(gr_D_sys[i], markerStyle[i], markerColor[i], markerSize[i]);
      setGraphLine(gr_D_sys[i], 1, markerColor[i], 2);
      drawSysError(gr_D_sys[i], 0.08, 0.02, 1, 1, 0);
    }
    setGraphMarker(gr_D[i], markerStyle[i], markerColor[i], markerSize[i]);
    setGraphLine(gr_D[i], 1, markerColor[i], 2);
    gr_D[i]->Draw("p");
  }
  
  TLegend *leg = new TLegend(0.62, 0.2, 0.89, 0.4);
  leg->SetLineColor(10);
  leg->SetTextSize(0.035);
  for(int i=0;i<NM;i++) 
    leg->AddEntry(gr_M[i], LegendName[i], "l");
  leg->Draw();

  leg = new TLegend(0.76, 0.74, 0.96, 0.94);
  leg->SetLineColor(10);
  leg->SetTextSize(0.035);
  leg->AddEntry(gr_D[2], "  STAR", "pl");
  leg->AddEntry(gr_D[4], "  ALICE", "pl");
  leg->AddEntry(gr_D[0], "  HypHI", "pl");
  leg->AddEntry(gr_D[1], "  E864", "pl");
  leg->Draw();

  
  drawHistBox(x1,x2,y1,y2);
  drawText(8.2, -0.07, "10");
  drawText(80, -0.07, "100");
  drawText(700, -0.07, "1000");
  
  c1->Update();
  c1->SaveAs("S3VsE.pdf");
  c1->SaveAs("S3VsE.png");

  
  TCanvas *c2 = new TCanvas("c2","c2",800,0,800,600);
  c2->cd();

  x1 = -50;
  x2 = 990;
  y1 = 0;
  y2 = 1.2;
  h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("#mu_{B} (MeV)");
  h0->GetXaxis()->SetLabelOffset(0.005);
  h0->SetYTitle("S_{3} #equiv ({}^{3}_{#Lambda}H/^{3}He)/(#Lambda/p)");
  h0->GetYaxis()->SetTitleOffset(1.1);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();

  //  drawColorBox(fMuB->Eval(19.6), y1, fMuB->Eval(3.0), y2, kGreen-4, 0.3);
  drawColorBox(fMuB->Eval(5.0), y1, fMuB->Eval(2.0), y2, 5, 0.5);

  for(int i=0;i<NM;i++) {
    gr_M_muB[i]->SetLineColor(lineColor[i]);
    gr_M_muB[i]->SetLineStyle(lineStyle[i]);
    gr_M_muB[i]->SetLineWidth(2);
    if(i==0) gr_M_muB[i]->SetLineWidth(1);
    gr_M_muB[i]->Draw("c");
  }

  for(int i=0;i<ND;i++) {
    if(sys[i]) {
      setGraphMarker(gr_D_muB_sys[i], markerStyle[i], markerColor[i], markerSize[i]);
      setGraphLine(gr_D_muB_sys[i], 1, markerColor[i], 2);
      drawSysError(gr_D_muB_sys[i], 12, 0.02, 1, 0, 0);
    }
    setGraphMarker(gr_D_muB[i], markerStyle[i], markerColor[i], markerSize[i]);
    setGraphLine(gr_D_muB[i], 1, markerColor[i], 2);
    gr_D_muB[i]->Draw("p");
  }
  
  leg = new TLegend(0.7, 0.74, 0.95, 0.94);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0.001);
  leg->SetLineColor(10);
  leg->SetTextSize(0.035);
  for(int i=0;i<NM;i++) 
    leg->AddEntry(gr_M_muB[i], LegendName[i], "l");
  leg->Draw();

  leg = new TLegend(0.55, 0.74, 0.7, 0.94);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0.001);
  leg->SetTextSize(0.035);
  leg->AddEntry(gr_D[2], "  STAR", "pl");
  leg->AddEntry(gr_D[4], "  ALICE", "pl");
  leg->AddEntry(gr_D[0], "  HypHI", "pl");
  leg->AddEntry(gr_D[1], "  E864", "pl");
  leg->Draw();

  drawHistBox(x1,x2,y1,y2);

  c2->Update();
  c2->SaveAs("S3VsMuB.pdf");
  c2->SaveAs("S3VsMuB.png");
}


