#include "draw.C+"
#include "style.C+"

void draw_3gev_h3l_dndy(){

  style();

  gROOT->LoadMacro("~/work/work/C/scaleData.C");

  TFile *infile;
  infile = new TFile("./hyper_rapidity_plot.root");
  const Double_t sc = 1e3;

  TGraphErrors*	h3l_syst_cent0010_tmp; h3l_syst_cent0010_tmp = (TGraphErrors*)infile->Get("h3l_syst_cent0010"); h3l_syst_cent0010_tmp->Print();
  TGraphErrors *h3l_syst_cent0010 = scale(h3l_syst_cent0010_tmp, sc); h3l_syst_cent0010->Print();
  TGraphErrors*	h3l_stat_cent0010_tmp; h3l_stat_cent0010_tmp = (TGraphErrors*)infile->Get("h3l_stat_cent0010"); TGraphErrors *h3l_stat_cent0010 = scale(h3l_stat_cent0010_tmp, sc);  
  TH1D*	h3l_coal_cent0010_tmp; h3l_coal_cent0010_tmp = (TH1D*)infile->Get("h3l_coal_cent0010"); h3l_coal_cent0010_tmp->Scale(sc); TGraph *h3l_coal_cent0010 = new TGraph(h3l_coal_cent0010_tmp);

  TGraphErrors*	h3l_syst_cent1050_tmp; h3l_syst_cent1050_tmp = (TGraphErrors*)infile->Get("h3l_syst_cent1050"); TGraphErrors *h3l_syst_cent1050 = scale(h3l_syst_cent1050_tmp, sc);
  TGraphErrors*	h3l_stat_cent1050_tmp; h3l_stat_cent1050_tmp = (TGraphErrors*)infile->Get("h3l_stat_cent1050"); TGraphErrors *h3l_stat_cent1050 = scale(h3l_stat_cent1050_tmp, sc);
  TH1D*	h3l_coal_cent1050_tmp; h3l_coal_cent1050_tmp = (TH1D*)infile->Get("h3l_coal_cent1050"); h3l_coal_cent1050_tmp->Scale(sc); TGraph *h3l_coal_cent1050 = new TGraph(h3l_coal_cent1050_tmp);
  

  TCanvas *x08 = new TCanvas("x08","x08",1000,800);
  // x08->cd()->SetRightMargin(0.05);
  // x08->cd()->SetTopMargin(0.04);
  x08->cd()->SetLeftMargin(0.15);
  // x08->cd()->SetBottomMargin(0.14);
  x08->cd();

  double x1 = -0.9;
  double x2 = 0.1;
  double y1 = 0.;
  double y2 = 5;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->SetXTitle("Rapidity y");
  h0->SetYTitle("B.R.#times dN/dy (#times 10^{-3})");
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetXaxis()->SetLabelSize(0.045);
  h0->GetYaxis()->SetTitleOffset(0.9);
  h0->Draw();

  
  drawSysBox(h3l_syst_cent0010, 0.02, 18);
  drawSysBox(h3l_syst_cent1050, 0.02, 18);

  h3l_coal_cent0010->SetLineWidth(2);  
  h3l_coal_cent0010->SetLineStyle(2);  
  h3l_coal_cent0010->Draw("c");
  h3l_coal_cent1050->SetLineWidth(2);  
  h3l_coal_cent1050->SetLineColor(2);  
  h3l_coal_cent1050->SetLineStyle(2);  
  h3l_coal_cent1050->Draw("c");

  h3l_stat_cent0010->SetLineWidth(2);
  h3l_stat_cent0010->Draw("psame");  
  h3l_stat_cent1050->SetMarkerColor(2);
  h3l_stat_cent1050->SetLineWidth(2);
  h3l_stat_cent1050->SetLineColor(2);
  h3l_stat_cent1050->Draw("psame");


  drawText(-0.84, 4.4, "Au + Au @ 3 GeV", 42, 0.05);
  
  TLegend *leg;
  leg = new TLegend(0.18, 0.70,0.4,0.84);
  leg->SetLineColor(10);
  leg->SetTextSize(0.045);
  leg->SetTextFont(42);
  leg->AddEntry(h3l_stat_cent0010," data","p");
  leg->AddEntry(h3l_coal_cent0010," JAM+Coal.","l");
  leg->Draw();

  drawText(-0.05, 4.4, "{}^{3}_{#Lambda}H", 52, 0.06);

  drawText(-0.07, 2.3, "0-10%", 52, 0.04);
  drawText(-0.07, 0.65, "10-50%", 52, 0.04, 0, 2);
  
  drawHistBox(x1,x2,y1,y2);

  x08->SaveAs("fig/h3l_dNdy_3GeV.pdf");
  x08->SaveAs("fig/h3l_dNdy_3GeV.png");





}
