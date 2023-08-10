#include "style.C+"
#include "draw.C+"

void draw_h3l_snn(){

  style();
  /*
  gStyle->SetTitleFont(62,"X");
  gStyle->SetTitleFont(62,"Y");
  gStyle->SetLabelFont(62,"X");
  gStyle->SetLabelFont(62,"Y");
  gStyle->SetTextFont(62);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);

  gStyle->SetLabelSize(0.0575,"X");
  gStyle->SetTitleSize(0.0575,"X");
  gStyle->SetLabelSize(0.0575,"Y");
  gStyle->SetTitleSize(0.0575,"Y");
  */

  double xxxxx  = 1600;
  double lrmargin = 120;
  double rrmargin = 10;
  double panelwid = 600;
  xxxxx = 1*panelwid+rrmargin+lrmargin;

  TCanvas *x0 = new TCanvas("x09","x09",850,xxxxx);
  x0->cd()->SetRightMargin(0.05);
  x0->cd()->SetTopMargin(0.04);
  x0->cd()->SetLeftMargin(0.1);
  x0->cd()->SetBottomMargin(0.14);
  x0->cd();

  TPad*x09 = new TPad("f5mass", "f5mass",0, 0,1,1);
  x09->SetLeftMargin(0.18);
  x09->SetRightMargin(0.02);
  x09->SetTopMargin(rrmargin/xxxxx);
  x09->SetBottomMargin(lrmargin/xxxxx);
  x09->SetFrameBorderMode(0);
  x09->SetFrameBorderMode(0);
  x09->Draw();
  x09->cd();


  //get data
  double sys_width = 0.11;

  double h3l_br = 0.25;
  double h4l_br = 0.5;

  TGraphErrors *alice_pbpb_2p76TeV; alice_pbpb_2p76TeV = new TGraphErrors();//from https://arxiv.org/pdf/1506.08453.pdf
  alice_pbpb_2p76TeV->SetMarkerStyle(20);
  alice_pbpb_2p76TeV->SetMarkerColorAlpha(kBlack,0);
  alice_pbpb_2p76TeV->SetLineColor(kBlack);
  alice_pbpb_2p76TeV->SetPoint(0,2760,3.86*1e-5/0.25);//0-10
  alice_pbpb_2p76TeV->SetPointError(0,0,0.77*1e-5/0.25);//0-10

  TGraphErrors *alice_pbpb_2p76TeV_sys; alice_pbpb_2p76TeV_sys = new TGraphErrors();
  alice_pbpb_2p76TeV_sys->SetMarkerStyle(20);
  alice_pbpb_2p76TeV_sys->SetMarkerColorAlpha(kBlack,0);
  alice_pbpb_2p76TeV_sys->SetLineColor(kBlack);
  alice_pbpb_2p76TeV_sys->SetFillColorAlpha(kBlack,0.15);
  alice_pbpb_2p76TeV_sys->SetPoint(0,2760,3.86*1e-5/0.25);//0-10
  alice_pbpb_2p76TeV_sys->SetPointError(0,2760*0.13,0.68*1e-5/0.25);//0-10

  TGraphErrors *alice_pbpb_2p76TeV_sys_clone; alice_pbpb_2p76TeV_sys_clone = (TGraphErrors*)alice_pbpb_2p76TeV_sys->Clone("alice_pbpb_2p76TeV_sys_clone");
  TGraphErrors *alice_pbpb_2p76TeV_clone; alice_pbpb_2p76TeV_clone = (TGraphErrors*)alice_pbpb_2p76TeV->Clone("alice_pbpb_2p76TeV_clone");

  double x_scale = 22; //x-axis contraction
  for(int i=0;i<1;i++){
    double xx,yy,exx,eyy;
    alice_pbpb_2p76TeV->GetPoint(i,xx,yy);
    eyy = alice_pbpb_2p76TeV->GetErrorY(i);
    alice_pbpb_2p76TeV_clone->SetPoint(i,xx/x_scale,yy);
    alice_pbpb_2p76TeV_clone->SetPointError(i,0,eyy);

    alice_pbpb_2p76TeV_sys->GetPoint(i,xx,yy);
    eyy = alice_pbpb_2p76TeV_sys_clone->GetErrorY(i);
    exx = alice_pbpb_2p76TeV_sys_clone->GetErrorX(i);
    alice_pbpb_2p76TeV_sys_clone->SetPoint(i,xx/x_scale,yy);
    alice_pbpb_2p76TeV_sys_clone->SetPointError(i,xx/x_scale*sys_width,eyy);
  }

  alice_pbpb_2p76TeV_clone->SetMarkerStyle(33);
  alice_pbpb_2p76TeV_clone->SetMarkerSize(3);
  alice_pbpb_2p76TeV_clone->SetMarkerColor(kAzure+7);
  alice_pbpb_2p76TeV_sys_clone->SetFillColorAlpha(kAzure+7,0.2);
  alice_pbpb_2p76TeV_sys_clone->SetLineColorAlpha(kAzure+7,0.);

  TGraphErrors *star_h3l_3GeV; star_h3l_3GeV = new TGraphErrors();
  TGraphErrors *star_h3l_3GeV_sys; star_h3l_3GeV_sys = new TGraphErrors();
  star_h3l_3GeV ->SetPoint(0, 3, 0.0126136);
  star_h3l_3GeV ->SetPointError(0, 0, 0.0016014);
  star_h3l_3GeV_sys->SetPoint(0, 3, 0.0126136);
  star_h3l_3GeV_sys ->SetPointError(0, 3*sys_width, 0.00402985);

  TGraphErrors *star_h4l_3GeV; star_h4l_3GeV = new TGraphErrors();
  TGraphErrors *star_h4l_3GeV_sys; star_h4l_3GeV_sys = new TGraphErrors();
  star_h4l_3GeV ->SetPoint(0, 3, 0.00565111);
  star_h4l_3GeV ->SetPointError(0, 0, 0.00024674);
  star_h4l_3GeV_sys->SetPoint(0, 3, 0.00565111);
  star_h4l_3GeV_sys ->SetPointError(0, 3*sys_width, 0.000938099);

  star_h3l_3GeV->SetLineColor(kBlack);
  star_h4l_3GeV->SetLineColor(kBlack);
  star_h3l_3GeV_sys->SetLineColor(kBlack);
  star_h4l_3GeV_sys->SetLineColor(kBlack);
  star_h3l_3GeV->SetMarkerColor(kBlack);
  star_h4l_3GeV->SetMarkerColor(kBlack);
  star_h3l_3GeV_sys->SetMarkerColor(kBlack);
  star_h4l_3GeV_sys->SetMarkerColor(kBlack);
  star_h3l_3GeV_sys->SetFillColorAlpha(kGray,1.);
  star_h4l_3GeV_sys->SetFillColorAlpha(kGray,1.);
  star_h3l_3GeV_sys->SetLineColorAlpha(kBlack,0.);
  star_h4l_3GeV_sys->SetLineColorAlpha(kBlack,0.);
  star_h3l_3GeV->SetMarkerStyle(20);
  star_h3l_3GeV_sys->SetMarkerStyle(20);
  star_h4l_3GeV->SetMarkerStyle(20);
  star_h4l_3GeV_sys->SetMarkerStyle(20);
  star_h3l_3GeV->SetMarkerSize(2.5);
  star_h3l_3GeV_sys->SetMarkerSize(2.5);
  star_h4l_3GeV->SetMarkerSize(2.5);
  star_h4l_3GeV_sys->SetMarkerSize(2.5);

  //may23
  // double dndy_3p9gev_h3l2b = 0.00329859/h3l_br;
  // double dndy_3p9gev_h3l2b_sys = 0.000785353/h3l_br;
  // double dndy_3p9gev_h3l2b_stat = 0.000700843/h3l_br;
  //june14 
  double dndy_3p9gev_h3l2b = 0.00372669/h3l_br;
  double dndy_3p9gev_h3l2b_sys = 0.00110904/h3l_br;
  double dndy_3p9gev_h3l2b_stat = 0.000608711/h3l_br;

  double dndy_4p5gev_h3l2b = 0.00250844/h3l_br;
  double dndy_4p5gev_h3l2b_sys = 0.000390693/h3l_br;
  double dndy_4p5gev_h3l2b_stat = 0.00058811/h3l_br;

  double dndy_27gev_h3l = 0.000227493;
  double dndy_27gev_stat = 2.71903e-05;
  double dndy_27gev_sys = 4.61068e-05;

  double dndy_19gev_h3l = 0.00043212;
  double dndy_19gev_stat = 3.29159e-05;
  double dndy_19gev_sys = 7.11475e-05;

  TGraphErrors *star_prelim;
  TGraphErrors *star_prelim_sys;
  star_prelim = new TGraphErrors();
  star_prelim_sys = new TGraphErrors();
  
  TGraphErrors *star_prelim_FXT;
  TGraphErrors *star_prelim_FXT_sys;
  star_prelim_FXT = new TGraphErrors();
  star_prelim_FXT_sys = new TGraphErrors();


  // FXT energies 
  int npoints_fxt=0;
  star_prelim_FXT->SetPoint(npoints_fxt,3.9,dndy_3p9gev_h3l2b);
  star_prelim_FXT_sys->SetPoint(npoints_fxt,3.9,dndy_3p9gev_h3l2b);
  star_prelim_FXT->SetPointError(npoints_fxt,0,dndy_3p9gev_h3l2b_stat);
  star_prelim_FXT_sys->SetPointError(npoints_fxt,3.9*sys_width,dndy_3p9gev_h3l2b_sys);
  npoints_fxt++;

  star_prelim_FXT->SetPoint(npoints_fxt,4.5,dndy_4p5gev_h3l2b);
  star_prelim_FXT_sys->SetPoint(npoints_fxt,4.5,dndy_4p5gev_h3l2b);
  star_prelim_FXT->SetPointError(npoints_fxt,0,dndy_4p5gev_h3l2b_stat);
  star_prelim_FXT_sys->SetPointError(npoints_fxt,4.5*sys_width,dndy_4p5gev_h3l2b_sys);
  npoints_fxt++;

  int npoints=0;
  // COL energies
  star_prelim->SetPoint(npoints,27,dndy_27gev_h3l);
  star_prelim_sys->SetPoint(npoints,27,dndy_27gev_h3l);
  star_prelim->SetPointError(npoints,0,dndy_27gev_stat);
  star_prelim_sys->SetPointError(npoints,27*sys_width,dndy_27gev_sys);
  npoints++;

  star_prelim->SetPoint(npoints,19.6,dndy_19gev_h3l);
  star_prelim_sys->SetPoint(npoints,19.6,dndy_19gev_h3l);
  star_prelim->SetPointError(npoints,0,dndy_19gev_stat);
  star_prelim_sys->SetPointError(npoints,19.6*sys_width,dndy_19gev_sys);
  npoints++;


  star_prelim_FXT_sys->SetLineColorAlpha(kWhite,0);
  star_prelim_FXT->SetMarkerSize(2.5);
  star_prelim_FXT_sys->SetMarkerSize(2.5);
  star_prelim_FXT_sys->SetLineColor(kRed-9);
  star_prelim_FXT->SetLineColor(kRed);
  star_prelim_FXT_sys->SetMarkerColor(kRed-9);
  star_prelim_FXT_sys->SetMarkerStyle(21);
  star_prelim_FXT_sys->SetFillColorAlpha(kRed-9,1);
  star_prelim_FXT->SetMarkerColor(kRed);
  star_prelim_FXT->SetMarkerStyle(21);


  star_prelim_sys->SetLineColorAlpha(kWhite,0);
  star_prelim->SetMarkerSize(2.5);
  star_prelim_sys->SetMarkerSize(2.5);

  star_prelim_sys->SetMarkerColor(kBlack);
  star_prelim_sys->SetMarkerStyle(21);
  star_prelim_sys->SetFillColorAlpha(kGray,1);
  star_prelim->SetMarkerColor(kBlack);
  star_prelim->SetMarkerStyle(21);

  TGraphErrors *qm2023_prelim;//placeholder
  qm2023_prelim = new TGraphErrors();
  qm2023_prelim->SetPoint(0,7.7, 0.00350599);//pT 0.1
  qm2023_prelim->SetPointError(0,0,0.000222685);
  qm2023_prelim->SetPoint(1,11.5, 0.00110522);//pT 0.1
  qm2023_prelim->SetPointError(1,0, 9.57651e-05);
  qm2023_prelim->SetPoint(2,14.6, 0.000525826);//pT 0.1
  qm2023_prelim->SetPointError(2,0,4.61124e-05);

  qm2023_prelim->SetLineColor(kBlack);
  qm2023_prelim->SetMarkerColor(kBlack);
  qm2023_prelim->SetMarkerStyle(21);
  qm2023_prelim->SetMarkerSize(2.5);
  //end get data

  //get theory

  double x[1000];
  double y[1000];
  double ymid[1000];
  double ymin[1000];
  double ymax[1000];

  TFile *_file1;
  _file1 = new TFile("hyper_model/anton_snn.root");

  TGraphAsymmErrors* tHlambda;
  TGraphAsymmErrors* fHlambda;

  tHlambda = (TGraphAsymmErrors*)_file1->Get("tHlambda");
  fHlambda = (TGraphAsymmErrors*)_file1->Get("fHlambda");

  tHlambda->RemovePoint(18);
  tHlambda->RemovePoint(17);
  tHlambda->RemovePoint(16);
  tHlambda->RemovePoint(15);
  tHlambda->RemovePoint(14);

  fHlambda->RemovePoint(18);
  fHlambda->RemovePoint(17);
  fHlambda->RemovePoint(16);
  fHlambda->RemovePoint(15);

  tHlambda->SetLineColor(kRed);
  fHlambda->SetLineColor(kRed);
  tHlambda->SetLineWidth(2);
  fHlambda->SetLineWidth(2);
  tHlambda->SetLineStyle(2);
  fHlambda->SetLineStyle(2);

  TGraph *sHlambda;
  sHlambda = new TGraph();
  for(int i=0;i<tHlambda->GetN();i++) {
  tHlambda->GetPoint(i,x[i],ymid[i]);
  ymax[i] = ymid[i]+tHlambda->GetErrorYhigh(i);;
  ymin[i] = ymid[i]-tHlambda->GetErrorYlow(i);;
  y[i] = ymid[i];
  }
  for(int i=0;i<tHlambda->GetN();i++) {
  sHlambda->SetPoint(i,x[i],ymax[i]);
  sHlambda->SetPoint(tHlambda->GetN()+i,x[tHlambda->GetN()-i-1],ymin[tHlambda->GetN()-i-1]);
  tHlambda->SetPointError(i,0,0,0,0);
  }
  sHlambda->SetFillColorAlpha(kRed,0.2);

  TGraph *qHlambda;
  qHlambda = new TGraph();
  for(int i=0;i<fHlambda->GetN();i++) {

    fHlambda->GetPoint(i,x[i],ymid[i]);
  ymax[i] = ymid[i]+fHlambda->GetErrorYhigh(i);;
  ymin[i] = ymid[i]-fHlambda->GetErrorYlow(i);;
  y[i] = ymid[i];
  }
  for(int i=0;i<fHlambda->GetN();i++) {
  qHlambda->SetPoint(i,x[i],ymax[i]);
  qHlambda->SetPoint(fHlambda->GetN()+i,x[fHlambda->GetN()-i-1],ymin[fHlambda->GetN()-i-1]);
  fHlambda->SetPointError(i,0,0,0,0);
  }
  qHlambda->SetFillColorAlpha(kRed,0.2);

  double horilinemargin = 0.16;

  TGraphErrors *thermal_central_h3l_alice_clone;
  thermal_central_h3l_alice_clone=new TGraphErrors();
  thermal_central_h3l_alice_clone->SetPoint(0,2760/x_scale*(1-horilinemargin)	,0.000150665	);
  thermal_central_h3l_alice_clone->SetPoint(1,2760/x_scale	,    0.000150665	);
  thermal_central_h3l_alice_clone->SetPoint(2,2760/x_scale*(1+horilinemargin)	,0.000150665	);

  thermal_central_h3l_alice_clone->SetLineColor(kRed);
  thermal_central_h3l_alice_clone->SetLineStyle(2);
  thermal_central_h3l_alice_clone->SetLineWidth(2);

  double urqmd_hybrid_2p76TeV = 8.35333E-05;
  double fist_2p76TeV = 0.000163;

  TGraphErrors *urqmd_hybrid_2p76TeV_clone; urqmd_hybrid_2p76TeV_clone = new TGraphErrors();
  urqmd_hybrid_2p76TeV_clone->SetPoint(0,2760/x_scale*(1-horilinemargin)	, urqmd_hybrid_2p76TeV	);
  urqmd_hybrid_2p76TeV_clone->SetPoint(1,2760/x_scale	, urqmd_hybrid_2p76TeV	);
  urqmd_hybrid_2p76TeV_clone->SetPoint(2,2760/x_scale*(1+horilinemargin)	, urqmd_hybrid_2p76TeV	);
  urqmd_hybrid_2p76TeV_clone->SetLineColor(kGreen+1);
  urqmd_hybrid_2p76TeV_clone->SetLineStyle(3);
  urqmd_hybrid_2p76TeV_clone->SetLineWidth(4);
  TGraphErrors *fist_2p76TeV_clone; fist_2p76TeV_clone = new TGraphErrors();
  fist_2p76TeV_clone->SetPoint(0, 2760/x_scale*(1-horilinemargin)	, fist_2p76TeV	);
  fist_2p76TeV_clone->SetPoint(1, 2760/x_scale	, fist_2p76TeV	);
  fist_2p76TeV_clone->SetPoint(2, 2760/x_scale*(1+horilinemargin)	, fist_2p76TeV	);
  fist_2p76TeV_clone->SetLineColor(kBlue);
  fist_2p76TeV_clone->SetLineStyle(7);
  fist_2p76TeV_clone->SetLineWidth(3);

  TFile *_file2;
  _file2 = new TFile("hyper_model/phqmd_yields_hypernuclei_auau_5centr.root");
  TDirectoryFile *h3;
  TDirectoryFile *h4;
  h3 = (TDirectoryFile*)_file2->Get("H3L+");
  h3->cd();
  TGraphErrors *phqmd_t80_h3l; phqmd_t80_h3l = (TGraphErrors*)h3->Get("phqmd_t80");
  h4 = (TDirectoryFile*)_file2->Get("H4L+");
  h3->cd();
  TGraphErrors *phqmd_t80_h4l; phqmd_t80_h4l = (TGraphErrors*)h4->Get("phqmd_t80");

  phqmd_t80_h3l->RemovePoint(phqmd_t80_h3l->GetN()-1);
  phqmd_t80_h3l->SetLineColor(kViolet+2);
  phqmd_t80_h3l->SetLineStyle(9);
  phqmd_t80_h4l->SetLineColor(kViolet+2);
  phqmd_t80_h4l->SetLineStyle(9);

  TGraphErrors *fist_h3l;
  fist_h3l = new TGraphErrors("hyper_model/thermal_fist_h3l_truncated.txt","%lg %lg");
  for(int i=0;i<7;i++){
    fist_h3l->RemovePoint(0);
  }
  fist_h3l->SetLineColor(kBlue);
  fist_h3l->SetLineStyle(7);
  fist_h3l->SetLineWidth(3);

  TGraphErrors *urqmd_coal_h3l;
  urqmd_coal_h3l = new TGraphErrors("hyper_model/urqmd_coal.txt","%lg %lg");
  urqmd_coal_h3l->SetLineColor(kGreen+2);
  urqmd_coal_h3l->SetLineWidth(3);
  //end get theory



  /////////////////////////
  ///Start ploting 
  /////////////////////////



  double cutoff_e = 80;
  double high_e = 300;
  double low_x=1.6;double high_x=high_e;double high_y=4e-1; double low_y=4e-5;

  TH1D *gdummy5;
  gdummy5 = new TH1D("gdummy5","",1,low_x,high_x);
  gdummy5 ->GetYaxis()->SetRangeUser(low_y,high_y);

  // gdummy5->GetYaxis()->SetTitleOffset(5.);
  gdummy5->GetYaxis()->SetTitle("{}^{3}_{#Lambda}H dN/dy (|y|<0.5)");
  gdummy5->GetXaxis()->SetTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  gdummy5->GetXaxis()->CenterTitle();
  gdummy5->GetXaxis()->SetLabelFont(42);
  gdummy5->GetXaxis()->SetLabelSize(0.0001);
  gdummy5->GetXaxis()->SetLabelOffset(999.);
  gdummy5->GetYaxis()->SetLabelFont(42);
  gdummy5->GetYaxis()->SetLabelSize(0.045);

  gdummy5->GetYaxis()->SetTitleFont(42);
  gdummy5->GetYaxis()->SetTitleSize(0.055);
  gdummy5->GetXaxis()->SetTitleFont(42);
  gdummy5->GetXaxis()->SetTitleSize(0.055);
  gdummy5->GetXaxis()->SetTitleOffset(1.2);
  gdummy5->GetYaxis()->SetTitleOffset(1.35);

  gdummy5 ->Draw();
  gPad->SetLogy();
  gPad->SetLogx();

  sHlambda->Draw("fsame");
  tHlambda->Draw("same");

  phqmd_t80_h3l->Draw("same");
  fist_h3l->Draw("same");
  urqmd_coal_h3l->Draw("same");

  thermal_central_h3l_alice_clone->Draw("same");
  urqmd_hybrid_2p76TeV_clone->Draw("same");
  fist_2p76TeV_clone->Draw("same");

  alice_pbpb_2p76TeV_sys_clone->Draw("p5same");
  alice_pbpb_2p76TeV_clone->Draw("psame");
  
  drawText(100,low_y*9,"Pb+Pb",42,0.035,0,kAzure+7);
  drawText(100,low_y*6,"2.76TeV",42,0.035,0,kAzure+7);

  star_h3l_3GeV_sys->Draw("p5same");
  star_h3l_3GeV->Draw("psame");

  //  star_prelim_sys->Draw("p5same");
  //  star_prelim->Draw("psame");

  //  star_prelim_FXT_sys->Draw("p5same");
  //  star_prelim_FXT->Draw("psame");


  //  qm2023_prelim->Draw("psame");

  const Int_t Label[3] = {3, 10, 30};
  const Double_t offset[3] = {0.95, 0.88, 0.88};
  for(int i=0;i<3;i++) {
    drawText(Label[i]*offset[i], low_y/1.6, Form("%d",Label[i]), 42, 0.045);
  }

  // urqmd_central_h3l_alice->Draw("same");
  // phqmd_central_h3l->Draw("same");
  // jam_central_h3l->Draw("same");

  // TPaveText *t_unc7 = new TPaveText(0.20,0.06,0.38,0.30,"NDCNB");
  // t_unc7->SetTextSize(60);
  // t_unc7->SetTextFont(63);
  // t_unc7->SetFillColorAlpha(kWhite,1);
  // t_unc7->SetTextColor(kBlack);
  // t_unc7->AddText("{}_{#Lambda}^{3}H");
  // t_unc7->Draw();
  // drawLatex( 0.7,0.9, "{}_{#Lambda}^{3}H",0.065);

  TLegend *legbes7;
  // legbes7 = new TLegend(0.22,0.98,0.56,0.8);
  legbes7 = new TLegend(0.64,0.96,0.95,0.84);
  legbes7->SetLineColor(10);
  legbes7->SetTextSize(0.04);
  legbes7->SetTextFont(42);
  legbes7->AddEntry(star_h3l_3GeV,"Au+Au (STAR)","p");
  //  legbes7->AddEntry(star_prelim,"Au+Au (Preliminary)","p");
  //  legbes7->AddEntry(star_prelim_FXT,"Au+Au (This analysis)","p");
  legbes7->AddEntry(alice_pbpb_2p76TeV_clone,"Pb+Pb (ALICE)","p");
  legbes7->Draw();


  TLegend *legbes9;
  legbes9 = new TLegend(0.64,0.84,0.97,0.54);
  legbes9->SetLineColor(10);
  legbes9->SetTextSize(0.035);
  legbes9->SetTextFont(42);
  legbes9->AddEntry(tHlambda,"Thermal(GSI)","l");
  legbes9->AddEntry(fist_h3l,"Thermal-FIST","l");
  legbes9->AddEntry(urqmd_coal_h3l,"Coal. (UrQMD)","l");
  legbes9->AddEntry(urqmd_hybrid_2p76TeV_clone,"Coal. (UrQMD-hy.)","l");
  legbes9->AddEntry(phqmd_t80_h3l,"PHQMD","l");
  legbes9->Draw();

  // qHlambda->Draw("fsame");
  // fHlambda->Draw("same");

  // phqmd_t80_h4l->Draw("same");

  // TPaveText *t_unc8 = new TPaveText(0.20,0.08+0.19,0.38,0.28+0.19,"NDCNB");
  // t_unc8->SetTextSize(60);
  // t_unc8->SetTextFont(63);
  // t_unc8->SetFillColorAlpha(kWhite,0);
  // t_unc8->SetTextColor(kBlack);
  // t_unc8->AddText("{}_{#Lambda}^{4}H");
  // t_unc8->Draw();
  //
  //
  // star_h4l_3GeV->SetMarkerColor(kBlack);
  // star_h4l_3GeV->SetLineColor(kBlack);
  // star_h4l_3GeV_sys->SetMarkerColor(kBlack);
  //
  // star_h4l_3GeV_sys->SetFillColorAlpha(kBlack,0.2);
  // star_h4l_3GeV ->Draw("psame");
  // star_h4l_3GeV_sys ->Draw("p5same");

  // x09->cd();
  // TPad *f7mass;
  // f7mass = new TPad("f7mass", "f7mass",0.35,0,0.75,0.11);
  //
  // f7mass->SetLeftMargin(0.);
  // f7mass->SetRightMargin(0.);
  // f7mass->SetTopMargin(0.0);
  // f7mass->SetBottomMargin(0.);
  // f7mass->SetFrameBorderMode(0);
  // f7mass->SetFrameBorderMode(0);
  // f7mass->SetFillColorAlpha(kWhite,0);
  // f7mass->Draw();
  // f7mass->cd();

  // TPaveText *t_3 = new TPaveText(0.6, .45,.6,.45,"NDCNB");
  // t_3->SetTextSize(55);
  // t_3->SetTextFont(43);
  // t_3->SetFillColorAlpha(kWhite,0);
  // t_3->AddText("#sqrt{s_{NN}} [GeV]");
  // t_3->Draw();
  // x09->cd();
  //
  // TPaveText *t_misc0 = new TPaveText(0.266,0.01,0.306,0.1905,"NBNDC");
  // t_misc0->SetTextSize(40);
  // t_misc0->SetTextFont(43);
  // t_misc0->SetFillColorAlpha(kYellow,0);
  // t_misc0->SetLineColor(kBlack);
  // t_misc0->AddText("3");
  // t_misc0->Draw();
  //
  // TPaveText *t_misc1 = new TPaveText(0.34+0.306,0.01,0.38+0.306,0.1905,"NBNDC");
  // t_misc1->SetTextSize(40);
  // t_misc1->SetTextFont(43);
  // t_misc1->SetFillColorAlpha(kYellow,0);
  // t_misc1->SetLineColor(kBlack);
  // t_misc1->AddText("30");
  // t_misc1->Draw();
  //
  //
  // TPaveText *t_misc2 = new TPaveText(0.34+0.12,0.01,0.38+0.12,0.1905,"NBNDC");
  // t_misc2->SetTextSize(40);
  // t_misc2->SetTextFont(43);
  // t_misc2->SetFillColorAlpha(kYellow,0);
  // t_misc2->SetLineColor(kBlack);
  // t_misc2->SetTextColor(kRed);
  // t_misc2->AddText("10");
  // // t_misc2->Draw();
  //
  // TPaveText *t_misc3 = new TPaveText(0.065,0.724,0.22,0.8,"NBNDC");
  // t_misc3->SetTextSize(40);
  // t_misc3->SetTextFont(43);
  // t_misc3->SetFillColorAlpha(kYellow,0);
  // t_misc3->SetLineColor(kBlack);
  // t_misc3->AddText("10^{#minus2}");
  // // t_misc3->Draw();
  //
  // TPaveText *t_misc4 = new TPaveText(0.065,0.724-0.173,0.22,0.8-0.173,"NBNDC");
  // t_misc4->SetTextSize(40);
  // t_misc4->SetTextFont(43);
  // t_misc4->SetFillColorAlpha(kYellow,0);
  // t_misc4->SetLineColor(kBlack);
  // t_misc4->AddText("10^{#minus4}");
  // // t_misc4->Draw();
  //
  // TPaveText *t_misc5 = new TPaveText(0.065,0.44,0.22,0.51,"NBNDC");
  // t_misc5->SetTextSize(40);
  // t_misc5->SetTextFont(43);
  // t_misc5->SetFillColorAlpha(kYellow,0);
  // t_misc5->SetLineColor(kBlack);
  // t_misc5->AddText("10^{#minus2}");
  // // t_misc5->Draw();
  //
  //
  // TPaveText *t_misc6 = new TPaveText(0.68,0.26,0.94,0.32,"NBNDC");
  // t_misc6->SetTextSize(20);
  // t_misc6->SetTextFont(53);
  // t_misc6->SetFillColorAlpha(kYellow,0);
  // t_misc6->SetLineColor(kBlack);
  // t_misc6->AddText("Assuming B.R.(^{3(4)}_{#Lambda}H             ");
  // t_misc6->AddText("#rightarrow^{3(4)}He + #pi^{-}) = 25%(50%)");
  // t_misc6->Draw("same");
  //
  //
  // TPaveText *t_misc7 = new TPaveText(0.84,0.64,0.97,0.68,"NBNDC");
  // t_misc7->SetTextSize(20);
  // t_misc7->SetTextFont(63);
  // t_misc7->SetFillColorAlpha(kYellow,0);
  // t_misc7->SetLineColor(kBlack);
  // t_misc7->SetTextColor(kAzure+7);
  // t_misc7->AddText("Pb+Pb");
  // t_misc7->AddText("2.76TeV");
  // t_misc7->Draw("same");
  //
  //
  //  drawText(0.23,0.26,"Assuming B.R.(^{3}_{#Lambda}H",42,0.03);
  //  drawText(0.23,0.21,"#rightarrow^{3}He + #pi^{-}) = 25%",42,0.03);
  //  drawText( 0.3, 0.5,"{}^{3}_{#Lambda}H",42,0.06 );
  // drawLatex( 0.3, 0.8,"{}^{3}_{#Lambda}H",0.06 );
  drawText( 2.2, high_y/2.5,"A + A  0-10%",22,0.055 );
  //  drawText( 0.22, 0.93,"STAR Preliminary",42, 0.05 );

  // drawLatex( 0.75, 0.87,"0-10%",0.05 );
  // drawLatex( 0.65, 0.93,"STAR Preliminary",0.05 );

  drawHistBox(low_x, high_x, low_y, high_y);
  
  //
  // TPad *b = new TPad("b","b",0.783,0.14205+0.0001,0.9711779,0.1892+0.0001);

  TPad *b = new TPad("b","b",0.7830189,0.1467,0.9716981,0.1938);
  b->SetBorderMode(0);
  b->Draw();
  b->cd();
  double arb_mar = 0.475-0.11;
  double arb_ymar = 0.7;
  drawLine(0.0,arb_mar,0.14,arb_mar,3);
  drawLine(0.24,arb_mar,1.,arb_mar,3);
  drawLine(0.14-0.03,arb_mar-arb_ymar,0.14+0.03,arb_mar+arb_ymar,2);
  drawLine(0.24-0.03,arb_mar-arb_ymar,0.24+0.03,arb_mar+arb_ymar,2);
  /*
  TLine *line = new TLine(0.0,arb_mar,0.14,arb_mar);
  line->SetLineWidth(3);
  line->Draw();
  line = new TLine(0.24,arb_mar,1.,arb_mar);
  line->SetLineWidth(3);
  line->Draw();
  line = new TLine(0.14-0.03,arb_mar-arb_ymar,0.14+0.03,arb_mar+arb_ymar);
  line->SetLineWidth(3);
  line->Draw();
  line = new TLine(0.24-0.03,arb_mar-arb_ymar,0.24+0.03,arb_mar+arb_ymar);
  line->SetLineWidth(3);
  line->Draw();
  */
  // drawline(100*0.95-2, low_y*0.9, 100*1.05-2, low_y*1.2, 1, 0, 1 );


  // drawline(low_x, low_y, high_x, low_y, 1, 0, 1 );
  // drawline(high_x, low_y, high_x, high_y, 1, 0, 1 );
  // drawline(high_x, high_y, low_x, high_y, 1, 0, 1 );
  // drawline(low_x, low_y, low_x, high_y, 1, 0, 1 );

  x0->cd();

  x0->Update();
  x0->SaveAs("fig/h3l_snn.pdf");
  x0->SaveAs("fig/h3l_snn.png");

}
