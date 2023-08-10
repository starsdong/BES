#include "draw.C+"
#include "style.C+"

void draw_3gev_dndy(){

  style();

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

  TFile *infile;
  infile = new TFile("./hyper_rapidity_plot.root");

  TH1D*	cent0010_dummy; cent0010_dummy = (TH1D*)infile->Get("cent0010_dummy");
  TGraphErrors*	h3l_syst_cent0010; h3l_syst_cent0010 = (TGraphErrors*)infile->Get("h3l_syst_cent0010");
  TGraphErrors*	h3l_stat_cent0010; h3l_stat_cent0010 = (TGraphErrors*)infile->Get("h3l_stat_cent0010");
  TGraphErrors*	h4l_syst_cent0010; h4l_syst_cent0010 = (TGraphErrors*)infile->Get("h4l_syst_cent0010");
  TGraphErrors*	h4l_stat_cent0010; h4l_stat_cent0010 = (TGraphErrors*)infile->Get("h4l_stat_cent0010");
  TH1D*	h3l_coal_cent0010; h3l_coal_cent0010 = (TH1D*)infile->Get("h3l_coal_cent0010");
  TH1D*	h4l_coal_cent0010; h4l_coal_cent0010 = (TH1D*)infile->Get("h4l_coal_cent0010");
  TH1D*	cent1050_dummy; cent1050_dummy = (TH1D*)infile->Get("cent1050_dummy");
  TGraphErrors*	h3l_syst_cent1050; h3l_syst_cent1050 = (TGraphErrors*)infile->Get("h3l_syst_cent1050");
  TGraphErrors*	h3l_stat_cent1050; h3l_stat_cent1050 = (TGraphErrors*)infile->Get("h3l_stat_cent1050");
  TGraphErrors*	h4l_syst_cent1050; h4l_syst_cent1050 = (TGraphErrors*)infile->Get("h4l_syst_cent1050");
  TGraphErrors*	h4l_stat_cent1050; h4l_stat_cent1050 = (TGraphErrors*)infile->Get("h4l_stat_cent1050");
  TH1D*	h3l_coal_cent1050; h3l_coal_cent1050 = (TH1D*)infile->Get("h3l_coal_cent1050");
  TH1D*	h4l_coal_cent1050; h4l_coal_cent1050 = (TH1D*)infile->Get("h4l_coal_cent1050");

  double xxxxx  = 1600;
  double lrmargin = 250;

  TCanvas *x08 = new TCanvas("x08","x08",xxxxx,800);
  x08->cd()->SetRightMargin(0.05);
  x08->cd()->SetTopMargin(0.04);
  x08->cd()->SetLeftMargin(0.25);
  x08->cd()->SetBottomMargin(0.14);
  x08->cd();

  TPad *f1mass;
  f1mass = new TPad("f1mass", "f1mass",0,0.,0.5,1);
  f1mass->SetTicky(1);
  f1mass->SetLeftMargin(lrmargin/xxxxx);
  f1mass->SetRightMargin(0.00);

  f1mass->SetTopMargin(0.02);
  f1mass->SetBottomMargin(0.16);
  f1mass->SetFrameBorderMode(0);
  f1mass->SetFrameBorderMode(0);
  f1mass->Draw();
  f1mass->cd();

  x08->cd();

  f1mass->Draw();
  f1mass->cd();
  cent0010_dummy->Draw();
  h3l_syst_cent0010->Draw("p5same");
  h4l_syst_cent0010->Draw("p5same");
  h3l_stat_cent0010->Draw("psame");
  h4l_stat_cent0010->Draw("psame");
  h3l_coal_cent0010->Draw("chistsame");
  h4l_coal_cent0010->Draw("chistsame");

  TLegend *legjam_header;
  legjam_header = new TLegend(0.02, 0.26,0.08,0.32);
  legjam_header->SetTextSize(0.0550);
  legjam_header->SetMargin(0.25);
  legjam_header->SetFillColorAlpha(kWhite,0);
  legjam_header->SetTextFont(42);
  legjam_header->AddEntry(h4l_coal_cent0010,"Coalesc. (JAM)","");

  TLegend *legjam;
  legjam = new TLegend(0.02, 0.17,0.42,0.29);
  legjam->SetTextSize(0.0550);
  legjam->SetMargin(0.25);
  legjam->SetFillColorAlpha(kWhite,0);
  legjam->SetTextFont(42);
  legjam->AddEntry(h3l_coal_cent0010,"{}^{3}_{#Lambda}H","l");

  TLegend *legjam2;
  legjam2 = new TLegend(0.02+0.20, 0.17,0.42+0.20,0.29);
  legjam2->SetTextSize(0.0550);
  legjam2->SetMargin(0.25);
  legjam2->SetFillColorAlpha(kWhite,0);
  legjam2->SetTextFont(42);
  legjam2->AddEntry(h4l_coal_cent0010,"{}^{4}_{#Lambda}H","l");

  TPaveText *g_0_0 = new TPaveText(0.12, 0.30,0.13,0.30,"NDCNB");
  g_0_0->SetTextSize(50);
  g_0_0->SetTextFont(43);
  g_0_0->SetFillColorAlpha(kWhite,0);
  g_0_0->AddText("1");
  g_0_0->Draw();

  double  dyy = 0.136;
  TPaveText *g_0_1 = new TPaveText(0.12, 0.30+dyy,0.13,0.30+dyy,"NDCNB");
  g_0_1->SetTextSize(50);
  g_0_1->SetTextFont(43);
  g_0_1->SetFillColorAlpha(kWhite,0);
  g_0_1->AddText("2");
  g_0_1->Draw();

  TPaveText *g_0_2 = new TPaveText(0.12, 0.30+2*dyy,0.13,0.30+2*dyy,"NDCNB");
  g_0_2->SetTextSize(50);
  g_0_2->SetTextFont(43);
  g_0_2->SetFillColorAlpha(kWhite,0);
  g_0_2->AddText("3");
  g_0_2->Draw();

  TPaveText *g_0_3 = new TPaveText(0.12, 0.30+3*dyy,0.13,0.30+3*dyy,"NDCNB");
  g_0_3->SetTextSize(50);
  g_0_3->SetTextFont(43);
  g_0_3->SetFillColorAlpha(kWhite,0);
  g_0_3->AddText("4");
  g_0_3->Draw();

  TPaveText *g_0_4 = new TPaveText(0.12, 0.30+4*dyy,0.13,0.30+4*dyy,"NDCNB");
  g_0_4->SetTextSize(50);
  g_0_4->SetTextFont(43);
  g_0_4->SetFillColorAlpha(kWhite,0);
  g_0_4->AddText("5");
  g_0_4->Draw();


  TPaveText *t_cent_0 = new TPaveText(0.26, 0.75+0.1,0.44,0.80+0.1,"NDCNB");
  t_cent_0->SetTextFont(63);
  t_cent_0->SetTextSize(55);
  t_cent_0->SetFillColorAlpha(kWhite,0);
  t_cent_0->AddText("(a) 0-10%");
  t_cent_0->Draw();

  TLegend *legcent_header;

  legcent_header = new TLegend(0.22-0.02, 0.64+0.1,0.32-0.02,0.70+0.1);
  legcent_header->SetTextSize(0.0550);
  legcent_header->SetMargin(0.1);
  legcent_header->SetFillColorAlpha(kWhite,0);
  legcent_header->SetTextFont(42);
  legcent_header->AddEntry(h3l_stat_cent0010,"Au+Au 3GeV","");

  TLegend *legcent;
  legcent = new TLegend(0.23-0.02, 0.54+0.10,0.53-0.02,0.66+0.10);
  legcent->SetTextSize(0.0550);
  legcent->SetMargin(0.1);
  legcent->SetFillColorAlpha(kWhite,0);
  legcent->SetTextFont(42);
  legcent->AddEntry(h3l_stat_cent0010,"{}^{3}_{#Lambda}H","p");

  TLegend *legcentc;
  legcentc = new TLegend(0.23+0.14-0.02, 0.54+0.10,0.53+0.14-0.02,0.66+0.10);
  legcentc->SetTextSize(0.0550);
  legcentc->SetMargin(0.1);
  legcentc->SetFillColorAlpha(kWhite,0);
  legcentc->SetTextFont(42);
  legcentc->AddEntry(h4l_stat_cent0010,"{}^{4}_{#Lambda}H","p");

  legcent_header->Draw();
  legcent->Draw();
  legcentc->Draw();

  x08->cd();

  TPad *f2mass;
  f2mass = new TPad("f2mass", "f2mass",0.5,0,1,1);
  f2mass->SetTicky(1);

  f2mass->SetLeftMargin(0.);
  f2mass->SetRightMargin(lrmargin/xxxxx);
  f2mass->SetTopMargin(0.02);
  f2mass->SetBottomMargin(0.16);
  f2mass->SetFrameBorderMode(0);
  f2mass->SetFrameBorderMode(0);
  f2mass->Draw();
  f2mass->cd();

  f2mass->Draw();
  f2mass->cd();
  cent1050_dummy->Draw();
  h3l_syst_cent1050->Draw("p5same");
  h4l_syst_cent1050->Draw("p5same");
  h3l_stat_cent1050->Draw("psame");
  h4l_stat_cent1050->Draw("psame");
  h3l_coal_cent1050->Draw("chistsame");
  h4l_coal_cent1050->Draw("chistsame");

  legjam_header->Draw();
  legjam->Draw();
  legjam2->Draw();

  TPaveText *t_cent_1 = new TPaveText(0.26+0.03, 0.75+0.1,0.90+0.03,0.80+0.1,"NDCNB");
  t_cent_1->SetTextSize(55);
  t_cent_1->SetTextFont(63);
  t_cent_1->SetFillColorAlpha(kWhite,0);
  t_cent_1->AddText("(b) 10-50%");
  t_cent_1->Draw();


  TPaveText *g_1_0 = new TPaveText(0.90, 0.30,0.91,0.30,"NDCNB");
  g_1_0->SetTextSize(50);
  g_1_0->SetTextFont(43);
  g_1_0->SetFillColorAlpha(kWhite,0);
  g_1_0->AddText("0.5");
  g_1_0->Draw();

  // dyy = 0.145;
  dyy = 0.136;

  TPaveText *g_1_1 = new TPaveText(0.90, 0.30+dyy,0.91,0.30+dyy,"NDCNB");
  g_1_1->SetTextSize(50);
  g_1_1->SetTextFont(43);
  g_1_1->SetFillColorAlpha(kWhite,0);
  g_1_1->AddText("1.0");
  g_1_1->Draw();

  TPaveText *g_1_2 = new TPaveText(0.90, 0.30+2*dyy,0.91,0.30+2*dyy,"NDCNB");
  g_1_2->SetTextSize(50);
  g_1_2->SetTextFont(43);
  g_1_2->SetFillColorAlpha(kWhite,0);
  g_1_2->AddText("1.5");
  g_1_2->Draw();


  TPaveText *g_1_3 = new TPaveText(0.90, 0.30+3*dyy,0.91,0.30+3*dyy,"NDCNB");
  g_1_3->SetTextSize(50);
  g_1_3->SetTextFont(43);
  g_1_3->SetFillColorAlpha(kWhite,0);
  g_1_3->AddText("2.0");
  g_1_3->Draw();

  TPaveText *g_1_4 = new TPaveText(0.90, 0.30+4*dyy,0.91,0.30+4*dyy,"NDCNB");
  g_1_4->SetTextSize(50);
  g_1_4->SetTextFont(43);
  g_1_4->SetFillColorAlpha(kWhite,0);
  g_1_4->AddText("2.5");
  g_1_4->Draw();

  x08->cd();
  TPad *f3mass;
  f3mass = new TPad("f3mass", "f3mass",0.45,0,0.55,0.08);


  f3mass->SetLeftMargin(0.);
  f3mass->SetRightMargin(0.);
  f3mass->SetTopMargin(0.0);
  f3mass->SetBottomMargin(0.);
  f3mass->SetFrameBorderMode(0);
  f3mass->SetFrameBorderMode(0);
  f3mass->Draw();
  f3mass->cd();


  x08->cd();
  TPaveText *t_2 = new TPaveText(0.26, .0,.75,.08,"NDCNB");
  t_2->SetTextSize(55);
  t_2->SetTextFont(43);
  t_2->SetFillColorAlpha(kWhite,0);
  t_2->AddText("Rapidity y");
  t_2->Draw();


  TPaveText *t_auaa = new TPaveText(0.26, 0.22,0.55,0.33,"NB");
  t_auaa->SetTextSize(55);
  t_auaa->SetTextFont(63);
  t_auaa->SetFillColorAlpha(kWhite,0);
  t_auaa->SetLineColor(kBlack);
  t_auaa->AddText("STAR");
  t_auaa->Draw();


  x08->SaveAs("fig/h3l_dNdy_3GeV.pdf");
  x08->SaveAs("fig/h3l_dNdy_3GeV.png");





}
