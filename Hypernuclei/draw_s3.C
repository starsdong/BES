#include "style.C+"
#include "draw.C+"

void draw_s3(){
  style();

  TCanvas *x05 = new TCanvas("x05","x05",800,800);
  x05->cd()->SetRightMargin(0.05);
  x05->cd()->SetTopMargin(0.02);
  x05->cd()->SetLeftMargin(0.16);
  x05->cd()->SetBottomMargin(0.16);
  x05->cd()->SetLogx();

  TGraphErrors *s3hybridurqmd;
  s3hybridurqmd = new TGraphErrors("xin_s3/s3hybridurqmd.txt","%lg %lg");
  TGraphErrors *s3stringmeltingampt;
  s3stringmeltingampt = new TGraphErrors("xin_s3/s3stringmeltingampt.txt","%lg %lg");
  TGraphErrors *s3defaultampt;
  s3defaultampt = new TGraphErrors("xin_s3/s3defaultampt.txt","%lg %lg");
  TGraphErrors *s3coalescencedcm;
  s3coalescencedcm = new TGraphErrors("xin_s3/s3coalescencedcm.txt","%lg %lg");
  TGraphErrors *s3thermalmodelgsi;
  s3thermalmodelgsi = new TGraphErrors("xin_s3/s3thermalmodelgsi.txt","%lg %lg");


  TH1D *s3dummy;
  s3dummy = new TH1D("s3dummy","",1000,1.5,6600);
  s3dummy->GetYaxis()->SetRangeUser(0,2.1);
  s3dummy->GetYaxis()->SetLabelFont(42);
  s3dummy->GetYaxis()->SetTitleFont(42);
  s3dummy->GetXaxis()->SetLabelFont(42);
  s3dummy->GetXaxis()->SetTitleFont(42);
  s3dummy->GetYaxis()->SetTitle("S_{3}");
  s3dummy->GetXaxis()->SetTitle("#sqrt{s_{NN}} [GeV]");
  s3dummy->GetXaxis()->SetTitleOffset(1.2);
  s3dummy->GetYaxis()->SetTitleOffset(1.1);
  s3dummy->Draw();

  s3hybridurqmd->SetLineColor(kGray+2);
  s3hybridurqmd->SetLineStyle(6);
  s3hybridurqmd->SetLineWidth(2);
  s3hybridurqmd->Draw("same");

  s3stringmeltingampt->SetLineColor(kOrange+7);
  s3stringmeltingampt->SetLineStyle(1);
  s3stringmeltingampt->SetLineWidth(2);
  s3stringmeltingampt->Draw("same");

  s3defaultampt->SetLineColor(kAzure+7);
  s3defaultampt->SetLineStyle(1);
  s3defaultampt->SetLineWidth(2);
  s3defaultampt->Draw("same");

  s3coalescencedcm->SetLineColor(kMagenta);
  s3coalescencedcm->SetLineStyle(5);
  s3coalescencedcm->SetLineWidth(2);
  s3coalescencedcm->Draw("same");

  s3thermalmodelgsi->SetLineColor(kBlue);
  s3thermalmodelgsi->SetLineStyle(2);
  s3thermalmodelgsi->SetLineWidth(2);
  s3thermalmodelgsi->Draw("same");

  double x_err_margin = 0.2;
  double alice_e = 2760;
  double alice_s3 = 0.60;
  double alice_s3_stat = 0.13;
  double alice_s3_syst = 0.21;

  TGraphErrors *s3alice;s3alice = new TGraphErrors();
  s3alice->SetPoint(0,alice_e,alice_s3);
  s3alice->SetPointError(0,0,alice_s3_stat);

  TGraphErrors *s3alice_sys;s3alice_sys = new TGraphErrors();
  s3alice_sys->SetPoint(0,alice_e,alice_s3);
  s3alice_sys->SetPointError(0,alice_e*x_err_margin,alice_s3_syst);

  s3alice->SetMarkerStyle(33);
  s3alice->SetMarkerSize(3);
  s3alice->SetMarkerColor(kBlack);
  s3alice->SetLineColor(kBlack);

  s3alice_sys->SetMarkerStyle(33);
  s3alice_sys->SetMarkerColor(kBlack);
  s3alice_sys->SetMarkerSize(3);
  s3alice_sys->SetLineColor(kBlack);
  s3alice_sys->SetFillColorAlpha(kGray+1,0.1);

  s3alice_sys->Draw("p5same");
  s3alice->Draw("psame");

  TGraphErrors *s3ags; s3ags = new TGraphErrors();
  double ags_s3 = 0.36;
  double ags_s3_stat = 0.26;
  s3ags->SetPoint(0,5,ags_s3);
  s3ags->SetPointError(0,0,ags_s3_stat);
  s3ags->SetMarkerStyle(21);
  s3ags->SetMarkerSize(2);
  s3ags->SetMarkerColor(kBlack);
  s3ags->SetLineColor(kBlack);
  s3ags->Draw("psame");

  TGraphErrors *s3star; s3star = new TGraphErrors();
  TGraphErrors *s3star_sys; s3star_sys = new TGraphErrors();
  // s3star->SetPoint(0,200,1.07895);
  // s3star->SetPoint(0,200,0.90);//star+phenix, refer to: http://cpc.ihep.ac.cn/fileZGWLC/journal/article/zgwlc/2020/11/PDF/CPC-2020-0211.pdf
  // s3star->SetPointError(0,0,0.218051);

  double h3l_he3_200 = 0.82;
  double h3l_he3_200_stat = 0.16;
  double h3l_he3_200_sys = 0.12;
  double lambda_200 = 16.7;
  double proton_200 = 18.4;
  double lambda_200_stat = 0.2;
  double lambda_200_sys = 1.1;
  double proton_200_sys = 2.6;
  double s3_200 = h3l_he3_200/lambda_200*proton_200;
  double s3_200_stat = s3_200*sqrt(h3l_he3_200_stat*h3l_he3_200_stat/h3l_he3_200/h3l_he3_200 + lambda_200_stat*lambda_200_stat/lambda_200/lambda_200);
  double s3_200_sys = s3_200*sqrt(h3l_he3_200_sys*h3l_he3_200_sys/h3l_he3_200/h3l_he3_200 + lambda_200_sys*lambda_200_sys/lambda_200/lambda_200 + proton_200_sys*proton_200_sys/proton_200/proton_200);
  double star_e = 200;

  s3star->SetPoint(0, star_e, s3_200   );
  s3star->SetPointError(0,0,  s3_200_stat  );
  s3star_sys->SetPoint(0,star_e, s3_200   );
  s3star_sys->SetPointError(0, star_e*x_err_margin,  s3_200_sys  );

  s3star->SetMarkerStyle(29);
  s3star->SetMarkerSize(4);
  s3star->SetMarkerColor(kBlack);
  s3star->SetLineColor(kBlack);

  s3star_sys->SetMarkerStyle(29);
  s3star_sys->SetMarkerSize(3);
  s3star_sys->SetMarkerColor(kBlack);
  s3star_sys->SetLineColor(kBlack);
  s3star_sys->SetFillColorAlpha(kGray+1,0.1);

  s3star->Draw("psame");
  s3star_sys->Draw("p5same");

  // double this_s3 = this_h3l/avg_h3_cent/this_lambda*lihui_proton;
  // double this_s3_00 = this_s3*sqrt(this_h3l_stat*this_h3l_stat/this_h3l/this_h3l + avg_h3_cent_stat*avg_h3_cent_stat/avg_h3_cent/avg_h3_cent + this_lambda_stat*this_lambda_stat/this_lambda/this_lambda + lihui_proton_stat*lihui_proton_stat/lihui_proton/lihui_proton) ;
  // double this_s3_01 = this_s3*sqrt(this_h3l_sys *this_h3l_sys /this_h3l/this_h3l + avg_h3_cent_sys *avg_h3_cent_sys /avg_h3_cent/avg_h3_cent + this_lambda_sys *this_lambda_sys /this_lambda/this_lambda + lihui_proton_sys *lihui_proton_sys /lihui_proton/lihui_proton) ;

  double this_e = 3.;
  double this_s3 = 0.230162;
  double this_s3_00 = 0.0325788;
  double this_s3_01 = 0.0743901;

  double this_27_e = 27.;
  double this_27_s3 = 0.559906;
  double this_27_s3_00 = 0.0962171;
  double this_27_s3_01 = 0.212674;

  TGraphErrors *this_s3_stat; this_s3_stat = new TGraphErrors();
  TGraphErrors *this_s3_sys; this_s3_sys = new TGraphErrors();

  this_s3_stat->SetPoint(0,this_e, this_s3);
  this_s3_stat->SetPointError(0,0,this_s3_00 );
  this_s3_sys->SetPoint(0,this_e, this_s3);
  this_s3_sys->SetPointError(0,this_e*x_err_margin, this_s3_01);

  this_s3_stat->SetPoint(1,this_27_e, this_27_s3);
  this_s3_stat->SetPointError(1,0,this_27_s3_00 );
  this_s3_sys->SetPoint(1,this_27_e, this_27_s3);
  this_s3_sys->SetPointError(1,this_27_e*x_err_margin, this_27_s3_01);

  this_s3_stat->SetMarkerColor(kBlack);
  this_s3_sys->SetMarkerColor(kBlack);
  this_s3_stat->SetLineColor(kBlack);
  this_s3_sys->SetLineColor(kBlack);

  this_s3_stat->SetMarkerStyle(20);
  this_s3_stat->SetMarkerSize(2);
  this_s3_sys->SetFillColorAlpha(kGray+1,0.1);

  this_s3_stat->Draw("psame");
  this_s3_sys->Draw("p5same");

  TLegend *legs3;
   legs3 = new TLegend(0.19, 0.67,0.58,0.97);
   legs3->SetTextSize(0.0300);
   legs3->SetMargin(0.12);
   legs3->SetFillColorAlpha(kWhite,0);
   legs3->SetTextFont(42);
   legs3->AddEntry(this_s3_stat,"Model","");
   legs3->AddEntry(s3defaultampt,"Default AMPT + Coal.","l");
   legs3->AddEntry(s3stringmeltingampt,"String Melting AMPT + Coal.","l");
   legs3->AddEntry(s3coalescencedcm,"Coalescence (DCM model)","l");
   legs3->AddEntry(s3thermalmodelgsi,"Thermal Model (GSI - Hei.)","l");
   legs3->AddEntry(s3hybridurqmd,"Hybrid URQMD","l");

   TLegend *legs3d;
    legs3d = new TLegend(0.62, 0.72,0.92,0.97);
    legs3d->SetTextSize(0.0300);
    legs3d->SetMargin(0.2);
    legs3d->SetFillColorAlpha(kWhite,0);
    legs3d->SetTextFont(42);
    legs3d->AddEntry(this_s3_stat,"Data","");
    legs3d->AddEntry(this_s3_stat,"Au+Au 0-10%","p");
    legs3d->AddEntry(s3ags,"E864 Au+Pt 0-10%","p");
    legs3d->AddEntry(s3star,"STAR Au+Au 0-80%","p");
    legs3d->AddEntry(s3alice,"ALICE Pb+Pb 0-10%","p");
    legs3->Draw();
   legs3d->Draw();


   TPaveText *t_cent2 = new TPaveText(0.48, 0.18,0.96,0.24,"NDCNB");
   t_cent2->SetTextSize(0.025);
   t_cent2->SetTextFont(52);
   t_cent2->SetFillColorAlpha(kWhite,0);
   t_cent2->AddText("Assuming B.R.({}^{3}_{#Lambda}H#rightarrow{}^{3}He + #pi^{-}) = 25%");
   t_cent2->Draw("same");



   double alice_e_ppb = 5020;
   double alice_s3_ppb = 0.44;
   double alice_s3_ppb_stat = 0.12;
   double alice_s3_ppb_syst = 0.10;

   TGraphErrors *s3alice_ppb;s3alice_ppb = new TGraphErrors();
   s3alice_ppb->SetPoint(0,alice_e_ppb,alice_s3_ppb);
   s3alice_ppb->SetPointError(0,0,alice_s3_ppb_stat);

   TGraphErrors *s3alice_ppb_sys; s3alice_ppb_sys = new TGraphErrors();
   s3alice_ppb_sys->SetPoint(0,alice_e_ppb,alice_s3_ppb);
   s3alice_ppb_sys->SetPointError(0,alice_e_ppb*x_err_margin,alice_s3_ppb_syst);

   s3alice_ppb->SetMarkerStyle(27);
   s3alice_ppb->SetMarkerSize(3);
   s3alice_ppb->SetMarkerColor(kRed);
   s3alice_ppb->SetLineColor(kRed);

   s3alice_ppb_sys->SetMarkerStyle(27);
   s3alice_ppb_sys->SetMarkerColor(kRed);
   s3alice_ppb_sys->SetMarkerSize(3);
   s3alice_ppb_sys->SetLineColor(kRed);
   s3alice_ppb_sys->SetFillColorAlpha(kRed,0.1);

   s3alice_ppb_sys->Draw("p5same");
   s3alice_ppb->Draw("psame");


   double hyphi_e = 2.7;
   double hyphi_s3 = 0.28;
   double hyphi_s3_stat = 0.14;


   TGraphErrors *s3hyphi;s3hyphi = new TGraphErrors();
   s3hyphi->SetPoint(0,hyphi_e,hyphi_s3);
   s3hyphi->SetPointError(0,0,hyphi_s3_stat);

   s3hyphi->SetMarkerStyle(28);
   s3hyphi->SetMarkerSize(3);
   s3hyphi->SetMarkerColor(kRed);
   s3hyphi->SetLineColor(kRed);

   s3hyphi->Draw("psame");
   this_s3_stat->Draw("psame");
   this_s3_sys->Draw("p5same");

   TLegend *legs4d;
    legs4d = new TLegend(0.62, 0.62,0.92,0.72);
    legs4d->SetTextSize(0.0300);
    legs4d->SetMargin(0.2);
    legs4d->SetFillColorAlpha(kWhite,0);
    legs4d->SetTextFont(42);
    legs4d->AddEntry(s3alice_ppb,"ALICE p+Pb 0-40%","p");
    legs4d->AddEntry(s3hyphi,"HypHI Li+C ","p");

   legs4d->Draw();


   x05->Update();
   x05->SaveAs("fig/s3.pdf");
   x05->SaveAs("fig/s3.png");

   cout << " == HypHI == " << endl;
   s3hyphi->Print();
   cout << " == AGS == " << endl;
   s3ags->Print();
   cout << " == STAR == " << endl;
   this_s3_stat->Print();
   this_s3_sys->Print();   
   s3star->Print();
   s3star_sys->Print();
   cout << " == ALICE == " << endl;
   s3alice->Print();
   s3alice_sys->Print();
   s3alice_ppb->Print();
   s3alice_ppb_sys->Print();

   s3hyphi->SetName("HypHI");
   s3ags->SetName("E864_AuPt_0_10");
   this_s3_stat->SetName("STAR_AuAu_0_10");
   this_s3_sys->SetName("STAR_AuAu_0_10_sys");
   s3star->SetName("STAR_AuAu_0_80");
   s3star_sys->SetName("STAR_AuAu_0_80_sys");
   s3alice->SetName("ALICE_PbPb_0_10");
   s3alice_sys->SetName("ALICE_PbPb_0_10_sys");  
   s3alice_ppb->SetName("ALICE_pPb_0_40");
   s3alice_ppb_sys->SetName("ALICE_pPb_0_40_sys");

   TFile *fout = new TFile("s3_data.root","recreate");
   s3hyphi->Write();
   s3ags->Write();
   this_s3_stat->Write();
   this_s3_sys->Write();
   s3star->Write();
   s3star_sys->Write();
   s3alice->Write();
   s3alice_sys->Write();
   s3alice_ppb->Write();
   s3alice_ppb_sys->Write();
   fout->Close();
   
}
