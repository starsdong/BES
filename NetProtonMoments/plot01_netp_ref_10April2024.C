//*--=======================================================
//*
//*       plot C4/C2 data different from model references
//*       and 70-80% perpherial data
//*       data from, Nov 1, 2020 // BES-II data plot against reference 
//*
//*--=======================================================
#include "style.C+"
#include "draw.C+"

void plot01_netp_ref_10April2024()
{
  style();

  //  STAR Data
  const Int_t NE = 9;
  const Double_t EE[NE] = {7.7, 11.5, 14.5, 19.6, 27, 39, 54.4, 62.4, 200};

  // C4/C2 diviations
  Double_t dd1[NE] = {0.62712, -0.466486, 0.996377, -2.36179, -3.06788, -1.30089, -2.45845, -0.742427, -0.396131};// HRG GCE, Skellman
  Double_t dd2[NE] = {0.696718, -0.213861, 1.14744, -1.7554, -2.32734, -0.264671, -1.39065, 0.0104701, 0.461737};// perpherial 70-80% data
  Double_t dd3[NE] = {1.07088, 0.11839, 1.73505, -1.42595, -2.62184, -0.58628, -0.76172, -0.124683, -0.0512915}; // UrQMD 5% central  
 
  // create TGraphErrors for data points
  TGraph *gr1;   // data with statistical errors
  TGraph *gr2; TGraph *gr22;   // data with statistical errors
  TGraph *gr3;   // data with statistical errors
  
  // Plotting  --============================================================
  TCanvas *c1 = new TCanvas("c1", "c1",0,0,1200,800);
  c1->Draw();
  
  TPad *p1 = new TPad("p1","",0.07,0.001,1.0,1.0,0,0);
  p1->SetLogx(); //p1->SetGrid();
  p1->SetLeftMargin(0.1);  p1->SetBottomMargin(0.16);  p1->SetTopMargin(0.005);  p1->SetRightMargin(0.1);  
  p1->Draw(); p1->cd();

  double x1 = 1.75;  double x2 = 275;
  double y1 = -3.75;  double y2 = 3.75;
  
  TH1D *d1 = new TH1D("d1","",1,x1,x2);
  d1->SetMinimum(y1);  d1->SetMaximum(y2);
  d1->GetXaxis()->SetNdivisions(408);  d1->GetXaxis()->SetTitle("");  d1->GetXaxis()->SetLabelOffset(100); d1->GetXaxis()->SetLabelSize(0.06);
  d1->GetYaxis()->SetNdivisions(210);  d1->GetYaxis()->SetTitle("");  d1->GetYaxis()->SetLabelOffset(100); d1->GetYaxis()->SetLabelSize(0.065);
  d1->Draw("c");
  
  //drawHistBox(x1,x2,y1,y2,2222);
  drawLine(x1,0,x2,0,1,9,1);

  TLatex*latex = new TLatex();
  latex->SetTextFont(42); latex->SetTextSize(0.055);  latex->SetTextAlign(22);

  //  latex->DrawLatex(2,y3,"2");
  double y3 = -4.05;
  latex->DrawLatex(2,y3,"2");  latex->DrawLatex(5,y3,"5");  latex->DrawLatex(10,y3,"10");  latex->DrawLatex(20,y3,"20");
  latex->DrawLatex(50,y3,"50");  latex->DrawLatex(100,y3,"100");  latex->DrawLatex(200,y3,"200");

  double x3=1.625;
  latex->SetTextAlign(32);
  latex->DrawLatex(x3,0.0,"0");  latex->DrawLatex(x3,1,"1");  latex->DrawLatex(x3,2,"2"); latex->DrawLatex(x3,3,"3");
  latex->DrawLatex(x3,-1,"-1");  latex->DrawLatex(x3,-2,"-2");  latex->DrawLatex(x3,-3,"-3"); 
  //  drawText(x3,1.0,"1.0",42,0.06,0);

  latex->SetTextAlign(20);  latex->SetTextSize(0.075); 
  latex->DrawLatex(22, -5.15, "Collision Energy #sqrt{s_{NN}} (GeV)");  

  latex->SetTextAngle(90);
  latex->DrawLatex(1.2, 0.0, "Deviations");
  drawText(285,3,"#copyrightCbR0604: Diviations STAR papers/neet-p/plot2_netp_model_01Nov2020.C",42,0.015,-90);

  drawText(3.25,-1.1,"Reference",42,0.035,0);
  drawLine(2.25,-1.65,3.75,-1.65,8,8,2);   //Draw_TGAE_Point_new_Symbol(2.95,-1.65, 0., 0., 0., 0., 20, kRed,2);
  drawText(4.0,-1.75,"HRG (GCE)",42,0.045,0);
  
//  Draw_TGAE_Point_new_Symbol(2.95,-2.35, 0., 0., 0., 0., 26,1,2);
  drawText(4.0,-2.475,"UrQMD (5%)",42,0.045,0);

//  Draw_TGAE_Point_new_Symbol(2.95,-3, 0., 0., 0., 0.,21,4,2);
  drawText(4.0,-3.1,"Data (70-80%)",42,0.045,0);


  drawText(30,3.0,"#kappa#sigma^{2}(data 5%) - Ref.",42,0.045,0);
  drawLine(28.0,2.875,155.0,2.875,1,1,1); 
  drawText(60,2.6,"#sigma_{total} ",42,0.045,0);
  drawText(10.5,2.775,"Deviation = ",42,0.045,0);
   
  // GCE C4/C2 --=========================================================================
  gr1 = new TGraph(NE, EE, dd1);
  gr1->SetMarkerStyle(20);   gr1->SetMarkerSize(2.25);  gr1->SetMarkerColor(2);   gr1->SetLineWidth(2); gr1->SetLineStyle(8);  gr1->SetLineColor(2);

  gr1->Draw("pl");
  // end of GCE

   // UrQMD C4/C2 --=========================================================================
  gr3 = new TGraph(NE, EE, dd3);
  gr3->SetMarkerStyle(26);   gr3->SetMarkerSize(3.0);  gr3->SetMarkerColor(1);   gr3->SetLineWidth(2); gr3->SetLineStyle(8);
  gr3->Draw("p");
  // end of UrQMD

   // data 70-80% C4/C2 --=========================================================================
  gr2 = new TGraph(NE, EE, dd2);
  gr2->SetMarkerStyle(25);   gr2->SetMarkerSize(3.0);  gr2->SetMarkerColor(4);   gr2->SetLineWidth(1); gr2->SetLineStyle(8);
  gr2->Draw("p");

  gr22 = new TGraph(NE, EE, dd2);
  gr22->SetMarkerSize(2.25); gr22->SetMarkerStyle(21);  gr22->SetMarkerColor(4);  
  gr22->Draw("psame");
  // end of data
  /*
  drawText(70,1.01,"(1)",42,0.085,0);
  drawText(169,1.01,"#sigma",42,0.085,0);
  latex->SetTextAlign(10); latex->SetTextSize(0.0785);latex->SetTextFont(32);
  latex->DrawLatex(135,1.01,"S");
 
  latex->SetTextAlign(10); latex->SetTextSize(0.055);
  latex->SetTextColor(41); latex->SetTextFont(22);
  latex->DrawLatex(2.55,0.15,"Au+Au Collisions");
  latex->SetTextColor(1);
  latex->DrawLatex(2.6,0.155,"Au+Au Collisions");
  latex->SetTextFont(132);  latex->SetTextSize(0.045);
  latex->DrawLatex(3.2,0.08,"Net-proton");
  latex->DrawLatex(3.2,0.005,"|y| < 0.5,  0.4 < p_{T}< 2.0 (GeV/c)");

  latex->SetTextColor(2);latex->SetTextSize(0.07);
  //latex->DrawLatex(2.5,0.175,"Au+Au Collisions");
  latex->SetTextColor(1);
  
  // legends --=================================================
  double y110=0.51;
  double y1101=0.49;
  drawLine(2.2,y110,3.4,y110,10,1,41);
  drawText(3.65,y1101,"UrQMD 0-5%",42,0.045,0);

  double y120=0.43;
  double y1201=0.41;
  drawLine(2.2,y120,3.4,y120,10,1,1);
  drawText(3.65,y1201,"HRG GCE",42,0.045,0);

  double y130=0.35;
  double y1301=0.33;
  drawLine(2.2,y130,3.4,y130,5,9,46);
  drawText(3.65,y1301,"HRG CE",42,0.045,0);
  */

  p1->Modified();  c1->Update();  c1->cd();
  //--=============================================================================

  // write out pdf file  --======================================================
  c1->SaveAs("plot01_netp_ref_10April2024.pdf");
}
