#include "draw.C+"
#include "style.C+"


void plotRatioOnePaper()
{
  style();

  //////////=====
  const int ntot_datapts=11;
  const int n_model = 4;
  const int plotflag[n_model] = {1, 1, 1, 1}; // UrQMD, HRG CE, Hydro EV, 70-80%

  ////data 0-5%
  double C42_ener[ntot_datapts]={7.7,9.2,11.5,14.6,17.3,19.6,27,39,54.4,62.4,200};
  double C42_val[ntot_datapts]={0.411909,0.540458,0.401178,0.430751,0.348682,0.339938,0.602046,0.739693,0.69715,0.792955,0.900669};
  double C42_stat[ntot_datapts]={0.246897,0.187029,0.134746,0.104085,0.119235,0.0728664,0.0950052,0.147006,0.053074,0.250823,0.208458};
  double C42_sys[ntot_datapts]={0.130209,0.068346,0.0647203,0.0403519,0.0831337,0.0395323,0.0304235,0.135754,0.0482269,0.1219,0.139363};
  double C42_tot[ntot_datapts]={0};
  
  ////UrQMD
  double cen_UQMD[ntot_datapts]={7.7,9.2,11.5,14.6,17.3,19.6,27,39,54.4,62.4,200};
  double urqmd_C42_05_data[ntot_datapts]={0.4753277, 0.4625844, 0.5253311, 0.594994, 0.6782202, 0.7327199, 0.903167,0.860031,0.83973130,0.829186,0.914062};
  double urqmd_C42_05_sterr[ntot_datapts]={0.06445722,0.06283556,0.06766454,0.0707438,0.05287615,0.05808925,0.063723,0.045722,0.055722,0.081655,0.072831};
  
  ////HRG CE
  double anar_cen[ntot_datapts]={7.7,9.2,11.5,14.6,17.3,19.6,27,39,54.4,62.4,200};
  double anar_C42[ntot_datapts]={0.40214780,0.44213425,0.5034468,0.58341961,0.63475107,0.67691625,0.73933998,0.76389,0.7953959,0.8117625,0.84005};
  
  ////Hydro
  double Hydro_cen[ntot_datapts]={7.7,9.2,11.5,14.6,17.3,19.6,27,39,54.4,62.4,200};
  double Hydro_C42[ntot_datapts]={0.34623,0.39264044,0.46380312,0.556624,0.58357922,0.605721,0.657514,0.758193,0.78623350,0.8008,0.838933};

  ////Hydro+EV
  double Hydro_cen_wtEV[ntot_datapts]={7.7,9.2,11.5,14.6,17.3,19.6,27,39,54.4,62.4,200};
  double Hydro_C42_wtEV[ntot_datapts]={0.239113,0.28855366,0.36436268,0.463244,0.49676714,0.524304,0.579964,0.693609,0.72515728,0.741546,0.785968};  

  ////70-80%
  double C42_val_cen80[ntot_datapts]={0.791157,0.724508,0.712454,0.731724,0.725412,0.733293,0.809095,0.79274,0.841178,0.79003,0.78476};
  double C42_stat_cen80[ntot_datapts]={0.00905602,0.00941409,0.00439071,0.00362469,0.00427153,0.00262129,0.0032776,0.0051496,0.00161249,0.0079563,0.0038362};
  double C42_sys_cen80[ntot_datapts]={0.00662053,0.00361412,0.00656009,0.00580586,0.00612719,0.00566322,0.0079558,0.010217,0.0106952,0.014472,0.011062};
  double C42_tot_cen80[ntot_datapts]={0};
  
  for(int y=0;y<ntot_datapts;y++)
    {
      C42_tot[y]=sqrt(C42_stat[y]*C42_stat[y]+C42_sys[y]*C42_sys[y]);
      cout << C42_tot[y] << endl;
      C42_tot_cen80[y]=sqrt(C42_stat_cen80[y]*C42_stat_cen80[y]+C42_sys_cen80[y]*C42_sys_cen80[y]);
      cout << C42_tot[y] << "\t" << C42_tot_cen80[y] << endl;
    }

  ////////////////////////////////////
  // rename arrays for simplicity
  ////////////////////////////////////
  double sNN[ntot_datapts];
  double C42_data[ntot_datapts];
  double C42_stat_data[ntot_datapts];
  double C42_sys_data[ntot_datapts];
  double C42_err_data[ntot_datapts];
  double C42_model[n_model][ntot_datapts];
  double C42_stat_model[n_model][ntot_datapts];
  double C42_sys_model[n_model][ntot_datapts];
  double C42_err_model[n_model][ntot_datapts];
  const Char_t *NameModel[n_model] = {"UrQMD (0-5%)", "HRG CE","Hydro EV","Data (70-80%)"};
  for(int i=0;i<ntot_datapts;i++) {
    sNN[i] = C42_ener[i];
    C42_data[i] = C42_val[i];
    C42_stat_data[i] = C42_stat[i];
    C42_sys_data[i] = C42_sys[i];
    C42_err_data[i] = C42_tot[i];
    
    for(int j=0;j<n_model;j++) {
      if(j==0) {
	C42_model[j][i] = urqmd_C42_05_data[i];
	C42_stat_model[j][i] = urqmd_C42_05_sterr[i];
	C42_sys_model[j][i] = 0.0;
	C42_err_model[j][i] = urqmd_C42_05_sterr[i];
      } else if(j==1) {
	C42_model[j][i] = anar_C42[i];
	C42_stat_model[j][i] = 0.0;
	C42_sys_model[j][i] = 0.0;
	C42_err_model[j][i] = 0.0;
      } else if(j==2) {
	C42_model[j][i] = Hydro_C42_wtEV[i];
	C42_stat_model[j][i] = 0.0;
	C42_sys_model[j][i] = 0.0;
	C42_err_model[j][i] = 0.0;
      } else {
	C42_model[j][i] = C42_val_cen80[i];
	C42_stat_model[j][i] = C42_stat_cen80[i];
	C42_sys_model[j][i] = C42_sys_cen80[i];
	C42_err_model[j][i] = C42_tot_cen80[i];
      }
    }
  }

  /////////////////
  // 3 GeV
  ////////////////
  const int ntot_fxt = 1;
  double C42_ener_fxt[ntot_fxt] = {3.0};
  double C42_fxt[ntot_fxt] = {-0.84574};
  double C42_stat_fxt[ntot_fxt] = {0.0862692};
  double C42_sys_fxt[ntot_fxt] = {0.821897};
  
  double C42_model_fxt[n_model][ntot_fxt] = {{-0.774123}, {0}, {0}, {0}};
  double C42_err_model_fxt[n_model][ntot_fxt] = {{0.0907071}, {0}, {0}, {0}};
  
  ///////////////////////////////////
  //////////===== chi2 test /////////
  ///////////////////////////////////
  
  const int n_test = 7;
  double chi2[n_model] = {0.0, 0.0, 0.0, 0.0};

  for(int im =0; im<n_model; im++) {
    cout << " Model\t" << NameModel[im] << endl;
    for(int i=0;i<n_test;i++) {
      double res = ( C42_data[i] - C42_model[im][i] )/sqrt(TMath::Power(C42_err_data[i], 2.));//+TMath::Power(C42_err_model[im][i], 2.));
      chi2[im] += res*res;
      //TMath::Power( C42_data[i]-C42_model[im][i], 2.0)/(TMath::Power(C42_err_data[i], 2.)+TMath::Power(C42_err_model[im][i], 2.));
      cout << "\t\t energy = " << sNN[i] << "\t residual = " << res << endl;
    }
    cout << " Model\t" << NameModel[im] << "\t" << setw(12) << chi2[im] << "\t p-value = " << TMath::Prob(chi2[im], n_test) << endl;
  }
  
  
  ////////////////////////////////////////////
  //////////===== likihood null test /////////
  ////////////////////////////////////////////


  //////////// making ratio plot
  const Double_t sc[n_model] = {-0.06, -0.03, 0., 0.03};
  double xe[n_model][ntot_datapts], r[n_model][ntot_datapts], re[n_model][ntot_datapts], res[n_model][ntot_datapts];
  double ener[n_model][ntot_datapts+ntot_fxt], r_unity[n_model][ntot_datapts+ntot_fxt], resm[n_model][ntot_datapts+ntot_fxt];
  for(int im=0;im<n_model;im++) {
    for(int i=0;i<ntot_datapts;i++) {
      xe[im][i] = C42_ener[i] * (1 + sc[im]);
      r[im][i] = C42_val[i] / C42_model[im][i];
      //      re[im][i] = fabs(C42_stat_Ref3X_besNEW[i] / C42_model[im][i]);
      re[im][i] = sqrt(pow(C42_stat[i]/C42_val[i], 2.0) + pow(C42_stat_model[im][i]/C42_model[im][i], 2.0)) * fabs(r[im][i]);
      res[im][i] = sqrt(pow(C42_sys[i]/C42_val[i], 2.0) + pow(C42_sys_model[im][i]/C42_model[im][i], 2.0)) * fabs(r[im][i]);
      // re[im][i] = sqrt(C42_stat[i]*C42_stat[i] + C42_stat_model[im][i]*C42_stat_model[im][i])  / C42_model[im][i];
      // res[im][i] = sqrt(C42_sys[i]*C42_sys[i] + C42_sys_model[im][i]*C42_sys_model[im][i])  / C42_model[im][i];

      ener[im][i+ntot_fxt] = C42_ener[i] * (1 + sc[im]);
      r_unity[im][i+ntot_fxt] = 1.0;
      resm[im][i+ntot_fxt] = fabs(C42_err_model[im][i] / C42_model[im][i]);  // model uncertainty
    }
  }
  TGraphErrors *gr_r[n_model], *gr_rs[n_model], *gr_rsm[n_model];
  for(int im=0;im<n_model;im++) {
    gr_r[im] = new TGraphErrors(ntot_datapts, xe[im], r[im], 0, re[im]);
    gr_r[im]->Print();
    gr_rs[im] = new TGraphErrors(ntot_datapts, xe[im], r[im], 0, res[im]);
  }

  double xe_fxt[n_model][ntot_fxt], r_fxt[n_model][ntot_fxt], re_fxt[n_model][ntot_fxt], res_fxt[n_model][ntot_fxt];
  for(int im=0;im<n_model;im++) {
    for(int i=0;i<ntot_fxt;i++) {
      if(fabs(C42_model_fxt[im][i])<1.e-4) continue;
      xe_fxt[im][i] = C42_ener_fxt[i];
      r_fxt[im][i] = C42_fxt[i] / C42_model_fxt[im][i];
      //      re_fxt[im][i] = fabs(C42_stat_fxt[i] / C42_model_fxt[im][i]);
      // re_fxt[im][i] = sqrt(C42_stat_fxt[i]*C42_stat_fxt[i] + C42_err_model_fxt[im][i]*C42_err_model_fxt[im][i]) / fabs(C42_model_fxt[im][i]);
      // res_fxt[im][i] = fabs(C42_sys_fxt[i] / C42_model_fxt[im][i]);
      re_fxt[im][i] = sqrt(pow(C42_stat_fxt[i]/C42_fxt[i], 2.0) + pow(C42_err_model_fxt[im][i]/C42_model_fxt[im][i], 2.0)) * r_fxt[im][i];
      res_fxt[im][i] = sqrt(pow(C42_sys_fxt[i]/C42_fxt[i], 2.0)) * r_fxt[im][i];

      ener[im][i] = C42_ener_fxt[i] * ( 1 + sc[im] );
      r_unity[im][i] = 1.0;      
      resm[im][i] = fabs(C42_err_model_fxt[im][i] / C42_model_fxt[im][i]);
    }
  }
  TGraphErrors *gr_r_fxt[n_model], *gr_rs_fxt[n_model];
  for(int im=0;im<n_model;im++) {
    gr_r_fxt[im] = new TGraphErrors(ntot_fxt, xe_fxt[im], r_fxt[im], 0, re_fxt[im]);
    gr_r_fxt[im]->Print();
    gr_rs_fxt[im] = new TGraphErrors(ntot_fxt, xe_fxt[im], r_fxt[im], 0, res_fxt[im]);
  }
  cout << " Model Uncertainties " << endl;
  for(int im=0;im<n_model;im++) {
    gr_rsm[im] = new TGraphErrors(ntot_fxt+ntot_datapts, ener[im], r_unity[im], 0, resm[im]);
    gr_rsm[im]->Print();
  }
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->SetLogx();
  //  c1->SetLogy();
  c1->Draw();

  double x1 = 1.8;
  double x2 = 300;
  double y1 = 0.0;
  double y2 = 2.2;
    
  TH1D *d0 = new TH1D("d0","",1, x1, x2);
  d0->SetMinimum(y1);
  d0->SetMaximum(y2);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  d0->GetXaxis()->SetLabelOffset(999.);
  d0->GetXaxis()->SetLabelSize(0.05);  
  d0->GetXaxis()->SetTitleOffset(1.2);  
  d0->GetXaxis()->SetTitleSize(0.065);  
  d0->GetYaxis()->SetNdivisions(205);  
  d0->GetYaxis()->SetTitle("[C_{4}/C_{2}]^{Data} / [C_{4}/C_{2}]^{Reference}");  
  //  d0->GetYaxis()->SetTitle("R_{42}^{data} / R_{42}^{model}");  
  d0->GetYaxis()->SetTitleOffset(1.1);  
  d0->GetYaxis()->SetTitleSize(0.068);  
  d0->GetYaxis()->SetLabelOffset(0.012);
  d0->GetYaxis()->SetLabelSize(0.05);
  d0->Draw("c");

  //  drawLine(x1, 0, x2, 0, 1, 9, 13);
  drawLine(x1, 1, x2, 1, 2, 9, 13);

  const Int_t markerColor[n_model] = {kBlue, kBlack, kBlack, kBlack};
  const Int_t markerStyle[n_model] = {25, 28, 26, 24};
  const Double_t markerSize[n_model] = {1.2, 1.4, 1.2, 1.6};
    
  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;
    cout << " Model " << im << endl;
    
    if(0) {
      setGraphFill(gr_rsm[im], 1001, kCyan-10, 1);
      gr_rsm[im]->Draw("e3");
    }
    
    drawSysBoxInRange(gr_rs[im], 0.03, 18, 1, y1, y2);
    setGraphMarker(gr_r[im], markerStyle[im], markerColor[im], markerSize[im]);
    setGraphLine(gr_r[im], 1, markerColor[im], 1);
    gr_r[im]->Draw("p");

    if(im==3) { // replot ratio to 70-80% data with solid marker
      TGraph *gr_d = (TGraph *)gr_r[im]->Clone(Form("Ratio_%d", im));
      gr_d->SetMarkerColor(kRed);
      gr_d->SetMarkerStyle(20);
      gr_d->SetMarkerSize(1.2);
      gr_d->Draw("p");
    }
    
    if(im==0) {
      drawSysBoxInRange(gr_rs_fxt[im], 0.03, 18, 1, y1, y2);
      setGraphMarker(gr_r_fxt[im], markerStyle[im], markerColor[im], markerSize[im]);
      setGraphLine(gr_r_fxt[im], 1, markerColor[im], 1);
      gr_r_fxt[im]->Draw("p");

      drawText(15, 1.95, "Au+Au Collisions at RHIC", 42, 0.05);
      drawText(13, 1.78, "0-5%, |y| < 0.5, 0.4 < p_{T} < 2.0 GeV/c", 42, 0.04);
    }
  }

  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;    
    gr_r[im]->Draw("p");
    gr_r_fxt[im]->Draw("p");
  }

  drawText(2.7, 1.1, "-0.5<y<0", 42, 0.03,90);
  
  drawText(72, 0.6, "Reference:", 42, 0.045);
  TLegend *leg = new TLegend(0.68, 0.2, 0.96, 0.38);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;
    leg->AddEntry(gr_r[im], NameModel[im], "p");
  }
  leg->Draw();

  const Int_t nt = 7;
  const Char_t *xlabel[nt] = {"2", "5", "10", "20", "50", "100", "200"};
  const Double_t xpos[nt] = {1.87, 4.7, 8.6, 17.5, 44, 82, 165};
  const Double_t ypos = -0.13;
  for(int it=0;it<nt;it++) {
    drawText(xpos[it], ypos, xlabel[it], 42, 0.05, 0);
  }
  const Int_t nl = 1;
  const Double_t xl[nl] = {57.58};
  const Double_t yl[nl] = {0.1184};
  TGraph *gr_l = new TGraph(nl, xl, yl);
  gr_l->SetMarkerColor(kRed);
  gr_l->SetMarkerStyle(20);
  gr_l->SetMarkerSize(1.2);
  gr_l->Draw("p");
  
  
  
  drawHistBox(x1, x2, y1, y2, 2222);
  c1->Update();
  c1->cd();
  
  c1->SaveAs("C42_ratio_one_paper.png");  
  c1->SaveAs("C42_ratio_one_paper.pdf");
  
}

