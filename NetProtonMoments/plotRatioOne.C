#include "draw.C+"
#include "style.C+"

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
    return -A1*xl*TMath::Exp(-0.5*TMath::Power(xl/s1,2.0)) + 1; // - A2*xl*TMath::Exp(-0.25*T
  } else {
    return -A1*(xl+A2*xl*xl)*TMath::Exp(-0.5*TMath::Power(xl/s2,2.0)) + 1;
  }
}


void plotRatioOne(const int conf = 1)  // conf:  1 - combine model uncertainty to data stat; 0 - separate model uncertainty
{
  style();

  //////////=====
  const int ntot_datapts=11;
  const int n_model = 4;
  const int plotflag[n_model] = {1, 1, 1, 1};
  
  ////data 0-5%
  double C42_ener_besNEW[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double C42_main_Ref3X_besNEW[ntot_datapts]={0.4119085,0.540458,0.4011784,0.430751,0.31512,0.339941,0.602046,0.739693,0.632837,0.792955,0.900669};
  double C42_stat_Ref3X_besNEW[ntot_datapts]={0.2468968,0.187029,0.1347458,0.104085,0.120367,0.0731331,0.0950052,0.147006,0.0553742,0.250823,0.2084582};
  double C42_sys_Ref3X_besNEW[ntot_datapts]={0.1302093,0.0680524,0.06472035,0.0403519,0.0880467,0.0389432,0.0304235,0.1357538,0.1387023,0.1219,0.1393631};
  double C42_tot_Ref3X_err_besNEW[ntot_datapts]={0};
  
  ////UrQMD
  double cen_UQMD[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double urqmd_C42_05_data[ntot_datapts]={0.4753277, 0.4625844, 0.5253311, 0.594994, 0.6782202, 0.7327199, 0.903167,0.860031,0.83973130,0.829186,0.914062};
  double urqmd_C42_05_sterr[ntot_datapts]={0.06445722,0.06283556,0.06766454,0.0707438,0.05287615,0.05808925,0.063723,0.045722,0.055722,0.081655,0.072831};
  
  ////HRG CE
  double anar_cen[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double anar_C42[ntot_datapts]={0.40214780,0.44213425,0.5034468,0.58341961,0.63475107,0.67691625,0.73933998,0.76389,0.7953959,0.8117625,0.84005};
  
  ////Hydro
  double Hydro_cen[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double Hydro_C42[ntot_datapts]={0.34623,0.39264044,0.46380312,0.556624,0.58357922,0.605721,0.657514,0.758193,0.78623350,0.8008,0.838933};

  ////Hydro+EV
  double Hydro_cen_wtEV[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double Hydro_C42_wtEV[ntot_datapts]={0.239113,0.28855366,0.36436268,0.463244,0.49676714,0.524304,0.579964,0.693609,0.72515728,0.741546,0.785968};  
  
  for(int y=0;y<ntot_datapts;y++)
    {
      C42_tot_Ref3X_err_besNEW[y]=sqrt(C42_stat_Ref3X_besNEW[y]*C42_stat_Ref3X_besNEW[y]+C42_sys_Ref3X_besNEW[y]*C42_sys_Ref3X_besNEW[y]);
      cout << C42_tot_Ref3X_err_besNEW[y] << endl;
    }

  ////////////////////////////////////
  // rename arrays for simplicity
  ////////////////////////////////////
  double sNN[ntot_datapts];
  double C42_data[ntot_datapts];
  double C42_err_data[ntot_datapts];
  double C42_model[n_model][ntot_datapts];
  double C42_err_model[n_model][ntot_datapts];
  const Char_t *NameModel[n_model] = {"UrQMD", "HRG-CE","Hydro", "Hydro+EV"};
  for(int i=0;i<ntot_datapts;i++) {
    sNN[i] = C42_ener_besNEW[i];
    C42_data[i] = C42_main_Ref3X_besNEW[i];
    C42_err_data[i] = C42_tot_Ref3X_err_besNEW[i];
    
    for(int j=0;j<n_model;j++) {
      if(j==0) {
	C42_model[j][i] = urqmd_C42_05_data[i];
	C42_err_model[j][i] = urqmd_C42_05_sterr[i];
      } else if(j==1) {
	C42_model[j][i] = anar_C42[i];
	C42_err_model[j][i] = 0.0;
      } else if(j==2) {
	C42_model[j][i] = Hydro_C42[i];
	C42_err_model[j][i] = 0.0;
      } else {
	C42_model[j][i] = Hydro_C42_wtEV[i];
	C42_err_model[j][i] = 0.0;
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
  
  double C42_model_fxt[n_model][ntot_fxt] = {{-0.774123}, {0}, {0}};
  double C42_err_model_fxt[n_model][ntot_fxt] = {{0.0907071}, {0}, {0}};
  
  ///////////////////////////////////
  //////////===== chi2 test /////////
  ///////////////////////////////////
  
  const int n_test = 7;
  double chi2[n_model] = {0.0, 0.0, 0.0};

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
  const Double_t sc[n_model] = {-0.07, -0.035, 0., 0.035};
  double xe[n_model][ntot_datapts], r[n_model][ntot_datapts], re[n_model][ntot_datapts], res[n_model][ntot_datapts];
  double ener[n_model][ntot_datapts+ntot_fxt], r_unity[n_model][ntot_datapts+ntot_fxt], resm[n_model][ntot_datapts+ntot_fxt];
  for(int im=0;im<n_model;im++) {
    for(int i=0;i<ntot_datapts;i++) {
      xe[im][i] = C42_ener_besNEW[i] * (1 + sc[im]);
      r[im][i] = C42_main_Ref3X_besNEW[i] / C42_model[im][i];
      //      re[im][i] = fabs(C42_stat_Ref3X_besNEW[i] / C42_model[im][i]);
      re[im][i] = sqrt(C42_stat_Ref3X_besNEW[i]*C42_stat_Ref3X_besNEW[i] + C42_err_model[im][i]*C42_err_model[im][i])  / C42_model[im][i];
      res[im][i] = fabs(C42_sys_Ref3X_besNEW[i] / C42_model[im][i]);

      ener[im][i+ntot_fxt] = C42_ener_besNEW[i] * (1 + sc[im]);
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
      re_fxt[im][i] = sqrt(C42_stat_fxt[i]*C42_stat_fxt[i] + C42_err_model_fxt[im][i]*C42_err_model_fxt[im][i]) / C42_model_fxt[im][i];
      res_fxt[im][i] = fabs(C42_sys_fxt[i] / C42_model_fxt[im][i]);

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
  d0->GetYaxis()->SetNdivisions(204);  
  d0->GetYaxis()->SetTitle("R_{42}^{data} / R_{42}^{model}");  
  d0->GetYaxis()->SetTitleOffset(1.0);  
  d0->GetYaxis()->SetTitleSize(0.065);  
  d0->GetYaxis()->SetLabelOffset(0.015);
  d0->GetYaxis()->SetLabelSize(0.05);
  d0->Draw("c");

  //  drawLine(x1, 0, x2, 0, 1, 9, 13);
  drawLine(x1, 1, x2, 1, 2, 9, 13);


  TF1 *fun = new TF1("fun",xgx,x1,x2,5);
  fun->SetParameters(4.5, 0.2, 11.5, 0.2, 0.18);
  // fun->SetLineWidth(4);
  // fun->SetLineStyle(9);
  // fun->SetLineColor(1);
  //  fun->Draw("same");

  TF1 *fun_p =  new TF1("fun_p",xgx,x1,fun->GetParameter(2),5);
  fun_p->SetParameters(4.5, 0.2, 11.5, 0.2, 0.18);
  //  fun_p->SetRange(x1, fun->GetParameter(2));
  fun_p->SetLineColor(kBlue);
  fun_p->SetLineWidth(8);
  fun_p->SetLineStyle(9);
  fun_p->Draw("same");

  TF1 *fun_n =  new TF1("fun_n",xgx,fun->GetParameter(2), x2,5);
  fun_n->SetParameters(4.5, 0.2, 11.5, 0.2, 0.18);
  //  fun_n->SetRange(fun->GetParameter(2), x2);
  fun_n->SetLineColor(kRed);
  fun_n->SetLineWidth(8);
  fun_n->SetLineStyle(9);
  fun_n->Draw("same");
  
  TFile *fout = new TFile("C42_Ratio_Func.root","recreate");
  fun->Write();
  fun_p->Write();
  fun_n->Write();
  fout->Close();
  

  TF1 *fitfun = new TF1("fitfun","-[0]*x*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2.0))",7., 200);
  fitfun->SetParameters(1., 11, 10);

  const Int_t markerColor[n_model] = {kBlack, kRed, kGreen+2, kGreen+2};
  const Int_t markerStyle[n_model] = {20, 21, 22, 23};
  // const Int_t markerColor[n_model] = {kBlack, kRed, kBlue};
  // const Int_t markerStyle[n_model] = {20, 21, 22};
    
  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;
    cout << " Model " << im << endl;
    
    if(0) {
      setGraphFill(gr_rsm[im], 1001, kCyan-10, 1);
      gr_rsm[im]->Draw("e3");
    }
    
    drawSysBoxInRange(gr_rs[im], 0.03, 18, 1, y1, y2);
    setGraphMarker(gr_r[im], markerStyle[im], markerColor[im], 1.2);
    setGraphLine(gr_r[im], 1, markerColor[im]);
    gr_r[im]->Draw("p");

    //gr_r[im]->Fit("fitfun","R");
    
    if(im==0) {
      drawSysBoxInRange(gr_rs_fxt[im], 0.03, 18, 1, y1, y2);
      setGraphMarker(gr_r_fxt[im], markerStyle[im], markerColor[im], 1.2);
      setGraphLine(gr_r_fxt[im], 1, markerColor[im]);
      gr_r_fxt[im]->Draw("p");

      drawText(50, 1.95, "Au+Au 0-5%", 42, 0.05);
      drawText(40, 1.8, "0.4<p_{T}<2.0 GeV/c, |y|<0.5", 42, 0.03);
    }
    
    
    if(0) {
      TLegend *leg = new TLegend(0.7, 0.66, 0.92, 0.84);
      leg->SetLineColor(10);
      leg->SetTextSize(0.04);
      leg->AddEntry(gr_r[im], "  STAR", "pl");
      leg->Draw();
    }
    
  }

  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;    
    gr_r[im]->Draw("p");
  }

  drawText(2.7, 1.1, "-0.5<y<0", 42, 0.03,90);
  
  TLegend *leg = new TLegend(0.74, 0.2, 0.96, 0.38);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  for(int im=n_model-1;im>=0;im--) {
    if(!plotflag[im]) continue;
    leg->AddEntry(gr_r[im], NameModel[im], "pl");
  }
  leg->Draw();

  const Int_t nt = 4;
  const Char_t *xlabel[nt] = {"3", "10", "30", "100"};
  const Double_t xpos[nt] = {2.8, 8.5, 27, 82};
  const Double_t ypos[nt] = {-0.15, -0.15, -0.15, -0.15};
  for(int it=0;it<nt;it++) {
    drawText(xpos[it], ypos[it], xlabel[it], 42, 0.05, 0);
  }
  
  drawHistBox(x1, x2, y1, y2);
  c1->Update();
  c1->cd();
  
  c1->SaveAs("C42_ratio_one.png");  
  c1->SaveAs("C42_ratio_one.pdf");
  
}

