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
    //    return -A1*xl*TMath::Exp(-0.5*TMath::Power(xl/s1,2.0)) + 1; // - A2*xl*TMath::Exp(-0.25*T
    return -A1*xl*TMath::Exp(-0.5*TMath::Power(xl/s1,2.0)); // - A2*xl*TMath::Exp(-0.25*T
  } else {
    return -A1*(xl+A2*xl*xl)*TMath::Exp(-0.5*TMath::Power(xl/s2,2.0));
  }
}

void plotDiffPHDGSI()
{
  style();

  TF1 *fun_muB = new TF1("fun_muB","[0]/(1+[1]*x)",2.0, 1000);
  fun_muB->SetParameters(1541, 0.359);

  const Int_t NEL = 5;
  const Double_t EL[NEL] = {3.0, 4.5, 7.7, 19.6, 200.};
  const Char_t* TextL[NEL] = {"3.0", "4.5", "7.7", "19.6", "200"};
  Double_t MuBL[NEL];
  for(int i=0;i<NEL;i++) {
      MuBL[i] = fun_muB->Eval(EL[i]);
  }
  
  //////////=====
  const int ntot_datapts=11;
  const int n_model = 4;
  const int plotflag[n_model] = {1, 1, 1, 0}; // UrQMD, HRG CE, Hydro EV, 70-80%

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
  double muB[ntot_datapts];
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
    muB[i] = fun_muB->Eval(sNN[i]);
    cout << " MuB = " << muB[i] << endl;
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
  // FXT energies
  ////////////////
  const int ntot_fxt = 4;
  double C42_ener_fxt[ntot_fxt] = {3.0, 3.2, 3.5, 3.9};
  double C42_fxt[ntot_fxt] = {-0.84574, -0.305, 0.155, 0.210};
  double C42_stat_fxt[ntot_fxt] = {0.0862692, 0.097, 0.162, 0.235};
  double C42_sys_fxt[ntot_fxt] = {0.821897, 0.183, 0.138, 0.166};
  double muB_fxt[ntot_fxt];
  for(int i=0;i<ntot_fxt;i++) {
    muB_fxt[i] = fun_muB->Eval(C42_ener_fxt[i]);
    cout << " muB fxt = " << muB_fxt[i] << endl;
  }
  
  double C42_model_fxt[n_model][ntot_fxt] = {{-0.774123, -0.135, 0.238, 0.431},
					     {0., 0., 0., 0.},
					     {0., 0., 0., 0.},
  					     {0., 0., 0., 0.}
  };
  double C42_err_model_fxt[n_model][ntot_fxt] = {{0.0907071, 0.065, 0.060, 0.064},
					     {0., 0., 0., 0.},
					     {0., 0., 0., 0.},
  					     {0., 0., 0., 0.}
  };
  
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
  const Double_t sc_muB[n_model] = {-10, -5, 0., 5};
  double xe[n_model][ntot_datapts], xmuB[n_model][ntot_datapts], diff[n_model][ntot_datapts], diffe[n_model][ntot_datapts], diffes[n_model][ntot_datapts];
  double ener[n_model][ntot_datapts+ntot_fxt], diff_zero[n_model][ntot_datapts+ntot_fxt], diffesm[n_model][ntot_datapts+ntot_fxt];
  for(int im=0;im<n_model;im++) {
    for(int i=0;i<ntot_datapts;i++) {
      xe[im][i] = C42_ener[i] * (1 + sc[im]);
      xmuB[im][i] = muB[i] + sc_muB[im];
      diff[im][i] = C42_val[i] - C42_model[im][i];
      diffe[im][i] = sqrt(pow(C42_stat[i], 2.0) + pow(C42_stat_model[im][i], 2.0));
      diffes[im][i] = sqrt(pow(C42_sys[i], 2.0) + pow(C42_sys_model[im][i], 2.0));

      ener[im][i+ntot_fxt] = C42_ener[i] * (1 + sc[im]);
      diff_zero[im][i+ntot_fxt] = 0.0;
      diffesm[im][i+ntot_fxt] = C42_err_model[im][i];  // model uncertainty
    }
  }
  TGraphErrors *gr_d[n_model], *gr_ds[n_model], *gr_dsm[n_model];
  TGraphErrors *gr_muB_d[n_model], *gr_muB_ds[n_model], *gr_muB_dsm[n_model];
  for(int im=0;im<n_model;im++) {
    gr_d[im] = new TGraphErrors(ntot_datapts, xe[im], diff[im], 0, diffe[im]);
    gr_d[im]->Print();
    gr_ds[im] = new TGraphErrors(ntot_datapts, xe[im], diff[im], 0, diffes[im]);
    
    gr_muB_d[im] = new TGraphErrors(ntot_datapts, xmuB[im], diff[im], 0, diffe[im]);
    gr_muB_d[im]->Print();
    gr_muB_ds[im] = new TGraphErrors(ntot_datapts, xmuB[im], diff[im], 0, diffes[im]);
  }

  double xe_fxt[n_model][ntot_fxt], xmuB_fxt[n_model][ntot_fxt], diff_fxt[n_model][ntot_fxt], diffe_fxt[n_model][ntot_fxt], diffes_fxt[n_model][ntot_fxt];
  for(int im=0;im<n_model;im++) {
    for(int i=0;i<ntot_fxt;i++) {
      if(fabs(C42_model_fxt[im][i])<1.e-4) continue;
      xe_fxt[im][i] = C42_ener_fxt[i];
      xmuB_fxt[im][i] = muB_fxt[i];
      diff_fxt[im][i] = C42_fxt[i] - C42_model_fxt[im][i];
      diffe_fxt[im][i] = sqrt(pow(C42_stat_fxt[i], 2.0) + pow(C42_err_model_fxt[im][i], 2.0));
      diffes_fxt[im][i] = C42_sys_fxt[i];

      ener[im][i] = C42_ener_fxt[i] * ( 1 + sc[im] );
      diff_zero[im][i] = 0.0;      
      diffesm[im][i] = C42_err_model_fxt[im][i];
    }
  }
  TGraphErrors *gr_d_fxt[n_model], *gr_ds_fxt[n_model];
  TGraphErrors *gr_muB_d_fxt[n_model], *gr_muB_ds_fxt[n_model];
  for(int im=0;im<n_model;im++) {
    gr_d_fxt[im] = new TGraphErrors(ntot_fxt, xe_fxt[im], diff_fxt[im], 0, diffe_fxt[im]);
    gr_d_fxt[im]->Print();
    gr_ds_fxt[im] = new TGraphErrors(ntot_fxt, xe_fxt[im], diff_fxt[im], 0, diffes_fxt[im]);

    gr_muB_d_fxt[im] = new TGraphErrors(ntot_fxt, xmuB_fxt[im], diff_fxt[im], 0, diffe_fxt[im]);
    gr_muB_d_fxt[im]->Print();
    gr_muB_ds_fxt[im] = new TGraphErrors(ntot_fxt, xmuB_fxt[im], diff_fxt[im], 0, diffes_fxt[im]);
  }
  cout << " Model Uncertainties " << endl;
  for(int im=0;im<n_model;im++) {
    gr_dsm[im] = new TGraphErrors(ntot_fxt+ntot_datapts, ener[im], diff_zero[im], 0, diffesm[im]);
    gr_dsm[im]->Print();
  }
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->SetLogx();
  //  c1->SetLogy();
  c1->Draw();

  double x1 = 1.8;
  double x2 = 300;
  double y1 = -0.65;
  double y2 = 0.59;
    
  TH1D *d0 = new TH1D("d0","",1, x1, x2);
  d0->SetMinimum(y1);
  d0->SetMaximum(y2);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  d0->GetXaxis()->SetLabelOffset(999.);
  d0->GetXaxis()->SetLabelSize(0.05);  
  d0->GetXaxis()->SetTitleOffset(1.2);  
  d0->GetXaxis()->SetTitleSize(0.065);  
  d0->GetYaxis()->SetNdivisions(405);  
  d0->GetYaxis()->SetTitle("[C_{4}/C_{2}]^{Data} #minus [C_{4}/C_{2}]^{Reference}");  
  //  d0->GetYaxis()->SetTitle("R_{42}^{data} / R_{42}^{model}");  
  d0->GetYaxis()->SetTitleOffset(1.1);  
  d0->GetYaxis()->SetTitleSize(0.068);  
  d0->GetYaxis()->SetLabelOffset(0.012);
  d0->GetYaxis()->SetLabelSize(0.05);
  d0->Draw("c");

  //  drawLine(x1, 0, x2, 0, 1, 9, 13);
  drawLine(x1, 0, x2, 0, 2, 9, 13);

  //  const int plotflag[n_model] = {1, 1, 1, 0}; // UrQMD, HRG CE, Hydro EV, 70-80%
  const Int_t markerColor[n_model] = {kBlack, kRed, kGreen+3, kBlack};
  const Int_t markerStyle[n_model] = {20, 21, 22, 24};
  const Double_t markerSize[n_model] = {1.2, 1.4, 1.2, 1.6};
    
  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;
    cout << " Model " << im << endl;
    
    if(0) {
      setGraphFill(gr_dsm[im], 1001, kCyan-10, 1);
      gr_dsm[im]->Draw("e3");
    }
    
    drawSysBoxInRange(gr_ds[im], 0.03, 18, 1, y1, y2);
    setGraphMarker(gr_d[im], markerStyle[im], markerColor[im], markerSize[im]);
    setGraphLine(gr_d[im], 1, markerColor[im], 2);
    gr_d[im]->Draw("p");

    if(im==3) { // replot ratio to 70-80% data with solid marker
      TGraph *gr_d3 = (TGraph *)gr_d[im]->Clone(Form("Diff_%d", im));
      gr_d3->SetMarkerColor(kRed);
      gr_d3->SetMarkerStyle(20);
      gr_d3->SetMarkerSize(1.2);
      gr_d3->Draw("p");
    }
    
    if(im==0) {
      drawSysBoxInRange(gr_ds_fxt[im], 0.03, 18, 1, y1, y2);
      setGraphMarker(gr_d_fxt[im], markerStyle[im], markerColor[im], markerSize[im]);
      setGraphLine(gr_d_fxt[im], 1, markerColor[im], 2);
      gr_d_fxt[im]->Draw("p");

      drawText(15, y1+(y2-y1)*0.90, "Au+Au Collisions at RHIC", 42, 0.05);
      drawText(13, y1+(y2-y1)*0.825, "0-5%, |y| < 0.5, 0.4 < p_{T} < 2.0 GeV/c", 42, 0.04);
    }
  }

  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;    
    gr_d[im]->Draw("p");
    gr_d_fxt[im]->Draw("p");
  }

  drawText(2.7, -0.2, "-0.5<y<0", 42, 0.03,90);
  
  drawText(72, y1+(y2-y1)*0.273, "Reference:", 42, 0.045);
  TLegend *leg = new TLegend(0.68, 0.2, 0.96, 0.38);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;
    leg->AddEntry(gr_d[im], NameModel[im], "pl");
  }
  leg->Draw();

  const Int_t nt = 7;
  const Char_t *xlabel[nt] = {"2", "5", "10", "20", "50", "100", "200"};
  const Double_t xpos[nt] = {1.87, 4.7, 8.6, 17.5, 44, 82, 165};
  const Double_t ypos = y1+(y2-y1)*(-0.059);
  for(int it=0;it<nt;it++) {
    drawText(xpos[it], ypos, xlabel[it], 42, 0.05, 0);
  }
  const Int_t nl = 1;
  const Double_t xl[nl] = {57.58};
  const Double_t yl[nl] = {y1+(y2-y1)*0.05382};
  TGraph *gr_l = new TGraph(nl, xl, yl);
  gr_l->SetMarkerColor(kRed);
  gr_l->SetMarkerStyle(20);
  gr_l->SetMarkerSize(1.2);
  //  gr_l->Draw("p");
  
  
  
  drawHistBox(x1, x2, y1, y2, 2222);
  c1->Update();
  c1->cd();
  
  c1->SaveAs("C42_diff_QM25.png");  
  c1->SaveAs("C42_diff_QM25.pdf");

  //////////////////////////////////////////////////
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  c2->SetTopMargin(0.06);
  //  c1->SetLogx();
  c2->SetTickx(0);
  c2->Draw();

  x1 = 0.0;
  x2 = 900;
    
  TH1D *d1 = new TH1D("d1","",1, x1, x2);
  d1->SetMinimum(y1);
  d1->SetMaximum(y2);
  d1->GetXaxis()->CenterTitle();
  d1->GetXaxis()->SetTitle("Baryon Chemical Potential #mu_{B} (MeV)");
  d1->GetXaxis()->SetLabelOffset(0.008);
  d1->GetXaxis()->SetLabelSize(0.05);  
  d1->GetXaxis()->SetTitleOffset(1.2);  
  d1->GetXaxis()->SetTitleSize(0.065);  
  d1->GetYaxis()->SetNdivisions(205);  
  d1->GetYaxis()->SetTitle("[C_{4}/C_{2}]^{Data} #minus [C_{4}/C_{2}]^{Reference}");  
  //  d1->GetYaxis()->SetTitle("R_{42}^{data} / R_{42}^{model}");  
  d1->GetYaxis()->SetTitleOffset(1.1);  
  d1->GetYaxis()->SetTitleSize(0.065);  
  d1->GetYaxis()->SetLabelOffset(0.01);
  d1->GetYaxis()->SetLabelSize(0.05);
  d1->Draw("c");

  //  drawLine(x1, 0, x2, 0, 1, 9, 13);
  drawLine(x1, 0, x2, 0, 2, 9, 13);

  TF1 *fun = new TF1("fun",xgx,x1,x2,5);
  //  fun->SetParameters(4.5, 0.2, 11.5, 0.2, 0.18);
  fun->SetParameters(3.0, 0.2, 11.5, 0.25, 0.18);
  // fun->SetLineWidth(4);
  // fun->SetLineStyle(9);
  // fun->SetLineColor(1);
  //  fun->Draw("same");

  TF1 *fun_p =  new TF1("fun_p",xgx,x1,fun->GetParameter(2),5);
  //  fun_p->SetParameters(4.5, 0.2, 11.5, 0.2, 0.18);
  fun_p->SetParameters(3.0, 0.2, 11.5, 0.25, 0.18);
  //  fun_p->SetRange(x1, fun->GetParameter(2));
  fun_p->SetLineColor(kBlue);
  fun_p->SetLineWidth(6);
  fun_p->SetLineStyle(9);
  //  fun_p->Draw("same");

  TF1 *fun_n =  new TF1("fun_n",xgx,fun->GetParameter(2), x2,5);
  //  fun_n->SetParameters(4.5, 0.2, 11.5, 0.2, 0.18);
  fun_n->SetParameters(3.0, 0.2, 11.5, 0.25, 0.18);
  //  fun_n->SetRange(fun->GetParameter(2), x2);
  fun_n->SetLineColor(kRed);
  fun_n->SetLineWidth(6);
  fun_n->SetLineStyle(9);
  //  fun_n->Draw("same");

  const Int_t np = 100;
  double ener_p[np], muB_p[np], y_p[np], y_p1[np];
  double ener_n[np], muB_n[np], y_n[np];
  for(int i=0;i<np;i++) {
    //    ener_p[i] = (11.5 - 2.0)/np*(i+0.1) + 2.0;
    ener_p[i] = (11.5 - 2.0)/np*(i+0.1) + 2.0;
    muB_p[i] = fun_muB->Eval(ener_p[i]);
    y_p[i] = fun_p->Eval(ener_p[i]);
    y_p1[i] = fun_p->Eval(ener_p[i]) * 0.5;
    
    //    ener_n[i] = (200 - 11.5)/np*(i+0.1) + 11.5;
    ener_n[i] = (200 - 11.5)/np*(i+0.1) + 11.5;
    muB_n[i] = fun_muB->Eval(ener_n[i]);
    y_n[i] = fun_p->Eval(ener_n[i]);
  }

  TGraph *gr_p = new TGraph(np, muB_p, y_p);
  gr_p->SetLineColor(kBlue);
  gr_p->SetLineWidth(6);
  gr_p->SetLineStyle(9);
  gr_p->Draw("c");
  
  TGraph *gr_p1 = new TGraph(np, muB_p, y_p1);
  gr_p1->SetLineColor(kGray);
  gr_p1->SetLineWidth(6);
  gr_p1->SetLineStyle(9);
  gr_p1->Draw("c");

  TGraph *gr_n = new TGraph(np, muB_n, y_n);
  gr_n->SetLineColor(kRed);
  gr_n->SetLineWidth(6);
  gr_n->SetLineStyle(9);
  gr_n->Draw("c");

  // energy ranges
  double energyR[6] = {27, 7.7, 4.9, 3.9, 3.0, 2.0};
  double muR[6];
  for(int i=0;i<6;i++) {
    muR[i] = fun_muB->Eval(energyR[i]);
  }


  TLine *ll[NEL];
  for(int i=0;i<NEL;i++) {
    ll[i] = new TLine(MuBL[i], y2*0.98, MuBL[i], y2);
    ll[i]->SetLineWidth(2);
    //    ll[i]->Draw("same");
    drawText(MuBL[i]-20, y2*1.02, TextL[i], 12, 0.04);
  }
  drawText(x1-80, y2*1.03, "#sqrt{s_{NN}}", 12, 0.04);
  drawText(x2-50, y2*1.02, "GeV", 12, 0.04);  
  
  TFile *fout = new TFile("C42_Diff_Func.root","recreate");
  fun->Write();
  fun_p->Write();
  fun_n->Write();
  fout->Close();
  

  TF1 *fitfun = new TF1("fitfun","-[0]*x*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2.0))",7., 200);
  fitfun->SetParameters(1., 11, 10);
    
  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;
    cout << " Model " << im << endl;
        
    drawSysBoxInRange(gr_muB_ds[im], 3, 18, 0, y1, y2);
    setGraphMarker(gr_muB_d[im], markerStyle[im], markerColor[im], 1.2);
    setGraphLine(gr_muB_d[im], 1, markerColor[im]);
    gr_muB_d[im]->Draw("p");

    //gr_r[im]->Fit("fitfun","R");
    
    if(im==0) {
      drawSysBoxInRange(gr_muB_ds_fxt[im], 3, 18, 0, y1, y2);
      setGraphMarker(gr_muB_d_fxt[im], markerStyle[im]+4, markerColor[im], 1.2);
      setGraphLine(gr_muB_d_fxt[im], 1, markerColor[im]);
      gr_muB_d_fxt[im]->Draw("p");

      drawText(50, 2.0, "Au+Au 0-5%", 42, 0.05);
      drawText(40, 1.85, "0.4<p_{T}<2.0 GeV/c, |y|<0.5", 42, 0.03);
    }
    
  }

  for(int im=0;im<n_model;im++) {
    if(!plotflag[im]) continue;    
    gr_muB_d[im]->Draw("p");
  }

  //    gr_n->Draw("c");


  drawText(770, -0.3, "-0.5<y<0", 42, 0.03,90);
  
  TLegend *leg1 = new TLegend(0.74, 0.77, 0.96, 0.93);
  leg1->SetLineColor(10);
  leg1->SetTextSize(0.035);
  for(int im=n_model-1;im>=0;im--) {
    if(!plotflag[im]) continue;
    leg1->AddEntry(gr_muB_d[im], NameModel[im], "pl");
  }
  leg1->Draw();

  drawArrow(fun_muB->Eval(27.0), y1+0.05, fun_muB->Eval(7.7), y1+0.05, 0.03, 30, 2, 1, 1, "<|>");
  drawText(fun_muB->Eval(17), y1+0.075, "BES-II COL", 52, 0.03);
  drawArrow(fun_muB->Eval(4.5), y1+0.05, fun_muB->Eval(3.0), y1+0.05, 0.03, 30, 2, 1, 1, "<|>");
  drawText(fun_muB->Eval(4.3), y1+0.075, "BES-II FXT", 52, 0.03);
  drawArrow(fun_muB->Eval(4.9), y1+0.16, fun_muB->Eval(2.0), y1+0.16, 0.03, 30, 2, 1, kGreen+3, "<|>");
  drawText(fun_muB->Eval(3.5), y1+0.185, "SIS100", 52, 0.03, 0, kGreen+3);
  drawArrow(fun_muB->Eval(8.0), y1+0.18, fun_muB->Eval(2.0), y1+0.18, 0.03, 30, 2, 1, kGreen+3, "<|>");
  drawText(fun_muB->Eval(3.5), y1+0.205, "SIS300", 52, 0.03, 0, kGreen+3);
  
  
  /*
  const Int_t nt = 4;
  const Char_t *xlabel[nt] = {"3", "10", "30", "100"};
  const Double_t xpos[nt] = {2.8, 8.5, 27, 82};
  const Double_t ypos[nt] = {-0.15, -0.15, -0.15, -0.15};
  for(int it=0;it<nt;it++) {
    drawText(xpos[it], ypos[it], xlabel[it], 42, 0.05, 0);
  }
  */
  drawHistBox(x1, x2, y1, y2);
  c2->Update();
  c2->cd();
  
  c2->SaveAs("C42_muB_diff_PHDGSI.png");  
  c2->SaveAs("C42_muB_diff_PHDGSI.pdf");
  
  
}

