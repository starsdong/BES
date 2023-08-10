#include "style.C+"
#include "draw.C+"
#include "findChi2Min.C"
#include "reset_error.C"

double ideogram0(double *x, double *par)
{
  const Int_t NP = 6;
  double func = 0;
  for(int i=0;i<NP;i++) {
    if(fabs(par[i*2])<0.1) continue; // skip zeros
    func += 1./par[i*2+1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[i*2])/par[i*2+1],2.0));
  }
  return func;
}

double ideogram(double *x, double *par)
{
  const Int_t NP = 6;
  double func = 0;
  for(int i=0;i<NP;i++) {
    if(fabs(par[i*3])<0.1) continue; // skip zeros
    double sig = (par[i*3+1]+par[i*3+2])/2.0;
    if(x[0]<par[i*3]) func += 1./sig/sig*TMath::Exp(-0.5*TMath::Power((x[0]-par[i*3])/par[i*3+1],2.0));
    else func += 1./sig/sig*TMath::Exp(-0.5*TMath::Power((x[0]-par[i*3])/par[i*3+2],2.0));
  }
  return func;
}

void BL_H3L(int config=0){
  style();

  const Int_t NP = 6;
  Int_t NC = NP;
  if(config) NC = NP-1;  // config=0:  all points   config=1:  exclude the last ALICE data
  const Double_t XMIN = -0.30;
  const Double_t XMAX = 0.70;
  
  const Double_t Tau_Lambda = 0.;
  const Int_t NP_HIS = 2; // number of history experimental data
  //BL [MeV]	Stat. (upper)	Stat. (lower)	Syst. (upper)	Syst. (lower)
  const double data_all[NP][5] = {
    {0.41,  0.12,  0.12,  0.0,   0.0}, // NPB1(67)
    {0.08,  0.07,  0.07,  0.0,   0.0}, // NPB4(68)
    {-0.13, 0.27,  0.27,  0.0,   0.0}, // PRD1(70)
    {0.27,  0.08,  0.08,  0.0,   0.0}, // NPB52(73)
    {0.41,  0.12,  0.12,  0.11,  0.11}, // STAR(20)
    {0.072, 0.063, 0.063, 0.036, 0.036}  // ALICE(22)
  };

  const Char_t Label[NP][20] = {
    "NPB1 1967",
    "NPB4 1968",
    "PRD1 1970",
    "NPB52 1973",
    "STAR 2020",
    "ALICE 2022"
  };

  // Calculate the average
  double xp[NP], yp[NP], eyl[NP], eyh[NP], ey[NP];
  for(int i=0;i<NP;i++) {
    xp[i] = i+1;
    yp[i] = data_all[i][0];
    eyh[i] = sqrt(data_all[i][1]*data_all[i][1] + data_all[i][3]*data_all[i][3]);
    eyl[i] = sqrt(data_all[i][2]*data_all[i][2] + data_all[i][4]*data_all[i][4]);
  }

  TGraphAsymmErrors *gr_data = new TGraphAsymmErrors(NC, xp, yp, 0, 0, eyl, eyh);
  TGraphAsymmErrors *gr_data_HI = new TGraphAsymmErrors(NP-NP_HIS, xp+NP_HIS, yp+NP_HIS, 0, 0, eyl+NP_HIS, eyh+NP_HIS);


  TGraphAsymmErrors *gr_data_y = new TGraphAsymmErrors(NC, yp, xp, eyl, eyh, 0, 0);
  
  // Test on Chi2 calculation
  const Int_t Nf = 200;
  double xf[Nf+1], chi2[Nf+1];
  for(int ip = 0;ip<Nf+1;ip++) {
    xf[ip] = XMIN + ip*(XMAX-XMIN)/Nf;
    chi2[ip] = 0;
    for(int i=0;i<NC;i++) {
      double err = reset_error(yp[i], eyl[i], eyh[i], xf[ip]);
      chi2[ip] += pow(fabs(xf[ip] - yp[i])/err, 2.0);
    }
  }

  //
  // binding energy plot, perform fit
  //  
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->SetLeftMargin(0.1);
  c1->SetBottomMargin(0.05);
  c1->cd();
  
  double x1 = 0;
  double x2 = NP+1;
  double y1 = XMIN;
  double y2 = XMAX;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetYaxis()->CenterTitle();
  h0->SetYTitle("{}^{3}_{#Lambda}H B_{#Lambda} (MeV)");
  h0->GetXaxis()->SetTitleOffset(999.);
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw();

  gr_data->SetMarkerSize(1.5);
  gr_data->SetMarkerStyle(20);
  gr_data->SetLineWidth(2);
  gr_data->Draw("p");

  drawLine(x1, Tau_Lambda, x2, Tau_Lambda, 2, 2, 1);
  
  TF1 *func = new TF1("func","pol0",x1,x2);
  func->SetParameter(0, 250);
  gr_data->Fit("func","R");

  double muAve = func->GetParameter(0);
  double sigAve = func->GetParError(0);
  double chi2min = func->GetChisquare();
  double muAve_old = 999.;
  double sigAve_old = 999.;
  double chi2min_old = 999.;
  while(fabs(chi2min_old - chi2min)>1.e-3) {
    muAve = func->GetParameter(0);
    sigAve = func->GetParError(0);
    chi2min = func->GetChisquare();
    for(int i=0;i<NP;i++) {
      ey[i] = reset_error(yp[i], eyl[i], eyh[i], muAve);
    }
    TGraphErrors *gr_data_reset = new TGraphErrors(NC, xp, yp, 0, ey);
    
    gr_data_reset->Fit("func","R");
    muAve_old = muAve;
    sigAve_old = sigAve;
    chi2min_old = chi2min;
    muAve = func->GetParameter(0);
    sigAve = func->GetParError(0);
    chi2min = func->GetChisquare();   
  }

  gr_data->SetName("data_raw");
  gr_data->Draw("p");
  
  for(int i=0;i<NP;i++) {
    ey[i] = reset_error(yp[i], eyl[i], eyh[i], muAve);
  }
  TGraphErrors *gr_data_reset = new TGraphErrors(NC, xp, yp, 0, ey);
  gr_data_reset->SetName("data_reset");

  double sigAve_wt = sigAve;
  double scale = sqrt(chi2min/(NC-1));
  sigAve_wt *= scale;
  cout << " +++++ Fit Chi2min/N = " << chi2min << "/" << NC << "\t tau = " << muAve << " +/- " << sigAve << endl;
  
  drawHistBox(x1,x2,y1,y2);
  
  c1->Update();
  c1->SaveAs(Form("fig/BL_H3L_data_%d.pdf", config));
  c1->SaveAs(Form("fig/BL_H3L_data_%d.png", config));

  /////////
  // remake binding energy figure, verticle with good lables
  /////////
  TCanvas *c11 = new TCanvas("c11","c11",600,800);
  c11->SetLeftMargin(0.05);
  c11->SetBottomMargin(0.14);
  c11->cd();
  
  x1 = XMIN;
  x2 = XMAX+0.22;
  y1 = 0;
  y2 = NP+3;
  TH1D *h01 = new TH1D("h01","",1,x1,x2);
  h01->GetXaxis()->CenterTitle();
  h01->SetXTitle("{}^{3}_{#Lambda}H B_{#Lambda} (MeV)");
  h01->GetXaxis()->SetTitleOffset(0.8);
  h01->GetYaxis()->SetNdivisions(101);
  h01->GetYaxis()->SetTitleOffset(999.);
  h01->GetYaxis()->SetLabelOffset(999.);
  h01->SetMaximum(y2);
  h01->SetMinimum(y1);
  h01->Draw();

  drawColorBox(muAve - sigAve_wt, 0, muAve + sigAve_wt, NC+1, 5, 0.5);

  gr_data_y->SetMarkerSize(1.5);
  gr_data_y->SetMarkerStyle(20);
  gr_data_y->SetLineWidth(2);
  gr_data_y->Draw("p");

  for(int i=0;i<NC;i++) {
    drawText(XMAX-0.04, i+0.9, Label[i], 42, 0.035);
  }
  
  drawText(0., y2*0.93, Form("Ave. = %4.2f #pm %4.2f (S=%3.1f)", muAve, sigAve_wt, scale), 42, 0.035);
  drawArrow(muAve, y2*0.91, muAve, y2*0.79, 0.025, 40, 2);

  drawLine(Tau_Lambda, y1, Tau_Lambda, NC+1, 2, 2, 1);
  // drawText(267, y2*0.89, "#Lambda", 42, 0.035, 0, 4);
  // drawArrow(Tau_Lambda, y2*0.90, Tau_Lambda, y2*0.84, 0.025, 40, 2, 2, 4);
  
  drawHistBox(x1,x2,y1,y2);
  
  c11->Update();
  c11->SaveAs(Form("fig/tau_H3L_lifetime_%d.pdf", config));
  c11->SaveAs(Form("fig/tau_H3L_lifetime_%d.png", config));


  TCanvas *c2 = new TCanvas("c2","c2",100, 100, 800,600);
  c2->SetLeftMargin(0.14);
  c2->cd();
  
  x1 = XMIN+0.2; //muAve - sigAve*5;
  x2 = XMAX-0.2; //muAve + sigAve*5;
  y1 = 8.; // floor(chi2min-2);
  y2 = 30; //floor(chi2min+18);
  TH1D *h1 = new TH1D("h1","",1,x1,x2);
  h1->GetXaxis()->CenterTitle();
  h1->SetXTitle("{}^{3}_{#Lambda}H B_{#Lambda} (MeV)");
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetTitleOffset(1.0);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->SetYTitle("#chi^{2}");
  h1->GetYaxis()->SetTitleOffset(0.9);
  h1->SetMaximum(y2);
  h1->SetMinimum(y1);
  h1->Draw();

  TGraph *gr_chi2 = new TGraph(Nf+1, xf, chi2);
  gr_chi2->SetName("chi2");
  gr_chi2->SetMarkerSize(1.0);
  gr_chi2->SetLineWidth(2);
  gr_chi2->Draw("pc");

  double chi2_min = 1.e99;
  double x_min = -1;
  double x_el = -1;
  double x_eh = -1;
  double x_e2l = -1; // 2sigma
  double x_e2h = -1;
  double x_e3l = -1; // 3sigma
  double x_e3h = -1;
  findChi2Min(gr_chi2, 1, &chi2_min, &x_min, &x_el, &x_eh);
  cout << " +++++ Scan Chi2min/N = " << chi2_min << "/" << NC << "\t tau = " << x_min << " -" << x_el << " +" << x_eh << endl;
  findChi2Min(gr_chi2, 2, &chi2_min, &x_min, &x_e2l, &x_e2h);
  cout << " +++++ Scan Chi2min/N = " << chi2_min << "/" << NC << "\t tau = " << x_min << " -" << x_e2l << " +" << x_e2h << endl;
  findChi2Min(gr_chi2, 3, &chi2_min, &x_min, &x_e3l, &x_e3h);
  cout << " +++++ Scan Chi2min/N = " << chi2_min << "/" << NC << "\t tau = " << x_min << " -" << x_e3l << " +" << x_e3h << endl;

  drawLine(x_min-x_el, chi2_min+1, x_min+x_eh, chi2_min+1);
  drawLine(x_min-x_e2l, chi2_min+4, x_min+x_e2h, chi2_min+4);
  drawLine(x_min-x_e3l, chi2_min+9, x_min+x_e3h, chi2_min+9);
  drawLine(x_min, y1, x_min, chi2_min, 2, 1, kGreen+4);
  drawLine(x_min-x_el, y1, x_min-x_el, chi2_min+1, 1, 2);
  drawLine(x_min+x_eh, y1, x_min+x_eh, chi2_min+1, 1, 2);
  drawLine(x_min-x_e2l, y1, x_min-x_e2l, chi2_min+4, 1, 2);
  drawLine(x_min+x_e2h, y1, x_min+x_e2h, chi2_min+4, 1, 2);
  drawLine(x_min-x_e3l, y1, x_min-x_e3l, chi2_min+9, 1, 2);
  drawLine(x_min+x_e3h, y1, x_min+x_e3h, chi2_min+9, 1, 2);

  
  drawHistBox(x1,x2,y1,y2);
  c2->Update();

  c2->SaveAs(Form("fig/BL_H3L_chi2scan_%d.pdf", config));
  c2->SaveAs(Form("fig/BL_H3L_chi2scan_%d.png", config));


  //////////////
  // ideogram
  //////////////
  TCanvas *c3 = new TCanvas("c3","c3",800, 800, 800,600);
  c3->SetLeftMargin(0.05);
  c3->SetRightMargin(0.3);
  c3->cd();
  

  double sigCut = 3.*sqrt(NC)*sigAve;
  cout << " sigCut = " << sigCut << " sigAve = " << sigAve << endl;
  TF1 *ideo = new TF1("ideogram", ideogram0, XMIN, XMAX, NP*2);
  double chi2_s[NP];
  double chi2_tot = 0.;
  for(int i=0;i<NP;i++) {
    double sig = reset_error(yp[i], eyl[i], eyh[i], muAve);
    ideo->SetParameter(i*2, yp[i]);
    ideo->SetParameter(i*2+1, sig);

    chi2_s[i] = 0;
    chi2_s[i] = pow(fabs(yp[i]-muAve)/sig, 2.0);
    //    if(yp[i]<muAve) chi2_s[i] = pow(fabs(yp[i]-muAve)/eyh[i], 2.0);
    //    else chi2_s[i] = pow(fabs(yp[i]-muAve)/eyl[i], 2.0);
    //    double sig = (eyl[i] + eyh[i])/2.0;
    cout << " sigma = " << sig << " sigCut = " << sigCut << " chi2 = " << chi2_s[i] << endl;
    if(sig>sigCut) {
      cout << " \t ==> sigma too large, throwing away data: " << Label[i] << endl;
      chi2_s[i] = 0.;
      ideo->SetParameter(i*2, 0.);
    }
  }
  if(config) {
    chi2_s[NP-1] = 0.;
    ideo->SetParameter((NP-1)*2, 0.);  // skip last ALICE data
  }
  for(int i=0;i<NP;i++) chi2_tot += chi2_s[i];

  x1 = XMIN;
  x2 = XMAX;
  y1 = 0;
  cout << "ideo maximum = " << ideo->GetMaximum() << endl;
  y2 = ideo->GetMaximum()/0.75;
  
  TH1D *h2 = new TH1D("h2","",1,x1,x2);
  h2->GetXaxis()->CenterTitle();
  h2->GetXaxis()->SetTitleOffset(1.0);
  h2->GetXaxis()->SetLabelOffset(0.01);
  h2->SetXTitle("{}^{3}_{#Lambda}H B_{#Lambda} (MeV)");
  h2->SetYTitle("");
  h2->GetYaxis()->SetNdivisions(101);
  h2->GetYaxis()->SetLabelOffset(999.);
  h2->GetYaxis()->SetTitleOffset(999.);
  h2->SetMaximum(y2);
  h2->SetMinimum(y1);
  h2->Draw();

  drawColorBox(muAve - sigAve_wt, 0, muAve + sigAve_wt, y2*0.8, 5, 0.5);
  cout << " +++++ Weighted Ave = " << muAve << " +/- " << sigAve_wt << endl;
  drawLine(0, y1, 0, y2*0.82, 2, 2, 1);

  ideo->SetLineWidth(2);
  ideo->SetTitle("");
  ideo->Draw("c same");

  double step = ideo->GetMaximum()/(NP*1.1);
  if(config) {
    drawText(XMAX+0.36, ideo->GetMaximum()*1.08, "#chi^{2}", 42, 0.035);
    drawLine(XMAX+0.1, ideo->GetMaximum()*1.03, XMAX+0.4, ideo->GetMaximum()*1.03, 2);
  } else {
    drawText(XMAX+0.35, ideo->GetMaximum()*1.1, "#chi^{2}", 42, 0.035);
    drawLine(XMAX+0.1, ideo->GetMaximum()*1.06, XMAX+0.4, ideo->GetMaximum()*1.06, 2);
  }
  int nn = 0;
  for(int i=0;i<NP;i++) {
    double x[1] = {yp[i]};
    double y[1] = {(i+1.5)*step};
    double yel[1] = {0.3*step};
    double yeh[1] = {0.3*step};
    double xel[1] = {eyl[i]};
    double xeh[1] = {eyh[i]};

    if(config && i==NP-1) continue;
    drawText(XMAX+0.1, y[0] - yel[0]*0.5, Label[i], 42, 0.035);
    TGraphAsymmErrors *gr_tmp = new TGraphAsymmErrors(1, x, y, xel, xeh, yel, yeh);
    gr_tmp->SetLineWidth(1);
    gr_tmp->SetMarkerSize(0.01);
    gr_tmp->Draw("p");

    if(fabs(chi2_s[i]>1.e-5)) {
      drawText(XMAX+0.35, y[0] - yel[0]*0.5, Form("%3.1f",chi2_s[i]), 42, 0.035);
      nn++;
    }
  }
  drawLine(XMAX+0.1, ideo->GetMaximum()*0.15, XMAX+0.4, ideo->GetMaximum()*0.15, 2);
  drawText(XMAX+0.33, ideo->GetMaximum()*0.06, Form("%4.1f",chi2_tot), 42, 0.035);
  drawText(XMAX+0.1, ideo->GetMaximum()*(-0.05), Form("C.L.             = %5.3f", TMath::Prob(chi2_tot, nn-1)), 42, 0.035);

  //  drawText(180, y2*0.93, Form("Weighted Ave = %5.1f - %3.1f + %3.1f", x_min, x_el*scale, x_eh*scale), 42, 0.03);
  drawText(0, y2*0.92, Form("Weighted Ave = %4.2f #pm %4.2f (S=%3.1f)", muAve, sigAve_wt, scale), 42, 0.035);
  drawArrow(muAve, y2*0.90, muAve, y2*0.82, 0.025, 40, 2);

  drawHistBox(x1,x2,y1,y2);
  c3->Update();
  c3->SaveAs(Form("fig/BL_H3L_ideogram_%d.pdf", config));
  c3->SaveAs(Form("fig/BL_H3L_ideogram_%d.png", config));


  TFile *fout = new TFile(Form("root/BL_data_chi2_%d.root",config),"recreate");
  fout->cd();
  gr_data->Write();
  gr_data_reset->Write();
  gr_chi2->Write();
  fout->Close();
}


