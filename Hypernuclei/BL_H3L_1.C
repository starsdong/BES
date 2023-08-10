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

void BL_H3L_1(int config=0){
  style();

  const Int_t NP = 6;
  Int_t NC = NP;
  if(config) NC = NP-1;  // config=0:  all points   config=1:  exclude the last ALICE data
  const Double_t XMIN = -0.30;
  const Double_t XMAX = 0.70;
  
  const Double_t Tau_Lambda = 0;
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
  // lifetime plot, perform fit
  //  
  TCanvas *c1 = new TCanvas("c1","c1",1000,600);
  c1->SetLeftMargin(0.1);
  c1->SetRightMargin(0.3);
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
  TCanvas *c11 = new TCanvas("c11","c11",1000,800);
  
  c11->SetLeftMargin(0.05);
  c11->SetBottomMargin(0.14);
  c11->SetRightMargin(0.3);
  c11->Draw();
  c11->cd();
  
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
  //  h2->GetXaxis()->SetNdivisions(104);
  h2->GetXaxis()->SetLabelOffset(0.008);
  h2->GetXaxis()->SetTitleOffset(1.1);
  h2->GetXaxis()->SetLabelSize(0.04);
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->SetXTitle("{}^{3}_{#Lambda}H B_{#Lambda} (MeV)");
  h2->SetYTitle("");
  h2->GetYaxis()->SetNdivisions(101);
  h2->GetYaxis()->SetLabelOffset(999.);
  h2->GetYaxis()->SetTitleOffset(999.);
  h2->SetMaximum(y2);
  h2->SetMinimum(y1);
  h2->Draw();

  drawColorBox(muAve - sigAve_wt, 0, muAve + sigAve_wt, y2*0.82, 5, 0.5);
  drawLine(0, y1, 0, y2*0.82, 2, 2, 1);
  cout << " +++++ Weighted Ave = " << muAve << " +/- " << sigAve_wt << endl;

  ideo->SetLineWidth(2);
  ideo->SetTitle("");
  ideo->Draw("c same");

  double step = ideo->GetMaximum()/(NP*1.1);
  if(config) {
    drawText(XMAX+0.37, ideo->GetMaximum()*1.08, "#chi^{2}", 42, 0.035);
    drawLine(XMAX+0.1, ideo->GetMaximum()*1.03, XMAX+0.4, ideo->GetMaximum()*1.03, 2);
  } else {
    drawText(XMAX+0.37, ideo->GetMaximum()*1.1, "#chi^{2}", 42, 0.035);
    drawLine(XMAX+0.1, ideo->GetMaximum()*1.06, XMAX+0.4, ideo->GetMaximum()*1.06, 2);
  }
  
  int nn = 0;
  for(int i=0;i<NP;i++) {
    double x[1] = {yp[i]};
    double y[1] = {(i+1.5)*step};
    double yel[1] = {0.03*step};
    double yeh[1] = {0.03*step};
    double xel[1] = {eyl[i]};
    double xeh[1] = {eyh[i]};

    if(config && i==NP-1) continue;
    drawText(XMAX+0.1, y[0] - yel[0], Label[i], 42, 0.035);
    TGraphAsymmErrors *gr_tmp = new TGraphAsymmErrors(1, x, y, xel, xeh, yel, yeh);
    gr_tmp->SetMarkerSize(1.5);
    gr_tmp->SetMarkerStyle(20);
    gr_tmp->SetLineWidth(2);
    gr_tmp->Draw("p");

    if(fabs(chi2_s[i]>1.e-5)) {
      drawText(XMAX+0.35, y[0] - yel[0], Form("%3.1f",chi2_s[i]), 42, 0.035);
      nn++;
    }
  }
  drawLine(XMAX+0.1, ideo->GetMaximum()*0.15, XMAX+0.4, ideo->GetMaximum()*0.15, 2);
  drawText(XMAX+0.33, ideo->GetMaximum()*0.06, Form("%4.1f",chi2_tot), 42, 0.035);
  drawText(XMAX+0.1, ideo->GetMaximum()*(-0.05), Form("C.L.          = %5.3f", TMath::Prob(chi2_tot, nn-1)), 42, 0.035);

  //  drawText(180, y2*0.93, Form("Weighted Ave = %5.1f - %3.1f + %3.1f", x_min, x_el*scale, x_eh*scale), 42, 0.03);
  drawText(0, y2*0.92, Form("Weighted Ave = %4.2f #pm %4.2f (S=%3.1f)", muAve, sigAve_wt, scale), 42, 0.035);
  drawArrow(muAve, y2*0.90, muAve, y2*0.82, 0.025, 40, 2);


  drawHistBox(x1,x2,y1,y2);
  
  c11->Update();
  c11->SaveAs(Form("fig/BL_H3L_1_%d.pdf", config));
  c11->SaveAs(Form("fig/BL_H3L_1_%d.png", config));


}


