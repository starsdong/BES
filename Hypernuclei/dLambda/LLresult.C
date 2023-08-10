#include "style.C+"
#include "draw.C+"
#define InvF0 1

void LLresult()
{
  style();

  const Int_t NS = 2;  // 2-states
  // Model p-Lambda f0-d0 parameters
  const Double_t f0_pL_m[NS] = {2.88, 1.67};
  const Double_t d0_pL_m[NS] = {2.88, 3.74};
  // Model d-Lambda f0-d0 parameters
  const Double_t f0_dL_m[NS] = {-16.3, 17.3};
  const Double_t d0_dL_m[NS] = {3.2, 3.6};

  const Int_t ND_pL = 1; // number of data points
  const Double_t f0_pL_d[ND_pL] = {2.20};
  const Double_t d0_pL_d[ND_pL] = {3.12};
  const Double_t f0e_pL_d[ND_pL] = {0.01};
  const Double_t d0e_pL_d[ND_pL] = {1.00};
  
  const Int_t ND_dL = 2; // number of data points
  const Double_t f0_dL_d[ND_dL] = {-23.9, 19.2};
  const Double_t d0_dL_d[ND_dL] = {3.2, 3.6};
  const Double_t f0e_dL_d[ND_dL] = {4.7, 2.4};
  const Double_t d0e_dL_d[ND_dL] = {0.01, 0.01};

  // plot inv_f0
  Double_t f0i_pL_m[NS], f0i_dL_m[NS];
  Double_t f0i_pL_d[ND_pL], f0ie_pL_d[ND_pL], f0i_dL_d[ND_dL], f0ie_dL_d[ND_dL];
  for(int i=0;i<NS;i++) {
    f0i_pL_m[i] = 1./f0_pL_m[i];
    f0i_dL_m[i] = 1./f0_dL_m[i];
  }
  for(int i=0;i<ND_pL;i++) {
    f0i_pL_d[i] = 1./f0_pL_d[i];
    f0ie_pL_d[i] = fabs(1./f0_pL_d[i]/f0_pL_d[i])*f0e_pL_d[i];
  }
  for(int i=0;i<ND_dL;i++) {
    f0i_dL_d[i] = 1./f0_dL_d[i];
    f0ie_dL_d[i] = fabs(1./f0_dL_d[i]/f0_dL_d[i])*f0e_dL_d[i];
  }
  
  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->Draw();

#ifdef InvF0
  double x1 = -0.2;
  double x2 = 0.8;
  double y1 = 0.;
  double y2 = 6.;
#else
  double x1 = -30;
  double x2 = 30;
  double y1 = 0.;
  double y2 = 6.;
#endif

  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
#ifdef InvF0
  h0->SetXTitle("1/f_{0} (fm^{-1})");
#else
  h0->SetXTitle("f_{0} (fm)");
#endif  
  h0->GetXaxis()->SetTitleOffset(0.9);
  h0->SetYTitle("d_{0} (fm)");
  h0->GetYaxis()->SetTitleOffset(0.8);
  h0->Draw();

#ifdef InvF0
  TGraph *gr_pL_m = new TGraph(NS, f0i_pL_m, d0_pL_m);
#else
  TGraph *gr_pL_m = new TGraph(NS, f0_pL_m, d0_pL_m);
#endif
  gr_pL_m->SetMarkerStyle(25);
  gr_pL_m->SetMarkerSize(2.0);
  gr_pL_m->Draw("p");

#ifdef InvF0
  TGraph *gr_dL_m = new TGraph(NS, f0i_dL_m, d0_dL_m);
#else
  TGraph *gr_dL_m = new TGraph(NS, f0_dL_m, d0_dL_m);
#endif
  gr_dL_m->SetMarkerStyle(26);
  gr_dL_m->SetMarkerSize(2.0);
  gr_dL_m->Draw("p");

#ifdef InvF0
  TGraphErrors *gr_pL_d = new TGraphErrors(ND_pL, f0i_pL_d, d0_pL_d, f0ie_pL_d, d0e_pL_d);
#else
  TGraphErrors *gr_pL_d = new TGraphErrors(ND_pL, f0_pL_d, d0_pL_d, f0e_pL_d, d0e_pL_d);
#endif
  gr_pL_d->SetMarkerStyle(21);
  gr_pL_d->SetMarkerSize(2.0);
  gr_pL_d->SetLineWidth(2);
  gr_pL_d->Draw("p");

#ifdef InvF0
  TGraphErrors *gr_dL_d = new TGraphErrors(ND_dL, f0i_dL_d, d0_dL_d, f0ie_dL_d, d0e_dL_d);
#else
  TGraphErrors *gr_dL_d = new TGraphErrors(ND_dL, f0_dL_d, d0_dL_d, f0e_dL_d, d0e_dL_d);
#endif
  gr_dL_d->SetMarkerStyle(22);
  gr_dL_d->SetMarkerSize(2.0);
  gr_dL_d->SetLineWidth(2);
  gr_dL_d->Draw("p");

  drawLine(0,y1,0,y2,1,2);

  
  double y_off = 0.3;
#ifdef InvF0
  double x_off = 0.01;
  drawText(f0i_pL_m[0]+x_off, d0_pL_m[0]+y_off, "p-#Lambda(s)", 42, 0.05, 90);
  drawText(f0i_pL_m[1]+x_off, d0_pL_m[1]+y_off, "p-#Lambda(t)", 42, 0.05, 90);
  drawText(f0i_dL_m[0]+x_off, d0_dL_m[0]+y_off, "d-#Lambda(D)", 42, 0.05, 90);
  drawText(f0i_dL_m[1]+x_off, d0_dL_m[1]+y_off, "d-#Lambda(Q)", 42, 0.05, 90);
#else
  double x_off = 1.;
  drawText(f0_pL_m[0]+x_off+1., d0_pL_m[0]+y_off, "p-#Lambda(s)", 42, 0.05, 75);
  drawText(f0_pL_m[1]+x_off-1., d0_pL_m[1]+y_off, "p-#Lambda(t)", 42, 0.05, 75);
  drawText(f0_dL_m[0]+x_off-3., d0_dL_m[0]+y_off, "d-#Lambda(D)", 42, 0.05, 75);
  drawText(f0_dL_m[1]+x_off+0., d0_dL_m[1]+y_off, "d-#Lambda(Q)", 42, 0.05, 75);
#endif

#ifdef InvF0
  drawText(0.4, 1.3, "Model Data", 42, 0.04);
  TLegend *leg = new TLegend(0.7, 0.18, 0.75, 0.3);
  TLegend *leg1 = new TLegend(0.75, 0.18, 0.95, 0.3);
#else
  drawText(7., 1.3, "Model Data", 42, 0.04);
  TLegend *leg = new TLegend(0.7, 0.18, 0.75, 0.3);
  TLegend *leg1 = new TLegend(0.75, 0.18, 0.95, 0.3);
#endif
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineColor(10);
  leg->SetTextSize(0.001);
  leg->AddEntry(gr_pL_m, "", "p");
  leg->AddEntry(gr_dL_m, "", "p");
  leg->Draw();
  
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetLineColor(10);
  leg1->SetTextSize(0.05);
  leg1->AddEntry(gr_pL_d, "  p-#Lambda", "p");
  leg1->AddEntry(gr_dL_d, "  d-#Lambda", "p");
  leg1->Draw();
  
  drawHistBox(x1,x2,y1,y2);

#ifdef InvF0
  c1->SaveAs("fig/dL_LL_inv.pdf");  
  c1->SaveAs("fig/dL_LL_inv.png");
#else
  c1->SaveAs("fig/dL_LL.pdf");  
  c1->SaveAs("fig/dL_LL.png");
#endif
  
}
