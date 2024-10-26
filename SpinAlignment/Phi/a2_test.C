#include "draw.C+"
#include "style.C+"

Double_t a2_m(Double_t *x, Double_t *par)
{
  const Double_t m_th = 0.99;
  double m = x[0];
  double fS = par[0]; // S/B ratio - averaged
  double mu = par[1];
  double sig = par[2];
  double a0 = par[3];
  double a0m = par[4];  // Bkgd a2 vs m:  a0 - a0m * m
  double a2S = par[5];   // Signal a2 vs m: constant

  double rhoS = TMath::Exp(-0.5*TMath::Power((m-mu)/sig, 2.0));
  double rhoB = TMath::Power(m-m_th, 0.5);
  double f = fS*rhoS/(rhoS+rhoB);
  double a2B = a0 + a0m*(m-mu);
  return f*a2S + (1-f)*a2B;
}


void a2_test()
{
  style();
  const Double_t MassPhi = 1.019;
  const Double_t Width = 0.002;

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Draw();

  TH1D *h0 = new TH1D("h0","",1, 0.99, 1.06);
  h0->SetMinimum(-0.08);
  h0->SetMaximum(0.0);
  h0->Draw("c");
  
  TF1 *f1 = new TF1("f1", a2_m, 0.99, 1.06, 6);
  f1->SetParameters(0.1, MassPhi, Width, -0.042, -1.3, -0.2);
  f1->Draw("same");
}
