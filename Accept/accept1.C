Double_t Perp(Double_t *x, Double_t *par)
{
  double y = x[0];
  double m = par[0];
  double eta = par[1];
  double a = TMath::CosH(y);
  double b = TMath::CosH(eta);
  return m*TMath::Sqrt((a*a-1)/(b*b-a*a));
}

void accept1()
{
  gROOT->Reset();

  const Int_t NE = 8;
  const Double_t E[NE] = {100, 70, 44, 31.2, 9.8, 7.3, 4.55, 3.85};
  const Double_t m = 0.9315;
  TLorentzVector targ(0,0,0,m);

  Double_t E_cms[NE];
  Double_t y0[NE];
  for(int i=0;i<NE;i++) {
    double p = sqrt(E[i]*E[i]-m*m);
    TLorentzVector proj(0,0,-p,E[i]);
    TLorentzVector s = targ + proj;
    E_cms[i] = s.Mag();
    y0[i] = s.Rapidity();
    cout << "proj E = " << E[i] << "\t sqrt(s) = " << E_cms[i] << "\t y0 = " << y0[i] << endl;
  }

  const Double_t zTar = 200.;
  const Double_t z[4] = {200., 0., -200., -200.};
  const Double_t r[4] = {200., 200., 200., 125.};
  Double_t eta_ref[4];
  for(int i=0;i<4;i++) {
    TVector3 a(r[i],0,z[i]-zTar);
    eta_ref[i] = a.Eta();
    cout << "eta_ref = " << eta_ref[i] << endl;    
  }


  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);

   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);

   // c1->SetGridx();
   // c1->SetGridy();

   c1->SetLeftMargin(0.13);
   c1->SetBottomMargin(0.13);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.02);

   double x1 = -2.4;
   double x2 = 0.1;
   double y1 = 0.0;
   double y2 = 3.0;
   
   TH1D *d0 = new TH1D("d0","",1,x1,x2);
   d0->SetMinimum(y1);
   d0->SetMaximum(y2);
   d0->GetXaxis()->SetNdivisions(208);
   d0->GetXaxis()->SetTitle("y_{lab}");
   d0->GetXaxis()->SetTitleOffset(0.9);
   d0->GetXaxis()->SetTitleSize(0.06);
   d0->GetXaxis()->SetLabelOffset(0.01);
   d0->GetXaxis()->SetLabelSize(0.045);
   d0->GetXaxis()->SetLabelFont(42);
   d0->GetXaxis()->SetTitleFont(42);
   d0->GetYaxis()->SetNdivisions(205);
   d0->GetYaxis()->SetTitle("p_{T} (GeV/c)");
   d0->GetYaxis()->SetTitleOffset(1.0);
   d0->GetYaxis()->SetTitleSize(0.06);
   d0->GetYaxis()->SetLabelOffset(0.005);
   d0->GetYaxis()->SetLabelSize(0.045);
   d0->GetYaxis()->SetLabelFont(42);
   d0->GetYaxis()->SetTitleFont(42);
   d0->Draw();


   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   TF1 *fun[4];
   for(int i=0;i<4;i++) {
     fun[i] = new TF1(Form("fun_%d",i),Perp,x1,x2,2);
     fun[i]->SetParameters(m,eta_ref[i]);
     fun[i]->SetRange(eta_ref[i],0.0);
     
     fun[i]->Draw("same");

     TLatex *tex = new TLatex(eta_ref[i]+0.0, 2.4, Form("#eta = %4.2f",eta_ref[i]));
     tex->SetTextFont(32);
     tex->SetTextSize(0.045);
     tex->SetTextAngle(90);
     tex->Draw("same");
   }
   TLine *l0 = new TLine(0,y1,0,y2);
   l0->SetLineColor(2);
   l0->SetLineWidth(2);
   l0->Draw("same");
   
   for(int i=0;i<NE;i++) {
     TLine *la = new TLine(y0[i], y1, y0[i], y1+0.5);
     la->SetLineWidth(2);
     la->SetLineStyle(2);
     la->Draw("same");

     TLatex *tex = new TLatex(y0[i]-0.05, y1+0.55+0.07*TMath::Power(-1.,i), Form("%3.1f",E_cms[i]));
     tex->SetTextFont(32);
     tex->SetTextSize(0.045);
     tex->Draw("same");
   }
   TLatex *tex = new TLatex(-2.3, 0.4, "#sqrt{s_{NN}}");
   tex->SetTextFont(32);
   tex->SetTextSize(0.05);
   tex->Draw("same");

   c1->Update();
   c1->SaveAs("fig/accept1.pdf");
   c1->SaveAs("fig/accept1.png");
}
