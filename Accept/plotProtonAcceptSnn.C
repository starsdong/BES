#include "style.C+"
#include "draw.C+"
Double_t Perp(Double_t *x, Double_t *par)
{
  double y0 = par[2];
  double y = - x[0] + y0;
  double m = par[0];
  double eta = par[1];
  double a = TMath::CosH(y);
  double b = TMath::CosH(eta);
  return m*TMath::Sqrt((a*a-1)/(b*b-a*a));
}

Double_t PerpCMS(Double_t *x, Double_t *par)
{
  double y = fabs(x[0]);
  double m = par[0];
  double eta = par[1];
  double a = TMath::CosH(y);
  double b = TMath::CosH(eta);
  return m*TMath::Sqrt((a*a-1)/(b*b-a*a));
}

void plotProtonAcceptSnn()
{
  style();

  const Double_t xmax = 3.2;  // x-axis maxim
  const Double_t ymax = 1.99;  // y-axis maxim
  const Double_t mu = 0.9315;  // nucleon mass
  const Double_t mp = 0.93827;
  const Int_t NE = 3;
  const Double_t E[NE] = {3.85, 9.8, 31.2};
  Double_t ybeam[NE];
  for(int i=0;i<NE;i++) {
    double mom = sqrt(E[i]*E[i] - mu*mu);
    TLorentzVector proj(0, 0, mom, E[i]);
    TLorentzVector targ(0.,0.,0.,mu);
    TLorentzVector s2 = proj + targ;
    ybeam[i] = s2.Rapidity();
    cout << E[i] << " y = " << ybeam[i] << " s1/2 = " << s2.Mag() << endl;
  }


  const Double_t zTar = 200.;
  const Double_t z[2] = {-200., -200.};
  const Double_t r[2] = {110., 70.};
  Double_t eta_ref[2];
  Double_t eta_ref_CMS[2];
  for(int i=0;i<2;i++) {
    TVector3 a(r[i],0,z[i]-zTar);
    eta_ref[i] = a.Eta();
    cout << "eta_ref = " << eta_ref[i] << endl;    
    TVector3 aCMS(r[i],0,z[i]);
    eta_ref_CMS[i] = aCMS.Eta();
    cout << "eta_ref_CMS = " << eta_ref_CMS[i] << endl;    
  }

  

  // Plotting
  TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,800);
  c1->SetLogz();
  c1->SetGridx();
  c1->SetGridy();
  c1->Draw();

  double x1 = -2.4;
  double x2 = 1.8;
  double y1 = 0.;
  double y2 = 2.4;

  TH1D *d0 = new TH1D("d0","",1,x1,x2);
  d0->SetMinimum(y1);
  d0->SetMaximum(y2);
  d0->GetXaxis()->SetNdivisions(408);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("y_{CMS}");
  d0->GetXaxis()->SetLabelOffset(0.01);
  d0->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  d0->Draw("c");
  drawHistBox(x1,x2,y1,y2);
  
  TF1 *fun[NE];
  const Int_t kColor[NE] = {kBlack, kRed, kBlue};
  for(int i=0;i<NE;i++) {
    fun[i] = new TF1(Form("fun_%d",i),Perp,-ybeam[i],x2,3);
    fun[i]->SetParameters(mp,eta_ref[1],-ybeam[i]);
    fun[i]->SetRange(-ybeam[i],-eta_ref[1]-ybeam[i]);
    fun[i]->SetLineColor(kColor[i]);
    fun[i]->Draw("same");

    drawLine(-ybeam[i], y1, -ybeam[i], 0.6, 2, 1, kColor[i]);

  }
  // TF1 *funCMS = new TF1("fun_CMS",PerpCMS,y0/2,-y0/2,2);
  // funCMS->SetParameters(mp,eta_ref_CMS[0]);
  // funCMS->SetRange(eta_ref_CMS[0],-eta_ref_CMS[0]);
  // funCMS->SetLineStyle(4);
  // funCMS->SetLineWidth(2);
  // funCMS->Draw("same");
  

  //  drawLine(-ybeam,y1,-ybeam,y2,2,1,1);
  // drawText(0.7, 0.8, "TPC dEdx", 42, 0.045, 70, kBlack);
  // drawText(0.68, 0.35, "iTPC + eTOF", 42, 0.045, 40, kBlue);

  drawLine(-0.5,0.4,0.5,0.4,5,2,1);
  drawLine(-0.5,2.0,0.5,2.0,5,2,1);
  drawLine(-0.5,0.4,-0.5,2.0,5,2,1);
  drawLine(0.5,0.4,0.5,2.0,5,2,1);

  
  drawText(1.5, 1.7, "3.0 GeV", 52, 0.045, 90, kBlack);
  drawText(1.05, 1.7, "4.5 GeV", 52, 0.045, 90, kRed);
  drawText(0.45, 1.7, "7.7 GeV", 52, 0.045, 90, kBlue);

  drawText(-1.3, 0.65, "Target", 32, 0.04, 0, kBlack);


  c1->SaveAs("fig/acceptProtonSnn.pdf");
  c1->SaveAs("fig/acceptProtonSnn.png");

}
