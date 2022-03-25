#include "style.C+"
#include "draw.C+"
Double_t Perp(Double_t *x, Double_t *par)
{
  const double y0 = -2.274;
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

void plotProtonAccept()
{
  style();

  const Double_t xmax = 3.2;  // x-axis maxim
  const Double_t ymax = 1.99;  // y-axis maxim
  const Double_t ybeam = 2.274; // beam rapidity
  const double y0 = -ybeam;
  const Double_t m = 0.9315;
  const Double_t mp = 0.93827;

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

  double x1 = -2.8;
  double x2 = 1.2;
  double y1 = 0.;
  double y2 = 2.4;

  TH1D *d0 = new TH1D("d0","",1,x1,x2);
  d0->SetMinimum(y1);
  d0->SetMaximum(y2);
  d0->GetXaxis()->SetNdivisions(104);
  d0->GetXaxis()->CenterTitle();
  d0->GetXaxis()->SetTitle("y_{CMS}");
  d0->GetXaxis()->SetLabelOffset(0.01);
  d0->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  d0->Draw("c");
  drawHistBox(x1,x2,y1,y2);
  
  TF1 *fun[2];
  const Int_t kColor[2] = {kBlack, kBlue};
  for(int i=0;i<2;i++) {
    if(i==1) continue;
    fun[i] = new TF1(Form("fun_%d",i),Perp,y0,x2,2);
    fun[i]->SetParameters(mp,eta_ref[i]);
    //    fun[i]->SetRange(y0,-eta_ref[i]+y0);
    fun[i]->SetRange(y0,-eta_ref[i]+y0);
    fun[i]->SetLineColor(kColor[i]);
    fun[i]->Draw("same");
  }
  TF1 *funCMS = new TF1("fun_CMS",PerpCMS,y0/2,-y0/2,2);
  funCMS->SetParameters(mp,eta_ref_CMS[0]);
  funCMS->SetRange(eta_ref_CMS[0],-eta_ref_CMS[0]);
  funCMS->SetLineStyle(4);
  funCMS->SetLineWidth(2);
  funCMS->Draw("same");
  

  drawLine(y0,y1,y0,y2,2,1,1);
  // drawText(0.7, 0.8, "TPC dEdx", 42, 0.045, 70, kBlack);
  // drawText(0.68, 0.35, "iTPC + eTOF", 42, 0.045, 40, kBlue);

  // drawLine(-0.5,0.4,0.5,0.4,3,2,1);
  // drawLine(-0.5,2.0,0.5,2.0,3,2,1);
  // drawLine(-0.5,0.4,-0.5,2.0,3,2,1);
  // drawLine(0.5,0.4,0.5,2.0,3,2,1);
  
  drawText(0.8, 0.15, " p", 52, 0.065, 0, kBlack);

  TPad *p6 = new TPad("p6","",0.33,0.91,0.80,0.979);
  p6->Draw();
  p6->cd();
  drawHistBox(0.001,0.999,0.001,0.999,3333);
  //  drawText(-0.7, 2.15, " Au+Au @ #sqrt{s_{NN}} = 3.0 GeV", 22, 0.045, 0, kBlack);
  drawText(0.1, 0.3, "Au+Au @ #sqrt{s_{NN}} = 9.2 GeV",42,0.6);
  p6->Modified();
  c1->Update();
  c1->cd();


  c1->SaveAs("fig/acceptProton.pdf");
  c1->SaveAs("fig/acceptProton.png");
}
