#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// a2 = < cos(2*(phi* - phi)) >  - phi*:  kaon momentum in phi rest frame, phi:  phi meson angle
// then calculate correction to rho_00:  delta_rho00 =  -4/3 * a2 * v2 
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Rapidity dependence - differential
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void a2_dy(const Double_t sigY = 9.9)
{
  style();
    
  TF1 *funA2 = new TF1("funA2","[1]*(1.+2.*[0]*x)",-1,1);
  funA2->SetParameter(0, 0.0);
  const Int_t NMAX = 200; // number of cos2phi bins

  const Int_t NY = 10;
  const Int_t i_edge[NY] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  //  const Double_t rap[NY] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const Double_t rap[NY] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
  const Int_t NV2 = 4;
  const Double_t v2[NV2] = {0.0, 0.1, 0.2, 0.3};
  const Double_t sc = -4./3.;  

  double_t a2[NV2][NY], a2e[NV2][NY];
  double_t drho[NV2][NY], drhoe[NV2][NY];
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->Draw();
  c1->Divide(NY,NV2);
  //  TFile *fin = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root",v2, sigY));
  TFile *fin[NV2];
  TGraphErrors *gr_e[NV2][NY];
  TH1D *fMc[NV2][NY], *fRc[NV2][NY];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    fin[iv2] = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root", v2[iv2], sigY));
    for(int i=0;i<NY;i++) {
      c1->cd(i+1+iv2*NY);
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hYCos2Phi"))->ProjectionY(Form("Mc_%d_%d",iv2, i),i_edge[i]-1,i_edge[i]); // differential
      fMc[iv2][i]->Rebin(20);
      fRc[iv2][i] = ((TH2D *)fin[iv2]->Get("hYCos2PhiRc3"))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]-1,i_edge[i]);
      fRc[iv2][i]->Rebin(20);

      double xx[NMAX];
      double eff[NMAX], err[NMAX];
      for(int j=0;j<fMc[iv2][i]->GetNbinsX();j++) {
	xx[j] = fMc[iv2][i]->GetBinCenter(j+1);
	double n = fMc[iv2][i]->GetBinContent(j+1);
	double m = fRc[iv2][i]->GetBinContent(j+1);
	eff_err(m, n, &eff[j], &err[j]);
      }
      gr_e[iv2][i] = new TGraphErrors(fMc[iv2][i]->GetNbinsX(), xx, eff, 0, err);
      gr_e[iv2][i]->SetName(Form("a2_eff_%d_%d", iv2, i));

      TH1D *h0 = new TH1D("h0","",1,-1,1);
      h0->SetMaximum(eff[0]/0.7);
      h0->GetXaxis()->CenterTitle();
      h0->GetXaxis()->SetTitle("< cos(2*(#phi* - #phi)) >");
      h0->GetYaxis()->SetTitle("Efficiency*Acceptance");  
      h0->Draw();

      gr_e[iv2][i]->Draw("p");
      gr_e[iv2][i]->Fit("funA2","R");
      a2[iv2][i] = funA2->GetParameter(0);
      a2e[iv2][i] = funA2->GetParError(0);

      drawText(-0.4, eff[0]/0.7*0.9, Form("v_{2} = %3.1f", v2[iv2]), 42, 0.08);
      drawText(-0.4, eff[0]/0.7*0.8, Form("|y| = %3.1f", rap[i]), 42, 0.08);
      drawText(-0.4, eff[0]*0.15, Form("a_{2} = %7.4f", a2[iv2][i]), 42, 0.08);
      

      drho[iv2][i] = sc * a2[iv2][i] * v2[iv2];
      drhoe[iv2][i] = sc * a2e[iv2][i] * v2[iv2];
      c1->Update();
    }
    //    fin[iv2]->Close();
  }
  c1->Update();
  c1->SaveAs(Form("fig/a2_dy_fit_Y_%3.1f.pdf", sigY));
  c1->SaveAs(Form("fig/a2_dy_fit_Y_%3.1f.png", sigY));  

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,1.0,1,-0.3, 0.05);
  h2->GetXaxis()->SetTitle("#phi-meson rapidity");
  h2->GetYaxis()->SetTitle("Eff. Slope a_{2}");
  h2->Draw();
  drawLine(0.0, 0.0, 1.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  const Int_t markerStyle[NV2] = {24, 20, 21, 22};
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NY, rap, a2[i], 0, a2e[i]);
    gr[i]->SetName(Form("a2_%d",i));
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }

  TLegend *leg = new TLegend(0.26, 0.24, 0.5, 0.54);
  leg->SetTextSize(0.05);
  for(int i=0;i<NV2;i++) {
    leg->AddEntry(gr[i], Form(" v_{2} = %3.1f", v2[i]));
  }
  leg->Draw();
  drawHistBox(0., 1.0, -0.3, 0.05);
  c2->Update();

  c2->SaveAs(Form("fig/a2_dy_Y_%3.1f.pdf", sigY));
  c2->SaveAs(Form("fig/a2_dy_Y_%3.1f.png", sigY));
  
  TCanvas *c3 = new TCanvas("c3","",800,0,800,600);
  c3->Draw();
  TH2D *h3 = new TH2D("h3","",1,0.0,1.0,1,-0.02, 0.1);
  h3->GetXaxis()->SetTitle("#phi-meson rapidity");
  h3->GetYaxis()->SetTitle("Corretion #equiv -4/3*a_{2}*v_{2}");
  h3->Draw();
  drawLine(0.0, 0.0, 1.0, 0.0, 2, 8, 1);

  TGraphErrors *gr_drho[NV2];
  for(int i=0;i<NV2;i++) {
    gr_drho[i] = new TGraphErrors(NY, rap, drho[i], 0, drhoe[i]);
    gr_drho[i]->SetName(Form("drho_%d",i));
    gr_drho[i]->SetMarkerStyle(markerStyle[i]);
    gr_drho[i]->SetMarkerSize(2.0);
    gr_drho[i]->SetLineWidth(2);
    gr_drho[i]->Draw("p");
  }
  leg = new TLegend(0.26, 0.64, 0.5, 0.94);
  leg->SetTextSize(0.05);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr[i], Form(" v_{2} = %3.1f", v2[i]));
  }
  leg->Draw();
  drawHistBox(0., 1.0, -0.02, 0.1);
  c3->Update();

  TFile *fout = new TFile(Form("root/drho_a2_dy_Y_%3.1f.root", sigY),"recreate");
  for(int i=0;i<NV2;i++) {
    gr[i]->Write();
    gr_drho[i]->Write();
    for(int j=0;j<NY;j++) {
      gr_e[i][j]->Write();
    }
  }
  fout->Close();

  c3->SaveAs(Form("fig/drho_a2_dy_Y_%3.1f.pdf", sigY));
  c3->SaveAs(Form("fig/drho_a2_dy_Y_%3.1f.png", sigY));

  /*
  TH1D *fMc = ((TH2D *)fin->Get("hPtCosTheta"))->ProjectionY("Mc",i1,i2);
  funCosTheta->SetParameter(1, fMc->GetMaximum());
  fMc->Fit("funCosTheta","R");
  c1->cd(2);
  TH1D *fRc = ((TH2D *)fin->Get("hPtCosThetaRc"))->ProjectionY("Rc",i1,i2);
  fRc->Fit("funCosTheta","R");
  c1->cd(3);
  TH1D *fRc1 = ((TH2D *)fin->Get("hPtCosThetaRc1"))->ProjectionY("Rc1",i1,i2);
  fRc1->Fit("funCosTheta","R");
  c1->cd(4);
  TH1D *fRc2 = ((TH2D *)fin->Get("hPtCosThetaRc2"))->ProjectionY("Rc2",i1,i2);
  fRc2->Fit("funCosTheta","R");
  c1->cd(5);
  TH1D *fRc3 = ((TH2D *)fin->Get("hPtCosThetaRc3"))->ProjectionY("Rc3",i1,i2);
  fRc3->Fit("funCosTheta","R");
    */

    /*
  c1->Update();
  c1->SaveAs(Form("fig/accept_Sergei_v2_%3.1f_Y_%3.1f.pdf",v2, sigY));
  c1->SaveAs(Form("fig/accept_Sergei_v2_%3.1f_Y_%3.1f.png",v2, sigY));
    */
}
