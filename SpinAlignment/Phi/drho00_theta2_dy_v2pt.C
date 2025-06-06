#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

/////////////////////////////////////////////////////
// delta_rho00 = 5/2 (<cos^2 theta*> - 1/3)
/////////////////////////////////////////////////////
// Rapidity dependence
/////////////////////////////////////////////////////

void drho00_theta2_dy_v2pt()
{
  style();
    
  const Int_t NMAX = 200; // number of cos2phi bins

  const Int_t NY = 10;
  const Int_t i_edge[NY] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  //  const Double_t rap[NY] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const Double_t rap[NY] = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95};
  const Int_t NV2 = 2;
  const Char_t *v2Name[NV2] = {"v2_0.0_Y_9.9", "v2pt"};
  const Char_t *v2Label[NV2] = {"v200", "v2pt"};

  double_t v2[NV2][NY], v2e[NV2][NY];
  double_t v2Rc[NV2][NY], v2eRc[NV2][NY];
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->Draw();
  c1->Divide(NY, NV2);
  TFile *fin[NV2];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    //    fin[iv2] = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root", v2Name[iv2], sigY));
    fin[iv2] = new TFile(Form("accept_Sergei_%s.root", v2Name[iv2]));
    TH1D *fMc[NV2][NY], *fRc[NV2][NY];
    TGraphErrors *gr_e[NV2][NY];
    for(int i=0;i<NY;i++) {
      c1->cd(i+1+iv2*NY);
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hYCosTheta2"))->ProjectionY(Form("Mc_%d_%d",iv2, i),i_edge[i]-1,i_edge[i]);
      fRc[iv2][i] = ((TH2D *)fin[iv2]->Get("hYCosTheta2Rc3"))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]-1,i_edge[i]);


      TH1D *h0 = new TH1D("h0","",1,0,1);
      h0->SetMaximum(fRc[iv2][i]->GetMaximum()/0.7);
      h0->SetMinimum(0);
      h0->Draw();

      fRc[iv2][i]->Draw("same");

      v2[iv2][i] = 2.5 * (fMc[iv2][i]->GetMean() - 1./3.);
      v2e[iv2][i] = 2.5 * fMc[iv2][i]->GetMeanError();
      v2Rc[iv2][i] = 2.5 * (fRc[iv2][i]->GetMean() - 1./3.);
      v2eRc[iv2][i] = 2.5 * fRc[iv2][i]->GetMeanError();
      c1->Update();
    }
    //    fin[iv2]->Close();
  }

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,1.0,1,-0.05,0.2);
  h2->GetXaxis()->SetTitle("#phi-meson rapidity");
  h2->GetYaxis()->SetTitle("#Delta#rho_{00} = 5/2 * (<cos^{2}#theta*> - 1/3) ");  
  h2->Draw();
  drawLine(0.0, 0.0, 1.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  const Int_t markerStyle[NV2] = {24, 20};
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NY, rap, v2[i], 0, v2e[i]);
    gr[i]->SetName(Form("drho00_CosTheta2_Mc_%d",i));
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }

  TGraphErrors *gr_rc[NV2];
  for(int i=0;i<NV2;i++) {
    gr_rc[i] = new TGraphErrors(NY, rap, v2Rc[i], 0, v2eRc[i]);
    gr_rc[i]->SetName(Form("drho00_CosTheta2_Rc_%d",i));
    gr_rc[i]->SetMarkerStyle(markerStyle[i]);
    gr_rc[i]->SetMarkerSize(2.0);
    gr_rc[i]->SetMarkerColor(2);
    gr_rc[i]->SetLineWidth(2);
    gr_rc[i]->SetLineColor(2);
    gr_rc[i]->Draw("p");
  }
  TLegend *leg = new TLegend(0.25, 0.64, 0.4, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.4, 0.64, 0.55, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_rc[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  drawText(0.12, 0.175, "MC");
  drawText(0.34, 0.175, "RC", 42, 0.05, 0, 2);
  drawHistBox(0., 1.0, -0.05, 0.2);

  c2->Update();

  c2->SaveAs(Form("fig/drho00_theta2_dy_v2pt.pdf"));
  c2->SaveAs(Form("fig/drho00_theta2_dy_v2pt.png"));
  
  TFile *fout = new TFile(Form("root/drho00_CosTheta2_dy_v2pt.root"),"recreate");
  for(int i=0;i<NV2;i++) {
    gr[i]->Write();
    gr_rc[i]->Write();
  }
  fout->Close();
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
