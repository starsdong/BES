#include "style.C+"
#include "draw.C+"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// rho00 extracted using conventional costheta*-bin method, fit costheta* distributions with funCosTheta
// Try to apply corrections delta_rho00 =  -4/3 * a2 * v2  extracted from a2.C macro
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Rapidity dependence
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void rho00_y_v2pt()
{
  style();
  TF1 *funCosTheta = new TF1("funCosTheta","[1]*((1-[0])+(3*[0]-1)*x*x)",-1,1);
  funCosTheta->SetParameter(0, 1./3.);

  const Int_t NY = 10;
  const Int_t i_edge[NY] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  const Double_t rap[NY] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const Double_t rap1[NY] = {0.12, 0.22, 0.32, 0.42, 0.52, 0.62, 0.72, 0.82, 0.92, 1.02};
  const Int_t NV2 = 2;
  const Char_t *v2Name[NV2] = {"v2_0.0_Y_9.9", "v2pt"};
  const Char_t *v2Label[NV2] = {"v200", "v2pt"};
  const Int_t NMAX = 200;

  double_t r[NV2][NY], re[NV2][NY];
  double_t rmc[NV2][NY], remc[NV2][NY];
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->Draw();
  c1->Divide(NY, NV2);
  TFile *fin[NV2];
  TH1D *fMc[NV2][NY], *fRc[NV2][NY];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    //    fin[iv2] = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root", v2Name[iv2], sigY));
    fin[iv2] = new TFile(Form("accept_Sergei_%s.root", v2Name[iv2]));
    for(int i=0;i<NY;i++) {
      c1->cd(i+1+iv2*NY);
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hYCosTheta"))->ProjectionY(Form("Mc_%d_%d",iv2,i),1,i_edge[i]);
      fMc[iv2][i]->Rebin(fMc[iv2][i]->GetNbinsX()/NMAX);
      fMc[iv2][i]->Draw();
      funCosTheta->SetParameter(1, fMc[iv2][i]->GetMaximum());
      fMc[iv2][i]->Fit("funCosTheta","R");
      rmc[iv2][i] = funCosTheta->GetParameter(0) - 1./3.;
      remc[iv2][i] = funCosTheta->GetParError(0);

      fRc[iv2][i] = ((TH2D *)fin[iv2]->Get("hYCosThetaRc3"))->ProjectionY(Form("Rc_%d_%d",iv2,i),1,i_edge[i]);
      fRc[iv2][i]->Rebin(fRc[iv2][i]->GetNbinsX()/NMAX);
      fRc[iv2][i]->Draw();
      funCosTheta->SetParameter(1, fRc[iv2][i]->GetMaximum());
      fRc[iv2][i]->Fit("funCosTheta","R");
      r[iv2][i] = funCosTheta->GetParameter(0) - 1./3.;
      re[iv2][i] = funCosTheta->GetParError(0);
      c1->Update();
    }
  }

  double r_corr[NV2][NY], re_corr[NV2][NY];
  TFile *fCorr = new TFile(Form("root/drho_a2_y_v2pt.root"));
  TGraphErrors *gr_drho[NV2];
  for(int i=0;i<NV2;i++) {
    gr_drho[i] = (TGraphErrors *)fCorr->Get(Form("drho_%d",i));
    for(int j=0;j<NY;j++) {
      r_corr[i][j] = r[i][j] - gr_drho[i]->GetY()[j];
      re_corr[i][j] = TMath::Sqrt(re[i][j]*re[i][j] + gr_drho[i]->GetEY()[j]*gr_drho[i]->GetEY()[j]);    
    }
  }

  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,1.2,1,-0.01,0.025);
  h2->GetXaxis()->SetTitle("#phi-meson rapidity");
  h2->GetYaxis()->SetTitle("#Delta#rho_{00} from cos#theta* fit");  
  h2->Draw();
  drawLine(0.0, 0.0, 1.2, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  TGraphErrors *gr_rc[NV2];
  TGraphErrors *gr_corr[NV2];
  const Int_t markerStyle[NV2] = {24, 20};

  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NY, rap, rmc[i], 0, remc[i]);
    gr[i]->SetName(Form("rho00_Mc_CosThetaBin_%d",i));
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);

    gr_rc[i] = new TGraphErrors(NY, rap, r[i], 0, re[i]);
    gr_rc[i]->SetName(Form("rho00_CosThetaBin_%d",i));
    gr_rc[i]->SetMarkerStyle(markerStyle[i]);
    gr_rc[i]->SetMarkerSize(2.0);
    gr_rc[i]->SetMarkerColor(2);
    gr_rc[i]->SetLineWidth(2);
    gr_rc[i]->SetLineColor(2);
    
    gr_corr[i] = new TGraphErrors(NY, rap1, r_corr[i], 0, re_corr[i]);
    gr_corr[i]->SetName(Form("rho00Corr_CosThetaBin_%d",i));
    gr_corr[i]->SetMarkerStyle(markerStyle[i]);
    gr_corr[i]->SetMarkerSize(2.0);
    gr_corr[i]->SetMarkerColor(4);
    gr_corr[i]->SetLineWidth(2);
    gr_corr[i]->SetLineColor(4);
    gr_corr[i]->Draw("p");
  }


  // gr[0]->SetFillColor(5);
  // gr[0]->Draw("e3");

  for(int i=0;i<NV2;i++) {
    gr[i]->Draw("p");
    gr_rc[i]->Draw("p");	
    gr_corr[i]->Draw("p");	
  }

  TLegend *leg = new TLegend(0.2, 0.64, 0.35, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.35, 0.64, 0.5, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_rc[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.5, 0.64, 0.65, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_corr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  drawText(0.09, 0.021, "MC");
  drawText(0.3, 0.021, "RC", 42, 0.05, 0, 2);
  drawText(0.52, 0.021, "RC corr.", 42, 0.05, 0, 4);
  drawHistBox(0., 1.2, -0.01, 0.025);

  
  c2->Update();

  c2->SaveAs(Form("fig/rho00_y_v2pt.pdf"));
  c2->SaveAs(Form("fig/rho00_y_v2pt.png"));
  
  TFile *fout = new TFile(Form("root/rho00_y_v2pt.root"),"recreate");
  for(int i=0;i<NV2;i++) {
    gr[i]->Write();
    gr_rc[i]->Write();
    gr_corr[i]->Write();
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
