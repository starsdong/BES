#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

//void rho00(const Double_t v2 = 0.1, const Double_t sigY = 5.0)
void drho00_v2pt(const Char_t *RcName = "Rc4")
{
  style();
    
  const Int_t NMAX = 200; // number of cos2phi bins

  // const Float_t i1 = 25;
  // const Float_t i2 = 48;  // 1.2 - 2.4 GeV/c
  const Int_t NPt = 12;
  const Int_t i_edge[NPt+1] = {12, 18, 24, 30, 36, 42, 48, 54, 60, 72, 84, 96, 108};
  const Double_t pT[NPt] = {0.75, 1.05, 1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.3, 3.9, 4.5, 5.1};
  const Double_t pT1[NPt] = {0.78, 1.08, 1.38, 1.68, 1.98, 2.28, 2.58, 2.88, 3.33, 3.93, 4.53, 5.13};
  // const Int_t NPt = 6;
  // const Int_t i_edge[NPt+1] = {12, 24, 36, 48, 60, 84, 108};
  // const Double_t pT[NPt] = {0.9, 1.5, 2.1, 2.7, 3.6, 4.8};
  // const Double_t pT1[NPt] = {0.95, 1.55, 2.15, 2.75, 3.65, 4.85};
  const Int_t NV2 = 2;
  const Char_t *v2Name[NV2] = {"v2_0.0_Y_9.9", "v2pt_Minv"};
  const Char_t *v2Label[NV2] = {"v200", "v2pt"};
  const Double_t sc = -4./3.;

  double_t v2[NV2][NPt], v2e[NV2][NPt];
  double_t v2Rc[NV2][NPt], v2eRc[NV2][NPt];
  TCanvas *c1 = new TCanvas("c1","",1600,900);
  c1->Draw();
  c1->Divide(6,3);
  //  TFile *fin = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root",v2, sigY));
  TFile *fin[NV2];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    fin[iv2] = new TFile(Form("accept_Sergei_%s.root", v2Name[iv2]));
    TH1D *fMc[NV2][NPt], *fRc[NV2][NPt];
    TGraphErrors *gr_e[NV2][NPt];
    for(int i=0;i<NPt;i++) {
      c1->cd(i+1+iv2*NPt);
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hPtCos2PhiRP"))->ProjectionY(Form("Mc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      if(iv2==0) {
	fRc[iv2][i] = ((TH2D *)fin[iv2]->Get(Form("hPtCos2PhiRPRc3")))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      } else {
	fRc[iv2][i] = ((TH2D *)fin[iv2]->Get(Form("hPtCos2PhiRP%s",RcName)))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      }

      TH1D *h0 = new TH1D("h0","",1,-1,1);
      h0->SetMaximum(fRc[iv2][i]->GetMaximum()/0.7);
      h0->SetMinimum(0);
      h0->Draw();

      fMc[iv2][i]->SetLineWidth(2);
      fMc[iv2][i]->Draw("esame");

      v2[iv2][i] = fMc[iv2][i]->GetMean() * sc;
      v2e[iv2][i] = fMc[iv2][i]->GetMeanError() * sc;
      v2Rc[iv2][i] = fRc[iv2][i]->GetMean() * sc;
      v2eRc[iv2][i] = fRc[iv2][i]->GetMeanError() * sc;
      c1->Update();
    }
    fin[iv2]->Close();
  }

  double r_corr[NV2][NPt], re_corr[NV2][NPt];
  TFile *fCorr = new TFile(Form("root/drho_a2_pT_v2pt_%s.root",RcName));
  TGraphErrors *gr_drho[NV2];
  for(int i=0;i<NV2;i++) {
    gr_drho[i] = (TGraphErrors *)fCorr->Get(Form("drho_%d",i));
    for(int j=0;j<NPt;j++) {
      r_corr[i][j] = v2Rc[i][j] - gr_drho[i]->GetY()[j];
      re_corr[i][j] = TMath::Sqrt(v2eRc[i][j]*v2eRc[i][j] + gr_drho[i]->GetEY()[j]*gr_drho[i]->GetEY()[j]);    
    }
  }
  
  
  TCanvas *c2 = new TCanvas("c2","", 800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,5.0,1,-0.01, 0.04);
  h2->GetXaxis()->SetTitle("#phi-meson p_{T} (GeV/c)");
  h2->GetYaxis()->SetTitle("#Delta#rho_{00} = -4/3 * <cos[2(#phi* - #Psi)]>");  
  h2->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  const Int_t markerStyle[NV2] = {24, 20};
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NPt, pT, v2[i], 0, v2e[i]);
    gr[i]->SetName(Form("drho00_Cos2PhiRP_Mc_%d",i));
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }

  TGraphErrors *gr_rc[NV2];
  for(int i=0;i<NV2;i++) {
    gr_rc[i] = new TGraphErrors(NPt, pT1, v2Rc[i], 0, v2eRc[i]);
    gr_rc[i]->SetName(Form("drho00_Cos2PhiRP_Rc_%d",i));
    gr_rc[i]->SetMarkerStyle(markerStyle[i]);
    gr_rc[i]->SetMarkerSize(2.0);
    gr_rc[i]->SetMarkerColor(2);
    gr_rc[i]->SetLineWidth(2);
    gr_rc[i]->SetLineColor(2);
    gr_rc[i]->Draw("p");
  }

  TGraphErrors *gr_corr[NV2];
  for(int i=0;i<NV2;i++) {
    gr_corr[i] = new TGraphErrors(NPt, pT1, r_corr[i], 0, re_corr[i]);
    gr_corr[i]->SetName(Form("drho00Corr_Cos2PhiRP_Rc_%d",i));
    gr_corr[i]->SetMarkerStyle(markerStyle[i]);
    gr_corr[i]->SetMarkerSize(2.0);
    gr_corr[i]->SetMarkerColor(4);
    gr_corr[i]->SetLineWidth(2);
    gr_corr[i]->SetLineColor(4);
    gr_corr[i]->Draw("p");
  }


  TLegend *leg = new TLegend(0.5, 0.64, 0.65, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.65, 0.64, 0.8, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_rc[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.8, 0.64, 0.95, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_corr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  drawText(2.3, 0.035, "MC");
  drawText(3.2, 0.035, "RC", 42, 0.05, 0, 2);
  drawText(4.0, 0.035, "RC corr.", 42, 0.05, 0, 4);
  drawHistBox(0., 5.0, -0.01, 0.04);

  c2->Update();

  c2->SaveAs(Form("fig/Cos2PhiRP_pT_v2pt_%s.pdf",RcName));
  c2->SaveAs(Form("fig/Cos2PhiRP_pT_v2pt_%s.png",RcName));

  TFile *fout = new TFile(Form("root/drho00_Cos2PhiRP_pT_v2pt_%s.root",RcName),"recreate");
  for(int i=0;i<NV2;i++) {
    gr[i]->Write();
    gr_rc[i]->Write();
    gr_corr[i]->Write();
  }
  fout->Close();
}
