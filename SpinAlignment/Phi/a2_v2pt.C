#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

//void rho00(const Double_t v2 = 0.1, const Double_t sigY = 5.0)
void a2_v2pt(const Char_t *RcName = "Rc4")
{
  style();
  const Double_t MassPhi = 1.01946;
  
  TF1 *funA2 = new TF1("funA2","[1]*(1.+2.*[0]*x)",-1,1);
  funA2->SetParameter(0, 0.0);
  const Int_t NMAX = 200; // number of cos2phi bins

  // const Float_t i1 = 25;
  // const Float_t i2 = 48;  // 1.2 - 2.4 GeV/c
  const Int_t NPt = 12;
  const Int_t i_edge[NPt+1] = {12, 18, 24, 30, 36, 42, 48, 54, 60, 72, 84, 96, 108};
  const Double_t pT[NPt] = {0.75, 1.05, 1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.3, 3.9, 4.5, 5.1};
  // const Int_t NPt = 6;
  // const Int_t i_edge[NPt+1] = {12, 24, 36, 48, 60, 84, 108};
  // const Double_t pT[NPt] = {0.9, 1.5, 2.1, 2.7, 3.6, 4.8};
  const Int_t NV2 = 2;
  const Char_t *v2Name[NV2] = {"v2_0.0_Y_9.9", "v2pt_Minv"};
  const Double_t sc = -4./3.;

  Double_t a2[NV2][NPt], a2e[NV2][NPt];
  Double_t drho[NV2][NPt], drhoe[NV2][NPt];
  Double_t v2[NV2][NPt];
  for(int i=0;i<NV2;i++) {
    for(int j=0;j<NPt;j++) {
      if(i==0) { v2[i][j] = 0.;}
      else {
	double mT = TMath::Sqrt(pT[j]*pT[j] + MassPhi*MassPhi);
	v2[i][j] = (mT - MassPhi) > 1.2 ? 0.2 :  0.2/1.2 * (mT - MassPhi);
      }
    }
  }
  
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
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hPtCos2Phi"))->ProjectionY(Form("Mc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      fMc[iv2][i]->Rebin(20);
      if(iv2==0) {
	fRc[iv2][i] = ((TH2D *)fin[iv2]->Get(Form("hPtCos2PhiRc3")))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      } else {
	fRc[iv2][i] = ((TH2D *)fin[iv2]->Get(Form("hPtCos2Phi%s",RcName)))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      }
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

      TH1D *h0 = new TH1D("h0","",1,-1,1);
      h0->SetMaximum(eff[0]/0.7);
      h0->Draw();

      gr_e[iv2][i]->Draw("p");
      gr_e[iv2][i]->Fit("funA2","R");
      a2[iv2][i] = funA2->GetParameter(0);
      a2e[iv2][i] = funA2->GetParError(0);

      drho[iv2][i] = sc * a2[iv2][i] * v2[iv2][i];
      drhoe[iv2][i] = sc * a2e[iv2][i] * v2[iv2][i];
      c1->Update();
    }
    fin[iv2]->Close();
  }

  TCanvas *c2 = new TCanvas("c2","");
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,5.0,1,-0.1, 0.01);
  h2->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  const Int_t markerStyle[NV2] = {24, 20};
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NPt, pT, a2[i], 0, a2e[i]);
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }
  c2->Update();

  c2->SaveAs(Form("fig/a2_pT_v2pt_%s.pdf",RcName));
  c2->SaveAs(Form("fig/a2_pT_v2pt_%s.png",RcName));
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->Draw();
  TH2D *h3 = new TH2D("h3","",1,0.0,5.0,1,-0.01, 0.04);
  h3->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  TGraphErrors *gr_drho[NV2];
  for(int i=0;i<NV2;i++) {
    gr_drho[i] = new TGraphErrors(NPt, pT, drho[i], 0, drhoe[i]);
    gr_drho[i]->SetName(Form("drho_%d",i));
    gr_drho[i]->SetMarkerStyle(markerStyle[i]);
    gr_drho[i]->SetMarkerSize(2.0);
    gr_drho[i]->SetLineWidth(2);
    gr_drho[i]->Draw("p");
  }
  c3->Update();

  TFile *fout = new TFile(Form("root/drho_a2_pT_v2pt_%s.root", RcName),"recreate");
  for(int i=0;i<NV2;i++) {
    gr_drho[i]->Write();
  }
  fout->Close();

  c3->SaveAs(Form("fig/drho_a2_pT_v2pt_%s.pdf", RcName));
  c3->SaveAs(Form("fig/drho_a2_pT_v2pt_%s.png", RcName));
}
