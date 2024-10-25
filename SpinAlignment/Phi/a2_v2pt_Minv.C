#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

//void rho00(const Double_t v2 = 0.1, const Double_t sigY = 5.0)
void a2_v2pt_Minv()
{
  style();
  const Double_t MassPhi = 1.01946;
  
  TF1 *funA2 = new TF1("funA2","[1]*(1.+2.*[0]*x)",-1,1);
  funA2->SetParameter(0, 0.0);
  const Int_t NMAX = 200; // number of cos2phi bins

  const Double_t Mmin = 0.98;
  const Double_t Mmax = 1.08;
  const Int_t NMmax = 200;
  const Double_t M_binW = (Mmax-Mmin)/NMmax;
  
  const Double_t M1 = 1.01;
  const Double_t M2 = 1.03;
  const Int_t NM = (M2-M1+1e-4)/M_binW; // 40 bin from 1.01 to 1.03
  const Int_t i_start = (M1-Mmin+1e-4)/M_binW + 1;
  
  Double_t Minv[NM];
  const Int_t NV2 = 1;
  const Char_t *v2Name[NV2] = {"v2pt_Minv"};
  const Double_t sc = -4./3.;

  Double_t a2[NV2][NM], a2e[NV2][NM];

  TCanvas *c1 = new TCanvas("c1","",1600,900);
  c1->Draw();
  c1->Divide(8, 5);
  //  TFile *fin = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root",v2, sigY));
  TFile *fin[NV2];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    fin[iv2] = new TFile(Form("accept_Sergei_%s.root", v2Name[iv2]));
    TH1D *fMc[NV2][NM], *fRc[NV2][NM];
    TGraphErrors *gr_e[NV2][NM];
    for(int i=0;i<NM;i++) {
      c1->cd(i+1+iv2*NM);
      Minv[i] = M1 + (M2-M1)/NM * (i+0.5);
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hMinvCos2Phi"))->ProjectionY(Form("Mc_%d_%d",iv2, i),i+i_start,i+1+i_start);
      fMc[iv2][i]->Rebin(20);
      fRc[iv2][i] = ((TH2D *)fin[iv2]->Get("hMinvCos2PhiRc3"))->ProjectionY(Form("Rc_%d_%d",iv2, i),i+i_start, i+1+i_start);
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

      c1->Update();
    }
    fin[iv2]->Close();
  }

  TCanvas *c2 = new TCanvas("c2","");
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,M1,M2,1,-0.1, 0.01);
  h2->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
  h2->GetYaxis()->SetTitle("a_{2}");
  h2->Draw();
  drawLine(M1, 0.0, M2, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  const Int_t markerStyle[NV2] = {20};
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NM, Minv, a2[i], 0, a2e[i]);
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }
  c2->Update();

  c2->SaveAs(Form("fig/a2_pT_v2pt_Minv.pdf"));
  c2->SaveAs(Form("fig/a2_pT_v2pt_Minv.png"));
  
}
