#include "style.C+"
#include "draw.C+"

//void rho00(const Double_t v2 = 0.1, const Double_t sigY = 5.0)
void rho00_AccCorr_v2pt(const Char_t *RcName = "Rc4")
{
  style();
  TF1 *funCosTheta = new TF1("funCosTheta","[1]*((1-[0])+(3*[0]-1)*x*x)",-1,1);
  funCosTheta->SetParameter(0, 1./3.);

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
  const Char_t *v2Label[NV2] = {"v200", "v2pt"};
  const Int_t NMAX = 200;
  
  double_t r0[NV2][NPt], r0e[NV2][NPt];
  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->Draw();
  c1->Divide(NPt, NV2);
  //  TFile *fin = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root",v2, sigY));
  TFile *fin[NV2];
  TH1D *fRc[NV2][NPt];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    fin[iv2] = new TFile(Form("accept_Sergei_%s.root", v2Name[iv2]));
    for(int i=0;i<NPt;i++) {
      c1->cd(i+1+iv2*NPt);
      if(iv2==0) {
	fRc[iv2][i] = ((TH2D *)fin[iv2]->Get("hPtCosThetaRc3"))->ProjectionY(Form("Rc_%d_%d",iv2,i),i_edge[i]+1,i_edge[i+1]);
      } else {
	fRc[iv2][i] = ((TH2D *)fin[iv2]->Get(Form("hPtCosTheta%s",RcName)))->ProjectionY(Form("Rc_%d_%d",iv2,i),i_edge[i]+1,i_edge[i+1]);
      }
      fRc[iv2][i]->Rebin(fRc[iv2][i]->GetNbinsX()/NMAX);
      fRc[iv2][i]->Draw();
      funCosTheta->SetParameter(1, fRc[iv2][i]->GetMaximum());
      fRc[iv2][i]->Fit("funCosTheta","R");
      r0[iv2][i] = funCosTheta->GetParameter(0) - 1./3.;
      r0e[iv2][i] = funCosTheta->GetParError(0);
      c1->Update();
    }
  }

  cout << "============================================================" << endl;
  cout << "============================================================" << endl;
  cout << "============================================================" << endl;
  cout << "============================================================" << endl;

  // Apply acceptance*efficiency correction v2=0 case
  double_t r[NV2][NPt], re[NV2][NPt];
  TGraphErrors *fRcCorr[NV2][NPt];
  double xx[NV2][NPt][NMAX], yCorr[NV2][NPt][NMAX], yeCorr[NV2][NPt][NMAX];
  for(int iv2 = 0; iv2<NV2; iv2++) {
    for(int i=0;i<NPt;i++) {
      funCosTheta->SetParameters(r0[0][i] + 1./3., 1.);  // v2=0 case
      for(int j=0;j<fRc[iv2][i]->GetNbinsX();j++) {
	xx[iv2][i][j] = fRc[iv2][i]->GetBinCenter(j+1);
	yCorr[iv2][i][j] = fRc[iv2][i]->GetBinContent(j+1)/funCosTheta->Eval(xx[iv2][i][j]);
	yeCorr[iv2][i][j] = TMath::Sqrt(fRc[iv2][i]->GetBinContent(j+1))/funCosTheta->Eval(xx[iv2][i][j]);
      }
      fRcCorr[iv2][i] = new TGraphErrors(NMAX, xx[iv2][i], yCorr[iv2][i], 0, yeCorr[iv2][i]);
      funCosTheta->SetParameters(1./3, yCorr[iv2][i][0]);
      fRcCorr[iv2][i]->Fit("funCosTheta","R");
      r[iv2][i] = funCosTheta->GetParameter(0) - 1./3.;
      re[iv2][i] = funCosTheta->GetParError(0);
    }
  }
  
  
  double r_corr[NV2][NPt], re_corr[NV2][NPt];
  TFile *fCorr = new TFile(Form("root/drho_a2_pT_v2pt_%s.root",RcName));
  TGraphErrors *gr_drho[NV2];
  for(int i=1;i<NV2;i++) {
    gr_drho[i] = (TGraphErrors *)fCorr->Get(Form("drho_%d",i));
    for(int j=0;j<NPt;j++) {
      r_corr[i][j] = r[i][j] - gr_drho[i]->GetY()[j];
      re_corr[i][j] = TMath::Sqrt(re[i][j]*re[i][j] + gr_drho[i]->GetEY()[j]*gr_drho[i]->GetEY()[j]);    
    }
  }
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,5.0,1,-0.02,0.04);
  h2->GetXaxis()->SetTitle("#phi-meson p_{T} (GeV/c)");
  h2->GetYaxis()->SetTitle("#Delta#rho_{00} (Accep. Corr.) from cos#theta* fit");  
  h2->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  TGraphErrors *gr_corr[NV2];
  const Int_t markerStyle[NV2] = {24, 20};
  
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NPt, pT, r[i], 0, re[i]);
    gr[i]->SetName(Form("rho00_AccCorr_CosThetaBin_%d",i));
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);

    
    gr_corr[i] = new TGraphErrors(NPt, pT, r_corr[i], 0, re_corr[i]);
    gr_corr[i]->SetName(Form("rho00Corr_AccCorr_CosThetaBin_%d",i));
    gr_corr[i]->SetMarkerStyle(markerStyle[i]);
    gr_corr[i]->SetMarkerSize(2.0);
    gr_corr[i]->SetMarkerColor(2);
    gr_corr[i]->SetLineWidth(2);
    gr_corr[i]->SetLineColor(2);
    gr_corr[i]->Draw("p");
  }

  for(int i=0;i<NV2;i++) {
    gr[i]->Draw("p");
    gr_corr[i]->Draw("p");	
  }

  TLegend *leg = new TLegend(0.65, 0.64, 0.8, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.8, 0.64, 0.95, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_corr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  drawText(3.2, 0.035, "RC", 42, 0.05, 0, 2);
  drawText(4.0, 0.035, "RC corr.", 42, 0.05, 0, 4);
  drawHistBox(0., 5.0, -0.02, 0.04);
  
  c2->Update();

  c2->SaveAs(Form("fig/rho00_AccCorr_pT_v2pt_%s.pdf",RcName));
  c2->SaveAs(Form("fig/rho00_AccCorr_pT_v2pt_%s.png",RcName));

  TFile *fout = new TFile(Form("root/rho00_AccCorr_v2pt_%s.root",RcName),"recreate");
  for(int i=0;i<NV2;i++) {
    gr[i]->Write();
    gr_corr[i]->Write();
  }
  fout->Close();
  
}
