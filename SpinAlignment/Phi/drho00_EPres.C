#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

/////////////////////////////////////////////////////
// delta_rho00 = -4/3 < cos(2*(phi* - Psi_RP)) >
/////////////////////////////////////////////////////
void drho00_EPres(const Double_t sigY = 5.0)
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
  // const Int_t NV2 = 4;
  // const Double_t v2Name[NV2] = {0.0, 0.1, 0.2, 0.3};
  const Int_t NV2 = 1;
  const Double_t v2Name[NV2] = {0.1};
  const Int_t NRes = 2;
  const Double_t resName[NRes] = {0.0, 0.5};
  const Double_t sc = -4./3.;

  TFile *fres = new TFile("EPRes.root");
  TGraph *gr_res = (TGraph *)fres->Get("EPRes");
  

  double_t v2[NV2][NRes][NPt], v2e[NV2][NRes][NPt];
  double_t v2Rc[NV2][NRes][NPt], v2eRc[NV2][NRes][NPt];
  double_t v2RcT[NV2][NRes][NPt], v2eRcT[NV2][NRes][NPt];  // EP res correction
  double_t v2RcCorr[NV2][NRes][NPt], v2eRcCorr[NV2][NRes][NPt];  // EP res correction
  TH1D *fMc[NV2][NRes][NPt], *fRc[NV2][NRes][NPt], *fRcT[NV2][NRes][NPt];
  TGraphErrors *gr_e[NV2][NRes][NPt];

  TCanvas *c1 = new TCanvas("c1","",2000,1000);
  c1->Draw();
  c1->Divide(NPt/2, NRes*2);
  //  TFile *fin = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root",v2, sigY));
  TFile *fin[NV2][NRes];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    for(int ires = 0; ires<NRes;ires++) {
      fin[iv2][ires] = new TFile(Form("accept_Sergei_v2_%3.1f_res_%3.1f_Y_%3.1f.root", v2Name[iv2], resName[ires], sigY));
      double EPRes = gr_res->Eval(resName[ires]);
      if(EPRes>1.0) EPRes = 1.0;
      cout << " Input Psi Res = " << resName[ires] << "  EPRes = " << EPRes << endl;
      for(int i=0;i<NPt;i++) {
	c1->cd(i+1+ires*NPt+iv2*NPt*NRes);
	fMc[iv2][ires][i] = ((TH2D *)fin[iv2][ires]->Get("hPtCos2PhiRP"))->ProjectionY(Form("Mc_%d_%d_%d",iv2, ires, i),i_edge[i]+1,i_edge[i+1]);
	fRc[iv2][ires][i] = ((TH2D *)fin[iv2][ires]->Get("hPtCos2PhiRPRc3"))->ProjectionY(Form("Rc_%d_%d_%d",iv2, ires, i),i_edge[i]+1,i_edge[i+1]);
	fRcT[iv2][ires][i] = ((TH2D *)fin[iv2][ires]->Get("hPtCos2PhiRPRc3T"))->ProjectionY(Form("RcT_%d_%d_%d",iv2, ires, i),i_edge[i]+1,i_edge[i+1]);
	
	
	TH1D *h0 = new TH1D("h0","",1,fRc[iv2][ires][i]->GetXaxis()->GetXmin(), fRc[iv2][ires][i]->GetXaxis()->GetXmax());
	h0->SetMaximum(fRc[iv2][ires][i]->GetMaximum()/0.7);
	h0->SetMinimum(0);
	h0->Draw();
	
	fMc[iv2][ires][i]->SetLineWidth(2);
	fMc[iv2][ires][i]->Draw("esame");
	
	v2[iv2][ires][i] = fMc[iv2][ires][i]->GetMean() * sc;
	v2e[iv2][ires][i] = fMc[iv2][ires][i]->GetMeanError() * sc;
	v2Rc[iv2][ires][i] = fRc[iv2][ires][i]->GetMean() * sc;
	v2eRc[iv2][ires][i] = fRc[iv2][ires][i]->GetMeanError() * sc;
	// v2RcCorr[iv2][ires][i] = v2Rc[iv2][ires][i]/EPRes;
	// v2eRcCorr[iv2][ires][i] = v2eRc[iv2][ires][i]/EPRes;
	v2RcT[iv2][ires][i] = fRcT[iv2][ires][i]->GetMean() * sc;
	v2eRcT[iv2][ires][i] = fRcT[iv2][ires][i]->GetMeanError() * sc;
	c1->Update();
      }
    }
    //    fin[iv2]->Close();
  }

  /*
  double r_corr[NV2][NPt], re_corr[NV2][NPt];
  TFile *fCorr = new TFile(Form("root/drho_a2_pT_Y_%3.1f.root", sigY));
  TGraphErrors *gr_drho[NV2];
  for(int i=0;i<NV2;i++) {
    gr_drho[i] = (TGraphErrors *)fCorr->Get(Form("drho_%d",i));
    for(int j=0;j<NPt;j++) {
      r_corr[i][j] = v2Rc[i][j] - gr_drho[i]->GetY()[j];
      re_corr[i][j] = TMath::Sqrt(v2eRc[i][j]*v2eRc[i][j] + gr_drho[i]->GetEY()[j]*gr_drho[i]->GetEY()[j]);    
    }
  }
  */

  TCanvas *c2 = new TCanvas("c2","", 800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,5.0,1,-0.01, 0.04);
  h2->GetXaxis()->SetTitle("#phi-meson p_{T} (GeV/c)");
  h2->GetYaxis()->SetTitle("#Delta#rho_{00} = -4/3 * <cos[2(#phi* - #Psi)]>");  
  h2->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2][NRes];
  const Int_t markerStyle[NRes] = {24, 20};
  for(int i=0;i<NV2;i++) {
    for(int j=0;j<NRes;j++) {
      gr[i][j] = new TGraphErrors(NPt, pT, v2[i][j], 0, v2e[i][j]);
      gr[i][j]->SetName(Form("drho00_Cos2PhiRP_Mc_%d_%d",i,j));
      gr[i][j]->SetMarkerStyle(markerStyle[j]);
      gr[i][j]->SetMarkerSize(2.0);
      gr[i][j]->SetLineWidth(2);
      gr[i][j]->Draw("p");
    }
  }

  //  return;
  TGraphErrors *gr_rc[NV2][NRes];
  for(int i=0;i<NV2;i++) {
    for(int j=0;j<NRes;j++) {
      gr_rc[i][j] = new TGraphErrors(NPt, pT, v2Rc[i][j], 0, v2eRc[i][j]);
      gr_rc[i][j]->SetName(Form("drho00_Cos2PhiRP_Rc_%d_%d",i,j));
      gr_rc[i][j]->SetMarkerStyle(markerStyle[j]);
      gr_rc[i][j]->SetMarkerSize(2.0);
      gr_rc[i][j]->SetMarkerColor(2);
      gr_rc[i][j]->SetLineWidth(2);
      gr_rc[i][j]->SetLineColor(2);
      gr_rc[i][j]->Draw("p");
    }
  }

  TGraphErrors *gr_rcT[NV2][NRes];
  for(int i=0;i<NV2;i++) {
    for(int j=0;j<NRes;j++) {
      gr_rcT[i][j] = new TGraphErrors(NPt, pT1, v2RcT[i][j], 0, v2eRcT[i][j]);
      gr_rcT[i][j]->SetName(Form("drho00_Cos2PhiRP_RcT_%d_%d",i,j));
      gr_rcT[i][j]->SetMarkerStyle(markerStyle[j]);
      gr_rcT[i][j]->SetMarkerSize(2.0);
      gr_rcT[i][j]->SetMarkerColor(4);
      gr_rcT[i][j]->SetLineWidth(2);
      gr_rcT[i][j]->SetLineColor(2);
      //      gr_rcT[i][j]->Draw("p");
    }
  }

  TGraphErrors *gr_corr[NV2][NRes];
  for(int i=0;i<NV2;i++) {
    for(int j=0;j<NRes;j++) {
      gr_corr[i][j] = new TGraphErrors(NPt, pT1, v2RcCorr[i][j], 0, v2eRcCorr[i][j]);
      gr_corr[i][j]->SetName(Form("drho00Corr_Cos2PhiRP_Rc_%d_%d",i,j));
      gr_corr[i][j]->SetMarkerStyle(markerStyle[j]);
      gr_corr[i][j]->SetMarkerSize(2.0);
      gr_corr[i][j]->SetMarkerColor(4);
      gr_corr[i][j]->SetLineWidth(2);
      gr_corr[i][j]->SetLineColor(4);
      //      gr_corr[i][j]->Draw("p");
    }
  }


  TLegend *leg = new TLegend(0.5, 0.70, 0.65, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int j=NRes-1;j>=0;j--) {
    leg->AddEntry(gr[0][j], Form("#sigma_{#Psi} = %3.1f", resName[j]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.65, 0.70, 0.8, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int j=NRes-1;j>=0;j--) {
    leg->AddEntry(gr_rc[0][j], Form("#sigma_{#Psi} = %3.1f", resName[j]), "pl");
  }
  leg->Draw();
  /*
  leg = new TLegend(0.8, 0.70, 0.95, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int j=NRes-1;j>=0;j--) {
    leg->AddEntry(gr_rcT[0][j], Form("#sigma_{#Psi} = %3.1f", resName[j]), "pl");
  }
  leg->Draw();
  */
  drawText(2.3, 0.035, "MC");
  drawText(3.2, 0.035, "RC", 42, 0.05, 0, 2);
  //  drawText(4.0, 0.035, "RCT", 42, 0.05, 0, 4);
  drawHistBox(0., 5.0, -0.01, 0.04);

  c2->Update();

  c2->SaveAs(Form("fig/Cos2PhiRP_pT_Res_Y_%3.1f.pdf", sigY));
  c2->SaveAs(Form("fig/Cos2PhiRP_pT_Res_Y_%3.1f.png", sigY));

  TFile *fout = new TFile(Form("root/drho00_Cos2PhiRP_pT_Res_Y_%3.1f.root", sigY),"recreate");
  for(int i=0;i<NV2;i++) {
    for(int j=0;j<NRes;j++) {
      gr[i][j]->Write();
      gr_rc[i][j]->Write();
      gr_rcT[i][j]->Write();
    }
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
