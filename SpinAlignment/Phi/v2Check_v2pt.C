#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

//void v2Check(const Double_t sigY = 9.9)
void v2Check_v2pt()
{
  style();
    
  const Int_t NMAX = 200; // number of cos2phi bins

  // const Float_t i1 = 25;
  // const Float_t i2 = 48;  // 1.2 - 2.4 GeV/c
  const Int_t NPt = 6;
  const Int_t i_edge[NPt+1] = {12, 24, 36, 48, 60, 84, 108};
  const Double_t pT[NPt] = {0.9, 1.5, 2.1, 2.7, 3.6, 4.8};
  const Double_t pT1[NPt] = {0.95, 1.55, 2.15, 2.75, 3.65, 4.85};
  const Int_t NV2 = 2;
  const Char_t *v2Name[NV2] = {"v2_0.0_Y_9.9", "v2pt"};
  const Char_t *v2Label[NV2] = {"v200", "v2pt"};

  double_t v2[NV2][NPt], v2e[NV2][NPt];
  double_t v2Rc[NV2][NPt], v2eRc[NV2][NPt];
  TCanvas *c1 = new TCanvas("c1","",1600,900);
  c1->Draw();
  c1->Divide(6,3);
  //  TFile *fin = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root",v2, sigY));
  TFile *fin[NV2];
  for(int iv2 = 0;iv2<NV2;iv2++) {
    fin[iv2] = new TFile(Form("accept_Sergei_%s.root", v2Name[iv2]));
    //    fin[iv2] = new TFile(Form("accept_Sergei_v2_%3.1f_Y_%3.1f.root", v2Name[iv2], sigY));
    TH1D *fMc[NV2][NPt], *fRc[NV2][NPt];
    TGraphErrors *gr_e[NV2][NPt];
    for(int i=0;i<NPt;i++) {
      c1->cd(i+1+iv2*NPt);
      fMc[iv2][i] = ((TH2D *)fin[iv2]->Get("hPtV2"))->ProjectionY(Form("Mc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);
      fRc[iv2][i] = ((TH2D *)fin[iv2]->Get("hPtV2Rc3"))->ProjectionY(Form("Rc_%d_%d",iv2, i),i_edge[i]+1,i_edge[i+1]);


      TH1D *h0 = new TH1D("h0","",1,-1,1);
      h0->SetMaximum(fRc[iv2][i]->GetMaximum()/0.7);
      h0->SetMinimum(0);
      h0->Draw();

      fRc[iv2][i]->Draw("same");

      v2[iv2][i] = fMc[iv2][i]->GetMean();
      v2e[iv2][i] = fMc[iv2][i]->GetMeanError();
      v2Rc[iv2][i] = fRc[iv2][i]->GetMean();
      v2eRc[iv2][i] = fRc[iv2][i]->GetMeanError();
      c1->Update();
    }
    fin[iv2]->Close();
  }
  c1->Update();
  

  TCanvas *c2 = new TCanvas("c2","");
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,5.0,1,-0.05, 0.3);
  h2->GetXaxis()->SetTitle("#phi-meson p_{T} (GeV/c)");
  h2->GetYaxis()->SetTitle("v_2");  
  h2->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  TGraphErrors *gr[NV2];
  const Int_t markerStyle[NV2] = {24, 20};
  for(int i=0;i<NV2;i++) {
    gr[i] = new TGraphErrors(NPt, pT, v2[i], 0, v2e[i]);
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetLineWidth(2);
    gr[i]->Draw("p");
  }

  TGraphErrors *gr_rc[NV2];
  for(int i=0;i<NV2;i++) {
    gr_rc[i] = new TGraphErrors(NPt, pT1, v2Rc[i], 0, v2eRc[i]);
    gr_rc[i]->SetMarkerStyle(markerStyle[i]);
    gr_rc[i]->SetMarkerSize(2.0);
    gr_rc[i]->SetMarkerColor(2);
    gr_rc[i]->SetLineWidth(2);
    gr_rc[i]->SetLineColor(2);
    gr_rc[i]->Draw("p");
  }


  TLegend *leg = new TLegend(0.55, 0.36, 0.75, 0.54);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  leg = new TLegend(0.75, 0.36, 0.95, 0.54);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=NV2-1;i>=0;i--) {
    leg->AddEntry(gr_rc[i], Form("v_{2} = %s", v2Label[i]), "pl");
  }
  leg->Draw();
  drawText(2.6, 0.12, "MC");
  drawText(4.0, 0.12, "RC", 42, 0.05, 0, 2);
  drawHistBox(0., 5.0, -0.05, 0.3);
  
  c2->Update();

  c2->SaveAs(Form("fig/v2_pT_v2pt.pdf"));
  c2->SaveAs(Form("fig/v2_pT_v2pt.png"));
  
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
