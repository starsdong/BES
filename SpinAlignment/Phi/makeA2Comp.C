#include "style.C+"
#include "draw.C+"
#include "/Users/starsdong/work/work/C/eff_err.C"

void makeA2Comp(const Int_t iV2 = 2, const Double_t sigY = 9.9)
{
  style();
  
  const Int_t NRC = 4;
  // Rc: |eta|<0.9;  Rc1: Rc+pT>0.2;  Rc2: Rc1+tpcEff;  Rc3: Rc2+tofEff
  const Char_t *RCName[NRC] = {"Rc","Rc1","Rc2","Rc3"}; //
  const Char_t *RCLabel[NRC] = {"Rc0: |#eta|<0.9", "Rc1: Rc0+pT>0.2", "Rc2: Rc1+tpcEff", "Rc3: Rc2+tofEff"};

  const Int_t NV2 = 4;
  const Double_t v2[NV2] = {0.0, 0.1, 0.2, 0.3};

  TGraphErrors *gr[NRC];
  for(int i=0;i<NRC;i++) {
    TFile *fin = new TFile(Form("root/drho_a2_Accept_%d_pT_Y_%3.1f.root", i, sigY));
    gr[i] = (TGraphErrors *)fin->Get(Form("a2_%d", iV2));
  }

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,5.0,1,-0.1, 0.01);
  h2->GetXaxis()->SetTitle("#phi-meson p_{T} (GeV/c)");
  h2->GetYaxis()->SetTitle("Eff. Slope a_{2}");
  h2->Draw();
  drawLine(0.0, 0.0, 5.0, 0.0, 2, 8, 1);

  const Int_t markerStyle[NRC] = {24, 20, 21, 22};
  for(int i=0;i<NRC;i++) {
    gr[i]->SetMarkerStyle(markerStyle[i]);
    gr[i]->SetMarkerSize(2.0);
    gr[i]->SetMarkerColor(i%2+1);
    gr[i]->SetLineWidth(2);
    gr[i]->SetLineColor(i%2+1);
    gr[i]->Draw("p");
  }

  
  TLegend *leg = new TLegend(0.66, 0.24, 0.9, 0.54);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  for(int i=0;i<NRC;i++) {
    leg->AddEntry(gr[i], RCLabel[i], "pl");
  }
  leg->Draw();
  drawHistBox(0., 5.0, -0.1, 0.01);
  c2->Update();

  c2->SaveAs(Form("fig/a2_Comp_pT_Y_%3.1f.pdf", sigY));
  c2->SaveAs(Form("fig/a2_Comp_pT_Y_%3.1f.png", sigY));
 
}
