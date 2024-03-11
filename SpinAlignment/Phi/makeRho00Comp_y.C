////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compare Rho00 extracted from various methods:
//      1:  fit to costheta* distribution using (1 - rho_00) + (3*rho_00 - 1) * costheta*^2
//      2:  delta_rho00 =  -4/3 < cos(2*(phi* - Psi_RP)) >
//      3:  delta_rho00 = 5/2 (<cos^2 theta*> - 1/3)
///////////////////////////////////////////////////////////////////////////////////////////////////////
#include "style.C+"
#include "draw.C+"

void makeRho00Comp_y(const Int_t Iv2 = 2)
{
  style();
  
  const Int_t NF = 3; // configurations
  const Char_t *FileName[NF] = {"drho00_Cos2PhiRP", "drho00_CosTheta2", "rho00"};
  const Int_t NV2 = 4;
  const Double_t v2Name[NV2] = {0.0, 0.1, 0.2, 0.3};
  const Char_t *HistName[NF] = {"drho00_Cos2PhiRP_Rc", "drho00_CosTheta2_Rc", "rho00_CosThetaBin"};
  const Char_t *LabelName[NF] = {"-4/3*<cos[2(#phi*-#Psi)]>", "5/2*(<cos^{2}#theta*>-1/3)", "cos#theta* fit"};
  const Int_t NY = 10;
  const Double_t rap[NY] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  
  TFile *fin[NF];
  TGraphErrors *gr[NF][NV2];
  for(int i=0;i<NF;i++) {
    fin[i] = new TFile(Form("root/%s_y_Y_9.9.root", FileName[i]));
    for(int j=0;j<NV2;j++) {
      gr[i][j] = (TGraphErrors *)fin[i]->Get(Form("%s_%d",HistName[i],j));
    }
  }

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Draw();
  TH2D *h2 = new TH2D("h2","",1,0.0,1.2,1,-0.01, 0.02);
  h2->GetXaxis()->SetTitle("#phi-meson rapidity");
  h2->GetYaxis()->SetTitle("#Delta#rho_{00} Raw");  
  h2->Draw();
  drawLine(0.0, 0.0, 1.2, 0.0, 2, 8, 1);

  const Int_t markerStyle[NV2] = {24, 21, 20, 22};
  const Int_t markerColor[NF] = {1, 2, 4};
  
  for(int i=0;i<NF;i++) {
    for(int j=0;j<NV2;j++) {
      if(j!=0 && j!=Iv2) continue;
      gr[i][j]->SetMarkerStyle(markerStyle[j]+i);
      gr[i][j]->SetMarkerColor(markerColor[i]);
      gr[i][j]->SetLineColor(markerColor[i]);      
      gr[i][j]->Draw("p");
    }
  }

  TLegend *leg = new TLegend(0.25, 0.64, 0.5, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=0;i<NF;i++) {
    leg->AddEntry(gr[i][0], LabelName[i], "pl");
  }
  leg->Draw();

  leg = new TLegend(0.3, 0.64, 0.6, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=0;i<NF;i++) {
    leg->AddEntry(gr[i][2], LabelName[i], "pl");
  }
  leg->Draw();

  drawText(0.05, 0.017, "v_{2} = 0    0.2", 42, 0.045);
  drawLine(0.19, 0.007, 0.19, 0.018, 1, 2);
  drawLine(0.1, 0.0165, 0.28, 0.0165, 1, 2);
  
  drawHistBox(0.0, 1.2, -0.01, 0.02);
  c2->Update();
  c2->SaveAs(Form("fig/drho00Comp_y_%d.png", Iv2));
  c2->SaveAs(Form("fig/drho00Comp_y_%d.pdf", Iv2));


  // Acceptance corrected for v2=0 case
  double r_corr[NF][NV2][NY], re_corr[NF][NV2][NY]; // Acceptance corrected (subtrack v2=0 case)
  TGraphErrors *gr_corr[NF][NV2];
  for(int i=0;i<NF;i++) {
    for(int j=0;j<NV2;j++) {
      for(int k=0;k<gr[i][j]->GetN();k++) {
	r_corr[i][j][k] = gr[i][j]->GetY()[k] - gr[i][0]->GetY()[k];
	re_corr[i][j][k] = TMath::Sqrt(gr[i][j]->GetEY()[k]*gr[i][j]->GetEY()[k]);// + gr[i][0]->GetEY()[k]*gr[i][0]->GetEY()[k]);
      }
      gr_corr[i][j] = new TGraphErrors(NY, rap, r_corr[i][j], 0, re_corr[i][j]);
      gr_corr[i][j]->SetMarkerStyle(markerStyle[j]+4+i);
      gr_corr[i][j]->SetMarkerColor(markerColor[i]);
      gr_corr[i][j]->SetMarkerSize(2.0);
      gr_corr[i][j]->SetLineColor(markerColor[i]);      
    }
  }
  
  TCanvas *c3 = new TCanvas("c3","",800,0,800,600);
  c3->Draw();
  TH2D *h3 = new TH2D("h3","",1,0.0,1.2,1,-0.005, 0.025);
  h3->GetXaxis()->SetTitle("#phi-meson rapidity");
  h3->GetYaxis()->SetTitle("#Delta#rho_{00}");  
  h3->Draw();
  drawLine(0.0, 0.0, 1.2, 0.0, 2, 8, 1);

  for(int i=0;i<NF;i++) {
    for(int j=0;j<NV2;j++) {
      if(j!=Iv2) continue;
      gr_corr[i][j]->Draw("p");
    }
  }

  // apply -4/3*a2*v2 correction
  double r_corr0[NF][NV2][NY], re_corr0[NF][NV2][NY];
  TFile *fCorr = new TFile(Form("root/drho_a2_y_Y_9.9.root"));
  TGraphErrors *gr_drho[NV2];
  for(int i=0;i<NV2;i++) {
    gr_drho[i] = (TGraphErrors *)fCorr->Get(Form("drho_%d",i));
  }
  for(int i=0;i<NF;i++) {
    for(int j=0;j<NV2;j++) {
      for(int k=0;k<NY;k++) {
	r_corr0[i][j][k] = r_corr[i][j][k] - gr_drho[j]->GetY()[k];
	re_corr0[i][j][k] = TMath::Sqrt(re_corr[i][j][k]*re_corr[i][j][k] + gr_drho[j]->GetEY()[k]*gr_drho[j]->GetEY()[k]);
      }
    }
  }

  TGraphErrors *gr_corr0[NF][NV2];
  for(int i=0;i<NF;i++) {
    for(int j=0;j<NV2;j++) {
      gr_corr0[i][j] = new TGraphErrors(NY, rap, r_corr0[i][j], 0, re_corr0[i][j]);
      gr_corr0[i][j]->SetMarkerStyle(markerStyle[j]+i);
      gr_corr0[i][j]->SetMarkerColor(markerColor[i]);
      gr_corr0[i][j]->SetMarkerSize(2.0);
      gr_corr0[i][j]->SetLineColor(markerColor[i]);
    }
  }
  for(int i=0;i<NF;i++) {
    for(int j=0;j<NV2;j++) {  
      if(j!=Iv2) continue;
      gr_corr0[i][j]->Draw("p");
    }
  }
  
  leg = new TLegend(0.25, 0.64, 0.5, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=0;i<NF;i++) {
    leg->AddEntry(gr[i][0], "         ", "pl");
  }
  leg->Draw();

  leg = new TLegend(0.3, 0.64, 0.6, 0.88);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=0;i<NF;i++) {
    leg->AddEntry(gr[i][2], LabelName[i], "pl");
  }
  leg->Draw();

  drawText(0.08, 0.022, "Obs.  Corr.", 42, 0.045);
  drawLine(0.19, 0.012, 0.19, 0.023, 1, 2);
  drawLine(0.1, 0.0215, 0.28, 0.0215, 1, 2);

  
  drawText(0.75, 0.022, "Obs. = Raw^{v_{2}=0.2} - Raw^{v_{2}=0}", 42, 0.035);
  drawText(0.75, 0.020, "Corr. = Obs. + 4/3a_{2}v_{2}", 42, 0.035);
  
  drawHistBox(0.0, 1.2, -0.005, 0.025);
  c3->Update();
  c3->SaveAs(Form("fig/drho00CorrComp_y_%d.png", Iv2));
  c3->SaveAs(Form("fig/drho00CorrComp_t_%d.pdf", Iv2));


}
