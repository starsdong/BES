#include "TRandom3.h"

void EPRes()
{
  const Int_t NRes = 10;
  const Int_t Nevt = 1e6;

  TRandom3 *gRandom = new TRandom3;

  TH2D *hEPRes2D[NRes];
  for(int i=0;i<NRes;i++) {
    hEPRes2D[i] = new TH2D(Form("EPRes2D_%d",i),"",180,-TMath::Pi(),TMath::Pi(),200,-1,1);
  }

  double PsiRes[NRes];
  double RPRes[NRes];
  for(int i=0;i<NRes;i++) {
    PsiRes[i] = 1.0/NRes * (i+0.5);

    for(int j=0;j<Nevt;j++) {
      double Psi_RP_0 = gRandom->Rndm()*TMath::Pi()*2.0;
      double Psi_RP = gRandom->Gaus(Psi_RP_0, PsiRes[i]);

      hEPRes2D[i]->Fill(Psi_RP-Psi_RP_0,TMath::Cos(2*(Psi_RP-Psi_RP_0)));
    }

    RPRes[i] = hEPRes2D[i]->GetMean(2);
  }

  TGraph *gr_res = new TGraph(NRes,PsiRes,RPRes);
  gr_res->SetName("EPRes");
  gr_res->SetMarkerSize(2.5);
  
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Draw();

  gr_res->Draw("APC");
  c1->Update();
  

  
  TFile *fout = new TFile("EPRes.root","recreate");
  for(int i=0;i<NRes;i++) {
    hEPRes2D[i]->Write();
  }
  gr_res->Write();
  fout->Close();
}
