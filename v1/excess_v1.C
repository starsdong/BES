#include "draw.C+"
#include "style.C+"

void excess_v1(const Int_t conf = 1)
// conf:  1 - w/ models, 0 - w/o models
{
  style();
  gROOT->LoadMacro("~/work/work/C/kinematic.C");

   //////////////////////////////////////////////////////////////////////////////////////////////////////BES1
   //https://drupal.star.bnl.gov/STAR/files/starpublications/271/data.html
  const Int_t NP1 = 8;
  const Double_t mu = 0.9315;
  const Double_t bes1_energy[NP1] = { 7.7        , 11.5        , 14.5       , 19.6       , 27         , 39         , 62.4       , 200         };
  Double_t yb1[NP1];  // beam rapidity
  const Double_t bes1val_pv1[NP1] = { 0.022534 , -0.000486069, -0.00684062	, -0.00690107, -0.00561899, -0.00415636, -0.00127671, -0.000601299};
  const Double_t bes1staerr_pv1[NP1] = { 0.000863638, 0.000618966 ,0.000843574, 0.000853584, 0.000836661, 0.000565664, 0.00168219 , 0.000738392 };
  const Double_t bes1syserr_pv1[NP1] = { 0.00447883 , 0.00225644  , 	0.0021    , 0.000933001, 0.0008     , 0.0001     , 0.0031039  , 0.00011     };

  const Double_t bes1val_pbarv1[NP1] = {-0.0453219, -0.0365653, -0.0380232	, -0.0311992, -0.026994, -0.0206087, -0.00598958, -0.00552446};
  const Double_t bes1staerr_pbarv1[NP1] = {0.011382, 0.00371859, 0.00473727, 0.00264533, 0.0019531, 0.00107113, 0.00281943 , 0.00040525 };
  const Double_t bes1syserr_pbarv1[NP1] = {0.00523765, 0.00325287,0.00425287, 0.00015040, 9.125e-05, 5e-05     , 0.004238   , 3.1e-05    };
   
  const Double_t bes1val_netpv1[NP1] = {0.0229, 0.0006, -0.0050, -0.0041, -0.0008, 0.0024, 0.005, 0.0113};

  Double_t bes1val_excess_p[NP1];
  Double_t bes1staerr_excess_p[NP1];
  Double_t bes1syserr_excess_p[NP1];
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Double_t r1[NP1]; // BES-1 ratio
  for(int i=0;i<NP1;i++) {
    r1[i] = (bes1val_pv1[i] - bes1val_netpv1[i])/(bes1val_pbarv1[i] - bes1val_netpv1[i]);
    yb1[i] = rapidity(bes1_energy[i]);
    bes1val_excess_p[i] = (bes1val_pv1[i] - bes1val_pbarv1[i]) / (1. - r1[i]) * yb1[i];
    bes1staerr_excess_p[i] = sqrt(bes1staerr_pv1[i]*bes1staerr_pv1[i] + bes1staerr_pbarv1[i]*bes1staerr_pbarv1[i]) / (1. - r1[i]) * yb1[i];
    bes1syserr_excess_p[i] = sqrt(bes1syserr_pv1[i]*bes1syserr_pv1[i] + bes1syserr_pbarv1[i]*bes1syserr_pbarv1[i]) / (1. - r1[i]) * yb1[i];    
    cout << " Energy: " << bes1_energy[i] << "\t pbar/p = " << r1[i] << "\t yb = " << yb1[i] << "\t Excess V1 = " << bes1val_excess_p[i] << " +/- " << bes1staerr_excess_p[i] << " +/- " << bes1syserr_excess_p[i] << endl;
  }
  TGraph *gr_pbarp = new TGraph(NP1, bes1_energy, r1);
  TGraphErrors *gr_bes1_stat = new TGraphErrors(NP1, bes1_energy, bes1val_excess_p, 0, bes1staerr_excess_p);
  gr_bes1_stat->RemovePoint(NP1-2);
  TGraphErrors *gr_bes1_sys = new TGraphErrors(NP1, bes1_energy, bes1val_excess_p, 0, bes1syserr_excess_p);
  gr_bes1_sys->RemovePoint(NP1-2);
  

   //////////////////////////////////////////////////////////////////////////////////////////////////////BES2
  const Int_t NP2 = 5;
  const Double_t bes2_energy[NP2] = {7.7, 14.6, 19.6, 27., 54.4};
  Double_t yb2[NP2];  // beam rapidity
  const Double_t bes2val_pv1[NP2] = {0.0164966  ,  -0.00578296, -0.00823163, 0., 0.};
  const Double_t bes2staerr_pv1[NP2] = {0.000201101 , 0.000161092, 0.000165627, 0., 0.};
  const Double_t bes2syserr_pv1[NP2] = {0.00000000 , 0.0000     , 0.000184515, 0., 0.};

  const Double_t bes2val_pbarv1[NP2] = {-0.0521239, -0.0331744, -0.0311547, 0., 0.};
  const Double_t bes2staerr_pbarv1[NP2] = {0.00253558, 0.000668851, 0.000499875, 0., 0.};
  const Double_t bes2syserr_pbarv1[NP2] = {0.00000000, 0.00000000, 0.000479769, 0., 0.};
      
  // Double_t bes2val_excess_p_2[NP2] =   {0.0688809   ,0.0290301 ,0.0257154    };
  // Double_t bes2staerr_excess_p_2[NP2] = {0.00508717  ,0.0013760,0.00105332};
  // Double_t bes2syserr_excess_p_2[NP2] = {0.0        ,0.0      ,0.000781936};

  Double_t bes2val_excess_p[NP2];
  Double_t bes2staerr_excess_p[NP2];
  Double_t bes2syserr_excess_p[NP2];

  TFile *fin = new TFile("DataProtonV1TwoComp_27nG4p4GeV.root");  // 27 and 54.4 GeV data from Sooraj
  TGraphErrors *gr_SR = (TGraphErrors *)fin->Get("gSlope_excessv1_c1");
  
  Double_t r2[NP2]; // BES-2 ratio
  Double_t bes2_energy_d[NP2];
  for(int i=0;i<NP2;i++) {
    bes2_energy_d[i] = bes2_energy[i]*1.05;
    r2[i] = gr_pbarp->Eval(bes2_energy[i]);
    yb2[i] = rapidity(bes2_energy[i]);

    if(i<=2) { // Data from Emmy
      bes2val_excess_p[i] = (bes2val_pv1[i] - bes2val_pbarv1[i]) / (1. - r2[i]) * yb2[i];
      bes2staerr_excess_p[i] = sqrt(bes2staerr_pv1[i]*bes2staerr_pv1[i] + bes2staerr_pbarv1[i]*bes2staerr_pbarv1[i]) / (1. - r2[i]) * yb2[i];
      bes2syserr_excess_p[i] = sqrt(bes2syserr_pv1[i]*bes2syserr_pv1[i] + bes2syserr_pbarv1[i]*bes2syserr_pbarv1[i]) / (1. - r2[i]) * yb2[i];
    } else { // Data from Sooraj
      bes2val_excess_p[i] = gr_SR->GetY()[i-3] * yb2[i];
      bes2staerr_excess_p[i] = gr_SR->GetEY()[i-3] * yb2[i];
      bes2syserr_excess_p[i] = 0;      
    }
    cout << " Energy: " << bes2_energy[i] << "\t pbar/p = " << r2[i] << "\t yb = " << yb2[i] << "\t Excess V1 = " << bes2val_excess_p[i] << " +/- " << bes2staerr_excess_p[i] << " +/- " << bes2syserr_excess_p[i] << endl;
      //    cout << " Energy: " << bes2_energy[i] << "\t pbar/p = " << r2[i] << "\t yb = " << yb2[i] << "\t Excess V1 = " << bes2val_excess_p_2[i] << " +/- " << bes2staerr_excess_p_2[i] << " +/- " << bes2syserr_excess_p_2[i] << endl;
      //    cout << "---" << endl;
  }
  TGraphErrors *gr_bes2_stat = new TGraphErrors(NP2, bes2_energy_d, bes2val_excess_p, 0, bes2staerr_excess_p);
  TGraphErrors *gr_bes2_sys = new TGraphErrors(NP2, bes2_energy_d, bes2val_excess_p, 0, bes2syserr_excess_p);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Models
  const Int_t NM = 2;
  const Int_t NEMax = 10;
  TFile *fin = new TFile("JAMProtonV1TwoComp.root");
  const Char_t *GraphName[NM] = {"gExcessSlope_casc_c1","gExcessSlope_mfrmf_c1"};
  const Char_t *LabelName[NM] = {"JAM Cascade","JAM Mean Field"};
  TGraph *gr_m[NM], *gr_m_tmp[NM];
  Double_t model_energy[NM][NEMax], model_v1[NM][NEMax];
  for(int i=0;i<NM;i++) {
    gr_m_tmp[i] = (TGraph *)fin->Get(GraphName[i]);
    for(int j=0;j<gr_m_tmp[i]->GetN();j++) {
      model_energy[i][j] = gr_m_tmp[i]->GetX()[j];
      double y_m = rapidity(model_energy[i][j]);
      model_v1[i][j] = gr_m_tmp[i]->GetY()[j] * y_m;
    }
    gr_m[i] = new TGraph(gr_m_tmp[i]->GetN(), model_energy[i], model_v1[i]);
    gr_m[i]->Print();
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd()->SetLogx();

  double x1 = 3.5;
  double x2 = 4e2;
  double y1 = 0.05;
  double y2 = 0.2;
  if(conf) y2 = 0.6;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->SetYTitle("dv_{1}/d(y/y_{beam})");
  h0->GetYaxis()->SetTitleOffset(1.1);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw("c");

  drawLine(x1*1.1, 0.078, x2*0.9, 0.078, 2, 8, 1);

  const Int_t m_lineColor[NM] = {2, 4};
  const Int_t m_lineStyle[NM] = {2, 3};
  for(int i=0;i<NM;i++) {
    setGraphLine(gr_m[i], m_lineStyle[i], m_lineColor[i]);
    if(conf) gr_m[i]->Draw("l");
  }

  setGraphMarker(gr_bes1_sys, 24, 1, 1.5);
  setGraphLine(gr_bes1_sys, 1, 1, 2);
  drawSysBox(gr_bes1_sys, 0.08, 18, 1);

  setGraphMarker(gr_bes1_stat, 24, 1, 1.5);
  setGraphLine(gr_bes1_stat, 1, 1, 2);
  gr_bes1_stat->Draw("p");
  
  setGraphMarker(gr_bes2_sys, 20, 2, 1.5);
  setGraphLine(gr_bes2_sys, 1, 2, 2);
  drawSysBox(gr_bes2_sys, 0.08, 16, 1);

  setGraphMarker(gr_bes2_stat, 20, 2, 1.5);
  setGraphLine(gr_bes2_stat, 1, 2, 2);
  gr_bes2_stat->Draw("p");

  //  gr_bes2_stat->Fit("pol0","R","",10.,100.);
  
  drawHistBox(x1,x2,y1,y2);
  drawText(9.0, -0.06*(y2-y1)+y1, "10");
  drawText(86, -0.06*(y2-y1)+y1, "100");

  drawText(5.0, y2*0.9, "Excess Protons", 32, 0.06);
  drawText(70, y2*0.9, "Au+Au 10-40%", 42, 0.05);
  
  leg = new TLegend(0.7, 0.68, 0.92, 0.84);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_bes2_stat, "  BES-II", "pl");
  leg->AddEntry(gr_bes1_stat, "  BES-I", "pl");
  leg->Draw();
  
  leg = new TLegend(0.7, 0.68-0.06*NM, 0.92, 0.68);
  leg->SetLineColor(10);
  leg->SetTextSize(0.035);
  for(int i=0;i<NM;i++) {
    leg->AddEntry(gr_m[i], LabelName[i], "l");
  }
  if(conf) leg->Draw();

  
  c1->Update();
  c1->SaveAs(Form("excess_v1_%d.pdf",conf));  
  c1->SaveAs(Form("excess_v1_%d.png",conf));
  
}
