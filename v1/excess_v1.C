#include "draw.C+"
#include "style.C+"
void excess_v1(const Int_t conf = 1)
// conf:  1 - w/ models, 0 - w/o models
{
  double protonMass = 0.93827208816;
  style();
  //gROOT->LoadMacro("~/work/work/C/kinematic.C");

   //////////////////////////////////////////////////////////////////////////////////////////////////////BES1 2018 paper
   //https://drupal.star.bnl.gov/STAR/files/starpublications/271/data.html
  const Int_t NP1 = 8;
  Double_t bes1_energy[NP1] = { 7.7        , 11.5        , 14.5       , 19.6       , 27         , 39         , 62.4       , 200         };
  Double_t yb1[NP1];  // beam rapidity
  const Double_t bes1val_pv1[NP1] = { 0.022534 , -0.000486069, -0.00684062	, -0.00690107, -0.00561899, -0.00415636, -0.00127671, -0.000601299};
  const Double_t bes1staerr_pv1[NP1] = { 0.000863638, 0.000618966 ,0.000843574, 0.000853584, 0.000836661, 0.000565664, 0.00168219 , 0.000738392 };
  const Double_t bes1syserr_pv1[NP1] = { 0.00447883 , 0.00225644  , 	0.0021    , 0.000933001, 0.0008     , 0.0001     , 0.0031039  , 0.00011     };

  const Double_t bes1val_pbarv1[NP1] = {-0.0453219, -0.0365653, -0.0380232	, -0.0311992, -0.026994, -0.0206087, -0.00598958, -0.00552446};
  const Double_t bes1staerr_pbarv1[NP1] = {0.011382, 0.00371859, 0.00473727, 0.00264533, 0.0019531, 0.00107113, 0.00281943 , 0.00040525 };
  const Double_t bes1syserr_pbarv1[NP1] = {0.00523765, 0.00325287,0.00425287, 0.00015040, 9.125e-05, 5e-05     , 0.004238   , 3.1e-05    };
   
  const Double_t bes1val_netpv1[NP1] = {0.0229, 0.0006, -0.0050, -0.0041, -0.0008, 0.0024, 0.005, 0.0113};


   ////////////////////////////BES1 2014 paper
  // const Int_t NP1 = 7;
  // const Double_t bes1_energy[NP1] =       { 7.7        , 11.5        , 19.6       , 27         , 39         , 62.4       , 200        };
  // Double_t yb1[NP1];  // beam rapidity
  // const Double_t bes1val_pv1[NP1] =       { 0.0125633  ,-0.00544871  ,-0.00594486 ,-0.00497679 ,-0.00342149 ,-0.00267291 ,-0.000821396};
  // const Double_t bes1staerr_pv1[NP1] =    { 0.00110027 ,0.00111279   ,0.000777508 ,0.000761942 ,0.000927623 ,0.00195178  ,0.000211236 };
  // const Double_t bes1syserr_pv1[NP1] =    { 0.0029     ,0.0021       ,0.0014      ,0.001       ,0.0012      ,0.00317921  ,0.00052713  };
  
  
  // const Double_t bes1val_pbarv1[NP1] =    {-0.0358057  , -0.0351478  ,-0.0293238  ,-0.023209   ,-0.0186672  ,-0.0165848   ,-0.0048911 };
  // const Double_t bes1staerr_pbarv1[NP1] = { 0.02164369 ,0.00847682   ,0.00473378  ,0.00356216  ,0.00354199  ,0.00297745   ,0.00025203 };
  // const Double_t bes1syserr_pbarv1[NP1] = {0.019       ,0.0071       ,0.004       ,0.003       ,0.003       ,0.006104     ,0.000286368};
   
  // const Double_t bes1val_netpv1[NP1] =    {0.0131475   , -0.003985273,-0.003697327, -0.000788  , 0.00162967 ,0.00762832   , 0.00852520};


  Double_t bes1val_excess_p[NP1];
  Double_t bes1staerr_excess_p[NP1];
  Double_t bes1syserr_excess_p[NP1];
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Double_t r1[NP1]; // BES-1 ratio
  for(int i=0;i<NP1;i++) {
    if(bes1_energy[i]<30.) bes1_energy[i] *= 0.95;
    r1[i] = (bes1val_pv1[i] - bes1val_netpv1[i])/(bes1val_pbarv1[i] - bes1val_netpv1[i]);
	
    //yb1[i] = rapidity(bes1_energy[i]);
	yb1[i] = TMath::ACosH(bes1_energy[i]/protonMass);
    bes1val_excess_p[i] = (bes1val_pv1[i] - bes1val_pbarv1[i]) / (1. - r1[i]) * yb1[i];
    bes1staerr_excess_p[i] = sqrt(bes1staerr_pv1[i]*bes1staerr_pv1[i] + bes1staerr_pbarv1[i]*bes1staerr_pbarv1[i]) / (1. - r1[i]) * yb1[i];
    bes1syserr_excess_p[i] = sqrt(bes1syserr_pv1[i]*bes1syserr_pv1[i] + bes1syserr_pbarv1[i]*bes1syserr_pbarv1[i]) / (1. - r1[i]) * yb1[i];    
    cout << " Energy: " << bes1_energy[i] << "\t pbar/p = " << r1[i] << "\t yb = " << yb1[i] << "\t Excess V1 = " << bes1val_excess_p[i]/yb1[i] << " +/- " << bes1staerr_excess_p[i] << " +/- " << bes1syserr_excess_p[i] << endl;
  }
  TGraph *gr_pbarp = new TGraph(NP1, bes1_energy, r1);
  TGraphErrors *gr_bes1_stat = new TGraphErrors(NP1, bes1_energy, bes1val_excess_p, 0, bes1staerr_excess_p);
  gr_bes1_stat->RemovePoint(NP1-2);
  TGraphErrors *gr_bes1_sys = new TGraphErrors(NP1, bes1_energy, bes1val_excess_p, 0, bes1syserr_excess_p);
  gr_bes1_sys->RemovePoint(NP1-2);
  

   //////////////////////////////////////////////////////////////////////////////////////////////////////BES2
  // const Int_t NP2 = 8;
  // const Int_t NP2_emmy = 6;
  // const Double_t bes2_energy[NP2] = {7.7,9.2,11.5, 14.6,17.3, 19.6, 27., 54.4};
  // Double_t yb2[NP2];  // beam rapidity
  // const Double_t bes2val_pv1[NP2] = {0.0164966  ,6.10325e-003 ,-2.12082e-003   ,  -0.00578296,-7.57647e-003, -0.00823163, 0., 0.};
  // const Double_t bes2staerr_pv1[NP2] = {0.000201101 ,1.32043e-004,1.84170e-004 , 0.000161092, 1.83521e-004, 0.000165627, 0., 0.};
  // const Double_t bes2syserr_pv1[NP2] = {0.00000000 ,0.0 ,0.0 , 0.0000     ,0.0, 0.000184515, 0., 0.};

  // const Double_t bes2val_pbarv1[NP2] = {-0.0521239, -4.25345e-002, -4.06065e-002 , -0.0331744,-3.16780e-002, -0.0311547, 0., 0.};
  // const Double_t bes2staerr_pbarv1[NP2] = {0.00253558, 1.14476e-003, 5.39328e-004   , 0.000668851,6.30026e-004, 0.000499875, 0., 0.};
  // const Double_t bes2syserr_pbarv1[NP2] = {0.00000000,0.0,0.0, 0.00000000,0.0, 0.000479769, 0., 0.};
      
  // Double_t bes2val_excess_p[NP2]={0.0688809   ,4.69272e-002, 3.97100e-002   ,2.72578e-002 ,2.62833e-002, 0.0257154  ,0.0,0.0  };
  // Double_t bes2staerr_excess_p[NP2]={0.00508717  ,1.25621e-003,  8.72499e-004 ,1.07130e-003,1.31250e-003,0.00105332,0.0,0.0};
  // Double_t bes2syserr_excess_p[NP2]={0.006559631     ,0            ,  0.0 , 0.000473043  ,0.0,0.000781936,0.0,0.0};
  
  // const Int_t NP2 = 6;
  // const Int_t NP2_emmy = 4;
  // const Double_t bes2_energy[NP2] = {7.7,9.2, 14.6, 19.6, 27., 54.4};
  // Double_t yb2[NP2];  // beam rapidity
  // const Double_t bes2val_pv1[NP2] = {1.37180e-002  ,1.24595e-003,  -0.00578296, -0.00823163,0.0,0.0};
  // const Double_t bes2staerr_pv1[NP2] = {2.01159e-004  ,9.35210e-005, 0.000161092, 0.000165627,0.0,0.0};
  // const Double_t bes2syserr_pv1[NP2] = {0.00000000 ,0.0, 0.0000     , 0.000184515,0.0,0.0};

  // const Double_t bes2val_pbarv1[NP2] = {-5.15323e-002,-4.23434e-002 , -0.0331744, -0.0311547,0.0,0.0};
  // const Double_t bes2staerr_pbarv1[NP2] = {1.85570e-003,7.50955e-004,  0.000668851, 0.000499875,0.0,0.0};
  // const Double_t bes2syserr_pbarv1[NP2] = {0.00000000,  0.0            ,0.00000000, 0.000479769,0.0,0.0};
      
  // Double_t bes2val_excess_p[NP2]={0.0665038   ,4.42068e-002,2.72578e-002 ,0.0257154 ,0.0,0.0};
  // Double_t bes2staerr_excess_p[NP2]={0.00275253 ,1.51354e-003,1.07130e-003,0.00105332,0.0,0.0};
  // Double_t bes2syserr_excess_p[NP2]={0.006559631,0.0        , 0.000470764  ,0.001054285,0.0,0.0};
  
  const Int_t NP2 = 8;
  const Int_t NP2_emmy = 7;
  Double_t bes2_energy[NP2] =       {7.7        ,9.2        ,11.5      ,14.6       ,17.3       ,19.6       ,27.0         ,54.4};
  Double_t yb2[NP2];  // beam rapidity
  Double_t bes2val_pv1[NP2] =       {0.0130227  ,0.00302743 ,-0.00510691,-0.00739349,-0.00980025,-0.00960186,-0.0064442 , 0.0};
  Double_t bes2staerr_pv1[NP2] =    {0.00016079 ,0.000141831,0.000132198,0.000126254,0.000175722,0.000166092,0.000147467, 0.0};
  Double_t bes2syserr_pv1[NP2] =    {4.59E-03   ,3.78E-03   ,1.79E-03   ,0.000542573,7.78E-04   ,0.00052021 ,5.51E-04   , 0.0};

  Double_t bes2val_pbarv1[NP2] =    {-0.0520613 ,-0.0482222 ,-0.0412019 ,-0.0330864 ,-0.0329383 ,-0.0314665 ,-0.0258247 , 0.0};
  Double_t bes2staerr_pbarv1[NP2] = {0.00200259 ,0.00121533 ,0.000761715,0.000520525,0.000598985,0.000463355,0.000340297, 0.0};
  Double_t bes2syserr_pbarv1[NP2] = {2.13E-03   ,1.33E-03   ,1.27E-03   ,0.000625344,7.13E-04   ,0.000556806,1.99E-03   , 0.0};
 

 
  // Double_t bes2val_excess_p[NP2]=         {0.0653572  ,0.0518391  ,0.0371621  ,0.0272647  ,0.0252723  ,0.0250366  ,0.0238357  , 0.0};// 7.7 antiproton pid tracks
  // Double_t bes2staerr_excess_p[NP2]=      {0.00404388 ,0.00248054 ,0.00159364 ,0.00113734 ,0.00136454 ,0.00112728 ,0.000911952, 0.0};// 7.7 antiproton pid tracks
  // Double_t bes2syserr_excess_p[NP2]=      {9.49E-04   ,7.04E-04   ,1.10E-03   ,0.000217329,8.80E-04   ,0.000760106, 6.88E-04  , 0.0};// 7.7 antiproton pid tracks
  Double_t bes2val_excess_p[NP2]=         {0.0665038  ,0.0518391  ,0.0371621 ,2.72578e-002,0.0252723  ,0.0257154  ,0.0238357  , 0.0};
  Double_t bes2staerr_excess_p[NP2]=      {0.00275253 ,0.00248054 ,0.00159364,1.07130e-003,0.00136454 ,0.00105332 ,0.000911952, 0.0};
  Double_t bes2syserr_excess_p[NP2]=      {0.006559631,7.04E-04   ,1.10E-03  ,0.000470764 ,8.80E-04   ,0.001054285,6.88E-04 , 0.0};
  
  
  // Double_t v1_difference[NP2];
  // Double_t v1_difference_error[NP2];
  
  

  TFile *fin = new TFile("DataProtonV1TwoComp_27nG4p4GeV.root");  // 27 and 54.4 GeV data from Sooraj
  TGraphErrors *gr_SR = (TGraphErrors *)fin->Get("gSlope_excessv1_c1");
  
  Double_t r2[NP2]; // BES-2 ratio
  Double_t bes2_energy_d[NP2];
  for(int i=0;i<NP2;i++) {
    bes2_energy_d[i] = bes2_energy[i];
    r2[i] = gr_pbarp->Eval(bes2_energy[i]);
    //yb2[i] = rapidity(bes2_energy[i]);
	yb2[i] = TMath::ACosH(bes2_energy[i]/protonMass);

	// v1_difference[i]=(bes2val_pv1[i]-bes2val_pbarv1[i])*yb2[i];
	// v1_difference_error[i] =(bes2staerr_pv1[i]+bes2staerr_pbarv1[i])*yb2[i];

    if(i<NP2_emmy) { // Data from Emmy
	  bes2val_pv1[i]=bes2val_pv1[i]*yb2[i];
	  bes2staerr_pv1[i]=bes2staerr_pv1[i]*yb2[i];
	  bes2syserr_pv1[i]=bes2syserr_pv1[i]*yb2[i];
	  
	  bes2val_pbarv1[i]=bes2val_pbarv1[i]*yb2[i];
	  bes2staerr_pbarv1[i]=bes2staerr_pbarv1[i]*yb2[i];
	  bes2syserr_pbarv1[i]=bes2syserr_pbarv1[i]*yb2[i];
	  
	  bes2val_excess_p[i]=bes2val_excess_p[i]*yb2[i];
	  bes2staerr_excess_p[i]=bes2staerr_excess_p[i]*yb2[i];
	  bes2syserr_excess_p[i]=bes2syserr_excess_p[i]*yb2[i];
    } else { // Data from Sooraj
      bes2val_excess_p[i] = gr_SR->GetY()[i-NP2_emmy] * yb2[i];
      bes2staerr_excess_p[i] = gr_SR->GetEY()[i-NP2_emmy] * yb2[i];
      bes2syserr_excess_p[i] = 0;      
    }
    cout << " Energy: " << bes2_energy[i] << "\t pbar/p = " << r2[i] << "\t yb = " << yb2[i] << "\t Excess V1 = " << bes2val_excess_p[i]/yb2[i] << " +/- " << bes2staerr_excess_p[i] << " +/- " << bes2syserr_excess_p[i] << endl;
      //    cout << " Energy: " << bes2_energy[i] << "\t pbar/p = " << r2[i] << "\t yb = " << yb2[i] << "\t Excess V1 = " << bes2val_excess_p_2[i] << " +/- " << bes2staerr_excess_p_2[i] << " +/- " << bes2syserr_excess_p_2[i] << endl;
      //    cout << "---" << endl;
  }
  // TGraphErrors *gr_bes2_p_stat = new TGraphErrors(NP2_emmy, bes2_energy_d, bes2val_excess_p, 0, bes2staerr_excess_p);
  // TGraphErrors *gr_bes2_p_sys = new TGraphErrors(NP2_emmy, bes2_energy_d, bes2val_excess_p, 0, bes2syserr_excess_p);
  
  TGraphErrors *gr_bes2_stat = new TGraphErrors(NP2_emmy, bes2_energy_d, bes2val_excess_p, 0, bes2staerr_excess_p);
  TGraphErrors *gr_bes2_sys = new TGraphErrors(NP2_emmy, bes2_energy_d, bes2val_excess_p, 0, bes2syserr_excess_p);
  
  
  // TGraphErrors *gr_bes2_diff = new TGraphErrors(NP2, bes2_energy_d,v1_difference,0, v1_difference_error);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Models - version before April 2024
  // const Int_t NM = 2;
  // const Int_t NEMax = 10;
  // TFile *fin1 = new TFile("JAMProtonV1TwoComp.root");
  // const Char_t *GraphName[NM] = {"gExcessSlope_casc_c1","gExcessSlope_mfrmf_c1"};
  // const Char_t *LabelName[NM] = {"JAM Cascade","JAM Mean Field"};
  // TGraph *gr_m[NM], *gr_m_tmp[NM];
  // Double_t model_energy[NM][NEMax], model_v1[NM][NEMax];
  // for(int i=0;i<NM;i++) {
  //   gr_m_tmp[i] = (TGraph *)fin1->Get(GraphName[i]);
  //   for(int j=0;j<gr_m_tmp[i]->GetN();j++) {
  //     model_energy[i][j] = gr_m_tmp[i]->GetX()[j];
  //     // double y_m = rapidity(model_energy[i][j]);
  // 	  double y_m = TMath::ACosH(model_energy[i][j]/protonMass);
  //     model_v1[i][j] = gr_m_tmp[i]->GetY()[j] * y_m;
  //   }
  //   gr_m[i] = new TGraph(gr_m_tmp[i]->GetN(), model_energy[i], model_v1[i]);
  //   gr_m[i]->Print();
  // }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Models - version CPOD 2024
  const Int_t NM = 4;
  const Int_t NEMax = 10;
  TFile *fin1 = new TFile("JAMTwoComponentBES2_eosscan.root");
  const Char_t *GraphName[NM] = {"momdep_hard_K380MeV", "momdep_soft_K210MeV", "hard_K380MeV", "soft_K210MeV"};
  const Char_t *LabelName[NM] = {"p-dep #kappa = 380 MeV", "p-dep #kappa = 210 MeV", "#kappa = 380 MeV", "#kappa = 210 MeV"};
  TGraphErrors *gr_m[NM], *gr_m_tmp[NM];
  Double_t model_energy[NM][NEMax], model_v1[NM][NEMax];
  for(int i=0;i<NM;i++) {
    // gr_m_tmp[i] = (TGraph *)fin1->Get(GraphName[i]);
    // for(int j=0;j<gr_m_tmp[i]->GetN();j++) {
    //   model_energy[i][j] = gr_m_tmp[i]->GetX()[j];
    //   // double y_m = rapidity(model_energy[i][j]);
    //   //	  double y_m = TMath::ACosH(model_energy[i][j]/protonMass);
    //   //      model_v1[i][j] = gr_m_tmp[i]->GetY()[j] * y_m;
    // }
    // gr_m[i] = new TGraph(gr_m_tmp[i]->GetN(), model_energy[i], model_v1[i]);
    gr_m[i] = (TGraphErrors *)fin1->Get(Form("g_proton_excessv1_slopevsyoverybeam_1040_%s", GraphName[i]));
    //    gr_m[i] = new TGraph(gr_m_tmp[i]->GetN(), gr_m_tmp[i]->GetX(), gr_m_tmp[i]->GetY());
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd()->SetLogx();
  c1->cd()->SetLogy();

  double x1 = 3.5;
  double x2 = 3e2;
  double y1 = 0.05;
  double y2 = 0.25;
  if(conf) y2 = 0.7;
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->GetXaxis()->CenterTitle();
  h0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  h0->GetXaxis()->SetLabelOffset(999.);
  h0->GetXaxis()->SetLabelSize(0.5);
  h0->GetXaxis()->SetTitleSize(0.065);
  h0->GetXaxis()->SetTitleOffset(1.2);
  h0->SetYTitle("dv_{1}^{excess-p}/d(y/y_{beam})");
  h0->GetYaxis()->SetTitleSize(0.065);
  h0->GetYaxis()->SetTitleOffset(1.1);
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetLabelSize(0.5);
  h0->GetYaxis()->SetLabelOffset(999.);
  h0->SetMaximum(y2);
  h0->SetMinimum(y1);
  h0->Draw("c");

  drawLine(x1*1.0, 0.095, x2*1.0, 0.095, 2, 8, 1);
  // drawLine(x1*1.1, 0.078, x2*0.9, 0.078, 2, 8, 1);
  // drawLine(x1*1.0, bes1val_excess_p[NP1-1], x2*1.0, bes1val_excess_p[NP1-1], 2, 8, 1);
  // drawLine(x1*1.0, bes1val_excess_p[7]+bes1staerr_excess_p[7], x2*1.0, bes1val_excess_p[7]+bes1staerr_excess_p[7], 2, 1, 1);
  // drawLine(x1*1.0, bes1val_excess_p[7]-bes1staerr_excess_p[7], x2*1.0, bes1val_excess_p[7]-bes1staerr_excess_p[7], 2, 1, 1);
  //  const Int_t m_lineColor[NM] = {kRed+2, kRed+2, kGreen+4, kGreen+4};
  const Int_t m_fillColor[NM] = {kRed, kRed, kGreen+1, kGreen+1};
  const Int_t m_lineStyle[NM] = {2, 3, 2, 3};
  const Int_t m_fillStyle[NM] = {3001, 3244, 3001, 3244};
  int include_model = 1;
  if(include_model)
  {
	  for(int i=0;i<NM;i++) {
	    setGraphLine(gr_m[i], m_lineStyle[i], m_fillColor[i], 2);
	    setGraphFill(gr_m[i], m_fillStyle[i], m_fillColor[i], 0.5);
	    if(conf) gr_m[i]->Draw("e3");
	  }
  }

  // setGraphMarker(gr_bes1_sys, 24, 1, 1.5);
  // setGraphLine(gr_bes1_sys, 1, 1, 2);
  // drawSysBox(gr_bes1_sys, 0.08, 18, 1);

  setGraphMarker(gr_bes1_stat, 24, 1, 1.5);
  setGraphLine(gr_bes1_stat, 1, 1, 2);
  gr_bes1_stat->Draw("p");
  
  setGraphMarker(gr_bes2_sys, 20, 2, 1.5);
  setGraphLine(gr_bes2_sys, 1, 2, 2);
  drawSysBox(gr_bes2_sys, 0.05, 16, 1);

  setGraphMarker(gr_bes2_stat, 20, 2, 1.8);
  setGraphLine(gr_bes2_stat, 1, 2, 2);
  gr_bes2_stat->Draw("p");
  
  
  // gr_bes2_diff->SetName("difference");gr_bes2_diff->SetTitle("difference");
  // gr_bes2_diff->Draw("p");

  //  gr_bes2_stat->Fit("pol0","R","",10.,100.);
  
  drawHistBox(x1,x2,y1,y2);
  // drawText(9.0, -0.06*(y2-y1)+y1, "10");
  // drawText(86, -0.06*(y2-y1)+y1, "100");
  drawText(9.0, 0.042, "10");
  drawText(28,  0.042, "30");
  drawText(86,  0.042, "100");
  drawText(x1*0.72, 0.095, "0.1");
  drawText(x1*0.72, 0.19, "0.2");
  drawText(x1*0.72, 0.475, "0.5");

  //  drawText(30, y2*0.83, "Excess Flow of Protons", 32, 0.06);
  drawText(50, y2*0.78, "Au+Au 10-40%", 42, 0.05);
  
  TLegend *leg = new TLegend(0.7, 0.75, 0.92, 0.88);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_bes2_stat, "  BES-II", "pl");
  leg->AddEntry(gr_bes1_stat, "  BES-I", "pl");
  leg->Draw();
  if(include_model)
  {
	  leg = new TLegend(0.68, 0.74-0.06*NM, 0.92, 0.74);
	  leg->SetLineColor(10);
	  leg->SetTextSize(0.035);
	  for(int i=0;i<NM;i++) {
		leg->AddEntry(gr_m[i], LabelName[i], "f");
	  }
	  if(conf) leg->Draw();
  }
  
  drawText(35, 0.4, "STAR");
  drawText(35, 0.22, "JAM2");

  
  c1->Update();
  c1->SaveAs(Form("excess_v1_%d.pdf",conf));  
  c1->SaveAs(Form("excess_v1_%d.png",conf));
  
}
