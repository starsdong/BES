#include "draw.C+"
#include "style.C"

void compError()
{
  style();

  //////////=====
  const int n_bes2 = 7;
  float C42_proj_err[n_bes2] = {0.218041,  0.168534,  0.1316233,  0.0938992,  0.1117,  0.0655,  0.0860338};

  const int ntot_datapts=11;
  ////data 0-5%
  float C42_ener_besNEW[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  float C42_main_Ref3X_besNEW[ntot_datapts]={0.4119085,0.540458,0.4011784,0.430751,0.31512,0.339941,0.602046,0.739693,0.632837,0.792955,0.900669};
  float C42_stat_Ref3X_besNEW[ntot_datapts]={0.2468968,0.187029,0.1347458,0.104085,0.120367,0.0731331,0.0950052,0.147006,0.0553742,0.250823,0.2084582};
  float C42_sys_Ref3X_besNEW[ntot_datapts]={0.1302093,0.0680524,0.06472035,0.0403519,0.0880467,0.0389432,0.0304235,0.1357538,0.1387023,0.1219,0.1393631};
  float C42_tot_Ref3X_err_besNEW[ntot_datapts]={0};
  
  
  for(int y=0;y<ntot_datapts;y++)
    {
      C42_tot_Ref3X_err_besNEW[y]=sqrt(C42_stat_Ref3X_besNEW[y]*C42_stat_Ref3X_besNEW[y]+C42_sys_Ref3X_besNEW[y]*C42_sys_Ref3X_besNEW[y]);
      cout << C42_tot_Ref3X_err_besNEW[y] << endl;
    }

  TGraphErrors *gr_data_bes2 = new TGraphErrors(n_bes2, C42_ener_besNEW, C42_main_Ref3X_besNEW, 0, C42_stat_Ref3X_besNEW);
  TGraphErrors *gr_proj_bes2 = new TGraphErrors(n_bes2, C42_ener_besNEW, C42_main_Ref3X_besNEW, 0, C42_proj_err);


  const int ntot_besI_ref3=9;
  float ener_BESIref3[ntot_besI_ref3] = {7.7*0.95,11.5*0.95,14.6*0.95,19.6*0.95,27,39,54.4,62.4,200};
  float CR42_NetP_BESIref3cen5_data[ntot_besI_ref3]={1.76697,0.695796,1.46876,0.140553,0.196254,0.739693,0.632837,0.792955,0.900669};
  float CR42_NetP_BESIref3cen5_stat[ntot_besI_ref3]={1.15082,0.595077,0.395972,0.328768,0.219415,0.147006,0.0553742,0.250823,0.2084582};
  float CR42_NetP_BESIref3cen5_sys[ntot_besI_ref3]={0.413944,0.266722,0.2540532,0.1559867,0.1431581,0.1357538,0.1387023,0.1219,0.1393631};

  TGraphErrors *gr_data_bes1 = new TGraphErrors(ntot_besI_ref3, ener_BESIref3, CR42_NetP_BESIref3cen5_data, 0, CR42_NetP_BESIref3cen5_stat);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogx();

  TH1D *d0 = new TH1D("d0","",1,3, 300);
  d0->SetMinimum(-0.5);
  d0->SetMaximum(2.9);
  d0->GetXaxis()->CenterTitle();
  d0->SetXTitle("Collision Energy #sqrt{s_{NN}} (GeV)");
  d0->GetXaxis()->SetLabelOffset(0.000);
  d0->GetXaxis()->SetTitleOffset(1.1);  
  d0->SetYTitle("Net-proton C_{4}/C_{2}");
  d0->GetYaxis()->SetTitleOffset(1.);  
  d0->Draw("c");

  drawLine(3, 0, 300, 0, 1, 9, 16);
  drawLine(3, 1, 300, 1, 1, 9, 16);
  

  gr_proj_bes2->SetFillColor(kGreen);
  drawSysBox(gr_proj_bes2, 0.05, kGreen, 1);

  setGraphMarker(gr_data_bes1, 24);
  setGraphLine(gr_data_bes1);
  gr_data_bes1->Draw("p");

  setGraphMarker(gr_data_bes2);
  setGraphLine(gr_data_bes2);
  gr_data_bes2->Draw("p");

  drawHistBox(3, 300, -0.5, 2.9);


  drawText(60, 2.5, "Au+Au 0-5%");
  TLegend *leg = new TLegend(0.7, 0.66, 0.92, 0.84);
  leg->SetLineColor(10);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr_data_bes1, "  BES-I data", "pl");
  leg->AddEntry(gr_data_bes2, "  BES-II data", "pl");
  leg->AddEntry(gr_proj_bes2, "  BES-II proj.", "f");
  leg->Draw();
  
  c1->Update();
  c1->SaveAs("compError.png");  
  c1->SaveAs("compError.pdf");

}
