#include "draw.C+"
#include "style.C+"

void Fig14_knk1_Npart()
{
  style();

  ////////////////////////////////////////////////
  // Define constants
  ////////////////////////////////////////////////
  const Int_t NE = 9;  // All energies
  const Int_t index_54 = 6;  // 54.4 GeV index in energy array
  const Double_t Ene[NE] = {7.7, 11.5, 14.5, 19.6, 27, 39, 54.4, 62.4, 200};
  const Char_t *EneLabel[NE] = {"7.7", "11.5", "14.5", "19.6", "27", "39", "54.4", "62.4", "200"};
  const Char_t *EneDir[NE] = {"7", "11", "15", "19", "27", "39", "54", "62", "200"};  // directory names
  const Int_t NCum = 3; // 3 orders of kappa ratios
  const Int_t NP = 2; // number of particle categories: proton, anti-proton
  const Char_t *PName[NP] = {"Pro", "Apro"};
  const Char_t *PName_54[NP] = {"pro", "antipro"};
  const Int_t NCen = 9; // 9 centrality bins
  const double cen_offset[NP] = {3, 0}; // different rapidity offsets for different particles for plotting  
  
  ////////////////////////////////////////////////
  // Read in data: 54 GeV stored differently
  ////////////////////////////////////////////////
  TGraphErrors *gr_stat[NE][NCum][NP], *gr_sys[NE][NCum][NP];  // all central collisions
  double cen[NP][NCen];
  double cum[NE][NCum+1][NP][NCen], cum_e[NE][NCum+1][NP][NCen], cum_e_sys[NE][NCum+1][NP][NCen];   // name is cum, actually are kappa in this Fig.
  double cumR[NE][NCum][NP][NCen], cumR_e[NE][NCum][NP][NCen], cumR_e_sys[NE][NCum][NP][NCen]; // ratios to k1
  TFile *fin[NE];
  for(int i=0;i<NE;i++) {
    if(i!=index_54) {
      fin[i] = new TFile(Form("fig14_16_19/fig14_corrfun_centrality/%sGeV_centrality.root",EneDir[i]));
      for(int j=0;j<NP;j++) {
	//	fin[i][j] = new TFile(Form("rootfile_0517/final_corrfun_ycut/%sGeV/AuAu_sys_%s_Y0.5.root",EneDir[i],PName[j]));
	TGraphErrors *gr_stat_a[NCum], *gr_sys_a[NCum];
	for(int m=0;m<NCum;m++) {
	  gr_stat_a[m] = (TGraphErrors *)fin[i]->Get(Form("%s_K%d1_stat",PName[j],m+2));
	  gr_sys_a[m] = (TGraphErrors *)fin[i]->Get(Form("%s_K%d1_sys",PName[j],m+2));
	  
	  for(int k=0;k<NCen;k++) {
	    cen[j][k] = gr_stat_a[m]->GetX()[k] + cen_offset[j];
	    cumR[i][m][j][k] = gr_stat_a[m]->GetY()[k];
	    cumR_e[i][m][j][k] = gr_stat_a[m]->GetEY()[k];
	    cumR_e_sys[i][m][j][k] = gr_sys_a[m]->GetEY()[k];
	    
	    if(i==0&&m==NCum-1) { // 7.7 GeV C4 data points scaled down by 0.5
	      cumR[i][m][j][k] *= 0.5;
	      cumR_e[i][m][j][k] *= 0.5;
	      cumR_e_sys[i][m][j][k] *= 0.5;
	    }
	  } // end k->NCen
	} // end m->NCum
      } // end j->NP
      fin[i]->Close();

      for(int j=0;j<NP;j++) {
	for(int m=0;m<NCum;m++) {
	  gr_stat[i][m][j] = new TGraphErrors(NCen, cen[j], cumR[i][m][j], 0, cumR_e[i][m][j]);
	  gr_sys[i][m][j] = new TGraphErrors(NCen, cen[j], cumR[i][m][j], 0, cumR_e_sys[i][m][j]);
	  cout << "++ " << EneDir[i] << " " << PName[j] << Form(" k%d/k1",m+2) << endl;
	  gr_stat[i][m][j]->Print();
	}
      } // end j->NP
    } else { // 54 GeV data
      TFile *fin_54 = new TFile("rootfile_0517/54GeV/54GeV_CORRELATION_MAY17.root");
      for(int j=0;j<NP;j++) {
	TGraphErrors *gr_stat_b[NCum+1], *gr_sys_b[NCum+1];
	for(int m=0;m<NCum+1;m++) {
	  gr_stat_b[m] = (TGraphErrors *)fin_54->Get(Form("centrality_%s_correl_K%d_stat",PName_54[j], m+1));
	  gr_sys_b[m] = (TGraphErrors *)fin_54->Get(Form("centrality_%s_correl_K%d_sys",PName_54[j], m+1));
	  
	  for(int k=0;k<NCen;k++) {
	    cen[j][k] = gr_stat_b[m]->GetX()[k] + cen_offset[j];
	    cum[i][m][j][k] = gr_stat_b[m]->GetY()[k];
	    cum_e[i][m][j][k] = gr_stat_b[m]->GetEY()[k];
	    cum_e_sys[i][m][j][k] = gr_sys_b[m]->GetEY()[k];
	  }
	}
	for(int m=0;m<NCum;m++) {
	  for(int k=0;k<NCen;k++) {
	    cumR[i][m][j][k] = cum[i][m+1][j][k]/cum[i][0][j][k];
	    cumR_e[i][m][j][k] = fabs(cumR[i][m][j][k]) * sqrt(pow(cum_e[i][m+1][j][k]/cum[i][m+1][j][k], 2.0)+pow(cum_e[i][0][j][k]/cum[i][0][j][k], 2.0));
	    cumR_e_sys[i][m][j][k] = fabs(cumR[i][m][j][k]) * sqrt(pow(cum_e_sys[i][m+1][j][k]/cum[i][m+1][j][k], 2.0)+pow(cum_e_sys[i][0][j][k]/cum[i][0][j][k], 2.0));
	  }
	  
	  gr_stat[i][m][j] = new TGraphErrors(NCen, cen[j], cumR[i][m][j], 0, cumR_e[i][m][j]);
	  gr_sys[i][m][j] = new TGraphErrors(NCen, cen[j], cumR[i][m][j], 0, cumR_e_sys[i][m][j]);
	} // end m->NCum
	
      } // end j->NP    
    } // end if(i!=index_54)
  } // end i->NE
  
  ////////////////////////////////////////////////
  // Start Plotting
  ////////////////////////////////////////////////
  TCanvas *c1 = new TCanvas("c1","",1400,600);
  TPad *pad[NE][NCum];
  TH1D *h0[NE][NCum];
  const double pad_x1 = 0.06;
  const double pad_y1 = 0.10; // x,y pad offsets
  const double pad_x_gap = 0.; // gap between 27 GeV and 39 GeV (index 4-5)
  const double pad_x2 = 0.99;
  const double pad_y2 = 0.98;
  const int pad_gap_index = 5; // starting from index 5 (39 GeV)
  const double ymin[NCum] = {-0.11, -0.165, -0.95};
  const double ymax[NCum] = {0.022, 0.142, 1.5};
  const double xmin = -20;
  const double xmax = 380;
  const Int_t kColor[NP] = {kBlue, kBlack};
  const Int_t kStyle[NP] = {21, 34};
  const Double_t kSize[NP] = {0.8, 1.2};

  ////////////////////////////////////////////////
  // Draw main panels
  ////////////////////////////////////////////////
  for(int i=0;i<NE;i++) {
    for(int j=0;j<NCum;j++) {
      c1->cd();
      
      double x1 = pad_x1 + i*(pad_x2 - pad_x1 - pad_x_gap)/NE;
      double x2 = pad_x1 + (i+1)*(pad_x2 - pad_x1 - pad_x_gap)/NE;
      if(j>=pad_gap_index) {
	x1 += pad_x_gap;
	x2 += pad_x_gap;
      }
      double y1 = pad_y2 - (j+1)*(pad_y2 - pad_y1)/NCum;
      double y2 = pad_y2 - j*(pad_y2 - pad_y1)/NCum;
      pad[i][j] = new TPad(Form("pad_%d_%d",i,j),"",x1,y1,x2,y2);
      pad[i][j]->SetLeftMargin(0.001);
      pad[i][j]->SetRightMargin(0.001);
      pad[i][j]->SetTopMargin(0.001);
      pad[i][j]->SetBottomMargin(0.001);
      pad[i][j]->Draw();
      pad[i][j]->cd();

      h0[i][j] = new TH1D(Form("h0_%d_%d",i,j),"",1,xmin,xmax);
      h0[i][j]->SetMaximum(ymax[j]);
      h0[i][j]->SetMinimum(ymin[j]);
      h0[i][j]->GetXaxis()->SetNdivisions(504);
      h0[i][j]->GetXaxis()->SetTitleOffset(999.);
      h0[i][j]->GetXaxis()->SetLabelOffset(999.);
      h0[i][j]->GetYaxis()->SetNdivisions(503);
      h0[i][j]->GetYaxis()->SetTitleOffset(999.);
      h0[i][j]->GetYaxis()->SetLabelOffset(999.);
      h0[i][j]->Draw("c");
      drawLine(xmin,0,xmax,0,1,2,kBlack);

      for(int k=0;k<NP;k++) {
	cout << "++ " << EneDir[i] << " " << PName[k] << Form(" k%d/k1",j+2) << endl;
	gr_stat[i][j][k]->Print();
	gr_sys[i][j][k]->Print();

	drawSysError(gr_sys[i][j][k], 10, (ymax[j]-ymin[j])*0.02, kColor[k]);

	setGraphMarker(gr_stat[i][j][k], kStyle[k], kColor[k], kSize[k]);
	setGraphLine(gr_stat[i][j][k], 1, kColor[k], 1);
	gr_stat[i][j][k]->Draw("p");

      }

      /*
      TGraphErrors *gr_clone = (TGraphErrors *)gr_stat[i][j][NP-1]->Clone(Form("%s_clone",gr_stat[i][j][NP-1]->GetName()));
      // net-p replotting
      setGraphMarker(gr_clone, 24, kBlack, kSize[NP-1]+0.15);
      setGraphLine(gr_clone, 1, kRed, 1);
      gr_clone->Draw("p");
      */

      if(j==0) {
	drawText(100, (ymax[j]-ymin[j])*0.88+ymin[j], Form("%s GeV", EneLabel[i]), 42, 0.14);
      }

      if(i==0&&j==NCum-1) { // 7.7 GeV C4 legend
	drawText(20, 1.15, "#kappa_{4}/#kappa_{1} #times 0.5", 42, 0.14);
	drawText(20, 0.8, "7.7 GeV only", 42, 0.14);
      }
      
      pad[i][j]->Update();
    }
  }

  ////////////////////////////////////////////////
  // Draw axis labels and titles
  ////////////////////////////////////////////////
  c1->cd();
  TPad *pad_x = new TPad("pad_x", "", pad_x1, 0, pad_x2, pad_y1-0.003);
  pad_x->Draw();
  pad_x->cd();

  for(int i=0;i<NE;i++) {
    drawText(i*1./NE + 0.002, 0.75, "0", 42, 0.35);
    drawText(i*1./NE + 0.022, 0.75, "100", 42, 0.35);
    drawText(i*1./NE + 0.049, 0.75, "200", 42, 0.35);
    drawText(i*1./NE + 0.077, 0.75, "300", 42, 0.35);
  }
  drawText(0.33, 0.25, "Average Number of Participant Nucleons <N_{part}>", 42, 0.43);
  pad_x->Update();
  
  c1->cd();
  TPad *pad_y = new TPad("pad_y", "", 0, pad_y1, pad_x1-0.002, pad_y2);
  pad_y->Draw();
  pad_y->cd();

  const Int_t NYLabelMax = 5;
  const Int_t NYLabel[NCum] = {3, 3, 2};
  const Char_t *YLabel[NCum][NYLabelMax] = {{"#font[122]{-} 0.10","#font[122]{-} 0.05","   0.00","",""},
					    {"#font[122]{-} 0.1","   0.0","   0.1","",""},
					    {"   0","   1","",""}
  };
  const Double_t XOffset[NCum] = {0.01, 0.12, 0.27};
  const Double_t Offset[NCum] = {0.015, 0.077, 0.120};
  const Double_t Step[NCum] = {0.125, 0.093, 0.132}; 
  for(int i=0;i<NCum;i++) {
    for(int j=0;j<NYLabel[i];j++) {
      drawText(0.45+XOffset[i], (NCum-1-i)*1./NCum + j*Step[i] + Offset[i], YLabel[i][j], 42, 0.22);
    }
    drawText(0.35, (NCum-1-i+0.4)*1./NCum, Form("#kappa_{%d}/#kappa_{1}",i+2), 42, 0.3, 90);
  }
  pad_y->Update();

  ////////////////////////////////////////////////
  // Draw other legends and labels
  ////////////////////////////////////////////////
  c1->cd();
  TPad *pad_coll = new TPad("pad_coll","", 0.425, 0.405, 0.66, 0.495);
  pad_coll->Draw();
  pad_coll->cd();

  drawText(0.10, 0.65, "STAR Au+Au Collisions", 42, 0.48);
  drawText(0.02, 0.18, "0.4 < p_{T} < 2.0 (GeV/c), |y| < 0.5", 42, 0.45);
  pad_coll->Update();

  c1->cd();
  TPad *pad_leg = new TPad("pad_leg", "", 0.42, 0.29, 0.63, 0.35);
  pad_leg->Draw();
  pad_leg->cd();

  TLegend *leg[NP];
  const Int_t index_leg[NP] = {0, 1};
  const Char_t *name_leg[NP] = {" Proton", " Anti-proton"};
  const Double_t x1_leg[NP] = {0.1, 0.5};
  const Double_t x2_leg[NP] = {0.45, 0.85};
  TGraphErrors *gr_leg[NP];
  for(int i=0;i<NP;i++) {
    leg[i] = new TLegend(x1_leg[i], 0.1, x2_leg[i], 0.9);
    leg[i]->SetLineColor(10);
    leg[i]->SetTextSize(0.7);
    gr_leg[i] = (TGraphErrors *)gr_stat[0][0][index_leg[i]]->Clone();
    gr_leg[i]->SetMarkerSize(1.4);
    leg[i]->AddEntry(gr_leg[i], name_leg[i], "p");
    leg[i]->Draw();
  }
  
  pad_leg->Update();
  c1->cd();

  c1->Update();
  c1->SaveAs("fig/Fig14_knk1_Npart.pdf");
  c1->SaveAs("fig/Fig14_knk1_Npart.png");
  
  
}
