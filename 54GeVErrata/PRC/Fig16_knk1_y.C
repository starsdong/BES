#include "draw.C+"
#include "style.C+"

void Fig16_knk1_y()
{
  style();

  ////////////////////////////////////////////////
  // Define constants
  ////////////////////////////////////////////////
  const Int_t NE = 9;  // All energies
  const Int_t index_54 = 6;  // 54.4 GeV index in energy array
  const Double_t Ene[NE] = {7.7, 11.5, 14.5, 19.6, 27, 39, 54.4, 62.4, 200};
  const Char_t *EneLabel[NE] = {"7.7", "11.5", "14.5", "19.6", "27", "39", "54.4", "62.4", "200"};
  //  const Char_t *EneDir[NE] = {"7", "11", "15", "19", "27", "39", "54", "62", "200"};  // directory names
  const Char_t *EneDir[NE] = {"7.7", "11.5", "14.5", "19.6", "27", "39", "54", "62.4", "200"};  // directory names
  const Int_t NCum = 3; // 3 orders of kappa ratios
  const Int_t NP = 2; // number of particle categories: proton, anti-proton
  const Char_t *PName[NP] = {"Pro", "Apro"};
  //  const Char_t *PName_54[NP] = {"pro", "antipro"};
  const Char_t *PName_54[NP] = {"Pro", "Apro"};
  const Char_t *PName_His_54[NP] = {"Pro", "Pbar"};
  const Int_t NY = 5;  // rapidity bins
  const Double_t RAP[NY] = {0.1, 0.2, 0.3, 0.4, 0.5};
  const Int_t NCen = 9; // 9 centrality bins
  const double rap_offset[NP] = {0.008, 0.000}; // different rapidity offsets for different particles for plotting
  double rap[NP][NY];
  for(int i=0;i<NP;i++) {
    for(int j=0;j<NY;j++) {
      rap[i][j] = RAP[j] + rap_offset[i];
    }
  }
  
  
  ////////////////////////////////////////////////
  // Read in data: 54 GeV stored differently
  ////////////////////////////////////////////////
  TGraphErrors *gr_stat[NE][NCum][NP], *gr_sys[NE][NCum][NP];  // all central collisions
  double cum[NE][NCum+1][NP][NY], cum_e[NE][NCum+1][NP][NY], cum_e_sys[NE][NCum+1][NP][NY];   // name is cum, actually are kappa in this Fig.
  double cumR[NE][NCum][NP][NY], cumR_e[NE][NCum][NP][NY], cumR_e_sys[NE][NCum][NP][NY]; // ratios to k1
  TFile *fin[NE];
  for(int i=0;i<NE;i++) {
    if(i!=index_54) {
      fin[i] = new TFile(Form("fig14_16_19/fig16_corrfun_ycut/%sGeV_rapidity_scan.root",EneDir[i]));
      for(int j=0;j<NP;j++) {
	//	  fin[i][j][k] = new TFile(Form("rootfile_0517/final_corrfun_ycut/%sGeV/AuAu_sys_%s_Y%3.1f.root",EneDir[i],PName[j],RAP[k]));
	TGraphErrors *gr_stat_a[NCum], *gr_sys_a[NCum];
	for(int m=0;m<NCum;m++) {
	  gr_stat_a[m] = (TGraphErrors *)fin[i]->Get(Form("%s_K%d1_stat",PName[j],m+2));
	  gr_sys_a[m] = (TGraphErrors *)fin[i]->Get(Form("%s_K%d1_sys",PName[j],m+2));

	  for(int k=0;k<NY;k++) {
	    cumR[i][m][j][k] = gr_stat_a[m]->GetY()[k];
	    cumR_e[i][m][j][k] = gr_stat_a[m]->GetEY()[k];
	    cumR_e_sys[i][m][j][k] = gr_sys_a[m]->GetEY()[k];
	    //	    cout << EneDir[i] << " " << PName[j] << " " << rap[k] << Form(" C%d",m+1) << " " << cum[i][m][j][k] << endl;
	    
	    // cumR[i][m][j][k] = cum[i][m][j][k]/cum[i][0][j][k];
	    // cumR_e[i][m][j][k] = cumR[i][m][j][k] * sqrt(pow(cum_e[i][m][j][k]/cum[i][m][j][k], 2.0)+pow(cum_e[i][0][j][k]/cum[i][0][j][k], 2.0));
	    // cumR_e_sys[i][m][j][k] = cumR[i][m][j][k] * sqrt(pow(cum_e_sys[i][m][j][k]/cum[i][m][j][k], 2.0)+pow(cum_e_sys[i][0][j][k]/cum[i][0][j][k], 2.0));
	    
	    if(i==0&&m==NCum-1) { // 7.7 GeV C4 data points scaled down by 0.5
	      cumR[i][m][j][k] *= 0.5;
	      cumR_e[i][m][j][k] *= 0.5;
	      cumR_e_sys[i][m][j][k] *= 0.5;
	    }
	  } // end k->NY
	} // end m->NCum
      } // end j->NP
      fin[i]->Close();
	
      for(int j=0;j<NP;j++) {
	for(int m=0;m<NCum;m++) {
	  gr_stat[i][m][j] = new TGraphErrors(NY, rap[j], cumR[i][m][j], 0, cumR_e[i][m][j]);
	  gr_sys[i][m][j] = new TGraphErrors(NY, rap[j], cumR[i][m][j], 0, cumR_e_sys[i][m][j]);
	  cout << "++ " << EneDir[i] << " " << PName[j] << Form(" k%d/k1",m+2) << endl;
	  gr_stat[i][m][j]->Print();
	}
      } // end j->NP
    } else { // 54 GeV data
      //      TFile *fin_54 = new TFile("rootfile_0517/54GeV/54GeV_CORRELATION_MAY17.root");
      TFile *fin_54_stat[NY];
      TFile *fin_54_sys[NP][NY];
      
      for(int k=0;k<NY;k++) {
	fin_54_stat[k] = new TFile(Form("54GeV_data_Sep11/stat/stat.y0p%d.root",k+1));
	
	for(int j=0;j<NP;j++) {
	  fin_54_sys[j][k] = new TFile(Form("54GeV_data_Sep11/sys/%s/Sys_%s_y0p%d.root",PName_54[j],PName_His_54[j],k+1));
	  TGraphErrors *gr_stat_b[NCum], *gr_sys_b[NCum];
	  for(int m=0;m<NCum;m++) {
	    gr_stat_b[m] = (TGraphErrors *)fin_54_stat[k]->Get(Form("%s_k%d1",PName_His_54[j], m+2));
	    gr_sys_b[m] = (TGraphErrors *)fin_54_sys[j][k]->Get(Form("%s_k%d1_sys",PName_His_54[j], m+2));
	  
	    cumR[i][m][j][k] = gr_stat_b[m]->GetY()[0];
	    cumR_e[i][m][j][k] = gr_stat_b[m]->GetEY()[0];
	    cumR_e_sys[i][m][j][k] = gr_sys_b[m]->GetEY()[0];
	  } // end m->NCum
	  fin_54_sys[j][k]->Close();
	} // end j->NP
	fin_54_stat[k]->Close();  
      } // end k->NY

      for(int j=0;j<NP;j++) {
	for(int m=0;m<NCum;m++) { 	  
	  gr_stat[i][m][j] = new TGraphErrors(NY, rap[j], cumR[i][m][j], 0, cumR_e[i][m][j]);
	  gr_sys[i][m][j] = new TGraphErrors(NY, rap[j], cumR[i][m][j], 0, cumR_e_sys[i][m][j]);
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
  const double ymin[NCum] = {-0.11, -0.16, -1.1};
  const double ymax[NCum] = {0.031, 0.16, 1.5};
  const double xmin = 0.05;
  const double xmax = 0.55;
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
      h0[i][j]->GetXaxis()->SetNdivisions(404);
      h0[i][j]->GetXaxis()->SetTitleOffset(999.);
      h0[i][j]->GetXaxis()->SetLabelOffset(999.);
      h0[i][j]->GetYaxis()->SetNdivisions(505);
      h0[i][j]->GetYaxis()->SetTitleOffset(999.);
      h0[i][j]->GetYaxis()->SetLabelOffset(999.);
      h0[i][j]->Draw("c");
      drawLine(xmin,0,xmax,0,1,2,kBlack);

      for(int k=0;k<NP;k++) {
	cout << "++ " << EneDir[i] << " " << PName[k] << Form(" k%d/k1",j+2) << endl;
	gr_stat[i][j][k]->Print();

	drawSysError(gr_sys[i][j][k], 0.012, (ymax[j]-ymin[j])*0.02, kColor[k]);

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
	drawText(0.18, (ymax[j]-ymin[j])*0.88+ymin[j], Form("%s GeV", EneLabel[i]), 42, 0.14);
      }

      if(i==0&&j==NCum-1) { // 7.7 GeV C4 legend
	drawText(0.1, 1.15, "#kappa_{4}/#kappa_{1} #times 0.5", 42, 0.14);
	drawText(0.1, 0.8, "7.7 GeV only", 42, 0.14);
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
    drawText(i*1./NE + 0.024, 0.7, "0.2", 42, 0.35);
    drawText(i*1./NE + 0.069, 0.7, "0.4", 42, 0.35);
  }
  drawText(0.43, 0.25, "Rapidity Cut y_{max}", 42, 0.43);
  pad_x->Update();
  
  c1->cd();
  TPad *pad_y = new TPad("pad_y", "", 0, pad_y1, pad_x1-0.002, pad_y2);
  pad_y->Draw();
  pad_y->cd();

  const Int_t NYLabelMax = 5;
  const Int_t NYLabel[NCum] = {3, 3, 3};
  const Char_t *YLabel[NCum][NYLabelMax] = {{"#font[122]{-} 0.10","#font[122]{-} 0.05","   0.00","",""},
					    {"#font[122]{-} 0.1","   0.0","   0.1","",""},
					    {"#font[122]{-} 1","   0","   1","",""}
  };
  const Double_t XOffset[NCum] = {0.01, 0.12, 0.27};
  const Double_t Offset[NCum] = {0.018, 0.055, 0.005};
  const Double_t Step[NCum] = {0.115, 0.102, 0.125}; 
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
  TPad *pad_coll = new TPad("pad_coll","", 0.365, 0.58, 0.72, 0.67);
  pad_coll->Draw();
  pad_coll->cd();

  drawText(0.01, 0.6, "STAR Au+Au Collisions 0-5% most central", 42, 0.48);
  drawText(0.12, 0.18, "0.4 < p_{T} < 2.0 (GeV/c), |y| < y_{max}", 42, 0.45);
  pad_coll->Update();

  c1->cd();
  TPad *pad_leg = new TPad("pad_leg", "", 0.40, 0.28, 0.67, 0.34);
  pad_leg->Draw();
  pad_leg->cd();

  TLegend *leg[NP];
  const Int_t index_leg[NP] = {0, 1};
  const Char_t *name_leg[NP] = {" Proton", " Anti-proton"};
  const Double_t x1_leg[NP] = {0.1, 0.55};
  const Double_t x2_leg[NP] = {0.5, 0.9};
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
  c1->SaveAs("fig/Fig16_knk1_y.pdf");
  c1->SaveAs("fig/Fig16_knk1_y.png");
  
  
}
