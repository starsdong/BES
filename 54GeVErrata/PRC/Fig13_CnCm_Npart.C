#include "draw.C+"
#include "style.C+"

void Fig13_CnCm_Npart()
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
  const Int_t NCum = 3; // 3 orders of cumulant ratios
  const Char_t *CumName[NCum] = {"VM", "SD", "KV"};
  const Char_t *CumName_54[NCum] = {"R21", "R32", "R42"};
  const Int_t NP = 3; // number of particle categories: proton, anti-proton
  const Char_t *PName[NP] = {"Pro", "Apro", "Netp"};
  //  const Char_t *PName_54[NP] = {"pro", "antipro", "netp"};
  const Char_t *PName_54[NP] = {"Pro", "Pbar", "Netp"};
  const Char_t *PName_His_54[NP] = {"Pro", "Pbar", "Netp"};
  const Int_t NCen = 9; // 9 centrality bins
  const double cen_offset[NP] = {-3, 0, 3}; // different rapidity offsets for different particles for plotting  
  const int pad_gap_index = 5; // starting from index 5 (39 GeV)
  
  
  ////////////////////////////////////////////////
  // Read in data: 54 GeV stored differently
  ////////////////////////////////////////////////
  TGraphErrors *gr_stat[NE][NCum][NP], *gr_sys[NE][NCum][NP];  // all central collisions
  double cen[NP][NCen];
  double cum[NE][NCum+1][NP][NCen], cum_e[NE][NCum+1][NP][NCen], cum_e_sys[NE][NCum+1][NP][NCen];   // name is cum, actually are kappa in this Fig.
  double cumR[NE][NCum][NP][NCen], cumR_e[NE][NCum][NP][NCen], cumR_e_sys[NE][NCum][NP][NCen]; // ratios to k1
  TFile *fin[NE][NP];
  for(int i=0;i<NE;i++) {
    if(i!=index_54) {
      for(int j=0;j<NP;j++) {
	fin[i][j] = new TFile(Form("rootfile_0517/final_cum_ycut/%sGeV/AuAu_sys_%s_Y0.5.root",EneDir[i],PName[j]));
	TGraphErrors *gr_stat_a[NCum], *gr_sys_a[NCum];
	for(int m=0;m<NCum;m++) {
	  gr_stat_a[m] = (TGraphErrors *)fin[i][j]->Get(Form("%s_stat",CumName[m]));
	  gr_sys_a[m] = (TGraphErrors *)fin[i][j]->Get(Form("%s_sys",CumName[m]));

	  for(int k=0;k<NCen;k++) {
	    cen[j][k] = gr_stat_a[m]->GetX()[k] + cen_offset[j]; 
	    cumR[i][m][j][k] = gr_stat_a[m]->GetY()[k];
	    cumR_e[i][m][j][k] = gr_stat_a[m]->GetEY()[k];
	    cumR_e_sys[i][m][j][k] = gr_sys_a[m]->GetEY()[k];
	    //	    cout << EneDir[i] << " " << PName[j] << " " << cen[k] << Form(" C%d",m+1) << " " << cum[i][m][j][k] << endl;
	    if(i>=pad_gap_index && m==0) {
	      cumR[i][m][j][k] *= 2.0;
	      cumR_e[i][m][j][k] *= 2.0;
	      cumR_e_sys[i][m][j][k] *= 2.0;  // C2/C1 scaled up by 2 for making Ndivisions show up well (i>=gap)
	    }
	    if(i<pad_gap_index && m==1) {
	      cumR[i][m][j][k] *= 2.0;
	      cumR_e[i][m][j][k] *= 2.0;
	      cumR_e_sys[i][m][j][k] *= 2.0;  // C3/C2 scaled up by 2 for making Ndivisions show up well (i<gap)
	    }

	  } // end k->NCen
	} // end m->NCum
	fin[i][j]->Close();
	
	for(int m=0;m<NCum;m++) {
	  gr_stat[i][m][j] = new TGraphErrors(NCen, cen[j], cumR[i][m][j], 0, cumR_e[i][m][j]);
	  gr_sys[i][m][j] = new TGraphErrors(NCen, cen[j], cumR[i][m][j], 0, cumR_e_sys[i][m][j]);
	  cout << "++ " << EneDir[i] << " " << PName[j] << Form(" %s", CumName[m]) << endl;
	  gr_stat[i][m][j]->Print();
	}
      } // end j->NP
    } else { // 54 GeV data
      //      TFile *fin_54 = new TFile("rootfile_0517/54GeV/54GeV_CUMULANTS_MAY17.root");
      //      TFile *fin_54_stat = new TFile("54GeV_data_Sep11/stat/stat.y0p5.root");
      TFile *fin_54_stat = new TFile("54GeV_data_Sep18/eff_corrected/9bins/stat/stat.y0p5.root");
      TFile *fin_54_sys[NP];
      
      for(int j=0;j<NP;j++) {
	//	fin_54_sys[j] = new TFile(Form("54GeV_data_Sep11/sys/%s/Sys_%s_y0p5.root",PName_54[j],PName_His_54[j]));
	fin_54_sys[j] = new TFile(Form("54GeV_data_Sep18/eff_corrected/9bins/sys/%s/Sys_%s_y0p5.root",PName_54[j],PName_His_54[j]));
	TGraphErrors *gr_stat_b[NCum+1], *gr_sys_b[NCum+1];
	for(int m=0;m<NCum;m++) {
	  gr_stat_b[m] = (TGraphErrors *)fin_54_stat->Get(Form("%s_%s",PName_His_54[j], CumName_54[m]));
	  gr_sys_b[m] = (TGraphErrors *)fin_54_sys[j]->Get(Form("%s_%s_sys",PName_His_54[j], CumName_54[m]));
	  
	  for(int k=0;k<NCen;k++) {
	    cen[j][k] = gr_stat_b[m]->GetX()[k] + cen_offset[j];
	    cumR[i][m][j][k] = gr_stat_b[m]->GetY()[k];
	    cumR_e[i][m][j][k] = gr_stat_b[m]->GetEY()[k];
	    cumR_e_sys[i][m][j][k] = gr_sys_b[m]->GetEY()[k];

	    if(i>=pad_gap_index && m==0) {
	      cumR[i][m][j][k] *= 2.0;
	      cumR_e[i][m][j][k] *= 2.0;
	      cumR_e_sys[i][m][j][k] *= 2.0;  // C2/C1 scaled up by 2 for making Ndivisions show up well
	    }
	    if(i<pad_gap_index && m==1) {
	      cumR[i][m][j][k] *= 2.0;
	      cumR_e[i][m][j][k] *= 2.0;
	      cumR_e_sys[i][m][j][k] *= 2.0;  // C3/C2 scaled up by 2 for making Ndivisions show up well (i<gap)
	    }

	  } // end k->NCen
	  
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
  const double pad_x_gap = 0.03; // gap between 27 GeV and 39 GeV (index 4-5)
  const double pad_x2 = 0.99;
  const double pad_y2 = 0.98;
  const double ymin[NCum] = {0.78, 0.41*2, -0.6};
  const double ymax[NCum] = {1.88, 1.09*2, 3.4};  // C3/C2 for 39- GeV scaled up by x2 for ndiv to show up well
  const double ymin2[NCum] = {0.2*2, -0.03, 0.3};
  const double ymax2[NCum] = {9.8*2, 1.13, 1.32};  // C2/C1 for 39+ GeV scaled up by x2 for ndiv to show up well
  const double xmin = -20;
  const double xmax = 380;
  const Int_t kColor[NP] = {kBlue, kBlack, kRed};
  const Int_t kStyle[NP] = {21, 34, 20};
  const Double_t kSize[NP] = {0.6, 1.2, 0.7};
  const Int_t NDiv[NCum] = {504, 504, 404};
  const Int_t NDiv2[NCum] = {504, 505, 504};
  
  ////////////////////////////////////////////////
  // Draw main panels
  ////////////////////////////////////////////////
  for(int i=0;i<NE;i++) {
    for(int j=0;j<NCum;j++) {
      c1->cd();
      
      double x1 = pad_x1 + i*(pad_x2 - pad_x1 - pad_x_gap)/NE;
      double x2 = pad_x1 + (i+1)*(pad_x2 - pad_x1 - pad_x_gap)/NE;
      if(i>=pad_gap_index) {
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
      if(i<pad_gap_index) {
	h0[i][j]->SetMaximum(ymax[j]);
	h0[i][j]->SetMinimum(ymin[j]);
	h0[i][j]->GetYaxis()->SetNdivisions(NDiv[j]);
      } else {
	h0[i][j]->SetMaximum(ymax2[j]);
	h0[i][j]->SetMinimum(ymin2[j]);
	h0[i][j]->GetYaxis()->SetNdivisions(NDiv2[j]);
      }
      h0[i][j]->GetXaxis()->SetNdivisions(404);
      h0[i][j]->GetXaxis()->SetTitleOffset(999.);
      h0[i][j]->GetXaxis()->SetLabelOffset(999.);
      h0[i][j]->GetYaxis()->SetTitleOffset(999.);
      h0[i][j]->GetYaxis()->SetLabelOffset(999.);
      h0[i][j]->Draw("c");
      if(i>=pad_gap_index&&j==0 || i<pad_gap_index&&j==1) {
	drawLine(xmin,2,xmax,2,1,2,kBlack);
      } else {
	drawLine(xmin,1,xmax,1,1,2,kBlack);
      }

      for(int k=0;k<NP;k++) {
	cout << "++ " << EneDir[i] << " " << PName[k] << CumName[j] << endl;
	gr_stat[i][j][k]->Print();

	if(i<pad_gap_index) {
	  drawSysError(gr_sys[i][j][k], 10, (ymax[j]-ymin[j])*0.02, kColor[k]);
	} else {
	  drawSysError(gr_sys[i][j][k], 10, (ymax2[j]-ymin2[j])*0.02, kColor[k]);
	}

	setGraphMarker(gr_stat[i][j][k], kStyle[k], kColor[k], kSize[k]);
	setGraphLine(gr_stat[i][j][k], 1, kColor[k], 1);
	gr_stat[i][j][k]->Draw("p");

      }

      TGraphErrors *gr_clone = (TGraphErrors *)gr_stat[i][j][NP-1]->Clone(Form("%s_clone",gr_stat[i][j][NP-1]->GetName()));
      // net-p replotting
      setGraphMarker(gr_clone, 24, kBlack, kSize[NP-1]+0.15);
      setGraphLine(gr_clone, 1, kRed, 1);
      gr_clone->Draw("p");

      if(j==0) {
	if(i<pad_gap_index) {
	  drawText(100, (ymax[j]-ymin[j])*0.88+ymin[j], Form("%s GeV", EneLabel[i]), 42, 0.14);
	} else {
	  drawText(100, (ymax2[j]-ymin2[j])*0.88+ymin2[j], Form("%s GeV", EneLabel[i]), 42, 0.14);
	}
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
    double xx_1 = i*1.*(pad_x2-pad_x1-pad_x_gap)/(pad_x2-pad_x1)/NE + 0.002;
    double xx_2 = i*1.*(pad_x2-pad_x1-pad_x_gap)/(pad_x2-pad_x1)/NE + 0.021;
    double xx_3 = i*1.*(pad_x2-pad_x1-pad_x_gap)/(pad_x2-pad_x1)/NE + 0.048;
    double xx_4 = i*1.*(pad_x2-pad_x1-pad_x_gap)/(pad_x2-pad_x1)/NE + 0.075;
    if(i>=pad_gap_index) {
      xx_1 += pad_x_gap*(pad_x2-pad_x1)/(pad_x2-pad_x1-pad_x_gap)+0.001;
      xx_2 += pad_x_gap*(pad_x2-pad_x1)/(pad_x2-pad_x1-pad_x_gap)+0.001;
      xx_3 += pad_x_gap*(pad_x2-pad_x1)/(pad_x2-pad_x1-pad_x_gap)+0.001;
      xx_4 += pad_x_gap*(pad_x2-pad_x1)/(pad_x2-pad_x1-pad_x_gap)+0.001;
    }
    
    drawText(xx_1, 0.75, "0", 42, 0.35);
    drawText(xx_2, 0.75, "100", 42, 0.35);
    drawText(xx_3, 0.75, "200", 42, 0.35);
    drawText(xx_4, 0.75, "300", 42, 0.35);
  }
  drawText(0.33, 0.25, "Average Number of Participant Nucleons <N_{part}>", 42, 0.43);
  pad_x->Update();
  
  c1->cd();
  TPad *pad_y = new TPad("pad_y", "", 0, pad_y1, pad_x1-0.002, 1);
  pad_y->Draw();
  pad_y->cd();

  const Int_t NYLabelMax = 5;
  const Int_t NYLabel[NCum] = {2, 3, 2};
  const Char_t *YLabel[NCum][NYLabelMax] = {{"1.0","1.5","","",""},
					    {"0.50","0.75","1.00","",""},
					    {"0","2","","",""}
  };
  const Double_t XOffset[NCum] = {0.12, -0.01, 0.30};
  const Double_t Offset[NCum] = {0.04, 0.0255, 0.036};
  const Double_t Step[NCum] = {0.146, 0.1185, 0.163};
  const Char_t *YTitle[NCum] = {"C_{2}/C_{1}","C_{3}/C_{2}","C_{4}/C_{2}"}; 
  for(int i=0;i<NCum;i++) {
    for(int j=0;j<NYLabel[i];j++) {
      drawText(0.55+XOffset[i], (NCum-1-i)*1./NCum + j*Step[i] + Offset[i], YLabel[i][j], 42, 0.26);
    }
    drawText(0.4, (NCum-1-i+0.4)*1./NCum, YTitle[i], 42, 0.3, 90);
  }
  pad_y->Update();

  c1->cd();
  double x1_gap = pad_x1 + pad_gap_index*(pad_x2 - pad_x1 - pad_x_gap)/NE;
  double x2_gap = x1_gap + pad_x_gap;
  TPad *pad_yy = new TPad("pad_yy", "", x1_gap+0.002, pad_y1, x2_gap-0.002, 1);
  pad_yy->Draw();
  pad_yy->cd();

  const Int_t NYLabel2[NCum] = {3, 3, 2};
  const Char_t *YLabel2[NCum][NYLabelMax] = {{"2.5","5.0","7.5","",""},
					     {"0.0","0.5","1.0","",""},
					     {"0.5","1.0","","",""}
  };
  const Double_t XOffset2[NCum] = {0.12, 0.12, 0.12};
  const Double_t Offset2[NCum] = {0.052, -0.01, 0.053};
  const Double_t Step2[NCum] = {0.086, 0.141, 0.161};
  for(int i=0;i<NCum;i++) {
    for(int j=0;j<NYLabel2[i];j++) {
      drawText(0.15+XOffset2[i], (NCum-1-i)*1./NCum + j*Step2[i] + Offset2[i], YLabel2[i][j], 42, 0.58);
    }
  }
  pad_yy->Update();
  

  ////////////////////////////////////////////////
  // Draw other legends and labels
  ////////////////////////////////////////////////
  c1->cd();
  TPad *pad_coll = new TPad("pad_coll","", 0.07, 0.82, 0.33, 0.93);
  pad_coll->Draw();
  pad_coll->cd();

  drawText(0.08, 0.6, "STAR Au+Au Collisions", 42, 0.355);
  drawText(0.01, 0.18, "0.4 < p_{T} < 2.0 (GeV/c), |y| < 0.5", 42, 0.33);
  pad_coll->Update();

  c1->cd();
  TPad *pad_leg = new TPad("pad_leg", "", 0.60, 0.85, 0.90, 0.90);
  pad_leg->Draw();
  pad_leg->cd();

  TLegend *leg[NP];
  const Int_t index_leg[NP] = {2, 0, 1};
  const Char_t *name_leg[NP] = {" Net-Proton", " Proton", " Anti-proton"};
  const Double_t x1_leg[NP] = {0.05, 0.42, 0.65};
  const Double_t x2_leg[NP] = {0.4, 0.63, 0.95};
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

  double x_re[1] = {0.094};
  double y_re[1] = {0.5};
  TGraphErrors *gr_re = new TGraphErrors(1, x_re, y_re);
  gr_re->SetMarkerStyle(24);
  gr_re->SetMarkerSize(1.2);
  gr_re->Draw("p");

  
  pad_leg->Update();
  c1->cd();

  c1->Update();
  c1->SaveAs("fig/Fig13_CnCm_Npart.pdf");
  c1->SaveAs("fig/Fig13_CnCm_Npart.png");
  
  
}
