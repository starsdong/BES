#include "draw.C+"
#include "style.C+"

void Fig15_Cn_y()
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
  const Int_t NCum = 4; // 4 orders of cumulants
  const Int_t NP = 3; // number of particle categories: Net-p, proton, anti-proton
  const Char_t *PName[NP] = {"Pro", "Apro", "Netp"};
  const Char_t *PName_54[NP] = {"pro", "antipro", "netp"};
  const Int_t NY = 5;  // rapidity bins
  const Double_t RAP[NY] = {0.1, 0.2, 0.3, 0.4, 0.5};
  const Int_t NCen = 9; // 9 centrality bins
  const double rap_offset[NP] = {0., 0.008, 0.016}; // different rapidity offsets for different particles for plotting
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
  double cum[NE][NCum][NP][NY], cum_e[NE][NCum][NP][NY], cum_e_sys[NE][NCum][NP][NY];
  TFile *fin[NE][NP][NY];
  for(int i=0;i<NE;i++) {
    if(i!=index_54) {
      for(int j=0;j<NP;j++) {
	for(int k=0;k<NY;k++) {
	  fin[i][j][k] = new TFile(Form("rootfile_0517/final_cum_ycut/%sGeV/AuAu_sys_%s_Y%3.1f.root",EneDir[i],PName[j],RAP[k]));
	  TGraphErrors *gr_stat_tmp[NCum], *gr_sys_tmp[NCum];
	  for(int m=0;m<NCum;m++) {
	    gr_stat_tmp[m] = (TGraphErrors *)fin[i][j][k]->Get(Form("C%d_stat",m+1));
	    gr_sys_tmp[m] = (TGraphErrors *)fin[i][j][k]->Get(Form("C%d_sys",m+1));
	    cum[i][m][j][k] = gr_stat_tmp[m]->GetY()[0];
	    cum_e[i][m][j][k] = gr_stat_tmp[m]->GetEY()[0];
	    cum_e_sys[i][m][j][k] = gr_sys_tmp[m]->GetEY()[0];
	    //	    cout << EneDir[i] << " " << PName[j] << " " << rap[k] << Form(" C%d",m+1) << " " << cum[i][m][j][k] << endl;
	    if(i==0&&m==NCum-1) { // 7.7 GeV C4 data points scaled down by 0.5
	      cum[i][m][j][k] *= 0.5;
	      cum_e[i][m][j][k] *= 0.5;
	      cum_e_sys[i][m][j][k] *= 0.5;
	    }
	  } // end m->NCum
	  fin[i][j][k]->Close();
	} // end k->NY

	for(int m=0;m<NCum;m++) {
	  gr_stat[i][m][j] = new TGraphErrors(NY, rap[j], cum[i][m][j], 0, cum_e[i][m][j]);
	  gr_sys[i][m][j] = new TGraphErrors(NY, rap[j], cum[i][m][j], 0, cum_e_sys[i][m][j]);
	  //	  cout << "++ " << EneDir[i] << " " << PName[j] << Form(" C%d",m+1) << endl;
	  //	  gr_stat[i][m][j]->Print();
	}
      } // end j->NP
    } else { // 54 GeV data
      TFile *fin_54 = new TFile("rootfile_0517/54GeV/54GeV_CUMULANTS_MAY17.root");
      for(int j=0;j<NP;j++) {
	TGraphErrors *gr_stat_tmp[NCum], *gr_sys_tmp[NCum];
	for(int m=0;m<NCum;m++) {
	  gr_stat_tmp[m] = (TGraphErrors *)fin_54->Get(Form("rapidity_%s_C%d_stat",PName_54[j], m+1));
	  gr_sys_tmp[m] = (TGraphErrors *)fin_54->Get(Form("rapidity_%s_C%d_sys",PName_54[j], m+1));

	  gr_stat[i][m][j] = new TGraphErrors(gr_stat_tmp[m]->GetN(), rap[j], gr_stat_tmp[m]->GetY(), 0, gr_stat_tmp[m]->GetEY());
	  gr_sys[i][m][j] = new TGraphErrors(gr_sys_tmp[m]->GetN(), rap[j], gr_sys_tmp[m]->GetY(), 0, gr_sys_tmp[m]->GetEY());
	} // end m->NCum
	
      } // end j->NP    
    } // end if(i!=index_54)
  } // end i->NE

  ////////////////////////////////////////////////
  // Start Plotting
  ////////////////////////////////////////////////
  TCanvas *c1 = new TCanvas("c1","",1400,710);
  TPad *pad[NE][NCum];
  TH1D *h0[NE][NCum];
  const double pad_x1 = 0.06;
  const double pad_y1 = 0.09; // x,y pad offsets
  const double pad_x2 = 0.99;
  const double pad_y2 = 0.98;
  const double ymin[NE] = {-4.9, -4.9, -4, -16};
  const double ymax[NE] = {47, 44, 38, 63};
  const double xmin = 0.05;
  const double xmax = 0.55;
  const Int_t kColor[NP] = {kBlue, kBlack, kRed};
  const Int_t kStyle[NP] = {21, 34, 20};
  const Double_t kSize[NP] = {0.6, 1.1, 0.7};

  ////////////////////////////////////////////////
  // Draw main panels
  ////////////////////////////////////////////////
  for(int i=0;i<NE;i++) {
    for(int j=0;j<NCum;j++) {
      c1->cd();
      
      double x1 = pad_x1 + i*(pad_x2 - pad_x1)/NE;
      double x2 = pad_x1 + (i+1)*(pad_x2 - pad_x1)/NE;
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
      h0[i][j]->GetYaxis()->SetNdivisions(404);
      h0[i][j]->GetYaxis()->SetTitleOffset(999.);
      h0[i][j]->GetYaxis()->SetLabelOffset(999.);
      h0[i][j]->Draw("c");
      drawLine(xmin,0,xmax,0,1,2,kBlack);

      for(int k=0;k<NP;k++) {
	cout << "++ " << EneDir[i] << " " << PName[k] << Form(" C%d",j+1) << endl;
	gr_stat[i][j][k]->Print();

	drawSysError(gr_sys[i][j][k], 0.01, 0.7, kColor[k]);

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
	drawText(0.18, ymax[j]*0.84, Form("%s GeV", EneLabel[i]), 42, 0.14);
      }

      if(i==0&&j==NCum-1) { // 7.7 GeV C4 legend
	drawText(0.1, 48, "C_{4} #times 0.5", 42, 0.13);
	drawText(0.1, 38, "7.7 GeV only", 42, 0.13);
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
    drawText(i*1./NE + 0.024, 0.8, "0.2", 42, 0.25);
    drawText(i*1./NE + 0.069, 0.8, "0.4", 42, 0.25);
  }
  drawText(0.43, 0.25, "Rapidity Cut y_{max}", 42, 0.43);
  pad_x->Update();
  
  c1->cd();
  TPad *pad_y = new TPad("pad_y", "", 0, pad_y1, pad_x1-0.002, pad_y2);
  pad_y->Draw();
  pad_y->cd();

  const Int_t NYLabel[NCum] = {3, 3, 2, 4};
  const Char_t *YLabel[4] = {"  0","20","40","60"};
  const Double_t Offset[NCum] = {0.016, 0.016, 0.016, 0.038};
  const Double_t Step[NCum] = {0.095, 0.1, 0.12, 0.065}; 
  for(int i=0;i<NCum;i++) {
    for(int j=0;j<NYLabel[i];j++) {
      drawText(0.78, (3-i)*1./NCum + j*Step[i] + Offset[i], YLabel[j], 42, 0.2);
    }
    drawText(0.4, (3-i+0.4)*1./NCum, Form("C_{%d}",i+1), 42, 0.3, 90);
  }
  pad_y->Update();

  ////////////////////////////////////////////////
  // Draw other legends and labels
  ////////////////////////////////////////////////
  c1->cd();
  TPad *pad_coll = new TPad("pad_coll","", 0.43, 0.205, 0.66, 0.31);
  pad_coll->Draw();
  pad_coll->cd();

  drawText(0.05, 0.7, "STAR Au+Au Collisions", 42, 0.33);
  drawText(0.16, 0.4, "0-5% most central", 42, 0.30);
  drawText(0.05, 0.1, "0.4 < p_{T} < 2.0 (GeV/c), |y| < y_{max}", 42, 0.25);
  pad_coll->Update();

  c1->cd();
  TPad *pad_leg = new TPad("pad_leg", "", 0.38, 0.47, 0.67, 0.53);
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
    leg[i]->SetTextSize(0.5);
    gr_leg[i] = (TGraphErrors *)gr_stat[0][0][index_leg[i]]->Clone();
    gr_leg[i]->SetMarkerSize(1.3);
    leg[i]->AddEntry(gr_leg[i], name_leg[i], "p");
    leg[i]->Draw();
  }

  double x_re[1] = {0.094};
  double y_re[1] = {0.485};
  TGraphErrors *gr_re = new TGraphErrors(1, x_re, y_re);
  gr_re->SetMarkerStyle(24);
  gr_re->SetMarkerSize(1.2);
  gr_re->Draw("p");
  
  pad_leg->Update();
  c1->cd();

  c1->Update();
  c1->SaveAs("fig/Fig15_Cn_y.pdf");
  c1->SaveAs("fig/Fig15_Cn_y.png");
  
  
}
