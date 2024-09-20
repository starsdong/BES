#include "draw.C+"
#include "style.C+"

void makeComp(const Int_t iE = 0)
{
  style();

  const Int_t NPt = 5+1;  // 5 differential, 1 pT-ave
  const Double_t Pt[NPt] = {1.5, 2.1, 2.7, 3.6, 4.8, 0.6};  // 0.6 bin - pT-ave
  const Double_t PtEdge[NPt] = {1.2, 1.8, 2.4, 3.0, 4.2, 5.4};  //

  const Int_t NU = 5; // 3 results + 2 in pub
  const Int_t NE = 2; // 2 energy points
  const Char_t *NameUser[NE][NU] = {{"Baoshan-Run14", "CW-Run14", "Gavin-Run14", "Pub-Run14", "Pub-Run11"},
				    {"Baoshan-Run18", "CW-Run18", "Gavin-Run18", "Pub-Run18", "Pub-Run11"}
  };
  const Char_t *NameEnergy[NE] = {"200GeV@Run14", "27GeV@Run18"};
  const Char_t *LabelEnergy[NE] = {"200", "27"};

  const Double_t PtOffset[NU] = {-0.1, -0.05, 0., 0.05, 0.10};
  Double_t pT[NU][NPt];
  for(int i=0;i<NU;i++) {
    for(int j=0;j<NPt;j++) {
      pT[i][j] = Pt[j] + PtOffset[i];
    }
  }
  const Int_t NP[NE] = {5, 4};
  const Double_t rho[NE][NU][NPt] = {{{0.341654, 0.340816, 0.332985, 0.329835, 0.338988, 0.340398},  // Baoshan-200GeV
				      {0.339271, 0.338026, 0.337962, 0.333994, 0.331913, 0.33865},   // CW-200GeV
				      {0.342189, 0.340209, 0.334936, 0.333219, 0.339742, 0.340176}, // Gavin-200GeV
				      {0.336237, 0.333866, 0.329742, 0.327062, 0.331082, 0.3347}, // Pub-Run14
				      {0.340876, 0.341392, 0.334794, 0.331082, 0.332629, 0}}, // Pub-Run11
				     {{0.346184, 0.347417, 0.331666, 0.359454, 0.,       0.346015},  // Baoshan-27GeV
				      {0.34443,  0.34617,  0.33123,  0.34533,  0,        0.34430},   // CW-27GeV
				      {0.34536,  0.345596, 0.332992, 0.346875, 0,        0.344793},  // Gavin-27GeV
				      {0.3481,   0.3453,   0.3356,   0.3556,   0,        0.346832},  // Pub-Run18-27GeV
				      {0.342907, 0.344371, 0.360062, 0.332014, 0,        0.344039}}  // Pub-Run11-27GeV
  };
  
  const Double_t rho_e[NE][NU][NPt] = {{{0.000695697, 0.000910792, 0.00142713, 0.00211389, 0.00652381, 0.000513498}, // Baoshan-200GeV
					{0.000595797, 0.000793582, 0.00125364, 0.00189136, 0.00618794, 0.00044},     // CW-200GeV
					{0.000754189, 0.000979554, 0.00153717, 0.00227283, 0.00705258, 0.000539},   // Gavin-200GeV
					{0.000721649, 0.000927835, 0.0014433,  0.00206186, 0.00628866, 0.0005}, // Pub-Run14-200GeV
					{0.000721649, 0.00103093,  0.00154639, 0.00226804, 0.00701031, 0}}, // Pub-Run11-200GeV
				       {{0.00145959,  0.00253701,  0.00536096, 0.0119722,  0,          0.00123508},  // Baoshan-27GeV
					{0.00094,     0.00166,     0.00358,    0.00811,    0,          0.0008},      // CW-27GeV
					{0.00129392,  0.00224359,  0.00478058, 0.0108045,  0,          0.001086},    // Gavin-27GeV
					{0.0014,      0.0024,      0.005,      0.0113,     0,          0.00116909},  // Pub-Run18-27GeV
					{0.00264251,  0.00453337,  0.00951399, 0.0196433,  0,          0.00220591}}  // Pub-Run11-27GeV
  };

  TGraphErrors *gr[NE][NU];
  for(int i=0;i<NE;i++) {
    for(int j=0;j<NU;j++) {
      gr[i][j] = new TGraphErrors(NPt, pT[j], rho[i][j], 0, rho_e[i][j]);
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->Draw();

  const Double_t yy1[NE] = {0.324, 0.324};
  const Double_t yy2[NE] = {0.35, 0.36};
  double x1 = 0.0;
  double x2 = 6.0;
  double y1 = yy1[iE];
  double y2 = yy2[iE];
  TH1D *h0 = new TH1D("h0","",1,x1,x2);
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetTitle("Transverse Momentum p_{T} (GeV/c)");
  h0->GetYaxis()->SetNdivisions(105);
  h0->GetYaxis()->SetTitle("Uncorrected #phi-meson #rho_{00}");
  h0->Draw();

  drawLine(x1, 1./3, x2, 1./3, 1, 9, 1);
  //  drawLine(1.2, y1, 1.2, y2, 1, 2, 1);
  drawColorBox(x1, y1, 1.0, y2);
  drawLine(PtEdge[0], y1+0.04*(y2-y1), PtEdge[NPt-1], y1+0.04*(y2-y1), 1, 1, 2);
  for(int i=0;i<NPt;i++) {
    drawLine(PtEdge[i], y1+0.04*(y2-y1), PtEdge[i], y1+0.05*(y2-y1), 1, 1, 2);
  }
  
  const Int_t markerStyle[NU] = {21, 22, 23, 24, 20};
  const Int_t markerColor[NU] = {2, 6, 4, 1, 1};

  for(int j=0;j<NU;j++) {
    gr[iE][j]->SetMarkerSize(1.5);
    gr[iE][j]->SetMarkerStyle(markerStyle[j]);
    gr[iE][j]->SetMarkerColor(markerColor[j]);
    gr[iE][j]->SetLineColor(markerColor[j]);
    gr[iE][j]->SetLineWidth(2);
    gr[iE][j]->Draw("p");
  }

  drawText(Pt[NPt-1]-0.05, y1+0.04*(y2-y1), "p_{T}-average", 42, 0.03, 90);
  drawText(0.2, y2-0.1*(y2-y1), Form("Au+Au %s GeV, 20-60%%", LabelEnergy[iE]), 42, 0.05);

  drawHistBox(x1,x2,y1,y2);

  TLegend *leg = new TLegend(0.72, 0.68, 0.9, 0.94);
  leg->SetTextSize(0.035);
  leg->SetLineColor(10);
  for(int i=0;i<NU;i++) {
    leg->AddEntry(gr[iE][i], NameUser[iE][i], "pl");
  }
  leg->Draw();

  c1->Update();
  c1->SaveAs(Form("Comp_%sGeV.pdf",LabelEnergy[iE]));
  c1->SaveAs(Form("Comp_%sGeV.png",LabelEnergy[iE]));

}
