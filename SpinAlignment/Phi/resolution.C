Double_t powerlaw(Double_t *x, Double_t *par)
{
  double n = par[2];
  double p0 = par[1]*(n-3)/2.;
  double A = 2.*(n-1)*(n-2)/TMath::Pi()/(n-3)/(n-3)/par[1]/par[1]*par[0];
  if(x[0]<5.0) return A*x[0]*pow(1+x[0]/p0,-n);
  else return A*100.*x[0]*pow(1+x[0]/p0,-n);
}

void resolution(const double rho3 = 1.0, const int sc = 1, const Int_t Nevt = 10000000)  // rho3 = rho00*3
{

  const Double_t MassPhi = 1.01946;
  const Double_t GammaPhi = 0.00425;  // Width
  const Double_t SigmaPhi = 0.00425/2.355;  // Width
  const Double_t MassPi = 0.13957;
  const Double_t MassK = 0.49368;
      

  TF1 *funMom = new TF1("funMom",powerlaw,0.,10.,3);
  funMom->SetParameters(1, 1.0, 11);
  funMom->SetNpx(1000);

  TF1 *funCosTheta = new TF1("funCosTheta","(1-[0])+(3*[0]-1)*x*x",-1,1);
  funCosTheta->SetParameter(0, rho3/3.);

  TF1 *funBW = new TF1("funBW","1./(pow((x*x-[0]*[0]),2.0)+[0]*[0]*[1]*[1])",0.95,1.15);
  funBW->SetNpx(1000);
  funBW->SetParameters(MassPhi, GammaPhi);
  
  TRandom3 *gRandom = new TRandom3();

  //-------------------
  // Book histograms
  //-------------------
  TH2D *hPtY = new TH2D("hPtY","",200,0.,2.,200,-2.,2.);
  TH2D *hPtEta = new TH2D("hPtEta","",200,0.,2.,200,-2.,2.);
  TH2D *hPtCosTheta = new TH2D("hPtCosTheta","",200,0.,2.,200,-1.,1.);
  TH2D *hMass = new TH2D("hMass","",200,0.,2.,200,0.98,1.08);
  TH2D *hCosTheta2Mass = new TH2D("hCosTheta2Mass","",200,0.98,1.08,200,0,1.);
  
  //|eta|<1 pT>0.2
  TH2D *hPtYRc = new TH2D("hPtYRc","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRc = new TH2D("hPtEtaRc","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRc = new TH2D("hPtCosThetaRc","",200,0.,2.,200,-1.,1.);
  TH2D *hCosThetaDiff = new TH2D("hCosThetaDiff","",200,-1.,1.,200,-0.01,0.01);
  TH2D *hMassRc = new TH2D("hMassRc","",200,0.,2.,200,0.98,1.08);
  TH2D *hCosTheta2MassRc = new TH2D("hCosTheta2MassRc","",200,0.98,1.08,200,0,1.);

  // after smearing
  TH2D *hPtYRcRes = new TH2D("hPtYRcRes","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtEtaRcRes = new TH2D("hPtEtaRcRes","",200,0.,2.0,200,-2.,2.);
  TH2D *hPtCosThetaRcRes = new TH2D("hPtCosThetaRcRes","",200,0.,2.,200,-1.,1.);
  TH2D *hCosThetaResDiff = new TH2D("hCosThetaResDiff","",200,-1.,1.,200,-0.2,0.2);
  TH2D *hMassRcRes = new TH2D("hMassRcRes","",200,0.,2.,200,0.98,1.08);
  TH2D *hCosTheta2MassRcRes = new TH2D("hCosTheta2MassRcRes","",200,0.98,1.08,200,0,1.);

  TH2D *hPtRes = new TH2D("hPtRes","",200,0.,2.0,200,-0.1,0.1);
  TH2D *hThetaRes = new TH2D("hThetaRes","",200,0.,2.0,200,-0.05,0.05);
  TH2D *hPhiRes = new TH2D("hPhiRes","",200,0.,2.0,200,-0.05,0.05);
  TH2D *hMassRes = new TH2D("hMassRes","",200,0.,2.0,200,-0.01,0.01);
  

  TLorentzVector Phimom(0.,0.,0.,0.);
  TLorentzVector K1mom(0.,0.,0.,0.);   // Kaon momentum in the lab frame
  TLorentzVector K2mom(0.,0.,0.,0.);   // Kaon momentum in the lab frame
  TLorentzVector PhimomRc(0.,0.,0.,0.);   // Phi momentum Rc in the lab frame
  TLorentzVector K1momStar(0.,0.,0.,0.);   // Kaon momentum Rc in the Phi rest frame

  TLorentzVector K1momRes(0.,0.,0.,0.);   // Kaon momentum in the lab frame after smearing
  TLorentzVector K2momRes(0.,0.,0.,0.);   // Kaon momentum in the lab frame after smearing
  TLorentzVector PhimomRcRes(0.,0.,0.,0.);   // Phi momentum Rc in the lab frame after smearing
  TLorentzVector K1momStarRes(0.,0.,0.,0.);   // Kaon momentum Rc in the Phi rest frame after smearing

  for(int ievt = 0; ievt<Nevt; ievt++) {
    if(ievt%100000==0) cout << " processing " << ievt << "-th event ..." << endl;
    //    double MassMom = gRandom->Gaus(MassPhi, SigmaPhi);
    double MassMom;
    do {
      MassMom = funBW->GetRandom();
    } while (MassMom<MassK*2.+1e-5);
    double pT = gRandom->Rndm()*2.; //funMom->GetRandom();
    double y = gRandom->Rndm()*2. - 1. ; // rapidity
    /*
    double y;
    do {
      y = gRandom->Gaus(0, sigY);
    } while (fabs(y)>1.0);
    */
    double phi = gRandom->Rndm()*TMath::Pi()*2.;
    double mT = TMath::Sqrt(pT*pT+MassMom*MassMom);
    double pz = mT * TMath::SinH(y);
    double p = TMath::Sqrt(pT*pT+pz*pz);
    double eta = 0.5*TMath::Log((p+pz)/(p-pz));
    Phimom.SetPtEtaPhiM(pT, eta, phi, MassMom);

    hPtY->Fill(pT, y);
    hPtEta->Fill(pT, eta);
    hMass->Fill(pT, MassMom);

    double p_pri = TMath::Sqrt((MassMom+2*MassK)*(MassMom-2*MassK))/2.;
    double phi_pri = gRandom->Rndm()*TMath::Pi()*2.0;   // phi x-z
    //    double costheta_pri = gRandom->Rndm()*2. - 1.;
    double costheta_pri = funCosTheta->GetRandom();
    double py_pri = p_pri*costheta_pri;    // theta - angle between Y-axis (pol axis) and momentum
    double pxz_pri = TMath::Sqrt(p_pri*p_pri-py_pri*py_pri);
    double px_pri = pxz_pri*TMath::Cos(phi_pri);
    double pz_pri = pxz_pri*TMath::Sin(phi_pri);
    K1mom.SetXYZM(px_pri, py_pri, pz_pri, MassK);
    K2mom.SetXYZM(-px_pri, -py_pri, -pz_pri, MassK);
    hPtCosTheta->Fill(pT, costheta_pri);
    hCosTheta2Mass->Fill(MassMom, costheta_pri*costheta_pri);

    // transform to the lab frame
    K1mom.Boost(Phimom.BoostVector());
    K2mom.Boost(Phimom.BoostVector());

    PhimomRc = K1mom + K2mom;
    K1momStar = K1mom;
    K1momStar.Boost(-PhimomRc.BoostVector());
    double costheta_rc = K1momStar.Y()/K1momStar.Vect().Mag();
    hCosThetaDiff->Fill(costheta_pri, costheta_rc-costheta_pri);

    // Check the momentum
    if(0) {
      cout << " Phi momentum = " << Phimom.Px() << " " << Phimom.Py() << " " << Phimom.Pz() << " " << Phimom.E() << endl;
      cout << " K1 momentum = " << K1mom.Px() << " " << K1mom.Py() << " " << K1mom.Pz() << " " << K1mom.E() << endl;
      cout << " K2 momentum = " << K2mom.Px() << " " << K2mom.Py() << " " << K2mom.Pz() << " " << K2mom.E() << endl;
    }

    if(fabs(K1mom.Eta())<1. && fabs(K2mom.Eta())<1. && K1mom.Pt()>0.2 && K2mom.Pt()>0.2) {
      hPtYRc->Fill(pT, y);
      hPtEtaRc->Fill(pT, eta);
      hPtCosThetaRc->Fill(pT, costheta_rc);
      hMassRc->Fill(pT, PhimomRc.M());
      hCosTheta2MassRc->Fill(PhimomRc.M(), costheta_rc*costheta_rc);
    }

    // resolution smearing
    double pt1res = gRandom->Gaus(K1mom.Pt(), 0.02*sc); // 0.01
    double theta1res = gRandom->Gaus(K1mom.Theta(), 0.005*sc); //0.002);
    double phi1res = gRandom->Gaus(K1mom.Phi(), 0.01*sc); //0.005);
    TVector3 K1momVectRes(0.,0.,0);
    K1momVectRes.SetPtThetaPhi(pt1res, theta1res, phi1res);
    K1momRes.SetVectM(K1momVectRes, K1mom.M());
    
    double pt2res = gRandom->Gaus(K2mom.Pt(), 0.02*sc);  // 0.01
    double theta2res = gRandom->Gaus(K2mom.Theta(), 0.005*sc); //0.002);
    double phi2res = gRandom->Gaus(K2mom.Phi(), 0.01*sc); //0.005);
    TVector3 K2momVectRes(0.,0.,0);
    K2momVectRes.SetPtThetaPhi(pt2res, theta2res, phi2res);
    K2momRes.SetVectM(K2momVectRes, K2mom.M());

    PhimomRcRes = K1momRes + K2momRes;
    K1momStarRes = K1momRes;
    K1momStarRes.Boost(-PhimomRcRes.BoostVector());
    double costheta_rc_res = K1momStarRes.Y()/K1momStarRes.Vect().Mag();
    hCosThetaResDiff->Fill(costheta_pri, costheta_rc_res-costheta_pri);

    if(fabs(K1momRes.Eta())<1. && fabs(K2momRes.Eta())<1. && K1momRes.Pt()>0.2 && K2momRes.Pt()>0.2) {
      hPtYRcRes->Fill(PhimomRcRes.Pt(), PhimomRcRes.Rapidity());
      hPtEtaRcRes->Fill(PhimomRcRes.Pt(), PhimomRcRes.Eta());
      hPtCosThetaRcRes->Fill(PhimomRcRes.Pt(), costheta_rc_res);
      hMassRcRes->Fill(PhimomRcRes.Pt(), PhimomRcRes.M());
      hCosTheta2MassRcRes->Fill(PhimomRcRes.M(), costheta_rc_res*costheta_rc_res);

      hPtRes->Fill(pT, PhimomRcRes.Pt()-pT);
      hThetaRes->Fill(pT, PhimomRcRes.Theta()-Phimom.Theta());
      hPhiRes->Fill(pT, PhimomRcRes.Phi()-Phimom.Phi());
      hMassRes->Fill(pT, PhimomRcRes.M()-Phimom.M());
    }
    
  }

  TFile *fout = new TFile(Form("resolution_%3.1f_%d.root",rho3,sc),"recreate");
  hPtY->Write();
  hPtEta->Write();
  hPtCosTheta->Write();
  hPtYRc->Write();
  hPtEtaRc->Write();
  hPtCosThetaRc->Write();
  hCosThetaDiff->Write();

  hPtYRcRes->Write();
  hPtEtaRcRes->Write();
  hPtCosThetaRcRes->Write();
  hCosThetaResDiff->Write();

  hMass->Write();
  hMassRc->Write();
  hMassRcRes->Write();
  
  hPtRes->Write();
  hThetaRes->Write();
  hPhiRes->Write();
  hMassRes->Write();

  hCosTheta2Mass->Write();
  hCosTheta2MassRc->Write();
  hCosTheta2MassRcRes->Write();


  fout->Close();

}
