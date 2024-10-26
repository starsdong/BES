#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TLorentzVector.h"


Double_t mTExp(Double_t *x, Double_t *par)
{
  const Double_t MassPhi = 1.01946;
  return x[0]*TMath::Exp(-(x[0]-MassPhi)/par[0]);
}

Double_t eff(Double_t *x, Double_t *par)
{
  Double_t y = (TMath::Abs(x[0]-par[0])-par[2])/par[3];
  return par[1]/(1.+TMath::Exp(y))+par[4];
}


void accept_Sergei_v2pt_Minv(const Int_t Nevt = 100000000)
{

  const Double_t sigY = 9.9;
  
  const Double_t MassPhi = 1.01946;
  const Double_t SigmaPhi = 0.00182; // 1.82 MeV ~ Gamma = 4.2 MeV
  const Double_t MassPi = 0.13957;
  const Double_t MassK = 0.49368;

  // TF1 *fEff = new TF1("fun","[0]*(1-exp(-(x-[1])/[2]))",0,10);
  // fEff->SetParameters(0.659, 0.0686, 1.657);
  TGraph *gr_eff_tpc = new TGraph("eff_tpc_kaon.dat","%lg %lg");
  TGraph *gr_eff_tof = new TGraph("eff_tof_kaon.dat","%lg %lg");

  TF1 *fun_eff_phi = new TF1("fun_eff_phi","1+[0]*(1+cos(2*x))",-TMath::Pi()*4.,TMath::Pi()*4.);
  fun_eff_phi->SetParameter(0, -0.005); // x = phi(K) - Psi_2 / Fig. 26 in AN
  
  TF1 *fPhi = new TF1("fPhi","1. + 2.*[0]*cos(2.*x)", 0, TMath::Pi()*2.0);
  //  fPhi->SetParameter(0, 0.1);
  fPhi->SetParameter(0, 0.1); // preset, consider v2/nq = k * ET/nq

  TF1 *fmT = new TF1("fmT",mTExp,MassPhi,10.,3);
  fmT->SetParameter(0, 0.36);
  fmT->SetNpx(1000);

  TF1 *funCosTheta = new TF1("funCosTheta","(1-[0])+(3*[0]-1)*x*x",-1,1);
  funCosTheta->SetParameter(0, 1./3.);
  
  TRandom3 *gRandom = new TRandom3();

  //-------------------
  // Book histograms
  //-------------------
  TH2D *hPtY = new TH2D("hPtY","",200,0.,10.,200,-2.,2.);
  TH2D *hPtEta = new TH2D("hPtEta","",200,0.,10.,200,-2.,2.);
  TH2D *hPtV2 = new TH2D("hPtV2","",200,0.,10.,200,-1.,1.);  
  TH2D *hPtCosTheta = new TH2D("hPtCosTheta","",200,0.,10.,2000,-1.,1.);
  TH2D *hCosThetaDiff = new TH2D("hCosThetaDiff","",200,-1.,1.,200,-0.01,0.01);
  TH2D *hPtCos2Phi = new TH2D("hPtCos2Phi","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRP = new TH2D("hPtCos2PhiRP","",200,0.,10.,2000,-1.,1.); // delta_rho00 = -4/3 * <cos(2*(phi*-Psi_RP))> - alternative calculation
  TH2D *hPtCosTheta2 = new TH2D("hPtCosTheta2","",200,0.,10.,1000,0.,1.);
  TH2D *hMinvCos2Phi = new TH2D("hMinvCos2Phi","",200,0.98,1.08,2000,-1.,1.);
  
  // rapidity dependence, 1.2<pT<5.4
  TH2D *hYCosTheta = new TH2D("hYCosTheta","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2Phi = new TH2D("hYCos2Phi","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRP = new TH2D("hYCos2PhiRP","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCosTheta2 = new TH2D("hYCosTheta2","",20,0.,1.,1000,0.,1.);
  TH2D *hYV2 = new TH2D("hYV2","",20,0.,1.,200,-1.,1.);
  
  //|eta|<0.9
  TH2D *hPtYRc = new TH2D("hPtYRc","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtEtaRc = new TH2D("hPtEtaRc","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtV2Rc = new TH2D("hPtV2Rc","",200,0.,10.,200,-1.,1.);
  TH2D *hPtCosThetaRc = new TH2D("hPtCosThetaRc","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRc = new TH2D("hPtCos2PhiRc","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRPRc = new TH2D("hPtCos2PhiRPRc","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCosTheta2Rc = new TH2D("hPtCosTheta2Rc","",200,0.,10.,1000,0.,1.);
  TH2D *hMinvCos2PhiRc = new TH2D("hMinvCos2PhiRc","",200,0.98,1.08,2000,-1.,1.);

  TH2D *hYCosThetaRc = new TH2D("hYCosThetaRc","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRc = new TH2D("hYCos2PhiRc","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRPRc = new TH2D("hYCos2PhiRPRc","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCosTheta2Rc = new TH2D("hYCosTheta2Rc","",20,0.,1.,1000,0.,1.);
  TH2D *hYV2Rc = new TH2D("hYV2Rc","",20,0.,1.,200,-1.,1.);

  //|eta|<0.9 pT>0.2
  TH2D *hPtYRc1 = new TH2D("hPtYRc1","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtEtaRc1 = new TH2D("hPtEtaRc1","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtV2Rc1 = new TH2D("hPtV2Rc1","",200,0.,10.,200,-1.,1.);
  TH2D *hPtCosThetaRc1 = new TH2D("hPtCosThetaRc1","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRc1 = new TH2D("hPtCos2PhiRc1","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRPRc1 = new TH2D("hPtCos2PhiRPRc1","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCosTheta2Rc1 = new TH2D("hPtCosTheta2Rc1","",200,0.,10.,1000,0.,1.);
  TH2D *hMinvCos2PhiRc1 = new TH2D("hMinvCos2PhiRc1","",200,0.98,1.08,2000,-1.,1.);

  TH2D *hYCosThetaRc1 = new TH2D("hYCosThetaRc1","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRc1 = new TH2D("hYCos2PhiRc1","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRPRc1 = new TH2D("hYCos2PhiRPRc1","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCosTheta2Rc1 = new TH2D("hYCosTheta2Rc1","",20,0.,1.,1000,0.,1.);
  TH2D *hYV2Rc1 = new TH2D("hYV2Rc1","",20,0.,1.,200,-1.,1.);

  TH2D *hPtPhiK1 = new TH2D("hPtPhiK1","",200,0.,10.0,200,0.,10.0);
  TH2D *hPtPhiK2 = new TH2D("hPtPhiK2","",200,0.,10.0,200,0.,10.0);

  //|eta|<0.9 pT>0.2, - TPC efficiency
  TH2D *hPtYRc2 = new TH2D("hPtYRc2","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtEtaRc2 = new TH2D("hPtEtaRc2","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtV2Rc2 = new TH2D("hPtV2Rc2","",200,0.,10.,200,-1.,1.);
  TH2D *hPtCosThetaRc2 = new TH2D("hPtCosThetaRc2","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRc2 = new TH2D("hPtCos2PhiRc2","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRPRc2 = new TH2D("hPtCos2PhiRPRc2","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCosTheta2Rc2 = new TH2D("hPtCosTheta2Rc2","",200,0.,10.,1000,0.,1.);
  TH2D *hMinvCos2PhiRc2 = new TH2D("hMinvCos2PhiRc2","",200,0.98,1.08,2000,-1.,1.);
  
  TH2D *hYCosThetaRc2 = new TH2D("hYCosThetaRc2","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRc2 = new TH2D("hYCos2PhiRc2","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRPRc2 = new TH2D("hYCos2PhiRPRc2","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCosTheta2Rc2 = new TH2D("hYCosTheta2Rc2","",20,0.,1.,1000,0.,1.);
  TH2D *hYV2Rc2 = new TH2D("hYV2Rc2","",20,0.,1.,200,-1.,1.);

  //|eta|<0.9 pT>0.2, - TPC & TOF effficiency
  TH2D *hPtYRc3 = new TH2D("hPtYRc3","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtEtaRc3 = new TH2D("hPtEtaRc3","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtV2Rc3 = new TH2D("hPtV2Rc3","",200,0.,10.,200,-1.,1.);
  TH2D *hPtCosThetaRc3 = new TH2D("hPtCosThetaRc3","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRc3 = new TH2D("hPtCos2PhiRc3","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRPRc3 = new TH2D("hPtCos2PhiRPRc3","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCosTheta2Rc3 = new TH2D("hPtCosTheta2Rc3","",200,0.,10.,1000,0.,1.);
  TH2D *hMinvCos2PhiRc3 = new TH2D("hMinvCos2PhiRc3","",200,0.98,1.08,2000,-1.,1.);
  
  TH2D *hYCosThetaRc3 = new TH2D("hYCosThetaRc3","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRc3 = new TH2D("hYCos2PhiRc3","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRPRc3 = new TH2D("hYCos2PhiRPRc3","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCosTheta2Rc3 = new TH2D("hYCosTheta2Rc3","",20,0.,1.,1000,0.,1.);
  TH2D *hYV2Rc3 = new TH2D("hYV2Rc3","",20,0.,1.,200,-1.,1.);

  //|eta|<0.9 pT>0.2, - TBD
  TH2D *hPtYRc4 = new TH2D("hPtYRc4","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtEtaRc4 = new TH2D("hPtEtaRc4","",200,0.,10.0,200,-2.,2.);
  TH2D *hPtV2Rc4 = new TH2D("hPtV2Rc4","",200,0.,10.,200,-1.,1.);
  TH2D *hPtCosThetaRc4 = new TH2D("hPtCosThetaRc4","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRc4 = new TH2D("hPtCos2PhiRc4","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCos2PhiRPRc4 = new TH2D("hPtCos2PhiRPRc4","",200,0.,10.,2000,-1.,1.);
  TH2D *hPtCosTheta2Rc4 = new TH2D("hPtCosTheta2Rc4","",200,0.,10.,1000,0.,1.);
  TH2D *hMinvCos2PhiRc4 = new TH2D("hMinvCos2PhiRc4","",200,0.98,1.08,2000,-1.,1.);

  TH2D *hYCosThetaRc4 = new TH2D("hYCosThetaRc4","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRc4 = new TH2D("hYCos2PhiRc4","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCos2PhiRPRc4 = new TH2D("hYCos2PhiRPRc4","",20,0.,1.,2000,-1.,1.);
  TH2D *hYCosTheta2Rc4 = new TH2D("hYCosTheta2Rc4","",20,0.,1.,1000,0.,1.);
  TH2D *hYV2Rc4 = new TH2D("hYV2Rc4","",20,0.,1.,200,-1.,1.);

  TLorentzVector Phimom(0.,0.,0.,0.);
  TLorentzVector K1mom(0.,0.,0.,0.);   // Kaon momentum in the lab frame
  TLorentzVector K2mom(0.,0.,0.,0.);   // Kaon momentum in the lab frame
  TLorentzVector PhimomRc(0.,0.,0.,0.);   // Phi momentum Rc in the lab frame
  TLorentzVector K1momStar(0.,0.,0.,0.);   // Kaon momentum Rc in the Phi rest frame

  for(int ievt = 0; ievt<Nevt; ievt++) {
    if(ievt%100000==0) cout << " processing " << ievt << "-th event ..." << endl;
    double MassMom = gRandom->Gaus(MassPhi, SigmaPhi);
    //    double pT = gRandom->Rndm()*2.; //funMom->GetRandom();
    double y;
    if(sigY>9.) y = gRandom->Rndm()*4. - 2. ; // input 9.9 - flat rapidity
    else {
      do {
	y = gRandom->Gaus(0, sigY);
      } while (fabs(y)>2.0);
    }
    //    double phi = gRandom->Rndm()*TMath::Pi()*2.;
    //    double mT = TMath::Sqrt(pT*pT+MassMom*MassMom);
    double mT = fmT->GetRandom();
    double pT = TMath::Sqrt(mT*mT - MassMom*MassMom);
    if(pT<0.) continue;
    double v2 = (mT - MassMom) > 1.2 ? 0.2 :  0.2/1.2 * (mT - MassMom);

    double Psi_RP = gRandom->Rndm()*TMath::Pi()*2.0;

    fPhi->SetParameter(0, v2);
    double phi = fPhi->GetRandom() + Psi_RP;
    if(phi>TMath::Pi()*2.) phi -= TMath::Pi()*2.;
    
    double pz = mT * TMath::SinH(y);
    double p = TMath::Sqrt(pT*pT+pz*pz);
    double eta = 0.5*TMath::Log((p+pz)/(p-pz));
    Phimom.SetPtEtaPhiM(pT, eta, phi, MassMom);

    hPtY->Fill(pT, y);
    hPtEta->Fill(pT, eta);
    hPtV2->Fill(pT, TMath::Cos(2.*(phi - Psi_RP)));

    double p_pri = TMath::Sqrt((MassMom+2*MassK)*(MassMom-2*MassK))/2.;
    double phi_pri = gRandom->Rndm()*TMath::Pi()*2.0;   // phi x-z
    double costheta_pri = gRandom->Rndm()*2. - 1.;
    double py_pri = p_pri*costheta_pri;    // theta - angle between Y-axis (pol axis) and momentum
    double pxz_pri = TMath::Sqrt(p_pri*p_pri-py_pri*py_pri);
    double px_pri = pxz_pri*TMath::Cos(phi_pri);
    double pz_pri = pxz_pri*TMath::Sin(phi_pri);
    K1mom.SetXYZM(px_pri, py_pri, pz_pri, MassK);
    K2mom.SetXYZM(-px_pri, -py_pri, -pz_pri, MassK);
    double phi_star = K1mom.Phi();
    double cos2phi = TMath::Cos(2.*(phi_star - phi));

    // transform to the lab frame
    K1mom.Boost(Phimom.BoostVector());
    K2mom.Boost(Phimom.BoostVector());

    PhimomRc = K1mom + K2mom;
    K1momStar = K1mom;
    K1momStar.Boost(-PhimomRc.BoostVector());
    //    double costheta_rc = K1momStar.Y()/K1momStar.Vect().Mag();
    TVector3 L_unit(TMath::Cos(Psi_RP + TMath::Pi()*0.5), TMath::Sin(Psi_RP + TMath::Pi()*0.5), 0);
    double costheta_rc = K1momStar.Vect().Dot(L_unit)/K1momStar.Vect().Mag();
    if(fabs(y)<1.0) {
      hPtCosTheta->Fill(pT, costheta_rc);
      hPtCos2Phi->Fill(pT, cos2phi);
      hPtV2->Fill(pT, TMath::Cos(2.*(phi - Psi_RP)));
      hPtCos2PhiRP->Fill(pT, TMath::Cos(2.*(K1mom.Phi() - Psi_RP)));
      hPtCosTheta2->Fill(pT, costheta_rc*costheta_rc);

      hCosThetaDiff->Fill(costheta_pri, costheta_rc-costheta_pri);

      if(pT>1.2 && pT<5.4) {
	hYCosTheta->Fill(fabs(y), costheta_rc);
	hYCos2Phi->Fill(fabs(y), cos2phi);
	hYCos2PhiRP->Fill(fabs(y), TMath::Cos(2.*(phi_star - Psi_RP)));
	hYCosTheta2->Fill(fabs(y), costheta_rc*costheta_rc);
	hYV2->Fill(y, TMath::Cos(2.*(phi - Psi_RP)));
	hMinvCos2Phi->Fill(MassMom, cos2phi);
      }
    }

    // Check the momentum
    if(0) {
      cout << " Phi momentum = " << Phimom.Px() << " " << Phimom.Py() << " " << Phimom.Pz() << " " << Phimom.E() << endl;
      cout << " K1 momentum = " << K1mom.Px() << " " << K1mom.Py() << " " << K1mom.Pz() << " " << K1mom.E() << endl;
      cout << " K2 momentum = " << K2mom.Px() << " " << K2mom.Py() << " " << K2mom.Pz() << " " << K2mom.E() << endl;
    }

    if(fabs(K1mom.Eta())<0.9 && fabs(K2mom.Eta())<0.9) {
      hPtYRc->Fill(pT, y);
      hPtEtaRc->Fill(pT, eta);
      if(fabs(y)<1.0) {
	hPtCosThetaRc->Fill(pT, costheta_rc);
	hPtCos2PhiRc->Fill(pT, cos2phi);
	hPtV2Rc->Fill(pT, TMath::Cos(2.*(PhimomRc.Phi() - Psi_RP)));
	hPtCos2PhiRPRc->Fill(pT, TMath::Cos(2.*(phi_star - Psi_RP)));
	hPtCosTheta2Rc->Fill(pT, costheta_rc*costheta_rc);

	if(pT>1.2 && pT<5.4) {
	  hYCosThetaRc->Fill(fabs(y), costheta_rc);
	  hYCos2PhiRc->Fill(fabs(y), cos2phi);
	  hYCos2PhiRPRc->Fill(fabs(y), TMath::Cos(2.*(phi_star - Psi_RP)));
	  hYCosTheta2Rc->Fill(fabs(y), costheta_rc*costheta_rc);
	  hYV2Rc->Fill(y, TMath::Cos(2.*(phi - Psi_RP)));
	  hMinvCos2PhiRc->Fill(MassMom, cos2phi);
	}
      }

      if(K1mom.Pt()>0.2 && K2mom.Pt()>0.2) {
	hPtYRc1->Fill(pT, y);
	hPtEtaRc1->Fill(pT, eta);
	if(fabs(y)<1.0) {
	  hPtCosThetaRc1->Fill(pT, costheta_rc);
	  hPtCos2PhiRc1->Fill(pT, cos2phi);
	  hPtV2Rc1->Fill(pT, TMath::Cos(2.*(PhimomRc.Phi() - Psi_RP)));
	  hPtCos2PhiRPRc1->Fill(pT, TMath::Cos(2.*(phi_star - Psi_RP)));
	  hPtCosTheta2Rc1->Fill(pT, costheta_rc*costheta_rc);

	  if(pT>1.2 && pT<5.4) {
	    hYCosThetaRc1->Fill(fabs(y), costheta_rc);
	    hYCos2PhiRc1->Fill(fabs(y), cos2phi);
	    hYCos2PhiRPRc1->Fill(fabs(y), TMath::Cos(2.*(phi_star - Psi_RP)));
	    hYCosTheta2Rc1->Fill(fabs(y), costheta_rc*costheta_rc);
	    hYV2Rc1->Fill(y, TMath::Cos(2.*(phi - Psi_RP)));
	    hMinvCos2PhiRc1->Fill(MassMom, cos2phi);
	  }
	}
	
	hPtPhiK1->Fill(pT, K1mom.Pt());
	hPtPhiK2->Fill(pT, K2mom.Pt());

	if(gRandom->Rndm()<gr_eff_tpc->Eval(K1mom.Pt()) &&
	   gRandom->Rndm()<gr_eff_tpc->Eval(K2mom.Pt())) {
          hPtYRc2->Fill(pT, y);
	  hPtEtaRc2->Fill(pT, eta);
	  if(fabs(y)<1.0) {
	    hPtCosThetaRc2->Fill(pT, costheta_rc);
	    hPtCos2PhiRc2->Fill(pT, cos2phi);
	    hPtV2Rc2->Fill(pT, TMath::Cos(2.*(PhimomRc.Phi() - Psi_RP)));
	    hPtCos2PhiRPRc2->Fill(pT, TMath::Cos(2.*(phi_star - Psi_RP)));
	    hPtCosTheta2Rc2->Fill(pT, costheta_rc*costheta_rc);

	    if(pT>1.2 && pT<5.4) {
	      hYCosThetaRc2->Fill(fabs(y), costheta_rc);
	      hYCos2PhiRc2->Fill(fabs(y), cos2phi);
	      hYCos2PhiRPRc2->Fill(fabs(y), TMath::Cos(2.*(phi_star - Psi_RP)));
	      hYCosTheta2Rc2->Fill(fabs(y), costheta_rc*costheta_rc);
	      hYV2Rc2->Fill(y, TMath::Cos(2.*(phi - Psi_RP)));
	      hMinvCos2PhiRc2->Fill(MassMom, cos2phi);
	    }
	  }
	  
	  if(gRandom->Rndm()<gr_eff_tof->Eval(K1mom.Pt()) &&
	     gRandom->Rndm()<gr_eff_tof->Eval(K2mom.Pt())) {
	    hPtYRc3->Fill(pT, y);
	    hPtEtaRc3->Fill(pT, eta);
	    if(fabs(y)<1.0) {
	      hPtCosThetaRc3->Fill(pT, costheta_rc);
	      hPtCos2PhiRc3->Fill(pT, cos2phi);
	      hPtV2Rc3->Fill(pT, TMath::Cos(2.*(PhimomRc.Phi() - Psi_RP)));
	      hPtCos2PhiRPRc3->Fill(pT, TMath::Cos(2.*(phi_star - Psi_RP)));
	      hPtCosTheta2Rc3->Fill(pT, costheta_rc*costheta_rc);

	      if(pT>1.2 && pT<5.4) {
		hYCosThetaRc3->Fill(fabs(y), costheta_rc);
		hYCos2PhiRc3->Fill(fabs(y), cos2phi);
		hYCos2PhiRPRc3->Fill(fabs(y), TMath::Cos(2.*(phi_star - Psi_RP)));
		hYCosTheta2Rc3->Fill(fabs(y), costheta_rc*costheta_rc);
		hYV2Rc3->Fill(y, TMath::Cos(2.*(phi - Psi_RP)));
		hMinvCos2PhiRc3->Fill(MassMom, cos2phi);
	      }
	    }

	    // phi-Psi-dependent efficiency
	    double dphi1 = K1mom.Phi() - Psi_RP;
	    double dphi2 = K2mom.Phi() - Psi_RP;
	    if(gRandom->Rndm()<fun_eff_phi->Eval(dphi1) &&
	       gRandom->Rndm()<fun_eff_phi->Eval(dphi2)) {
	      hPtYRc4->Fill(pT, y);
	      hPtEtaRc4->Fill(pT, eta);
	      if(fabs(y)<1.0) {
		hPtCosThetaRc4->Fill(pT, costheta_rc);
		hPtCos2PhiRc4->Fill(pT, cos2phi);
		hPtV2Rc4->Fill(pT, TMath::Cos(2.*(PhimomRc.Phi() - Psi_RP)));
		hPtCos2PhiRPRc4->Fill(pT, TMath::Cos(2.*(phi_star - Psi_RP)));
		hPtCosTheta2Rc4->Fill(pT, costheta_rc*costheta_rc);
		
		if(pT>1.2 && pT<5.4) {
		  hYCosThetaRc4->Fill(fabs(y), costheta_rc);
		  hYCos2PhiRc4->Fill(fabs(y), cos2phi);
		  hYCos2PhiRPRc4->Fill(fabs(y), TMath::Cos(2.*(phi_star - Psi_RP)));
		  hYCosTheta2Rc4->Fill(fabs(y), costheta_rc*costheta_rc);
		  hYV2Rc4->Fill(y, TMath::Cos(2.*(phi - Psi_RP)));
		  hMinvCos2PhiRc4->Fill(MassMom, cos2phi);
		}
		
	      }
	      
	    }
	  } // end if	  

	}

      }
    }

  }

  TFile *fout = new TFile(Form("accept_Sergei_v2pt_Minv.root"),"recreate");
  hPtY->Write();
  hPtEta->Write();
  hPtCosTheta->Write();
  hPtCos2Phi->Write();
  hPtV2->Write();
  hPtCos2PhiRP->Write();
  hPtCosTheta2->Write();
  hYCosTheta->Write();
  hYCos2Phi->Write();
  hYCos2PhiRP->Write();
  hYCosTheta2->Write();
  hYV2->Write();
  
  hPtYRc->Write();
  hPtEtaRc->Write();
  hPtCosThetaRc->Write();
  hCosThetaDiff->Write();
  hPtCos2PhiRc->Write();
  hPtV2Rc->Write();
  hPtCos2PhiRPRc->Write();
  hPtCosTheta2Rc->Write();
  hYCosThetaRc->Write();
  hYCos2PhiRc->Write();
  hYCos2PhiRPRc->Write();
  hYCosTheta2Rc->Write();
  hYV2Rc->Write();
  
  hPtYRc1->Write();
  hPtEtaRc1->Write();
  hPtCosThetaRc1->Write();
  hPtCos2PhiRc1->Write();
  hPtV2Rc1->Write();
  hPtCos2PhiRPRc1->Write();
  hPtCosTheta2Rc1->Write();
  hYCosThetaRc1->Write();
  hYCos2PhiRc1->Write();
  hYCos2PhiRPRc1->Write();
  hYCosTheta2Rc1->Write();
  hYV2Rc1->Write();
  
  hPtYRc2->Write();
  hPtEtaRc2->Write();
  hPtCosThetaRc2->Write();
  hPtCos2PhiRc2->Write();
  hPtV2Rc2->Write();
  hPtCos2PhiRPRc2->Write();
  hPtCosTheta2Rc2->Write();
  hYCosThetaRc2->Write();
  hYCos2PhiRc2->Write();
  hYCos2PhiRPRc2->Write();
  hYCosTheta2Rc2->Write();
  hYV2Rc2->Write();

  hPtYRc3->Write();
  hPtEtaRc3->Write();
  hPtCosThetaRc3->Write();
  hPtCos2PhiRc3->Write();
  hPtV2Rc3->Write();
  hPtCos2PhiRPRc3->Write();
  hPtCosTheta2Rc3->Write();
  hYCosThetaRc3->Write();
  hYCos2PhiRc3->Write();
  hYCos2PhiRPRc3->Write();
  hYCosTheta2Rc3->Write();
  hYV2Rc3->Write();

  hPtYRc4->Write();
  hPtEtaRc4->Write();
  hPtCosThetaRc4->Write();
  hPtCos2PhiRc4->Write();
  hPtV2Rc4->Write();
  hPtCos2PhiRPRc4->Write();
  hPtCosTheta2Rc4->Write();
  hYCosThetaRc4->Write();
  hYCos2PhiRc4->Write();
  hYCos2PhiRPRc4->Write();
  hYCosTheta2Rc4->Write();
  hYV2Rc4->Write();

  hPtPhiK1->Write();
  hPtPhiK2->Write();


  hMinvCos2Phi->Write();
  hMinvCos2PhiRc->Write();
  hMinvCos2PhiRc1->Write();
  hMinvCos2PhiRc2->Write();
  hMinvCos2PhiRc3->Write();
  hMinvCos2PhiRc4->Write();
  
  fout->Close();

}
