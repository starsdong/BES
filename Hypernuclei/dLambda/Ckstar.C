#include "style.C+"
#include "draw.C+"
//#define __EFF__  // allow efficiency for doublet?
#define __FIT__
//#define __CONTOUR__

const double hbarc = 197.327;  // MeV * fm
const double sqrtPi = TMath::Sqrt( TMath::Pi() );
const Int_t NC = 3;
const Int_t NP = 40;  // read in 60 points
const Double_t KMAX = NP*5.0;  // K* max
#ifdef __EFF__
const Int_t NPar = 8;  // R0, R1, R2, f00, d00, f01, d01, eff
#else
const Int_t NPar = 7;  // R0, R1, R2, f00, d00, f01, d01,
#endif
const Char_t *Name[NC] = {"0to10","10to20","20to60"};
const Char_t *LabelName[NC] = {"0-10%","10-20%","20-60%"};
TGraphErrors *gr_data[NC];  //
TF1 *cf_fun[NC];

complex<double> fkstar(const double f0, const double d0, const double kstar)   // f0, d0 in fm,  kstar in MeV
{
  complex<double> val(1./f0 + 0.5*d0*kstar*kstar/hbarc/hbarc, -kstar/hbarc);
  return pow(val, -1.);
}

const double F1(const double z)
{
  //  return 0.5*sqrtPi*exp(-z*z)/z*ROOT::Math::erf(z);
    const Int_t N = 200;
    const double dt = z/N;
    
    double sum = 0;
    for(int i=0;i<N;i++) {
      double t = (i+0.5)*z/N;
      sum += exp(t*t-z*z)/z * dt;
    }
    return sum;
}

const double F2(const double z)
{
  return (1. - exp(-z*z)) / z;
}

const double LL(const double kstar, const double R, const double f0, const double d0)
{
  double sum = 0.;
  
  // Calculate fs
  complex<double> fsval = fkstar(f0, d0, kstar);
  double fsval_re = fsval.real();
  double fsval_im = fsval.imag();
  
  // Calculate first term in sum
  double fsval_mag2overR2 = (fsval_re * fsval_re + fsval_im * fsval_im) / (R*R);
  sum += 0.5 * fsval_mag2overR2 * ( 1. - d0/(2.*sqrtPi*R) );
  
  // Calculate second term
  double z = 2. * (kstar/hbarc) * R;
  
  sum += 2. * fsval_re * F1(z) / (sqrtPi * R);
  
  // Calculate third term
  sum += -1.0*fsval_im * F2(z) / R;
  
  return sum;
}


double CF_dl_0(double *x, double *par)
{
  return 1 + 1./3 * LL(x[0], par[0], par[1], par[2]);
}

double CF_dl_0_eff(double *x, double *par)
{
  double eff = par[3];
  return 1 + eff/(eff+2.) * LL(x[0], par[0], par[1], par[2]);
}

double CF_dl_1(double *x, double *par)
{
  return 1 + 2./3 * LL(x[0], par[0], par[1], par[2]);
}

double CF_dl(double *x, double *par)
{
  return 1 + 1./3 * LL(x[0], par[0], par[1], par[2]) + 2./3 * LL(x[0], par[0], par[3], par[4]);
}

double CF_dl_eff(double *x, double *par)
{
  double eff = par[5];
  return 1 + eff/(eff+2.) * LL(x[0], par[0], par[1], par[2]) + 2./(eff+2.) * LL(x[0], par[0], par[3], par[4]);
}

void draw()
{
  style();

  //  const Double_t par0[2][2] = {{-16.8, 2.3}, {-16.3, 3.2}};
  //  const Double_t par1[3][2] = {{17.3, 3.6}, {10.8, 3.8}, {7.6, 3.6}};
  const Double_t par0[2][2] = {{-16.8, 2.3}, {-40, 2.3}};
  const Double_t par1[3][2] = {{17.3, 3.6}, {17.3, 3.6}, {17.3, 3.6}};
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->Draw();
  
  TH1D *hist = new TH1D("hist","",1,0,200);
  hist->SetMinimum(0.);
  hist->SetMaximum(20.);
  hist->Draw();
  
  TF1 *cf[6];
  TF1 *cf_0[6];  // doublet only
  TF1 *cf_1[6];  // quartet only
  const Double_t RG = 2.5;
  ofstream outData;
  for(int i=0;i<6;i++) {
    //    if(i%3!=0) continue;
    outData.open(Form("model/dl_model_R_2.5_D_%d_Q_%d.txt", i/3, i%3));

    cf[i] = new TF1(Form("func_CF_dl_%d",i),CF_dl,0,KMAX,5);
    cf[i]->SetParameters(RG, par0[i/3][0], par0[i/3][1], par1[i%3][0], par1[i%3][1]);
    cf[i]->SetLineWidth(2);
    cf[i]->SetLineColor(i%3+1);
    //    cf[i]->SetLineColor(1);
    cf[i]->SetLineStyle(i/3+1);
    cf[i]->Draw("same");

    cf_0[i] = new TF1(Form("func_CF_dl_0_%d",i),CF_dl_0,0,KMAX,3);
    cf_0[i]->SetParameters(RG, par0[i/3][0], par0[i/3][1]);
    cf_0[i]->SetLineWidth(2);
    cf_0[i]->SetLineColor(i%3+1);
    //    cf_0[i]->SetLineColor(2);
    cf_0[i]->SetLineStyle(i/3+1);
    //    cf_0[i]->Draw("same");

    
    cf_1[i] = new TF1(Form("func_CF_dl_1_%d",i),CF_dl_1,0,KMAX,3);
    cf_1[i]->SetParameters(RG, par1[i%3][0], par1[i%3][1]);
    cf_1[i]->SetLineWidth(2);
    cf_1[i]->SetLineColor(i%3+1);
    //    cf_1[i]->SetLineColor(3);
    cf_1[i]->SetLineStyle(i/3+1);
    //    cf_1[i]->Draw("same");

    for(int j=0;j<200;j++) {
      double kstar = (j+0.5);
      std::cout << " kstar = " << kstar << " cf_0 = " << cf_0[i]->Eval(kstar) << "\t cf_1 = " << cf_1[i]->Eval(kstar) << "\t cf = " << cf[i]->Eval(kstar) << "\t cf_01 = " << cf_0[i]->Eval(kstar) + cf_1[i]->Eval(kstar) - 1.0 << std::endl;
      outData << setw(10) << kstar << setw(12) << cf_0[i]->Eval(kstar) << setw(12) << cf_1[i]->Eval(kstar) << setw(12) << cf[i]->Eval(kstar) << std::endl;
    }
    outData.close();
  }
   
}


void fcn(int &npar, double *gin, double &f, double *par,int iflag ) { 
  double sum=0;
  for(int ic=0;ic<NC;ic++) {
    //    if(ic==0) continue;
    cf_fun[ic]->SetParameter(0, par[ic]);
    cf_fun[ic]->SetParameter(1, par[3]);
    cf_fun[ic]->SetParameter(2, par[4]);
    cf_fun[ic]->SetParameter(3, par[5]);
    cf_fun[ic]->SetParameter(4, par[6]);
#ifdef __EFF__
    cf_fun[ic]->SetParameter(5, par[7]);
#endif
    for(int i=0;i<gr_data[ic]->GetN();i++){
      double kstar = gr_data[ic]->GetX()[i];
      
      sum+=pow((gr_data[ic]->GetY()[i]-cf_fun[ic]->Eval(kstar))/gr_data[ic]->GetEY()[i],2.0);
    }
  }
  f=sum;
  //  cout<<" chi2/ndf="<<sum/(gr_data[0]->GetN()*NC-NPar)<<endl;
}

double calChi2(double *par)
{
  double sum=0;
  for(int ic=0;ic<NC;ic++) {
    //    if(ic==0) continue;
    cf_fun[ic]->SetParameter(0, par[ic]);
    cf_fun[ic]->SetParameter(1, par[3]);
    cf_fun[ic]->SetParameter(2, par[4]);
    cf_fun[ic]->SetParameter(3, par[5]);
    cf_fun[ic]->SetParameter(4, par[6]);
#ifdef __EFF__
    cf_fun[ic]->SetParameter(5, par[7]);
#endif
    for(int i=0;i<gr_data[ic]->GetN();i++){
      double kstar = gr_data[ic]->GetX()[i];
      
      sum+=pow((gr_data[ic]->GetY()[i]-cf_fun[ic]->Eval(kstar))/gr_data[ic]->GetEY()[i],2.0);
    }
  }
  return sum;
}

void Ckstar()
{
  style();
  
  TGraphErrors *gr_tmp[NC];
  for(int i=0;i<NC;i++) {
    gr_tmp[i] = new TGraphErrors(Form("data/kstar_ppc_%s.txt",Name[i]),"%lg %lg %lg");
    double x_tmp[NP];
    for(int j=0;j<NP;j++) {
      x_tmp[j] = gr_tmp[i]->GetX()[j]*1e3; // GeV->MeV
    }
    gr_data[i] = new TGraphErrors(NP, x_tmp, gr_tmp[i]->GetY(), 0, gr_tmp[i]->GetEY());
  }

  TF1 *cf_fun_0[NC];  // doublet contribution only
  TF1 *cf_fun_1[NC];  // quartet contribution only
  for(int ic=0;ic<NC;ic++) {
#ifdef __EFF__
    cf_fun[ic] = new TF1(Form("cf_fun_%d",ic),CF_dl_eff,0,KMAX,6);
    cf_fun_0[ic] = new TF1(Form("cf_fun_0_%d",ic), CF_dl_0_eff,0,KMAX,4);
    cf_fun_1[ic] = new TF1(Form("cf_fun_1_%d",ic), CF_dl_1,0,KMAX,3);  // no difference between w/ and w/o efficiency
#else
    cf_fun[ic] = new TF1(Form("cf_fun_%d",ic),CF_dl,0,KMAX,5);
    cf_fun_0[ic] = new TF1(Form("cf_fun_0_%d",ic), CF_dl_0,0,KMAX,3);
    cf_fun_1[ic] = new TF1(Form("cf_fun_1_%d",ic), CF_dl_1,0,KMAX,3);
#endif
  }

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,900);
  c1->SetLeftMargin(0.1);
  c1->Divide(1,3);
  c1->Draw();
  
  //  double var[NPar] = {2.8, 2.6, 2.5, -16., 2.3, 10., 3.6};
#ifdef __EFF__
  double var[NPar] = {3.0, 2.9, 2.4, -16., 2.3, 17.3, 3.6, 1.0};
  double verr[NPar] = {0.3, 0.3, 0.3, 0.3, 0.1, 0.3, 0.1, 0.05};
  double var_min[NPar] = {2.0, 2.0, 2.0, -50., 2.3, 0., 3.6, 0.0};
  double var_max[NPar] = {4.0, 4.0, 4.0, -0., 2.3, 50., 3.6, 1.0};
#else
#ifdef __FIT__
  //  double var[NPar] = {3.0, 2.9, 2.4, -10., 2.3, 10, 3.6};
  double var[NPar] = {3.0, 2.9, 2.5, -20., 3.2, 10, 3.6};
  double verr[NPar] = {0.3, 0.3, 0.3, 0.3, 0.1, 0.3, 0.1};
  double var_min[NPar] = {2.5, 2.5, 2.0, -50., 1.0, 0., 2.0};
  double var_max[NPar] = {3.5, 3.5, 3.0, -0., 4.0, 20., 5.0};
#else
  // double var[NPar] = {2.96646, 2.8894, 2.44925, -5.68461, 2.3, 13.2108, 3.6};
  // double verr[NPar] = {0.3, 0.3, 0.3, 0.3, 0.1, 0.3, 0.1};
  // double var_min[NPar] = {2.96646, 2.8894, 2.44925, -5.68461, 2.3, 13.2108, 3.6};
  // double var_max[NPar] = {2.96646, 2.8894, 2.44925, -5.68461, 2.3, 13.2108, 3.6};
  double var[NPar] = {2.79479, 2.73791, 2.19929, -37.9044, 2.3, 6.5963, 3.6};
  double verr[NPar] = {0.3, 0.3, 0.3, 0.3, 0.1, 0.3, 0.1};
#endif  
#endif
  
  double err_mat[NPar][NPar];

#ifdef __FIT__
  TMinuit *gMinuit = new TMinuit(NPar);
  gMinuit->SetFCN(fcn);
#endif

  double arglist[10];
#ifdef __EFF__
  Char_t PARM_NAMES[NPar][255]={"R0","R1","R2","f0_0","d0_0","f0_1","d0_1", "eff"};
#else
  Char_t PARM_NAMES[NPar][255]={"R0","R1","R2","f0_0","d0_0","f0_1","d0_1"};
#endif
  
  int ivarbl,ierflg;
  double bnd1, bnd2;
  double minchi2=500;
  double a1=100;
  double a2=50;
  double eff = 1;
  double p1 = 0;
  int npars = NPar;
  int matrix = 0;
#ifdef __FIT__
  for(int i=0;i<NPar;i++) {
    gMinuit->mnparm( i, PARM_NAMES[i], var[i], verr[i], var_min[i], var_max[i], ivarbl);
    gMinuit->FixParameter(4);
    gMinuit->FixParameter(6);
  }

  arglist[0] = 0.0001;
  gMinuit->mnexcm("SET EPSmachine", arglist, 2, ierflg );

  arglist[0] = 100000;
  arglist[1] = 0.0001;
  gMinuit->mnexcm("MINOS", arglist, 2, ierflg );
  //  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg );
  cout<<"********************************************************************************************************"<<endl;
  cout<<""<<endl;
  cout<<"Here are the results:"<<endl;
  cout<<""<<endl;
  gMinuit->mnstat(minchi2,a2, eff,npars, npars, matrix);
  cout<<"chi2 / NDF = "<<minchi2 << " / " << (gr_data[0]->GetN()*3-NPar+2)<<endl;
  cout<<""<<endl;

  for(int i=0;i<NPar;i++) {
    gMinuit->GetParameter(i, var[i], verr[i]);
    std::cout << PARM_NAMES[i] << "\t" << var[i] << " +/- " << verr[i] << std::endl;
  }

#ifdef __CONTOUR__
  TGraph *gr_err[3];
  for(int i=0;i<3;i++) {
    //    if(i!=0) continue;
    gMinuit->SetErrorDef((i+1)*(i+1));
    gr_err[i] = (TGraph *)gMinuit->Contour(100, 3, 5);
  }
#endif
  
#endif  
  
  for(int ip=0;ip<NC;ip++) {
    c1->cd(ip+1);
  
    TH1D *hist = new TH1D("hist","",1,0,KMAX);
    hist->SetMinimum(0.);
    hist->SetMaximum(5.);
    hist->GetXaxis()->SetLabelOffset(0.015);
    hist->GetXaxis()->SetTitleOffset(1.0);
    hist->GetXaxis()->SetTitle("k* (MeV/c)");
    hist->GetYaxis()->SetTitle("d-#Lambda Correlation");
    hist->Draw();
    
    gr_data[ip]->Draw("p");
  
    cf_fun[ip]->SetParameter(0, var[ip]);
    cf_fun[ip]->SetParameter(1, var[3]);
    cf_fun[ip]->SetParameter(2, var[4]);
    cf_fun[ip]->SetParameter(3, var[5]);
    cf_fun[ip]->SetParameter(4, var[6]);
#ifdef __EFF__
    cf_fun[ip]->SetParameter(5, var[7]);
#endif
    cf_fun[ip]->SetLineStyle(1);
    cf_fun[ip]->SetLineWidth(3);
    cf_fun[ip]->Draw("c same");

    cf_fun_0[ip]->SetParameter(0, var[ip]);
    cf_fun_0[ip]->SetParameter(1, var[3]);
    cf_fun_0[ip]->SetParameter(2, var[4]);
#ifdef __EFF__
    cf_fun_0[ip]->SetParameter(3, var[7]);
#endif
    cf_fun_0[ip]->SetLineStyle(2);
    cf_fun_0[ip]->SetLineWidth(3);
    cf_fun_0[ip]->Draw("c same");
    gr_data[ip]->SetMarkerSize(1.0);
    gr_data[ip]->SetMarkerStyle(20);
    gr_data[ip]->Draw("p");

    drawText(70, 12, LabelName[ip], 42, 0.08);
    drawText(60, 9, Form("R_{G} = %4.2f #pm %4.2f fm", var[ip], verr[ip]), 42, 0.065);
  }
  double chi2 = calChi2(var);
  cout<<"chi2 / NDF = "<<chi2 << " / " << (gr_data[0]->GetN()*3-NPar+2)<<endl;

  c1->Update();
  c1->SaveAs("CF_dl.pdf");
  c1->SaveAs("CF_dl.png");


#ifdef __CONTOUR__
  TCanvas *c2 = new TCanvas("c2","c2",600,0,600,600);
  c1->SetBottomMargin(0.1);
  c1->SetLeftMargin(0.1);
  c1->Draw();

  TH1D *hC = new TH1D("hC","", 1, -100, 0);
  hC->SetMinimum(0);
  hC->SetMaximum(20);
  hC->GetXaxis()->SetLabelOffset(0.01);
  hC->GetXaxis()->SetTitle("f_{D} (fm)");
  hC->GetYaxis()->SetTitle("f_{Q} (fm)");
  hC->Draw();

  for(int i=0;i<3;i++) {
    //    if(i!=0) continue;
    gr_err[i]->SetLineColor(i+1);
    gr_err[i]->SetLineWidth(2);
    gr_err[i]->Draw("c");
  }

  drawHistBox(-100, 0, 0, 20);
#endif

}
