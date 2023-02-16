#include "style.C+"
#include "draw.C+"

const double hbarc = 197.327;  // MeV * fm
const double sqrtPi = TMath::Sqrt( TMath::Pi() );
const Int_t NC = 3;
const Int_t NP = 20;  // read in 40 points
const Int_t NPar = 7;  // R0, R1, R2, f00, d00, f01, d01
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
  return 0.5*sqrtPi*exp(-z*z)/z*ROOT::Math::erf(z);
  /*
    const Int_t N = 100;
    const double dt = z/N;
    
    double sum = 0;
    for(int i=0;i<N;i++) {
    double t = (i+0.5)*z/N;
    sum += exp(t*t-z*z)/z * dt;
    }
    return sum;
  */
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
//const double kstar, const double R, const double f00, const double d00, const double f01, const double d01)
{
  return 1 + 1./3*LL(x[0], par[0], par[1], par[2]);
}

double CF_dl(double *x, double *par)
//const double kstar, const double R, const double f00, const double d00, const double f01, const double d01)
{
  return 1 + 1./3*LL(x[0], par[0], par[1], par[2]) + 2./3*LL(x[0], par[0], par[3], par[4]);
}


void draw()
{
  style();

  const Double_t par0[2][2] = {{-16.8, 2.3}, {-16.3, 3.2}};
  const Double_t par1[3][2] = {{17.3, 3.6}, {10.8, 3.8}, {7.6, 3.6}};
  
  TCanvas *c1 = new TCanvas("c1","c1",0,0,800,600);
  c1->Draw();
  
  TH1D *hist = new TH1D("hist","",1,0,200);
  hist->SetMinimum(0.);
  hist->SetMaximum(20.);
  hist->Draw();
  
  TF1 *cf[6];
  for(int i=0;i<6;i++) {
    cf[i] = new TF1(Form("func_CF_dl_%d",i),CF_dl,0,200,5);
    cf[i]->SetParameters(2.5, par0[i/3][0], par0[i/3][1], par1[i%3][0], par1[i%3][1]);
    cf[i]->SetLineWidth(2);
    cf[i]->SetLineColor(i%3+1);
    cf[i]->SetLineStyle(i/3+1);
    cf[i]->Draw("same");
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
    for(int i=0;i<gr_data[ic]->GetN();i++){
      double kstar = gr_data[ic]->GetX()[i];
      
      sum+=pow((gr_data[ic]->GetY()[i]-cf_fun[ic]->Eval(kstar))/gr_data[ic]->GetEY()[i],2.0);
    }
  }
  f=sum;
  //  cout<<" chi2/ndf="<<sum/(gr_data[0]->GetN()*NC-NPar)<<endl;
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
  for(int ic=0;ic<NC;ic++) {
    cf_fun[ic] = new TF1(Form("cf_fun_%d",ic),CF_dl,0,200,5);
    cf_fun_0[ic] = new TF1(Form("cf_fun_0_%d",ic), CF_dl_0,0,200,3);
  }

  TCanvas *c1 = new TCanvas("c1","c1",0,0,600,900);
  c1->Divide(1,3);
  c1->Draw();
  
  //  double var[NPar] = {2.8, 2.6, 2.5, -16., 2.3, 10., 3.6};
  double var[NPar] = {2.8, 2.6, 2.5, -16., 2.3, 10., 3.6};
  double verr[NPar] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.01};
  double var_min[NPar] = {2.0, 2.0, 2.0, -100., 2.3, 0., 3.6};
  double var_max[NPar] = {4.0, 4.0, 4.0, -0., 2.3, 20., 3.6};
  
  double err_mat[NPar][NPar];
  
  TMinuit *gMinuit = new TMinuit(NPar);
  gMinuit->SetFCN(fcn);

  double arglist[10];
  Char_t PARM_NAMES[NPar][255]={"R0","R1","R2","f0_0","d0_0","f0_1","d0_1"};
  
  double PARM_START[2]={1,1};
  double PARM_STEP[2]={0.1,0.1};  
  int ivarbl,ierflg;
  double bnd1, bnd2;
  double minchi2=100;
  double a1=100;
  double a2=50;
  double eff = 1;
  double p1 = 0;
  int npars = NPar;
  int matrix = 0;
  for(int i=0;i<NPar;i++) {
    gMinuit->mnparm( i, PARM_NAMES[i], var[i], verr[i], var_min[i], var_max[i], ivarbl);
  }
  
  arglist[0] = 0.001;            // do at least 1000 function calls
  gMinuit->mnexcm("SET EPSmachine", arglist, 2, ierflg );
  
  arglist[0] = 10000;
  arglist[1] = 0.001;
  gMinuit->mnexcm("MINOS", arglist, 2, ierflg );
  cout<<"********************************************************************************************************"<<endl;
  cout<<""<<endl;
  cout<<""<<endl;
  cout<<"Here are the results:"<<endl;
  cout<<""<<endl;
  gMinuit->mnstat(minchi2,a2, eff,npars, npars, matrix);
  cout<<"chi2/NDF="<<minchi2/(gr_data[0]->GetN()*3-NPar+2)<<endl;
  cout<<""<<endl;

  for(int i=0;i<NPar;i++) {
    gMinuit->GetParameter(i, var[i], verr[i]);
  }

  for(int ip=0;ip<NC;ip++) {
    c1->cd(ip+1);
  
    TH1D *hist = new TH1D("hist","",1,0,100);
    hist->SetMinimum(0.);
    hist->SetMaximum(15.);
    hist->GetXaxis()->SetTitle("k* (MeV/c)");
    hist->GetYaxis()->SetTitle("d-#Lambda Correlation");
    hist->Draw();
    
    gr_data[ip]->Draw("p");
  
    cf_fun[ip]->SetParameter(0, var[ip]);
    cf_fun[ip]->SetParameter(1, var[3]);
    cf_fun[ip]->SetParameter(2, var[4]);
    cf_fun[ip]->SetParameter(3, var[5]);
    cf_fun[ip]->SetParameter(4, var[6]);
    cf_fun[ip]->SetLineStyle(1);
    cf_fun[ip]->SetLineWidth(3);
    cf_fun[ip]->Draw("c same");

    cf_fun_0[ip]->SetParameter(0, var[ip]);
    cf_fun_0[ip]->SetParameter(1, var[3]);
    cf_fun_0[ip]->SetParameter(2, var[4]);
    cf_fun_0[ip]->SetLineStyle(2);
    cf_fun_0[ip]->SetLineWidth(3);
    cf_fun_0[ip]->Draw("c same");
    //    gr_data[ip]->SetMarkerSize(1.5);
    //    gr_data[ip]->SetMarkerStyle(20);
    gr_data[ip]->Draw("p");

    drawText(70, 12, LabelName[ip], 42, 0.08);
    drawText(60, 9, Form("R_{G} = %3.1f #pm %3.1f fm", var[ip], verr[ip]), 42, 0.07);
  }

  c1->Update();
  c1->SaveAs("CF_dl.pdf");
  c1->SaveAs("CF_dl.png");
  
}
