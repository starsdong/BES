#include "draw.C+"
#include "style.C"

void statTest()
{
  style();

  //////////=====
  const int ntot_datapts=11;

  ////data 0-5%
  float C42_ener_besNEW[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  float C42_main_Ref3X_besNEW[ntot_datapts]={0.4119085,0.540458,0.4011784,0.430751,0.31512,0.339941,0.602046,0.739693,0.632837,0.792955,0.900669};
  float C42_stat_Ref3X_besNEW[ntot_datapts]={0.2468968,0.187029,0.1347458,0.104085,0.120367,0.0731331,0.0950052,0.147006,0.0553742,0.250823,0.2084582};
  float C42_sys_Ref3X_besNEW[ntot_datapts]={0.1302093,0.0680524,0.06472035,0.0403519,0.0880467,0.0389432,0.0304235,0.1357538,0.1387023,0.1219,0.1393631};
  float C42_tot_Ref3X_err_besNEW[ntot_datapts]={0};
  
  ////UrQMD
  double cen_UQMD[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double urqmd_C42_05_data[ntot_datapts]={0.4753277, 0.4625844, 0.5253311, 0.594994, 0.6782202, 0.7327199, 0.903167,0.860031,0.83973130,0.829186,0.914062};
  double urqmd_C42_05_sterr[ntot_datapts]={0.06445722,0.06283556,0.06766454,0.0707438,0.05287615,0.05808925,0.063723,0.045722,0.055722,0.081655,0.072831};
  
  ////HRG CE
  double anar_cen[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double anar_C42[ntot_datapts]={0.40214780,0.44213425,0.5034468,0.58341961,0.63475107,0.67691625,0.73933998,0.76389,0.7953959,0.8117625,0.84005};
  
  ////Hydro
  double Hydro_cen[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
  double Hydro_C42[ntot_datapts]={0.34623,0.39264044,0.46380312,0.556624,0.58357922,0.605721,0.657514,0.758193,0.78623350,0.8008,0.838933};
  
  for(int y=0;y<ntot_datapts;y++)
    {
      C42_tot_Ref3X_err_besNEW[y]=sqrt(C42_stat_Ref3X_besNEW[y]*C42_stat_Ref3X_besNEW[y]+C42_sys_Ref3X_besNEW[y]*C42_sys_Ref3X_besNEW[y]);
      cout << C42_tot_Ref3X_err_besNEW[y] << endl;
    }

  ////////////////////////////////////
  // rename arrays for simplicity
  ////////////////////////////////////
  double sNN[ntot_datapts];
  double C42_data[ntot_datapts];
  double C42_err_data[ntot_datapts];
  const int n_model = 3;
  double C42_model[n_model][ntot_datapts];
  double C42_err_model[n_model][ntot_datapts];
  const Char_t *NameModel[n_model] = {"UrQMD", "HRG-CE", "Hydro"};
  for(int i=0;i<ntot_datapts;i++) {
    sNN[i] = C42_ener_besNEW[i];
    C42_data[i] = C42_main_Ref3X_besNEW[i];
    C42_err_data[i] = C42_tot_Ref3X_err_besNEW[i];
    
    for(int j=0;j<n_model;j++) {
      if(j==0) {
	C42_model[j][i] = urqmd_C42_05_data[i];
	C42_err_model[j][i] = urqmd_C42_05_sterr[i];
      } else if(j==1) {
	C42_model[j][i] = anar_C42[i];
	C42_err_model[j][i] = 0.0;
      } else {
	C42_model[j][i] = Hydro_C42[i];
	C42_err_model[j][i] = 0.0;
      }
    }
  }

  ///////////////////////////////////
  //////////===== chi2 test /////////
  ///////////////////////////////////

  const int n_test = 7;
  double chi2[n_model] = {0.0, 0.0, 0.0};

  for(int im =0; im<n_model; im++) {
    cout << " Model\t" << NameModel[im] << endl;
    for(int i=0;i<n_test;i++) {
      double res = ( C42_data[i] - C42_model[im][i] )/sqrt(TMath::Power(C42_err_data[i], 2.));//+TMath::Power(C42_err_model[im][i], 2.));
      chi2[im] += res*res;
      //TMath::Power( C42_data[i]-C42_model[im][i], 2.0)/(TMath::Power(C42_err_data[i], 2.)+TMath::Power(C42_err_model[im][i], 2.));
      cout << "\t\t energy = " << sNN[i] << "\t residual = " << res << endl;
    }
    cout << " Model\t" << NameModel[im] << "\t" << setw(12) << chi2[im] << "\t p-value = " << TMath::Prob(chi2[im], n_test) << endl;
  }
  
  
  ////////////////////////////////////////////
  //////////===== likihood null test /////////
  ////////////////////////////////////////////

  
}
