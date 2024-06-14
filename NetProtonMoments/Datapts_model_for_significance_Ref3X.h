#ifndef __Datapts_arry__
#define __Datapts_arry__
#include <TF1.h>
using namespace std;


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
}


//
//cout<<"==============================="<<endl<<endl;
//cout<<"======BEST data BES I n II NEW combined==============="<<endl;
//cout<<"========================================================"<<endl;
//
//cout<<"==============================="<<endl<<endl;
//
//for(int y=0;y<ntot_datapts;y++)
//{
//    cout<<"::deviation BESIInew wt UrQMD_Ref3X::ener::"<<cen_UQMD[y]<<"::devn::"<<(C42_main_Ref3X_besNEW[y]-urqmd_C42_05_data[y])/sqrt(C42_tot_Ref3X_err_besNEW[y]*C42_tot_Ref3X_err_besNEW[y]+urqmd_C42_05_sterr[y]*urqmd_C42_05_sterr[y])<<endl;
//}
//
//cout<<"==============================="<<endl<<endl;
//
//for(int y=0;y<ntot_datapts;y++)
//{
//    cout<<"::deviation BESIInew wt HRG CE::ener::"<<anar_cen[y]<<"::devn::"<<(C42_main_Ref3X_besNEW[y]-anar_C42[y])/C42_tot_Ref3X_err_besNEW[y]<<endl;
//}
//
//cout<<"==============================="<<endl<<endl;
//
//for(int y=0;y<ntot_datapts;y++)
//{
//    cout<<"::deviation BESIInew wt hydro::ener::"<<Hydro_cen[y]<<"::devn::"<<(C42_main_Ref3X_besNEW[y]-Hydro_C42[y])/C42_tot_Ref3X_err_besNEW[y]<<endl;
//}
//
//cout<<"==============================="<<endl<<endl;





#endif
