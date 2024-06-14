#ifndef __Datapts_arry__
#define __Datapts_arry__
#include <TF1.h>
using namespace std;


//////////=====
const int ntot_datapts=11;


////data 0-5%
float C42_ener_besNEW[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
float C42_main_Ref3X_besNEW[ntot_datapts]={0.5405526,0.528129,0.4365745,0.421259,0.409057, 0.415694,0.602046,0.739693,0.632837,0.792955,0.900669};
float C42_stat_Ref3X_besNEW[ntot_datapts]={0.2612975,0.185524,0.1340996,0.105752,0.122237, 0.0726043,0.0950052,0.147006,0.0553742,0.250823,0.2084582};
float C42_sys_Ref3X_besNEW[ntot_datapts]={0.1524927,0.0541957,0.05197891,0.0680574,0.0957171, 0.0290416,0.0304235,0.1357538,0.1387023,0.1219,0.1393631};
float C42_tot_Ref3X_err_besNEW[ntot_datapts]={0};

////UrQMD
double cen_UQMD[ntot_datapts]={7.7,9.2,11.5,14.5,17.3,19.6,27,39,54.4,62.4,200};
double urqmd_C42_05_data[ntot_datapts]={0.5878511, 0.5694708, 0.612383, 0.718947, 0.7735332, 0.8779027, 0.903167,0.860031,0.83973130,0.829186,0.914062};
double urqmd_C42_05_sterr[ntot_datapts]={0.06807203,0.06603911,0.07039043,0.0735658,0.05434042,0.05982107,0.063723,0.045722,0.055722,0.081655,0.072831};

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
