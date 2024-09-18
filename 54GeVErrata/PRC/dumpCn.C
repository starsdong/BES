void dumpCn()
{
  const Int_t NCum = 4;
  const Int_t NP = 3;
  const Char_t *PName[NP] = {"Pro","Pbar","Netp"};
  const Int_t index4Print[NP] = {2, 0, 1};

  TGraphErrors *gr_stat[NP][NCum];
  TGraphErrors *gr_sys[NP][NCum];
  
  //  TFile *fin_stat = new TFile("54GeV_data_4Nu/stat.y0p5.root");
  TFile *fin_stat = new TFile("54GeV_data_Sep18/eff_corrected/9bins/stat/stat.y0p5.root");
  for(int i=0;i<NP;i++) {
    for(int j=0;j<NCum;j++) {
      gr_stat[i][j] = (TGraphErrors *)fin_stat->Get(Form("%s_C%d",PName[i],j+1));
    }
  }
  fin_stat->Close();

  TFile *fin_sys[NP];
  for(int i=0;i<NP;i++) {
    fin_sys[i] = new TFile(Form("54GeV_data_Sep18/eff_corrected/9bins/sys/%s/Sys_%s_y0p5.root",PName[i],PName[i]));
    for(int j=0;j<NCum;j++) {
      gr_sys[i][j] = (TGraphErrors *)fin_sys[i]->Get(Form("%s_C%d_sys",PName[i],j+1));
    }
    fin_sys[i]->Close();
  }

  for(int index=0;index<NP;index++) {
    int i = index4Print[index];
    cout << endl;
    cout << "====" << PName[i] << "====" << endl;
    for(int j=0;j<NCum;j++) {
      cout << "\t +++ C_" << j+1 << " +++ " << endl;
      cout << "{";
      for(int k=0;k<gr_stat[i][j]->GetN();k++) {
        cout << gr_stat[i][j]->GetY()[k];
	if(k!=gr_stat[i][j]->GetN()-1) cout << ", ";
      }
      cout << "}," << endl;
      cout << endl;
      cout << "\t +++ C_" << j+1 << " stat err +++ " << endl;
      cout << "{";
      for(int k=0;k<gr_stat[i][j]->GetN();k++) {
	cout << gr_stat[i][j]->GetEY()[k];
	if(k!=gr_stat[i][j]->GetN()-1) cout << ", ";
      }
      cout << "}," << endl;
      cout << endl;
      cout << "\t +++ C_" << j+1 << " sys err +++ " << endl;
      cout << "{";
      for(int k=0;k<gr_stat[i][j]->GetN();k++) {
	cout << gr_sys[i][j]->GetEY()[k];
	if(k!=gr_stat[i][j]->GetN()-1) cout << ", ";
      }
      cout << "}," << endl;
      cout << "-----------------------------" << endl;
    }
  }
  
}
