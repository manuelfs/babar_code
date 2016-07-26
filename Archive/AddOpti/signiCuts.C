//---------------------------------------------------------------------------------
// Description:
//      Third stage of the mode optimization
//      Macro that calculates "babar_code/txt/TableRun" + Run + "_optimization.txt"
//      with the significance and purity for different cuts of the BReco mode
//      purity
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/05/23 manuelf -- Adaptation from old SelectModes.C
//---------------------------------------------------------------------------------

TString round(double n, double d);

void signiCuts(TString truthName="babar_code/txt/Newtruth.txt", TString yieldName="babar_code/txt/Newyield.txt"){

  TString optimB0Name = "babar_code/txt/Newoptim.txt";

  fstream truthFile;
  truthFile.open(truthName,fstream::in);
  fstream yieldFile;
  yieldFile.open(yieldName,fstream::in);
  fstream optimB0File;
  optimB0File.open(optimB0Name,fstream::out);

  int SeedMode, YMode, type[90][54][6], S, N, Dss, Cr, Co;
  double purity[90][54], Tmatched[90][54][6];
  TString rub, pur;
  for(int i=0; i<87; i++){
    if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
    if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
    if(i==73) i = 80;if(i==85) i = 86;
    for(int j=1; j<54; j++){
      purity[i][j] = 0;
      for(int k=0; k<6; k++){
	type[i][j][k] = 0;
	Tmatched[i][j][k] = 0;
      }
    }
  }
  while(yieldFile){
    yieldFile>>SeedMode>>YMode>>rub>>S>>N>>Dss>>Cr>>Co>>pur;
    type[SeedMode][YMode][0] = S;
    type[SeedMode][YMode][1] = N;
    type[SeedMode][YMode][2] = Dss;
    type[SeedMode][YMode][3] = Cr;
    type[SeedMode][YMode][4] = Co;
    if(pur=="-") purity[SeedMode][YMode] = 0;
    else purity[SeedMode][YMode] = (double)pur.Atof();
  }
  int Track, Neutral;
  double Tpuri, Tave, Npuri, Nave;
  int S_B0, N_B0, Dss_B0, Cr_B0, Co_B0;
  while(truthFile){
    truthFile>>SeedMode>>YMode>>rub>>Track>>Tpuri>>Tave>>Neutral>>Npuri>>Nave;
    Tmatched[SeedMode][YMode][0] = (double)Track;
    Tmatched[SeedMode][YMode][1] = Tpuri;
    if(Track<3) Tmatched[SeedMode][YMode][1] = 0;
    Tmatched[SeedMode][YMode][2] = Tave;
    Tmatched[SeedMode][YMode][3] = (double)Neutral;
    Tmatched[SeedMode][YMode][4] = Npuri;
    Tmatched[SeedMode][YMode][5] = Nave;
  }
  double Ncuts[] = {10, 1, 0.5};
  for(int Ncut = 0; Ncut<3; Ncut++){
    for(int Ccut = 0; Ccut<100; Ccut += 1){
      S_B0 = N_B0 = Dss_B0 = Cr_B0 = Co_B0 = 0;
      for(int i=0; i<87; i++){
	if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
	if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
	if(i==73) i = 80;if(i==85) i = 86;

	for(int j=1; j<54; j++){
	  if (Tmatched[i][j][1]*100 >= Ccut  && Tmatched[i][j][5] <= Ncuts[Ncut]){
	    S_B0 += type[i][j][0];
	    N_B0 += type[i][j][1];
	    Dss_B0 += type[i][j][2];
	    Cr_B0 += type[i][j][3];
	    Co_B0 += type[i][j][4];
	  }
	}
      }
      double deno_B0 = (double)(S_B0+Dss_B0+Cr_B0+Co_B0);
      optimB0File<<round((double)Ccut,100.)<<"  "<<round((double)Ncuts[Ncut],1.)<<" |\t"<<round((double)S_B0,deno_B0);
      optimB0File<<"\t"<<round((double)(S_B0*S_B0),deno_B0)<<"\t";
      optimB0File<<S_B0<<"\t"<<N_B0<<"\t"<<Dss_B0<<"\t"<<Cr_B0<<"\t"<<Co_B0<<endl;
      if(Ccut%30 == 0)
	cout<<"Percentage of tracks truthmatched "<<Ccut<<" done"<<endl;
    }
    optimB0File<<endl;
    if(Ncut>3) Ncut=3;
  }

  cout<<"Yields for each cut stored in "<<optimB0Name<<endl;
}

TString round(double n, double d){
  if(d==0) return " - ";
  double b = ((int)(n/d *100+.5));
  b/=100.;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  return result;
}

