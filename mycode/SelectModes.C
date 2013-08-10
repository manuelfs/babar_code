void SelectModes(){

  fstream puricut;
  puricut.open("mycode/Table_purities_cuts.txt",fstream::in);
  fstream puri;
  puri.open("mycode/Table_purities.txt",fstream::in);

  int SeedMode, YMode, type[90][54][5], EmptyMode[90][54], S, N, Dss, Cr, Co;
  double purity[90][54];
  TString rub, pur;
  for(int i=0; i<87; i++){
    if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
    if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
    if(i==73) i = 80;if(i==85) i = 86;
    for(int j=1; j<54; j++){
      purity[i][j] = 0;
      EmptyMode[i][j] = 1;
      for(int k=0; k<5; k++)
	type[i][j][k] = 0;
    }
  }
  while(puricut){
    puricut>>SeedMode>>YMode>>rub>>S>>N>>Dss>>Cr>>Co>>pur;
    type[SeedMode][YMode][0] = S;
    type[SeedMode][YMode][1] = N;
    type[SeedMode][YMode][2] = Dss;
    type[SeedMode][YMode][3] = Cr;
    type[SeedMode][YMode][4] = Co;
    if(pur=="-") purity[SeedMode][YMode] = 0;
    else purity[SeedMode][YMode] = (double)pur.Atof();
  }

  fstream Empty;
  Empty.open("mycode/Empty.txt",fstream::out);
  for(int i=0; i<87; i++){
    int n = 1;
    if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
    if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
    if(i==73) i = 80;if(i==85) i = 86;
    for(int j=1; j<54; j++){
      double den = (double)(type[i][j][2]+type[i][j][4]);
      double sig = (double)type[i][j][0];
      bool keep = true;
      if(den==0&&type[i][j][0]>0) 
	keep = false;
      else if(den && (sig/(sig+den)>=.2))
	keep = false;
      if(keep) {
	Empty<<10000+i*100+j<<" ";
	if (!(n % 12)) Empty<<endl;
	n++;
      }
    }
    if (((n-1) % 12)) Empty<<endl;
    Empty<<endl;
  }
}


