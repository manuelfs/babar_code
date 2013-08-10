TString round(double n, double d);

void SelectModes(TString Run="24", TString Bden="All"){

  fstream puricut;
  puricut.open("modes/txt/TableRun"+Run+"_purities_cuts.txt",fstream::in);
  fstream puri;
  puri.open("modes/txt/TableRun"+Run+"_purities.txt",fstream::in);
  fstream truth;
  truth.open("modes/txt/TableRun24_truthmatched.txt",fstream::in);

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
  int Track, Neutral;
  double Tpuri, Tave, Npuri, Nave;
  while(truth){
    truth>>SeedMode>>YMode>>rub>>Track>>Tpuri>>Tave>>Neutral>>Npuri>>Nave;
    Tmatched[SeedMode][YMode][0] = (double)Track;
    Tmatched[SeedMode][YMode][1] = Tpuri;
    Tmatched[SeedMode][YMode][2] = Tave;
    Tmatched[SeedMode][YMode][3] = (double)Neutral;
    Tmatched[SeedMode][YMode][4] = Npuri;
    Tmatched[SeedMode][YMode][5] = Nave;
  }
  fstream optimi;
  TString fname = "modes/txt/TableRun"; fname+=Run;fname+="_optimization.txt";
  if(Bden=="Dss" || Bden=="Cr" || Bden=="Co") {
    fname="modes/txt/"; fname+=Bden; fname+=Run; fname+=".txt";
  }
  optimi.open(fname,fstream::out);
//   for(int Ncut = 900; Ncut>0; Ncut -= 100){
//     if(Ncut==800) Ncut = 200;
  for(int Ncut = 400; Ncut>50; Ncut -= 10){
    for(int Ccut = 400; Ccut>50; Ccut -= 10){
//   for(int Ncut = 0; Ncut<100; Ncut += 5){
//     for(int Ccut = 0; Ccut<100; Ccut += 5){
      S = N = Dss = Cr = Co = 0;
      for(int i=0; i<87; i++){
	if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
	if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
	if(i==73) i = 80;if(i==85) i = 86;
	int fcut = Ccut;
	if(i>=20&&i<40 || i>=60) fcut = Ncut;
	for(int j=1; j<54; j++){
	  if (Tmatched[i][j][2]*100<=fcut){
	    S += type[i][j][0];
	    N += type[i][j][1];
	    Dss += type[i][j][2];
	    Cr += type[i][j][3];
	    Co += type[i][j][4];
	  }
	}
      }
      double deno = (double)(S+Dss+Cr+Co);
      if(Bden=="Dss") deno = (double)(S+Dss);
      if(Bden=="Cr") deno = (double)(S+Cr);
      if(Bden=="Co") deno = (double)(S+Co);
      optimi<<round((double)Ccut,100.)<<"  "<<round((double)Ncut,100.)<<" |\t"<<round((double)S,deno);
      optimi<<"\t"<<round((double)(S*S),deno)<<"\t";
      optimi<<S<<"\t"<<N<<"\t"<<Dss<<"\t"<<Cr<<"\t"<<Co<<endl;
    }
    optimi<<endl;
    cout<<"Ncut "<<Ncut<<" done"<<endl;
  }

//   TString tcut = "55";
//   fstream Empty;
//   Empty.open("modes/txt/BSemiExclAdd_Empty_Truth"+tcut+".txt",fstream::out);
//   for(int i=0; i<87; i++){
//     int n = 1;
//     if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
//     if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
//     if(i==73) i = 80;if(i==85) i = 86;
//     for(int j=1; j<54; j++){
//       bool keep = false;
//       if(Tmatched[i][j][1]*100. < tcut.Atof()) 
// 	keep = true;
//       if(keep) {
// 	Empty<<10000+i*100+j<<" ";
// 	if (!(n % 12)) Empty<<endl;
// 	n++;
//       }
//     }
//     if (((n-1) % 12)) Empty<<endl;
//     Empty<<endl;
//   }
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

