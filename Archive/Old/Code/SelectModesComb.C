TString round(double n, double d);

void SelectModesComb(TString Run="24"){

  fstream puricut;
  puricut.open("modes/txt/TableRun"+Run+"_comb_cuts.txt",fstream::in);
  fstream puri;
  puri.open("modes/txt/TableRun"+Run+"_comb.txt",fstream::in);
  fstream truth;
  truth.open("modes/txt/TableRun24_truthmatched.txt",fstream::in);

  int SeedMode, YMode, type[90][54][6], Ds, ll, De, Dmu, Kl, els;
  double purity[90][54], Tmatched[90][54][6];
  TString rub;
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
    puricut>>SeedMode>>YMode>>rub>>Ds>>ll>>De>>Dmu>>Kl>>els;
    type[SeedMode][YMode][0] = Ds;
    type[SeedMode][YMode][1] = ll;
    type[SeedMode][YMode][2] = De;
    type[SeedMode][YMode][3] = Dmu;
    type[SeedMode][YMode][4] = Kl;
    type[SeedMode][YMode][5] = els;
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
  optimi.open("modes/txt/TableRun"+Run+"_opti_comb.txt",fstream::out);
  for(int Ncut = 900; Ncut>0; Ncut -= 100){
    if(Ncut==800) Ncut = 200;
    //for(int Ncut = 0; Ncut<100; Ncut += 50){
    for(int cut = 0; cut<100; cut += 1){
      Ds = ll = De = Dmu = Kl = els = 0;
      for(int i=0; i<87; i++){
	if(i==2) i = 10; if(i==19) i = 20; if(i==27) i = 30;
	if(i==49) i = 50;if(i==59) i = 60; if(i==62) i = 71;
	if(i==73) i = 80;if(i==85) i = 86;
	for(int j=1; j<54; j++){
	  if (Tmatched[i][j][1]*100>=cut && Tmatched[i][j][2]*100<=Ncut){
	    Ds += type[i][j][0];
	    ll += type[i][j][1];
	    De += type[i][j][2];
	    Dmu += type[i][j][3];
	    Kl += type[i][j][4];
	    els += type[i][j][5];
	  }
	}
      }
      optimi<<round((double)cut,100.)<<"  "<<round((double)Ncut,100.)<<" |\t"<<Ds+ll+De+Dmu+Kl+els;
      optimi<<"\t"<<Ds<<"\t"<<ll<<"\t"<<De<<"\t"<<Dmu<<"\t"<<Kl<<"\t"<<els<<endl;
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

