TString round(double n, int e, double d=1.);

void DssEff(TString Dss="", int expo = 4){
  TString semName = "AWG82/ntuples/small/FitDss_RunAll.root";
  TChain sem("ntp1");
  sem.Add(semName);
  TString Dss3Name = "AWG82/ntuples/small/FitDpipi_RunAll.root";
  TChain Dss3("ntp1");
  Dss3.Add(Dss3Name);
  TString texType = "\\bf{Signal }$\\boldsymbol{";
  if(Dss=="dss") {
    texType = "$\\boldsymbol{D^{(*)}\\pi^0\\,\\,";
    //expo = 5;
  }
  texType+="(10^{\\mbox{-}"; texType+=expo; texType+="})}$";

  double dss = 10282000+9912000,dsstau = 3354000,dss3 = 5124000+4940000, effi[6][4];
  double bfD1 = 56,bfD2 = 37,bfD0 = 20,bfD1p = 37,bfNon = 10+20+30+60;
  double bftotDss = bfD1+bfD2+bfD0+bfD1p+bfNon, bftotDssTau = bfD1+bfD2+bfD0+bfD1p;
  double bftotDss3 = bfD1+bfD1p;

  double genDss[] = {dss*bfD1/bftotDss,dss*bfD2/bftotDss,dss*bfD0/bftotDss,
		     dss*bfD1p/bftotDss,dss3*bfD1/bftotDss3,dss3*bfD1p/bftotDss3,
		     dsstau*bfD1/bftotDssTau,dsstau*bfD2/bftotDssTau,dsstau*bfD0/bftotDssTau,
		     dsstau*bfD1p/bftotDssTau};
  TString row[]={"$D_1(\\to D^*\\pi)\\ell\\nu$","$D_2^*(\\to D^{(*)}\\pi)\\ell\\nu$",
		 "$D_0^*(\\to D\\pi)\\ell\\nu$", "$D_1^{\'}(\\to D^*\\pi)\\ell\\nu$",
		 "$D_1(\\to D\\pi\\pi)\\ell\\nu$","$D_1^{\'}(\\to D^{(*)}\\pi\\pi)\\ell\\nu$",
		 "$D_1\\tau\\nu$","$D_2^*\\tau\\nu$","$D_0^*\\tau\\nu$","$D_1^{\'}\\tau\\nu$"};

  int lund[] = {10423, 425, 10421, 20423, 10423, 20423};
  TString texName = "babar_code/Tables/tex";texName+=Dss;texName+="DssEff.txt";
  fstream tex;
  tex.open(texName,fstream::out);
  TString txName = "babar_code/Tables/tx";txName+=Dss;txName+="DssEff.txt";
  fstream tx;
  tx.open(txName,fstream::out);
  tex<<"\\begin{tabular}{|l||c|c|c|c|}"<<endl<<"\\hline"<<endl;
  tex<<texType<<" & $D^0$ & $D^{*0}$ & $D^+$ & $D^{*+}$ \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;

  for(int i=0; i<4; i++){
    tex<<row[i];
    cout<<row[i]<<":\t";
    for(int j=1; j<5; j++){ 
      TString cut="MCTaumode==-1&&candType=="; cut+=j+4*(Dss=="dss"); 
      cut+="&&(abs(MCD)=="; cut+=lund[i]; cut+="||abs(MCD)=="; cut+=lund[i]-10; cut+=")";
      TCut tot = "1"; tot+=cut;
      double n = sem.GetEntries(tot);
      double N = genDss[i];
      effi[i][j-1] = n/N;
      double error = sqrt(n*N*N-n*n*N)/N/N;
      tex<<" & "<<round(n*pow(10,expo),1,N)<<" $\\pm$ "<<round(error*pow(10,expo),1);
      tx<<round(n*pow(10,expo),2,N)<<"\t"<<round(error*pow(10,expo),2)<<"\t";
      cout<<round(n*pow(10,expo),2,N)<<" +- "<<round(error*pow(10,expo),2)<<"\t";
    }
    tex<<" \\\\"<<endl;
    tx<<endl;
    cout<<endl;
  }
  tx<<endl;
  tex<<"\\hline"<<endl;
  for(int i=4; i<6; i++){
    tex<<row[i];
    cout<<row[i]<<":\t";
    for(int j=1; j<5; j++){ 
      TString cut="MCTaumode==-1&&candType=="; cut+=j+4*(Dss=="dss"); cut+="&&(abs(MCD)=="; cut+=lund[i];
      cut+="||abs(MCD)=="; cut+=lund[i]-10; cut+=")";
      TCut tot = "1"; tot+=cut;
      double n = Dss3.GetEntries(tot);
      double N = genDss[i];
      effi[i][j-1] = n/N;
      double error = sqrt(n*N*N-n*n*N)/N/N;
      tex<<" & "<<round(n*pow(10,expo),1,N)<<" $\\pm$ "<<round(error*pow(10,expo),1);
      cout<<round(n*pow(10,expo),2,N)<<" +- "<<round(error*pow(10,expo),2)<<"\t";
    }
    tex<<" \\\\"<<endl;
    cout<<endl;
  }
  tex<<"\\hline"<<endl;

  cout<<endl<<"Dss:\t";
  for(int chan=0;chan<4; chan++){
    double total = 0;
    for(int i=0; i<4; i++) total += effi[i][chan];
    cout<<round(total/4*pow(10,expo),2)<<", ";
  }
  cout<<endl<<"Dpipi:\t";
  for(int chan=0;chan<4; chan++){
    double total = 0;
    for(int i=4; i<6; i++) total += effi[i][chan];
    cout<<round(total/2*pow(10,expo),2)<<", ";
  }
  cout<<endl<<endl;

//   for(int i=6; i<10; i++){
//     tex<<row[i];
//     cout<<row[i]<<":\t\t\t";
//     for(int j=1; j<5; j++){ 
//       TString cut="MCTaumode>-1&&candType=="; cut+=j; cut+="&&abs(MCD)==";
//       if(j<3) cut+=lund[i-6];
//       else cut+=lund[i-6]-10;
//       TCut tot = selCut; tot+=cut;
//       double n = sem.GetEntries(tot);
//       double N = genDss[i]/(0.1778+0.1731); // Normalized with BF(tau->lnunu)
//       double error = sqrt(n*N*N-n*n*N)/N/N;
//       tex<<" & "<<round(n*pow(10,expo),1,N)<<" $\\pm$ "<<round(error*pow(10,expo),1);
//       tx<<round(n*pow(10,expo),2,N)<<"\t"<<round(error*pow(10,expo),2)<<"\t";
//       cout<<round(n*pow(10,expo),2,N)<<" +- "<<round(error*pow(10,expo),2)<<"\t";
//     }
//     tex<<" \\\\"<<endl;
//     cout<<endl;
//   }

  tex<<"\\hline \\end{tabular}"<<endl;
  tex<<"\\,\\,"<<endl<<endl;
  cout<<texName<<" done"<<endl;
}


TString round(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}


