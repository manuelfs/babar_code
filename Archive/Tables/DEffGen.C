#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void DEffGen(TString Dss="", int expo = 5){
  TString semName = "AWG82/ntuples/small/Add_Aug20Rest_Run123456.root";
  TChain sem("ntp1");
  sem.Add(semName);
  TCut selCut = MvaAll;
  TString texType = "\\bf{Signal }$\\boldsymbol{";
  if(Dss=="dss") {
    selCut=dssMvaAll;
    texType = "$\\boldsymbol{D^{(*)}\\pi^0\\,\\,";
    //expo = 6;
  }
  texType+="(10^{\\mbox{-}"; texType+=expo; texType+="})}$";

  double rawBs[] = {36968000., 103498000., 50556000., 167994000., 244322000., 68148000.,
		    37200000., 103124000., 49766000., 167466000., 244812000., 68016000.};
  double Bs[] = {0,0};
  for(int i=0; i<6; i++){
    Bs[0] += rawBs[i]*2/3.;
    Bs[1] += rawBs[i+6]*2/3.;
  }
  double BF[] = {2.07, 5.7, 2.24, 6.17};
  double BFtau[] = {0.7, 1.6, 0.7, 1.6};
  double genB0[4], genBp[4];
  for(int i=0; i<2; i++){
    genB0[i] = Bs[0]*2*BF[i]/100.*2;    // Times 2 because of e/mu
    genBp[i] = Bs[1]*2*BF[i+2]/100.*2;  // The other 2 comes because there are 2 B's
  }
  for(int i=2; i<4; i++){
    genB0[i] = Bs[0]*2*BFtau[i-2]/100.;
    genBp[i] = Bs[1]*2*BFtau[i]/100.;
  }
  TString row[]={"$D\\ell\\nu$","$D^*\\ell\\nu$","$D\\tau\\nu$","$D^*\\tau\\nu$"};

  int lund[] = {421,423};
  TString texName = "babar_code/Tables/tex";texName+=Dss;texName+="DEffGen.txt";
  fstream tex;
  tex.open(texName,fstream::out);
  TString txName = "babar_code/Tables/tx";txName+=Dss;txName+="DEffGen.txt";
  fstream tx;
  tx.open(txName,fstream::out);
  tex<<"\\begin{tabular}{|l||c|c|c|c|}"<<endl<<"\\hline"<<endl;
  tex<<texType<<" & $D^0$ & $D^{*0}$ & $D^+$ & $D^{*+}$ \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  double dummy = -1;
  for(int i=0; i<4; i++){
    if(i==2) tex<<"\\hline"<<endl;
    tex<<row[i];
    cout<<row[i]<<":\t";
    int ifix = i;
    for(int j=1; j<5; j++){ 
      TString cut="MCTaumode==-1&&MCType<13&&candType=="; 
      if(i>1) {
	cut="MCTaumode>-1&&MCType<13&&candType=="; 
	ifix = i-2;
      }
      cut+=j; cut+="&&abs(MCD)==";
      if(j<3) cut+=lund[ifix];
      else cut+=lund[ifix]-10;
      TCut tot = selCut; tot+=cut;
      double n = sem.GetEntries(tot); dummy = n;
      double N = genB0[i];
      if(j<3) N = genBp[i];
      double error = sqrt(n*N*N-n*n*N)/N/N;
      tex<<" & "<<round(n*pow(10,expo),1,N)<<" $\\pm$ "<<round(error*pow(10,expo),1);
      tx<<round(n*pow(10,expo),2,N)<<"\t"<<round(error*pow(10,expo),2)<<"\t";
      cout<<round(n*pow(10,expo),2,N)<<" +- "<<round(error*pow(10,expo),2)<<"\t";
    }
    tex<<" \\\\"<<endl;
    tx<<endl;
    cout<<endl;
  }

  tex<<"\\hline \\end{tabular}"<<endl;
  tex<<"\\,\\,"<<endl<<endl;
  cout<<texName<<" done"<<endl;
}


TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10,e)+0.5);
  b /= pow(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".")) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}


