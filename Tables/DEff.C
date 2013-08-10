#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void DEff(TString Dss="", int expo = 4){
  TString semName = "AWG82/ntuples/small/*Full_RunAll.root";
  TChain sem("ntp1");
  sem.Add(semName);
  TCut selCut = MvaAll;
  TString texType = "\\bf{Signal }$\\boldsymbol{";
  if(Dss=="dss") {
    selCut=dssMvaAll;
    texType = "$\\boldsymbol{D^{(*)}\\pi^0\\,\\,";
    expo = 5;
  }
  texType+="(10^{\\mbox{-}"; texType+=expo; texType+="})}$";
  //                D0l      D*0l     D+l      D*+l     D0tau    D*0tau   D+tau    D*+tau         
  double genBs[] = {4474000, 4477000, 4492000, 4475000, 4492000, 4488000, 4486000, 4492000+4486000};
  TString row[]={"$D\\ell\\nu$","$D^*\\ell\\nu$","$D\\tau\\nu$","$D^*\\tau\\nu$"};

  int lund[] = {421,423};
  TString texName = "babar_code/Tables/tex"; texName+=Dss;texName+="DEff.txt";
  fstream tex;
  tex.open(texName,fstream::out);
  TString txName = "babar_code/Tables/tx";txName+=Dss;txName+="DEff.txt";
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
      double N = genBs[i+(j>2)*2+(i>1)*2];
      if(i>1) N = N/(0.1778+0.1731); // Normalized with BF(tau->lnunu)
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


