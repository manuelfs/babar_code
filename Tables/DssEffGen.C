#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void DssEffGen(TString Dss="", int expo = 4){
  TString semName = "AWG82/ntuples/small/Add_Aug20Rest_Run123456.root";
  TChain sem("ntp1");
  sem.Add(semName);
  TCut selCut = MvaAll;
  TString texType = "\\bf{Signal }$\\boldsymbol{";
  if(Dss=="dss") {
    selCut=dssMvaAll;
    texType = "$\\boldsymbol{D^{(*)}\\pi^0\\,\\,";
    //expo = 5;
  }
  texType+="(10^{\\mbox{-}"; texType+=expo; texType+="})}$";

  double rawBs[] = {36968000., 103498000., 50556000., 167994000., 244322000., 68148000.,
		    37200000., 103124000., 49766000., 167466000., 244812000., 68016000.};
  double Bs[] = {0,0};
  for(int i=0; i<6; i++){
    Bs[0] += rawBs[i]*2/3.;
    Bs[1] += rawBs[i+6]*2/3.;
  }
  double bfD1[] = {52,56},bfD2[] = {23,30},bfD0[] = {45,49},bfD1p[] = {83,90},bfNon[] = {40,38};
  double bfD1tau[] = {13,13},bfD2tau[] = {20,20},bfD0tau[] = {13,13},bfD1ptau[] = {20,20};
  double BF[] = {52, 23, 45, 83, 56, 30, 49, 90};
  double BFtau[] = {13, 20, 13, 20, 13, 20, 13, 20};
  double genB0[8], genBp[8];
  for(int i=0; i<4; i++){
    genB0[i] = Bs[0]*BF[i]/10000.*2;
    genBp[i] = Bs[1]*BF[i+4]/10000.*2;
  }
  for(int i=4; i<8; i++){
    genB0[i] = Bs[0]*BFtau[i-4]/10000.;
    genBp[i] = Bs[1]*BFtau[i]/10000.;
  }
  TString row[]={"$D_1(\\to D^*\\pi)\\ell\\nu$","$D_2^*(\\to D^{(*)}\\pi)\\ell\\nu$",
		 "$D_0^*(\\to D\\pi)\\ell\\nu$", "$D_1^{\'}(\\to D^*\\pi)\\ell\\nu$",
		 "$D_1\\tau\\nu$","$D_2^*\\tau\\nu$","$D_0^*\\tau\\nu$","$D_1^{\'}\\tau\\nu$"};

  int lund[] = {10423, 425, 10421, 20423};
  TString texName = "babar_code/Tables/tex";texName+=Dss;texName+="DssEffGen.txt";
  fstream tex;
  tex.open(texName,fstream::out);
  TString txName = "babar_code/Tables/tx";txName+=Dss;txName+="DssEffGen.txt";
  fstream tx;
  tx.open(txName,fstream::out);
  tex<<"\\begin{tabular}{|l||c|c|c|c|}"<<endl<<"\\hline"<<endl;
  tex<<texType<<" & $D^0$ & $D^{*0}$ & $D^+$ & $D^{*+}$ \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  double dummy = -1;
  for(int i=0; i<8; i++){
    if(i==4) tex<<"\\hline"<<endl;
    tex<<row[i];
    cout<<row[i]<<":\t";
    int ifix = i;
    for(int j=1; j<5; j++){ 
      TString cut="MCTaumode==-1&&candType=="; 
      if(i>3) {
	cut="MCTaumode>-1&&candType=="; 
	ifix = i-4;
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


