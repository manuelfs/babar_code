TString round(double n, int e, double d=1.);

void ratiosGenDss(TString Dss="", int expo=0, int expTex=2){
  TString texType = "\\bf{Signal }$\\boldsymbol{";
  if(Dss=="dss") {
    texType = "$\\boldsymbol{D^{(*)}\\pi^0\\,\\,";
  }
  texType+="(10^{\\mbox{-}"; texType+=expTex; texType+="})}$";
  TString row[]={"$D_1(\\to D^*\\pi)\\ell\\nu$","$D_2^*(\\to D^{(*)}\\pi)\\ell\\nu$",
		 "$D_0^*(\\to D\\pi)\\ell\\nu$", "$D_1^{\'}(\\to D^*\\pi)\\ell\\nu$",
		 "$D_1\\tau\\nu$","$D_2^*\\tau\\nu$","$D_0^*\\tau\\nu$","$D_1^{\'}\\tau\\nu$"};

  int lund[] = {421,423};
  TString texName = "babar_code/Tables/texRatioGenDss"; texName+=Dss;texName+="DEff.txt";
  fstream tex;
  tex.open(texName,fstream::out);
  TString txName = "babar_code/Tables/tx";txName+=Dss;txName+="DssEff.txt";
  fstream txCoc;
  txCoc.open(txName,fstream::in);
  TString txGenName = "babar_code/Tables/tx";txGenName+=Dss;txGenName+="DssEffGen.txt";
  fstream txGen;
  txGen.open(txGenName,fstream::in);
  tex<<"\\begin{tabular}{|l||c|c|c|c|}"<<endl<<"\\hline"<<endl;
  tex<<texType<<" & $D^0$ & $D^{*0}$ & $D^+$ & $D^{*+}$ \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  double coc, coce, gen, gene;
  for(int i=0; i<8; i++){
    if(i==4) tex<<"\\hline"<<endl;
    tex<<row[i];
    cout<<row[i]<<":\t";
    int ifix = i;
    for(int j=0; j<4; j++){ 
      txCoc>>coc>>coce;
      txGen>>gen>>gene;
      double err = -1;
      if(coc!=0) err = sqrt(pow(gene/coc,2)+pow(gen*coce/coc/coc,2));
      tex<<" & "<<round(gen*pow(10,expo),1,coc)<<" $\\pm$ "<<round(err*pow(10,expo),1);
      cout<<round(gen*pow(10,expo),2,coc)<<" +- "<<round(err*pow(10,expo),2)<<"\t";
    }
    tex<<" \\\\"<<endl;
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


