TString round(double n, int e, double d=1.);

void table(TString filename) {
  fstream File;
  File.open(filename,fstream::in);
  filename.Remove(0,filename.Last('/')+1);

  TString texName = "babar_code/Reweight/Tables/tex"; texName += filename;
  fstream tex;
  tex.open(texName,fstream::out);

  double yields[14][4];
  for(int i=0; i<14*4; i++){
    File >> yields[i%14][i/14];
  }
  TString modes[] = {"D^0", "D^{*0}", "D^+", "D^{*+}", "D_1^0", "D_2^{*0}",
		     "D_0^{*0}", "D_1^{\'0}", "D_1^+", "D_2^{*+}",
		     "D_0^{*+}", "D_1^{\'+}", "NR^0", "NR^+"};

  tex<<"\\begin{tabular}{|l||rr|rr|}"<<endl<<"\\hline"<<endl;
  tex<<" & D^0 & D^{*0} & D^+ & D^{*+} \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;

  for(int i=0; i<14; i++) {
    tex<<modes[i];
    for(int j=0; j<4; j++) {
      tex<<" & "<<round(yields[i][j],0);
    }
    tex<<" \\\\"<<endl;
    if(i%4 == 3) tex<<"\\hline"<<endl;
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


