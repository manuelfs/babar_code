// Macro that collapses a bit table into a presentation friendly one
// type 1 is for a table with only the electron channel, 2 for muons
// and 3 collapses D0 with D*0, and D+ with D*+

TString round(double n, double d);

void TableCollapse(TString tname="", int type=3){
  if(tname==""){
    cout<<"I require the name of the table to collapse. Thanks"<<endl;
    return;
  }
  if(type==3)
    cout<<"Collapsing tables into D(*)0 and D(*)+"<<endl;
  if(type==1)
    cout<<"Collapsing tables for the electron mode"<<endl;
  if(type==2)
    cout<<"Collapsing tables for the muon mode"<<endl;

  fstream Other;
  Other.open(tname,fstream::in);

  char rubbish[100];
  int MYields[28][4];
  for(int j=0; j<4; j++){
    for(int i=0; i<8; i++){
      Other>>rubbish;
    }
    for(int i=0; i<7; i++){
      if(i<5) Other>>rubbish;
      int row = 7*j+i;
      Other>>rubbish>>MYields[row][0]>>MYields[row][1];
      Other>>rubbish>>MYields[row][2]>>MYields[row][3]>>rubbish;
    }
  }

  int S[5][4];
  for(int j=0; j<4; j++){
    for(int i=0; i<5; i++){
      S[i][j] = 0;
    }
  }
  int c1=0, c2=1, lag = 0;

  // Electrons
  if(type==1){
    c1 = 0; c2 = 1;
  } else if(type==2){
    c1 = 2; c2 = 3;
  }
  for(int k=0; k<2; k++){
    if(type<3 && k==1) break;
    if(k==1){
      c1 = 2; c2 = 3;
    }
    for(int i=0; i<4; i++){
      if(type==3 && i>1) lag = 2;
      else lag=0;
      int row = 7*i;
      S[0][0+lag] += MYields[row][c1];
      S[0][1+lag] += MYields[row][c2];
      if(i%2==1){
	S[0][0+lag] += MYields[row+1][c1];
	S[0][1+lag] += MYields[row+1][c2];
      } else {
	S[4][0+lag] += MYields[row+1][c1];
	S[4][1+lag] += MYields[row+1][c2];
      }
      S[1][0+lag] += MYields[row+2][c1]+MYields[row+3][c1];
      S[1][1+lag] += MYields[row+2][c2]+MYields[row+3][c2];
      S[2][0+lag] += MYields[row+4][c1];
      S[2][1+lag] += MYields[row+4][c2];
      S[3][0+lag] += MYields[row+5][c1];
      S[3][1+lag] += MYields[row+5][c2];
      S[4][0+lag] += MYields[row+6][c1];
      S[4][1+lag] += MYields[row+6][c2];
    }
  }
  TString basename(tname);
  Int_t slashpos = basename.Last('/');
  TString directory;
  if (slashpos>=0) {
    basename.Remove(0,slashpos+1);      
  } else {
    cout<<"Wrong name"<<endl;
    return;
  }
  
  TString tag = "elec", teg = "e";
  if(type==2) {
    tag = "muon";
    teg = "\\mu";
  } else if(type == 3){
    tag = "cha";
    teg = "D^{(*)0}";
  }
  TString tabName = "mycode/collapsed/"; tabName+=tag; tabName+=basename;
  fstream Table;
  Table.open(tabName,fstream::out);
  fstream tex;
  tex.open("mycode/collapsed/"+tag+"tex"+basename,fstream::out);

  lag = 0;
  TString TaNames[] ={"Signal", "Norm.", "D**", "Cross.", "Comb."};
  for(int k=0; k<2; k++){
    if(type<3 && k==1) break;
    if(k==1){
      lag = 2;
      tag = "Neutral";
      teg = "D^{(*)+}";
    }
    Table<<tag<<"\t  New\t  Base     Ratio "<<"\t"<<endl;
    Table<<"===================================="<<endl;
    tex<<"\\begin{tabular}{|l||cc|c|}"<<endl<<"\\hline"<<endl;
    tex<<"$\\mathbf{"<<teg<<"}$ \\textbf{mode} & New & Base & Ratio \\\\"<<endl;
    tex<<"\\hline \\hline"<<endl;
    for(int i=0; i<5; i++){
      Table<<TaNames[i]<<"\t  "<<S[i][0+lag]<<"\t   "<<S[i][1+lag]<<"\t    ";
      Table<<round((double)S[i][0+lag],(double)S[i][1+lag])<<endl;
      tex<<TaNames[i]<<" & "<<S[i][0+lag]<<" & "<<S[i][1+lag]<<" & ";
      tex<<round((double)S[i][0+lag],(double)S[i][1+lag])<<" \\\\"<<endl;
    }
    Table<<endl<<endl;
    tex<<"\\hline \\end{tabular}"<<endl;
    tex<<"\\,\\,"<<endl<<endl;
  }
  cout << tabName<<" stored"<<endl;
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

