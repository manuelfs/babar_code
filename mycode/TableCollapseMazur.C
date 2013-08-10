// Macro that collapses a bit table into a presentation friendly one

TString round(double n, double d);

void TableCollapseMazur(TString tname=""){
  if(tname==""){
    cout<<"I require the name of the table to collapse. Thanks"<<endl;
    return;
  }

  fstream Other;
  Other.open(tname,fstream::in);

  char rubbish[100];
  int MYields[28][4];
  for(int j=0; j<4; j++){
    for(int i=0; i<5; i++){
      Other>>rubbish;
    }
    for(int i=0; i<7; i++){
      if(i<5) Other>>rubbish;
      int row = 7*j+i;
      Other>>rubbish>>MYields[row][0]>>MYields[row][1]>>rubbish;
    }
  }

  int S[5][4];
  for(int j=0; j<4; j++){
    for(int i=0; i<5; i++){
      S[i][j] = 0;
    }
  }
  int c1=0, c2=1, lag = 0;

  for(int k=0; k<1; k++){
    for(int i=0; i<4; i++){
      if(i>1) lag = 2;
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
  
  TString tag = "Mcha", teg = "D^{(*)0}";
  TString tabName = "mycode/collapsed/"; tabName+=tag; tabName+=basename;
  fstream Table;
  Table.open(tabName,fstream::out);
  fstream tex;
  tex.open("mycode/collapsed/"+tag+"tex"+basename,fstream::out);

  double nMazur = 690.1;
  double nGen = 894.6; // 2/3 of runs 1-6 BSemiExclAdd
  lag = 0;
  TString TaNames[] ={"Signal", "Norm.", "D**", "Cross.", "Comb."};
  for(int k=0; k<2; k++){
    if(k==1){
      lag = 2;
      tag = "Neutral";
      teg = "D^{(*)+}";
    }
    Table<<tag<<"\t  New\t  Base     Ratio "<<"\t"<<endl;
    Table<<"===================================="<<endl;
    tex<<"\\begin{tabular}{|l||cc|c|}"<<endl<<"\\hline"<<endl;
    tex<<"$\\mathbf{"<<teg<<"}$  & New & M14 & Ratio \\\\"<<endl;
    tex<<"\\hline \\hline"<<endl;
    for(int i=0; i<5; i++){
      Table<<TaNames[i]<<"\t  "<<S[i][0+lag]<<"\t   "<<S[i][1+lag]<<"\t    ";
      Table<<round((double)S[i][0+lag]*nMazur/nGen,(double)S[i][1+lag])<<endl;
      tex<<TaNames[i]<<" & "<<S[i][0+lag]<<" & "<<S[i][1+lag]<<" & ";
      tex<<round((double)S[i][0+lag]*nMazur/nGen,(double)S[i][1+lag])<<" \\\\"<<endl;
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

