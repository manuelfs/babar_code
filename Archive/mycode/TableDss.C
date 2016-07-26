TString round(double n, double d);
void toTable(int cand, TChain &ch, int *DYields, int MYields[][2], TString *TaNames, TString *TeNames, 
	     TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag, bool doTruth);

void TableDss(TString folder, TString otherTag = "NN", TString Run="24", TString tag="All", bool doTruth=true){

  TCut cosT = "abs(candCosT)<.9";
  TCut mpi0 = "mpi0>.125&&mpi0<.145&&bestepi0==1";
  TCut dssee = "eextrapi0<.5&&ppi0>.4";
  TCut dssacc = "mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&pmisspi0>.2";
  TCut Mva = "(candMvaDssComb>0.18&&candMvaDssDl>0.1)";
//   TCut Mva = "(candMvaDssComb>-0.24&&candMvaDssDl>-0.2&&candType==1)||(candMvaDssComb>-0.09&&candMvaDssDl>-0.09&&candType==2)||(candMvaDssComb>-0.13&&candMvaDssDl>-0.24&&candType==3)||(candMvaDssComb>-0.09&&candMvaDssDl>-0.2&&candType==4)";
//   TCut Mva2 = "(candMvaDssComb>0.18&&candMvaDssDl>0.1)";
//   TCut Mva1 = "(candMvaDssDl>0.15&&candType==1)||(candMvaDssDl>0.22&&candType==2)||(candMvaDssDl>0.35&&candType==3)||(candMvaDssDl>0.19&&candType==4)";
  TCut basic = dssacc;

  if(tag=="All")
    basic = dssacc+mpi0+dssee;
  else if(tag=="afterBestB")
    basic = dssacc;
  else if(tag=="Mva")
    basic = dssacc+Mva;
  else if(tag=="Mva2")
    basic = dssacc+Mva2;
  else if(tag=="cosT")
    basic = dssacc+Mva2+cosT;

  // Loading files
  TString isAdd;
  int nfiles = 0;
  if(folder.Contains("Add"))
    isAdd="Add";
  else
    isAdd="noAdd";

  TChain ch("ntp1");
  //TString filename = "AWG82/ntuples/small/"; filename += folder; filename += ".root";
  //TString filename = "AWG82/ntuples/Add_BDT_KM/MergedRest2/MergedRest12345.root";
  TString filename = "AWG82/ntuples/Add_Final/MergedRest1/MergedRest123456.root";
  ch.Add(filename);
  cout<<"Added "<<filename<<" with entries: "<<ch.GetEntries()<<endl;   

  TString Mdir = "mycode/DssTables/"; Mdir += folder;
  gSystem->mkdir(Mdir);
  TString sTruth = "noTru"; if(doTruth) sTruth="";
  TString taFile = "mycode/DssTables/table"; taFile+=Run; taFile+=sTruth; taFile+="_";
  taFile+=folder; taFile+="_"; taFile+=otherTag; taFile+=tag; taFile+=".txt";
  fstream Table;
  Table.open(taFile,fstream::out);
  fstream tex;
  tex.open("mycode/DssTables/tex"+Run+sTruth+"_"+folder+"_"+otherTag+tag+".txt",fstream::out);
  fstream Other;
  Other.open("mycode/DssTables/table"+Run+sTruth+"_"+otherTag+".txt",fstream::in);
  //Other.open("mycode/DssTables/table_"+otherTag+".txt",fstream::in);
  //basic += "candBntCha==0&&candDntCha==0";
  char rubbish[100];


  int DYields[] = {133,8,2848,99,104,69,0};
  int MYields[7][2];
  for(int i=0; i<7; i++){
    int total = 4;
    if(i==0) total = 10;
    if(i>4) total = 3;
    for(int j=0; j<total; j++)
      Other>>rubbish;
    Other>>MYields[i][0];
    for(int j=0; j<2; j++)
      Other>>rubbish;
    Other>>MYields[i][1];
  }
  TString TaNames[]={"D*0 pi0","D0 pi0","Oth D**","D*0 l","D0 l","Cross","Comb.","D*0"};
  TString TeNames[]={"D^{**} \\to \\text{D}^{*0}\\pi^{0} & ","D^{**} \\to \\text{D}^{0}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*0} \\,l \\,\\nu & ","D^{0} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{*0}"};
  TCut cuts[]={"MCType>12&&MCDssmode==2","MCType>12&&MCDssmode==1","MCType>12&&MCDssmode>2",
	       "(MCType==2||MCType==4||MCType==6)","(MCType==1||MCType==3||MCType==5)","(MCType>6&&MCType<13)","MCType==0"};
  toTable(2, ch, DYields, MYields, TaNames, TeNames, cuts, Table, tex, basic, folder, otherTag, doTruth);


  int DYields2[] = {110,155,3875,1186,346,163,0};
  for(int i=0; i<7; i++){
    int total = 4;
    if(i==0) total = 12;
    if(i>4) total = 3;
    for(int j=0; j<total; j++)
      Other>>rubbish;
    Other>>MYields[i][0];
    for(int j=0; j<2; j++)
      Other>>rubbish;
    Other>>MYields[i][1];
  }
  TString TaNames2[]={"D0 pi0","D*0 pi0","Oth D**","D*0 l","D0 l","Cross","Comb.","D0"};
  TString TeNames2[]={"D^{**} \\to \\text{D}^{0}\\pi^{0} & ","D^{**} \\to \\text{D}^{*0}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*0} \\,l \\,\\nu & ","D^{0} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{0}"};
  TCut cuts2[]={"MCType>12&&MCDssmode==1","MCType>12&&MCDssmode==2","MCType>12&&MCDssmode>2",
	       "(MCType==2||MCType==4||MCType==6)","(MCType==1||MCType==3||MCType==5)","(MCType>6&&MCType<13)","MCType==0"};
  toTable(1, ch, DYields2, MYields, TaNames2, TeNames2, cuts2, Table, tex, basic, folder, otherTag, doTruth);


  int DYields3[] = {87,1,1866,9,92,65,0};
  for(int i=0; i<7; i++){
    int total = 4;
    if(i==0) total = 12;
    if(i>4) total = 3;
    for(int j=0; j<total; j++)
      Other>>rubbish;
    Other>>MYields[i][0];
    for(int j=0; j<2; j++)
      Other>>rubbish;
    Other>>MYields[i][1];
  }
  TString TaNames3[]={"D*+ pi0","D+ pi0","Oth D**","D*+ l","D+ l","Cross","Comb.","D*+"};
  TString TeNames3[]={"D^{**} \\to \\text{D}^{*+}\\pi^{0} & ","D^{**} \\to \\text{D}^{+}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*+} \\,l \\,\\nu & ","D^{+} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{*+}"};
  TCut cuts3[]={"MCType>12&&MCDssmode==4","MCType>12&&MCDssmode==3","MCType>12&&(MCDssmode<3||MCDssmode>4)",
	       "(MCType==8||MCType==10||MCType==12)","(MCType==7||MCType==9||MCType==11)","(MCType>0&&MCType<7)","MCType==0"};
  toTable(4, ch, DYields3, MYields, TaNames3, TeNames3, cuts3, Table, tex, basic, folder, otherTag, doTruth);



  int DYields4[] = {46,14,448,401,93,63,0};
  for(int i=0; i<7; i++){
    int total = 4;
    if(i==0) total = 12;
    if(i>4) total = 3;
    for(int j=0; j<total; j++)
      Other>>rubbish;
    Other>>MYields[i][0];
    for(int j=0; j<2; j++)
      Other>>rubbish;
    Other>>MYields[i][1];
  }
  TString TaNames4[]={"D+ pi0","D*+ pi0","Oth D**","D*+ l","D+ l","Cross","Comb.","D+"};
  TString TeNames4[]={"D^{**} \\to \\text{D}^{+}\\pi^{0} & ","D^{**} \\to \\text{D}^{*+}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*+} \\,l \\,\\nu & ","D^{+} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{+}"};
  TCut cuts4[]={"MCType>12&&MCDssmode==3","MCType>12&&MCDssmode==4","MCType>12&&(MCDssmode<3||MCDssmode>4)",
	       "(MCType==8||MCType==10||MCType==12)","(MCType==7||MCType==9||MCType==11)","(MCType>0&&MCType<7)","MCType==0"};
  toTable(3, ch, DYields4, MYields, TaNames4, TeNames4, cuts4, Table, tex, basic, folder, otherTag, doTruth);

  gSystem->CopyFile("mycode/DssTables/table"+Run+sTruth+"_"+folder+"_"+otherTag+tag+".txt",
		   "mycode/DssTables/table"+Run+sTruth+"_"+folder+".txt",kTRUE);
  cout<<taFile<<" written"<<endl;
}

//=============================================================================================================

void toTable(int cand, TChain &ch, int *DYields, int MYields[][2], TString *TaNames, TString *TeNames, 
	     TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag, bool doTruth){
  TCut eCut = "candIsMu==0";
  TCut muCut = "candIsMu==1";
  TCut truthMatch = "candLepTru==1";
  TString candCut = "candType=="; candCut += cand;
  Table<<" "<<TaNames[7]<<"\t  "<<folder<<"_e\t  Base_e   Ratio_e "<<"\t"<<folder;
  Table<<"_mu\tBase_mu   Ratio_mu "<<endl;
  Table<<"===================================================================="<<endl;
  tex<<"\\begin{tabular}{|l||cc|c|cc|c|}"<<endl<<"\\hline"<<endl;
  tex<<"$\\mathbf{"<<TeNames[7]<<"}$ & "<<folder<<"_e & Base_e & R_e & ";
  tex<<folder<<"_{\\mu} & Base_{\\mu} & R_{\\mu} \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  basic += candCut;
  TCut tCut = basic; tCut += eCut;
  int bkg_e = ch.GetEntries(tCut);
  tCut = basic; tCut += muCut;
  int bkg_mu = ch.GetEntries(tCut);
  for(int i=0; i<7; i++){
    TCut finalCut = basic; if(doTruth)finalCut += truthMatch; finalCut += cuts[i];
    int entries=0;
    if (i < 6){ 
      TCut feCut = finalCut; feCut += eCut;
       entries = ch.GetEntries(feCut);
       bkg_e -= entries;
    } else entries = bkg_e;
    Table<<TaNames[i]<<"\t   "<<entries<<"\t   "<<MYields[i][0];
    Table<<"\t    "<<round((double)entries,(double)MYields[i][0]);
    tex<<TeNames[i]<<entries<<" & "<<MYields[i][0];
    tex<<" & "<<round((double)entries,(double)MYields[i][0])<<" & ";
    if (i!= 6){ 
      TCut fmuCut = finalCut; fmuCut += muCut;
      entries = ch.GetEntries(fmuCut);
      bkg_mu -= entries;
    }else entries = bkg_mu;
    Table<<"   \t "<<entries<<"\t "<<MYields[i][1];
    Table<<"  \t    "<<round((double)entries,(double)MYields[i][1])<<endl;
    tex<<entries<<" & "<<MYields[i][1];
    tex<<" & "<<round((double)entries,(double)MYields[i][1])<<" \\\\"<<endl;
  }
  Table<<endl;
  tex<<"\\hline \\end{tabular}"<<endl;
  tex<<"\\,\\,"<<endl<<endl;
  cout<<TaNames[i]<<" done"<<endl;
  return;

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

