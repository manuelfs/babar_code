#include "mycode/cuts.hh"

TString round(double n, double d);
void toTable(int cand, TChain &ch, int MYields[][2], TString *TaNames, TString *TeNames, 
	     TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString Run, 
	     bool doTruth, TChain &uds, TChain &ccbar);

void TableDssconti(TString folder, TString otherTag = "NN", TString Run="24", TString tag="All", bool doTruth=true){

  if(tag.Contains("All"))
    basic = dss+mpi0+dssee;
  else if(tag=="afterBestB")
    basic = dss;
  else if(tag.Contains("Mva"))
    basic = dss+dssMva;
  if(tag.Contains("cosT"))
    basic += cosT;

  // Loading files
  TString isAdd;
  int nfiles = 0;
  if(folder.Contains("Add"))
    isAdd="Add";
  else
    isAdd="noAdd";

  TChain ch("ntp1");
  TString filename = "AWG82/ntuples/small/"; filename += folder; filename+="Rest_Run"; filename+=Run; filename += ".root";
  ch.Add(filename);
  cout<<"Added "<<filename<<" with entries: "<<ch.GetEntries()<<endl;   
  TChain uds("ntp1");
  filename = "AWG82/ntuples/small/uds30_Run12356.root";
  uds.Add(filename);
  cout<<"Added "<<filename<<" with entries: "<<uds.GetEntries()<<endl;   
  TChain ccbar("ntp1");
  filename = "AWG82/ntuples/small/ccbar30_Run12356.root";
  ccbar.Add(filename);
  cout<<"Added "<<filename<<" with entries: "<<ccbar.GetEntries()<<endl;   

  TString Mdir = "mycode/DssTables/"; Mdir += folder;
  gSystem->mkdir(Mdir);
  TString sTruth = "noTru"; if(doTruth) sTruth="";
  TString taFile = "mycode/DssTables/contable"; taFile+=Run; taFile+=sTruth; taFile+="_";
  taFile+=folder; taFile+="_"; taFile+=otherTag; taFile+=tag; taFile+=".txt";
  fstream Table;
  Table.open(taFile,fstream::out);
  fstream tex;
  tex.open("mycode/DssTables/contex"+Run+sTruth+"_"+folder+"_"+otherTag+tag+".txt",fstream::out);
  fstream Other;
  Other.open("mycode/DssTables/contable"+Run+sTruth+"_"+otherTag+".txt",fstream::in);
  //Other.open("mycode/DssTables/table_"+otherTag+".txt",fstream::in);
  //basic += "candBntCha==0&&candDntCha==0";
  char rubbish[100];

  otherTag = Run;
  int MYields[9][2];
  for(int i=0; i<9; i++){
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
  TString TaNames[]={"D*0 pi0","D0 pi0","Oth D**","D*0 l","D0 l","Cross","Comb.","uds", "ccbar","D*0"};
  TString TeNames[]={"D^{**} \\to \\text{D}^{*0}\\pi^{0} & ","D^{**} \\to \\text{D}^{0}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*0} \\,l \\,\\nu & ","D^{0} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","$uds$ & ",
		     "$c \\bar{c}$ & ","D^{*0}"};
  TCut cuts[]={"MCType>12&&MCDssmode==2","MCType>12&&MCDssmode==1","MCType>12&&MCDssmode>2",
	       "(MCType==2||MCType==4||MCType==6)","(MCType==1||MCType==3||MCType==5)","(MCType>6&&MCType<13)","MCType==0","",""};
  toTable(2, ch, MYields, TaNames, TeNames, cuts, Table, tex, basic, folder, otherTag, doTruth, uds, ccbar);

  for(int i=0; i<9; i++){
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
  TString TaNames2[]={"D0 pi0","D*0 pi0","Oth D**","D*0 l","D0 l","Cross","Comb.","uds", "ccbar","D0"};
  TString TeNames2[]={"D^{**} \\to \\text{D}^{0}\\pi^{0} & ","D^{**} \\to \\text{D}^{*0}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*0} \\,l \\,\\nu & ","D^{0} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","$uds$ & ",
		      "$c \\bar{c}$ & ","D^{0}"};
  TCut cuts2[]={"MCType>12&&MCDssmode==1","MCType>12&&MCDssmode==2","MCType>12&&MCDssmode>2",
	       "(MCType==2||MCType==4||MCType==6)","(MCType==1||MCType==3||MCType==5)","(MCType>6&&MCType<13)","MCType==0","",""};
  toTable(1, ch, MYields, TaNames2, TeNames2, cuts2, Table, tex, basic, folder, otherTag, doTruth, uds, ccbar);


  for(int i=0; i<9; i++){
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
  TString TaNames3[]={"D*+ pi0","D+ pi0","Oth D**","D*+ l","D+ l","Cross","Comb.","uds", "ccbar","D*+"};
  TString TeNames3[]={"D^{**} \\to \\text{D}^{*+}\\pi^{0} & ","D^{**} \\to \\text{D}^{+}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*+} \\,l \\,\\nu & ","D^{+} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","$uds$ & ",
		      "$c \\bar{c}$ & ","D^{*+}"};
  TCut cuts3[]={"MCType>12&&MCDssmode==4","MCType>12&&MCDssmode==3","MCType>12&&(MCDssmode<3||MCDssmode>4)",
	       "(MCType==8||MCType==10||MCType==12)","(MCType==7||MCType==9||MCType==11)","(MCType>0&&MCType<7)","MCType==0","",""};
  toTable(4, ch, MYields, TaNames3, TeNames3, cuts3, Table, tex, basic, folder, otherTag, doTruth, uds, ccbar);



  for(int i=0; i<9; i++){
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
  TString TaNames4[]={"D+ pi0","D*+ pi0","Oth D**","D*+ l","D+ l","Cross","Comb.","uds", "ccbar","D+"};
  TString TeNames4[]={"D^{**} \\to \\text{D}^{+}\\pi^{0} & ","D^{**} \\to \\text{D}^{*+}\\pi^{0} & ","Other D^{**} & ",
		     "D^{*+} \\,l \\,\\nu & ","D^{+} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","$uds$ & ",
		      "$c \\bar{c}$ & ","D^{+}"};
  TCut cuts4[]={"MCType>12&&MCDssmode==3","MCType>12&&MCDssmode==4","MCType>12&&(MCDssmode<3||MCDssmode>4)",
	       "(MCType==8||MCType==10||MCType==12)","(MCType==7||MCType==9||MCType==11)","(MCType>0&&MCType<7)","MCType==0","",""};
  toTable(3, ch, MYields, TaNames4, TeNames4, cuts4, Table, tex, basic, folder, otherTag, doTruth, uds, ccbar);

  gSystem->CopyFile(taFile, "mycode/DssTables/contable"+Run+sTruth+"_"+folder+".txt",kTRUE);
  cout<<taFile<<" written"<<endl;
}

//=============================================================================================================

void toTable(int cand, TChain &ch, int MYields[][2], TString *TaNames, TString *TeNames, 
	     TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString Run, 
	     bool doTruth, TChain &uds, TChain &ccbar){

  double nccbar[] = {58900000., 168844000., 83974000., 0., 366758000., 104778000.}; // Run 4 is 252830000
  double nuds[] = {47180000., 130858000., 66892000., 0., 317846000., 84414000}; // Run 4 is 213380000
  double MCBs[] = {36968000+37200000., 103498000+103124000., 50556000+49766000., 167994000+167466000., 244322000+244812000., 68148000+68016000.};
  double realBs = 0, totuds = 0, totccbar = 0;
  for(int i=0; i<6; i++) {
    totuds += nuds[i];
    totccbar += nccbar[i];
    TString irun = ""; irun += i+1;
    if(Run.Contains(irun)){
      realBs += MCBs[i]*2/3;
    }
  }
  double wuds = realBs/totuds*2.09/1.05; 
  double wccbar = realBs/totccbar*1.3/1.05; 
  cout<<"Weight uds: "<<wuds<<"; Weight ccbar: "<<wccbar<<endl;

  TCut eCut = "candIsMu==0";
  TCut muCut = "candIsMu==1";
  TCut truthMatch = "candLepTru==1";
  TString candCut = "candType=="; candCut += cand;
  Table<<" "<<TaNames[9]<<"\t  "<<folder<<"_e\t  Base_e   Ratio_e "<<"\t"<<folder;
  Table<<"_mu\tBase_mu   Ratio_mu "<<endl;
  Table<<"===================================================================="<<endl;
  tex<<"\\begin{tabular}{|l||cc|c|cc|c|}"<<endl<<"\\hline"<<endl;
  tex<<"$\\mathbf{"<<TeNames[9]<<"}$ & "<<folder<<"_e & Base_e & R_e & ";
  tex<<folder<<"_{\\mu} & Base_{\\mu} & R_{\\mu} \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  basic += candCut;
  TCut tCut = basic; tCut += eCut;
  int bkg_e = ch.GetEntries(tCut);
  tCut = basic; tCut += muCut;
  int bkg_mu = ch.GetEntries(tCut);
  for(int i=0; i<9; i++){
    if(i==7) tex<<"\\hline"<<endl;
    TCut finalCut = basic; if(doTruth && i<7)finalCut += truthMatch; finalCut += cuts[i];
    int entries=0;
    TCut feCut = finalCut; feCut += eCut;
    if (i < 6){ 
       entries = ch.GetEntries(feCut);
       bkg_e -= entries;
    } else if(i==6) entries = bkg_e;
    else if(i==7) {
      entries = wuds * uds.GetEntries(feCut) + 0.5;
    } else if(i==8) entries = wccbar*ccbar.GetEntries(feCut)+.5;
    Table<<TaNames[i]<<"\t   "<<entries<<"\t   "<<MYields[i][0];
    Table<<"\t    "<<round((double)entries,(double)MYields[i][0]);
    tex<<TeNames[i]<<entries<<" & "<<MYields[i][0];
    tex<<" & "<<round((double)entries,(double)MYields[i][0])<<" & ";
    TCut fmuCut = finalCut; fmuCut += muCut;
    if (i< 6){ 
      entries = ch.GetEntries(fmuCut);
      bkg_mu -= entries;
    } else if(i==6) entries = bkg_mu;
    else if(i==7) entries = wuds*uds.GetEntries(fmuCut)+.5;
    else if(i==8) entries = wccbar*ccbar.GetEntries(fmuCut)+.5;
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

