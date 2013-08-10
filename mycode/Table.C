#include "mycode/cuts.hh"

TString round(double n, double d);
void toTable(int cand, TChain &ch, int *DYields, int MYields[][2], TString *TaNames, TString *TeNames, 
	     TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag, bool doTruth);

void Table(TString folder, TString otherTag = "NN", TString Run="24", TString tag="AllM2", bool doTruth=true){

  TCut M2 = "candM2>1";
  if(tag.Contains("All"))
    basic += ee;
  else if(tag=="afterBestB")
    basic = "";
  else if(tag.Contains("Mva"))
    basic += Mva;
  if(tag.Contains("M2"))
     basic += M2;

  int nFiles = 0;
  TChain ch("ntp1");
  TString directory = "AWG82/ntuples/";
  directory += folder;
  TString filename;
  TString unifile = "AWG82/ntuples/Add_Final/MergedRest1/MergedRest123456.root";
  //TString unifile = "AWG82/ntuples/Add_BDT_KM/MergedRest1/MergedRest123456.root";
  cout<<"Using "<<unifile<<endl;
  ch.Add(unifile);
//   for(int i=1; i<7; i++){
//     TString irun = ""; irun += i;
//     if(Run.Contains(irun)){
//       TString dirs = "/MergedRest/*Run"; dirs += irun; dirs+="*";
//       filename=directory;
//       filename += dirs;
//       nFiles += ch.Add(filename);
//       cout<<"Added "<<nFiles<<" files from "<<filename<<" with entries: "<<ch.GetEntries()<<endl;   
//     }
//   }

  TString Mdir = "mycode/tables/"; Mdir += folder;
  gSystem->mkdir(Mdir);
  TString sTruth = "noTru"; if(doTruth) sTruth="";
  TString taFile = "mycode/tables/table"; taFile+=Run; taFile+=sTruth; taFile+="_";
  taFile+=folder; taFile+="_"; taFile+=otherTag; taFile+=tag; taFile+=".txt";
  fstream Table;
  Table.open(taFile,fstream::out);
  fstream tex;
  tex.open("mycode/tables/tex"+Run+sTruth+"_"+folder+"_"+otherTag+tag+".txt",fstream::out);
  fstream Other;
  Other.open("mycode/tables/table"+Run+sTruth+"_"+otherTag+".txt",fstream::in);
  //Other.open("mycode/tables/table_"+otherTag+".txt",fstream::in);
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
  TString TaNames[]={"D*0 Tau","D0 Tau","D*0 l","D0 l","D** l","Cross","Comb.","D*0"};
  TString TeNames[]={"D^{*0} \\,\\tau \\,\\nu & ","D^{0} \\,\\tau \\,\\nu & ","D^{*0} \\,l \\,\\nu & ",
		     "D^{0} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{*0}"};
  TCut cuts[]={"MCType==6","MCType==5","(MCType==2||MCType==4)","(MCType==1||MCType==3)",
	       "(MCType==13||MCType==14)","(MCType>6&&MCType<13)","MCType==0"};
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
  TString TaNames2[]={"D0 Tau","D*0 Tau","D*0 l","D0 l","D** l","Cross","Comb.","D0"};
  TString TeNames2[]={"D^{0} \\,\\tau \\,\\nu & ","D^{*0} \\,\\tau \\,\\nu & ","D^{*0} \\,l \\,\\nu & ",
		     "D^{0} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{0}"};
  TCut cuts2[]={"MCType==5","MCType==6","(MCType==2||MCType==4)","(MCType==1||MCType==3)",
	       "(MCType==13||MCType==14)","(MCType>6&&MCType<13)","MCType==0"};
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
  TString TaNames3[]={"D*+ Tau","D+ Tau","D*+ l","D+ l","D** l","Cross","Comb.","D*+"};
  TString TeNames3[]={"D^{*+} \\,\\tau \\,\\nu & ","D^{+} \\,\\tau \\,\\nu & ","D^{*+} \\,l \\,\\nu & ",
		     "D^{+} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{*+}"};
  TCut cuts3[]={"MCType==12","MCType==11","(MCType==8||MCType==10)","(MCType==7||MCType==9)",
	       "(MCType==13||MCType==14)","(MCType>0&&MCType<7)","MCType==0"};
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
  TString TaNames4[]={"D+ Tau","D*+ Tau","D*+ l","D+ l","D** l","Cross","Comb.","D+"};
  TString TeNames4[]={"D^{+} \\,\\tau \\,\\nu & ","D^{*+} \\,\\tau \\,\\nu & ","D^{*+} \\,l \\,\\nu & ",
		     "D^{+} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{+}"};
  TCut cuts4[]={"MCType==11","MCType==12","(MCType==8||MCType==10)","(MCType==7||MCType==9)",
	       "(MCType==13||MCType==14)","(MCType>0&&MCType<7)","MCType==0"};
  toTable(3, ch, DYields4, MYields, TaNames4, TeNames4, cuts4, Table, tex, basic, folder, otherTag, doTruth);

  gSystem->CopyFile("mycode/tables/table"+Run+sTruth+"_"+folder+"_"+otherTag+tag+".txt",
		   "mycode/tables/table"+Run+sTruth+"_"+folder+".txt",kTRUE);
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

