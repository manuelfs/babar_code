#include "mycode/cuts.hh"

TString round(double n, double d);
void toTable(int cand, TChain &ch, int *DYields, int *MYields, TString *TaNames, 
	     TString *TeNames, TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag);

void Table_Mazur(TString folder, TString Run="123456", TString tag="All"){

  TCut M2 = "candM2>1";
  if(tag.Contains("All"))
    basic += ee;
  else if(tag=="afterBestB")
    basic = "";
  else if(tag.Contains("Mva"))
    basic += Mva;
  if(tag.Contains("M2"))
     basic += M2;

  // Loading files
  int nFiles = 0;
  TChain ch("ntp1");
  TString directory = "AWG82/ntuples/";
  directory += folder;
  TString filename = "AWG82/ntuples/small/Add_Truth30Rest_Run123456.root";
  cout<<"Using "<<filename<<endl;
  ch.Add(filename);
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
  TString taFile = "mycode/Mazurtables/tableM_"; taFile += folder; taFile += tag; taFile += ".txt";
  fstream Table;
  Table.open(taFile,fstream::out);
  fstream tex;
  tex.open("mycode/Mazurtables/texM_"+folder+tag+".txt",fstream::out);

  TString otherTag = "M14";

  int DYields[] = {133,8,2848,99,104,69,0};
  int MYields[] = {135,16,2798,72,97,35,44};
  TString TaNames[]={"D*0 Tau","D0 Tau","D*0 l","D0 l","D** l","Cross","Comb.","D*0"};
  TString TeNames[]={"D^{*0} \\,\\tau \\,\\nu & ","D^{0} \\,\\tau \\,\\nu & ","D^{*0} \\,l \\,\\nu & ",
		     "D^{0} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{*0}"};
  TCut cuts[]={"MCType==6","MCType==5","(MCType==2||MCType==4)","(MCType==1||MCType==3)",
	       "(MCType==13||MCType==14)","((MCType>6&&MCType<13))","(MCType==0)"};
  toTable(2, ch, DYields, MYields, TaNames, TeNames, cuts, Table, tex, basic, folder, otherTag);


  int DYields2[] = {110,155,3875,1186,346,163,0};
  int MYields2[] = {142,153,3010,1143,270,93,94};
  TString TaNames2[]={"D0 Tau","D*0 Tau","D*0 l","D0 l","D** l","Cross","Comb.","D0"};
  TString TeNames2[]={"D^{0} \\,\\tau \\,\\nu & ","D^{*0} \\,\\tau \\,\\nu & ","D^{*0} \\,l \\,\\nu & ",
		     "D^{0} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{0}"};
  TCut cuts2[]={"MCType==5","MCType==6","(MCType==2||MCType==4)","(MCType==1||MCType==3)",
	       "(MCType==13||MCType==14)","((MCType>6&&MCType<13))","(MCType==0)"};
  toTable(1, ch, DYields2, MYields2, TaNames2, TeNames2, cuts2, Table, tex, basic, folder, otherTag);


  int DYields3[] = {87,1,1866,9,92,65,0};
  int MYields3[] = {92,4,1695,5,50,2,27};
  TString TaNames3[]={"D*+ Tau","D+ Tau","D*+ l","D+ l","D** l","Cross","Comb.","D*+"};
  TString TeNames3[]={"D^{*+} \\,\\tau \\,\\nu & ","D^{+} \\,\\tau \\,\\nu & ","D^{*+} \\,l \\,\\nu & ",
		     "D^{+} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{*+}"};
  TCut cuts3[]={"MCType==12","MCType==11","(MCType==8||MCType==10)","(MCType==7||MCType==9)",
	       "(MCType==13||MCType==14)","(MCType>0&&MCType<7)","(MCType==0)"};
  toTable(4, ch, DYields3, MYields3, TaNames3, TeNames3, cuts3, Table, tex, basic, folder, otherTag);



  int DYields4[] = {46,14,448,401,93,63,0};
  int MYields4[] = {43,22,380,427,34,9,31};
  TString TaNames4[]={"D+ Tau","D*+ Tau","D*+ l","D+ l","D** l","Cross","Comb.","D+"};
  TString TeNames4[]={"D^{+} \\,\\tau \\,\\nu & ","D^{*+} \\,\\tau \\,\\nu & ","D^{*+} \\,l \\,\\nu & ",
		     "D^{+} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Crossfeed & ","Comb. & ","D^{+}"};
  TCut cuts4[]={"MCType==11","MCType==12","(MCType==8||MCType==10)","(MCType==7||MCType==9)",
	       "(MCType==13||MCType==14)","(MCType>0&&MCType<7)","(MCType==0)"};
  toTable(3, ch, DYields4, MYields4, TaNames4, TeNames4, cuts4, Table, tex, basic, folder, otherTag);

  cout<<taFile<<" written"<<endl;

}

//====================================================================================================================

void toTable(int cand, TChain &ch, int *DYields, int *MYields, TString *TaNames, 
	     TString *TeNames, TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag){

  double nMazur = 690.1;
  double nGen = 894.6; // 2/3 of runs 1-6 BSemiExclAdd
  //double nGen = 716.6; // runs 1-4
  cout<<"Doing "<<TaNames[7]<<" mode"<<endl;
  TString candCut = "candType=="; candCut += cand;
  Table<<" "<<TaNames[7]<<"\t  "<<folder<<"\t "<<otherTag<<"\t  Ratio"<<endl;
  Table<<"============================================"<<endl;
  tex<<"\\begin{tabular}{|l||cc|c|}"<<endl<<"\\hline"<<endl;
  tex<<"$\\mathbf{"<<TeNames[7]<<"}$ & New & M14 & Ratio \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  for(int i=0; i<7; i++){
    TCut finalCut = basic; finalCut+=candCut; 
    finalCut += cuts[i];
    int entries = ch.GetEntries(finalCut);
    
    Table<<TaNames[i]<<"\t  "<<entries;
    Table<<"\t   "<<MYields[i]<<"\t   "<<round((double)(entries)*nMazur/nGen,(double)(MYields[i]))<<endl;
    tex<<TeNames[i]<<entries;
    tex<<" & "<<MYields[i]<<" & "<<round((double)(entries)*nMazur/nGen,(double)(MYields[i]))<<" \\\\"<<endl;
  }
  Table<<endl;
  tex<<"\\hline \\end{tabular}"<<endl;
  tex<<"\\,\\,"<<endl<<endl;
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

