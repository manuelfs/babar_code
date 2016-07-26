TString round(double n, double d);
void toTable(int cand, TChain &ch, int *DYields, int *MYields, int *M22Yields, TString *TaNames, 
	     TString *TeNames, TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag);

void Table_David_Mazur(TString folder, TString otherTag = "NN"){
  // Cuts not yet done candExtraTracks==0&&candThetaLep>.4&&candThetaLep<2.6&&candM2>-4&&candM2<12&&candPstarLep>0&&candPstarLep<2.4
  TCut MMiss = "candPMiss>.2";
  TCut deltaE = "abs(candDeltaE)<.072";
  TCut Q2 = "candQ2>4";
  TCut fit = "candFitProb>-.001";
  TCut mes = "candMES>5.27";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut kin1 = "(candType==1||candType==3)&&candPstarLep<2.31";
  TCut kin2 = "(candType==2||candType==4)&&candPstarLep<2.26";
  TCut kin = (kin1||kin2);
  TCut cleanD0 = "(candType<3||(candType==4&&candDstarType==1)) && (candDType<5||candDType==7||candDType==9)";
  TCut cleanDplus = "(candType==3||(candType==4&&candDstarType==2)) && (candDType<4||candDType>5)";
  TCut clean = (cleanD0 || cleanDplus);
  TCut basic = MMiss+Q2+ee;
  TString tag="All";
  if(tag=="All")
    basic = MMiss+Q2+ee;
  else if(tag=="afterBestB")
    basic = fit;
  else if(tag=="EExtra")
    basic = ee;
  else if(tag=="MMiss")
    basic = MMiss;
  else if(tag=="Q2")
    basic = Q2;
  else if(tag=="Clean")
    basic = clean+MMiss+Q2+ee;

  // Loading files
  TString isAdd;
  int nfiles = 0;
  if(folder.Contains("Add"))
    isAdd="Add";
  else
    isAdd="noAdd";
  TChain ch("ntp1");
  if(isAdd=="noAdd"){
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1235-BSemiExcl-Run2-R22d-v08/*");
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1235-BSemiExcl-Run4-R22d-v08/*");
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1237-BSemiExcl-Run2-R22d-v08/*");
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1237-BSemiExcl-Run4-R22d-v08/*");
  } else {
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1235-BSemiExclAdd-Run2-R22f-v01/*");
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1235-BSemiExclAdd-Run4-R22f-v01/*");
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1237-BSemiExclAdd-Run2-R22f-v01/*");
    nfiles += ch.Add("AWG82/ntuples/"+folder+"/SP-1237-BSemiExclAdd-Run4-R22f-v01/*");
  }
  cout << "Added "<<nfiles<<" files"<<endl;
  fstream Table;
  Table.open("mycode/tables/tableMicroMD_"+folder+".txt",fstream::out);
  fstream tex;
  tex.open("mycode/tables/textableMicroMD_"+folder+".txt",fstream::out);

  otherTag = "M14";

  int DYields[] = {133,8,2848,99,104,69,0};
  int MYields[] = {131,16,2798,72,97,79};
  int M22Yields[] = {119,4,2271,69,71,64};
  TString TaNames[]={"D*0 Tau","D0 Tau","D*0 l","D0 l","D** l","Bkg","Cross","D*0"};
  TString TeNames[]={"D^{*0} \\,\\tau \\,\\nu & ","D^{0} \\,\\tau \\,\\nu & ","D^{*0} \\,l \\,\\nu & ",
		     "D^{0} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Bkg & ","Crossfeed & ","D^{*0}"};
  TCut cuts[]={"MCType==6","MCType==5","(MCType==2||MCType==4)","(MCType==1||MCType==3)",
	       "(MCType==13||MCType==14)","(MCType==0||(MCType>6&&MCType<13))","(MCType>6&&MCType<13)"};
  toTable(2, ch, DYields, MYields, M22Yields, TaNames, TeNames, cuts, Table, tex, basic, folder, otherTag);


  int DYields2[] = {110,155,3875,1186,346,163,0};
  int MYields2[] = {142,153,3010,1143,270,187};
  int M22Yields2[] = {75,101,2298,916,195,151};
  TString TaNames2[]={"D0 Tau","D*0 Tau","D*0 l","D0 l","D** l","Bkg","Cross","D0"};
  TString TeNames2[]={"D^{0} \\,\\tau \\,\\nu & ","D^{*0} \\,\\tau \\,\\nu & ","D^{*0} \\,l \\,\\nu & ",
		     "D^{0} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Bkg & ","Crossfeed & ","D^{0}"};
  TCut cuts2[]={"MCType==5","MCType==6","(MCType==2||MCType==4)","(MCType==1||MCType==3)",
	       "(MCType==13||MCType==14)","(MCType==0||(MCType>6&&MCType<13))","(MCType>6&&MCType<13)"};
  toTable(1, ch, DYields2, MYields2, M22Yields2, TaNames2, TeNames2, cuts2, Table, tex, basic, folder, otherTag);


  int DYields3[] = {87,1,1866,9,92,65,0};
  int MYields3[] = {92,4,1695,5,50,29};
  int M22Yields3[] = {54,2,1205,5,37,32};
  TString TaNames3[]={"D*+ Tau","D+ Tau","D*+ l","D+ l","D** l","Bkg","Cross","D*+"};
  TString TeNames3[]={"D^{*+} \\,\\tau \\,\\nu & ","D^{+} \\,\\tau \\,\\nu & ","D^{*+} \\,l \\,\\nu & ",
		     "D^{+} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Bkg & ","Crossfeed & ","D^{*+}"};
  TCut cuts3[]={"MCType==12","MCType==11","(MCType==8||MCType==10)","(MCType==7||MCType==9)",
	       "(MCType==13||MCType==14)","(MCType>=0&&MCType<7)","(MCType>0&&MCType<7)"};
  toTable(4, ch, DYields3, MYields3, M22Yields3, TaNames3, TeNames3, cuts3, Table, tex, basic, folder, otherTag);



  int DYields4[] = {46,14,448,401,93,63,0};
  int MYields4[] = {43,22,380,427,34,40};
  int M22Yields4[] = {29,15,357,296,62,46};
  TString TaNames4[]={"D+ Tau","D*+ Tau","D*+ l","D+ l","D** l","Bkg","Cross","D+"};
  TString TeNames4[]={"D^{+} \\,\\tau \\,\\nu & ","D^{*+} \\,\\tau \\,\\nu & ","D^{*+} \\,l \\,\\nu & ",
		     "D^{+} \\,l \\,\\nu & ","D^{**} \\,l \\,\\nu & ","Bkg & ","Crossfeed & ","D^{+}"};
  TCut cuts4[]={"MCType==11","MCType==12","(MCType==8||MCType==10)","(MCType==7||MCType==9)",
	       "(MCType==13||MCType==14)","(MCType>=0&&MCType<7)","(MCType>0&&MCType<7)"};
  toTable(3, ch, DYields4, MYields4, M22Yields4, TaNames4, TeNames4, cuts4, Table, tex, basic, folder, otherTag);


}

//====================================================================================================================

void toTable(int cand, TChain &ch, int *DYields, int *MYields, int *M22Yields, TString *TaNames, 
	     TString *TeNames, TCut *cuts, fstream &Table, fstream &tex, TCut basic, TString folder, TString otherTag){
  TString candCut = "candType=="; candCut += cand;
  Table<<" "<<TaNames[7]<<"\t  "<<folder<<"\t David\t  M/D\t "<<otherTag<<"\t  "<<folder<<"/"<<otherTag<<endl;
  Table<<"==================================================="<<endl;
  tex<<"\\begin{tabular}{|l||c|cc|cc|}"<<endl<<"\\hline"<<endl;
  tex<<"$\\mathbf{"<<TeNames[7]<<"}$ & R22_{M} & R22_{D} & M22/D22 & R14_{M} & M22/M14 \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  for(int i=0; i<6; i++){
    TCut finalCut = basic; finalCut+=candCut; 
    finalCut += cuts[i];
    int entries = M22Yields[i];//ch.Draw("MCType",finalCut);
    
    Table<<TaNames[i]<<"\t  "<<entries<<"\t "<<DYields[i]<<"\t  "<<round((double)entries,(double)DYields[i]);
    Table<<"\t   "<<MYields[i]<<"\t   "<<round((double)(entries)*711./535.,(double)(MYields[i]))<<endl;
    tex<<TeNames[i]<<entries<<" & "<<DYields[i]<<" & "<<round((double)(entries),(double)DYields[i]);
    tex<<" & "<<MYields[i]<<" & "<<round((double)(entries)*711./535.,(double)(MYields[i]))<<" \\\\"<<endl;
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

