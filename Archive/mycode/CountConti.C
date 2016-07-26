//---------------------------------------------------------------------------------
// Description:
//      Counts the uds and ccbar yield
//
// Author List:
//      Manuel Franco Sevilla                     Stanford University
//
// Revision History:
//      09/06/02 manuelf -- Creation
//---------------------------------------------------------------------------------

TString round(double n, int e, double d=1.);

void CountConti(TString tag="30") {

  double nccbar[] = {58900000., 168844000., 83974000., 0., 366758000., 104778000.}; // Run 4 is 252830000
  double nuds[] = {47180000., 130858000., 66892000., 0., 317846000., 84414000}; // Run 4 is 213380000
  double MCBs[] = {36968000+37200000., 103498000+103124000., 50556000+49766000., 167994000+167466000., 244322000+244812000., 68148000+68016000.};
  double realBs = 0, totuds = 0, totccbar = 0;
  for(int i=0; i<6; i++) {
    realBs += MCBs[i]*2/3;
    totuds += nuds[i];
    totccbar += nccbar[i];
  }
  double wuds = realBs/totuds*2.09/1.05; 
  double wccbar = realBs/totccbar*1.3/1.05; 
  cout<<"Weight uds: "<<wuds<<"; Weight ccbar: "<<wccbar<<endl;

  TChain uds("ntp1");
  TChain ccbar("ntp1");
  if(tag==""){
  uds.Add("AWG82/ntuples/uds"+tag+"/Merged/*Run1*");
  uds.Add("AWG82/ntuples/uds"+tag+"/Merged/*Run2*");
  uds.Add("AWG82/ntuples/uds"+tag+"/Merged/*Run3*");
  uds.Add("AWG82/ntuples/uds"+tag+"/Merged/*Run5*");
  uds.Add("AWG82/ntuples/uds"+tag+"/Merged/*Run6*");
  ccbar.Add("AWG82/ntuples/ccbar"+tag+"/Merged/*Run1*");
  ccbar.Add("AWG82/ntuples/ccbar"+tag+"/Merged/*Run2*");
  ccbar.Add("AWG82/ntuples/ccbar"+tag+"/Merged/*Run3*");
  ccbar.Add("AWG82/ntuples/ccbar"+tag+"/Merged/*Run5*");
  ccbar.Add("AWG82/ntuples/ccbar"+tag+"/Merged/*Run6*");
  }else{
   uds.Add("AWG82/ntuples/small/uds30_Run12356.root");
   ccbar.Add("AWG82/ntuples/small/ccbar30_Run12356.root");
  }

//   TChain udsDss("ntp1");
//   TChain ccbarDss("ntp1");
//   udsDss.Add("AWG82/ntuples/Dss/uds"+tag+"_Run126.root");
//   ccbarDss.Add("AWG82/ntuples/Dss/ccbar"+tag+"_Run126.root");

  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut Mva = "(candMvaDl>0.18&&candType==1)||(candMvaDl>0.22&&candType==2)||(candMvaDl>0.45&&candType==3)||(candMvaDl>0.19&&candType==4)";
//   TCut sigee = MMiss+Q2+ee+"candM2>1";
//   TCut sigMva = MMiss+Q2+Mva+"candM2>1";
  TCut sigee = MMiss+Q2+ee;
  TCut sigMva = MMiss+Q2+Mva;

  TCut cosT = "abs(candCosT)<.9";
  TCut dssacc = "mm2pi0>-4&&mm2pi0<12&&candPstarLep>0&&candPstarLep<2.4&&candMES>5.2&&candMES<5.3&&pmisspi0>.2&&bestepi0==1";
  TCut mpi0 = "mpi0>.125&&mpi0<.145";
  TCut dssee = "eextrapi0<.5&&ppi0>.4";
  TCut Mva = "(candMvaDssComb>-0.24&&candMvaDssDl>-0.2&&candType==1)||(candMvaDssComb>-0.09&&candMvaDssDl>-0.09&&candType==2)||(candMvaDssComb>-0.13&&candMvaDssDl>-0.24&&candType==3)||(candMvaDssComb>-0.09&&candMvaDssDl>-0.2&&candType==4)";
  TCut Mva2 = "(candMvaDssComb>0.18&&candMvaDssDl>0.1)";
  TCut dssee = mpi0+dssee+dssacc;
  TCut dssMva = Mva2+dssacc;

  for(int i=1; i<5; i++) {
    TString scut = "candType=="; scut+=i;
    TCut cut = scut;
    double uMva = wuds*(double)uds.GetEntries(sigMva+cut);
    double uee = wuds*(double)uds.GetEntries(sigee+cut);
    double cMva = wccbar*(double)ccbar.GetEntries(sigMva+cut);
    double cee = wccbar*(double)ccbar.GetEntries(sigee+cut);

    double udss = wuds*(double)uds.GetEntries(dssee+cut);
    double cdss = wccbar*(double)ccbar.GetEntries(dssee+cut);
    double uMvadss = wuds*(double)uds.GetEntries(dssMva+cut);
    double cMvadss = wccbar*(double)ccbar.GetEntries(dssMva+cut);
    cout<<i<<"\tMva \tEExtra \tRatio \tMvaDss \teeDss"<<endl;
    cout<<"uds\t"<< round(uMva,0) <<"\t"<< round(uee,0) <<"\t"<< round(uMva,2,uee);
    cout<<"\t"<< round(uMvadss,0) <<"\t"<< round(udss,0) <<"\t"<<  round(uMvadss,2,udss)<<endl;
    cout<<"ccbar\t"<< round(cMva,0) <<"\t"<< round(cee,0) <<"\t"<<  round(cMva,2,cee);
    cout<<"\t"<< round(cMvadss,0) <<"\t"<< round(cdss,0) <<"\t"<<  round(cMvadss,2,cdss)<<endl<<endl;
  }
  cout<<"======================================"<<endl<<endl;
  for(int i=1; i<3; i++) {
    TCut cut = "candType<3";
    if(i==2) cut = "candType>2";
    double uMva = wuds*(double)uds.GetEntries(sigMva+cut);
    double uee = wuds*(double)uds.GetEntries(sigee+cut);
    double cMva = wccbar*(double)ccbar.GetEntries(sigMva+cut);
    double cee = wccbar*(double)ccbar.GetEntries(sigee+cut);

    double udss = wuds*(double)uds.GetEntries(dssee+cut);
    double cdss = wccbar*(double)ccbar.GetEntries(dssee+cut);
    double uMvadss = wuds*(double)uds.GetEntries(dssMva+cut);
    double cMvadss = wccbar*(double)ccbar.GetEntries(dssMva+cut);
    cout<<i<<"\tMva \tEExtra \tRatio \tMvaDss \teeDss"<<endl;
    cout<<"uds\t"<< round(uMva,0) <<"\t"<< round(uee,0) <<"\t"<< round(uMva,2,uee);
    cout<<"\t"<< round(uMvadss,0) <<"\t"<< round(udss,0) <<"\t"<<  round(uMvadss,2,udss)<<endl;
    cout<<"ccbar\t"<< round(cMva,0) <<"\t"<< round(cee,0) <<"\t"<<  round(cMva,2,cee);
    cout<<"\t"<< round(cMvadss,0) <<"\t"<< round(cdss,0) <<"\t"<<  round(cMvadss,2,cdss)<<endl<<endl;
  }
}

TString round(double n, int e, double d){
  if(d==0) return " - ";
  double b = (int)(n/d*pow(10,e)+0.5);
  b /= pow(10,e);
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(b<10 && e==2){
    if(result.Length() == 1)
      result += ".00";
    if(result.Length() == 3)
      result += "0";
  }
  return result;
}

