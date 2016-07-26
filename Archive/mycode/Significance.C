TString round(double d);

void Significance(TString tag=""){
  // Loading files
  TChain ch("ntp1");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run4-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run4-R22d-v06/*");

  int entries = ch.GetEntries();
  TCut MMiss = "candPMiss>.2";
  TCut Q2 = "candQ2>4";
  TCut ee1 = "candType<3&&candEExtra<.2";
  TCut ee2 = "candType==3&&candEExtra<.15";
  TCut ee3 = "candType==4&&candEExtra<.3";
  TCut ee = (ee1||ee2||ee3);
  TCut basic = ee+MMiss+Q2;

  TString D0names[]={"K^- \\pi^+","K^-\\pi^+\\pi^0","K^-\\pi^+\\pi^+\\pi^-",
		     "K_S^0\\,\\pi^+\\pi^-","K_S^0\\,\\pi^+\\pi^-\\pi^0",
		     "K_S^0\\,K_S^0","K_S^0\\,\\pi^0","\\pi^+\\pi^-","K^+K^-"};
  TString Dplusnames[]={"K^-\\pi^+\\pi^+","K^-\\pi^+\\pi^+\\pi^0","K_S^0\\,\\pi^+",
			"K_S^0\\,\\pi^+\\pi^0","K_S^0\\,\\pi^+\\pi^+\\pi^-",
			"K^+K^-\\pi^+","K_S^0\\,K^+"};
  fstream tex;
  tex.open("mycode/tables/Significance"+tag+".txt",fstream::out);
  tex<<"\\begin{tabular}{|l|ccc||l|ccc|}"<<endl<<"\\hline"<<endl;
  tex<<"$\\mathbf{D^{(*)0}}$ & $S$ & B & $\\frac{S}{S+B}$(\\%) & ";
  tex<<"$\\mathbf{D^{(*)+}}$ & $S$ & B & $\\frac{S}{S+B}$(\\%)\\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;

  double S=0,B=0;
  for(int i=1; i<10;i++){
    TString tS = "(MCType==5||MCType==6)&&candType<3&&candDType=="; tS+=i;
    TString tB = "(MCType==0||MCType>6)&&candType<3&&candDType==";tB+=i;
    TCut cutS(tS), cutB(tB);
    S = ch.GetEntries(basic+cutS);
    B = ch.GetEntries(basic+cutB);
    int puri = S/(S+B)*100;
    tex<<D0names[i-1]<<" & "<<S<<" & "<<B<<" & "<<puri<<" & ";
    if(i<8){
      tS = "(MCType==11||MCType==12)&&(candType==3||(candType==4&&candDstarType==5))&&candDType=="; 
      tB = "(MCType<7||MCType>12)&&(candType==3||(candType==4&&candDstarType==5))&&candDType==";
      tS+=i; tB+=i;
      TCut cutS(tS), cutB(tB);
      S = ch.GetEntries(basic+cutS);
      B = ch.GetEntries(basic+cutB);
      puri = S/(S+B)*100;
      tex<<Dplusnames[i-1]<<" & "<<S<<" & "<<B<<" & "<<puri<<" \\\\"<<endl;  
    } else {
      tex<<" & & & \\\\"<<endl;
    }
    cout<<"DMode "<<i<<" has S "<<S<<", S/sqrt(S+B) "<<round(S/sqrt(S+B));
    cout<<" and S/(S+B) "<<round(S/(S+B))<<endl;
  }
  tex<<"\\hline \\end{tabular}"<<endl<<endl;
}

TString round(double d){
  double b = (int)(d *10+.5);
  TString sb="";
  if(!((int)b%100)){
    b = b/100.;
    sb += b; sb += ".0";
  }else if(!((int)b%10)){
    b = b/10.;
    sb += b; sb += ".0";
  }else {
    sb+=b/10;
  }
  return sb;
}

