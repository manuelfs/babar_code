#include "mycode/cuts.hh"
TString round(double n, int e, double d=1.);

void expectedEvents(TString Dss=""){
  TString semName = "AWG82/ntuples/small/Add_Aug20Rest_Run123456.root";
  TChain sem("ntp1");
  sem.Add(semName);
  TCut selCut = MvaAll;
  TString texType = "\\bf{Signal }";
  if(Dss=="dss") {
    selCut=dssMvaAll;
    texType = "$\\boldsymbol{D^{(*)}\\pi^0}$";
  }

  TString row[]={"$D\\ell\\nu$","$D^*\\ell\\nu$","$D\\tau\\nu$","$D^*\\tau\\nu$",
		 "$D^{**}\\ell\\nu$","$D^{**}\\tau\\nu$","$D_1^{\'}(\\to D^*\\pi)\\ell\\nu$",
		 "$D_1(\\to D\\pi\\pi)\\ell\\nu$"};

  double MCBs[] = {36968000+37200000., 103498000+103124000., 50556000+49766000., 
		   167994000+167466000., 244322000+244812000., 68148000+68016000.};
  double dataBs[] = {22389980.4, 67394307.5, 35569248.8, 110449802.7, 147190396.5, 84767412.6};
  double nMCB = 0, ndataBs = 0;
  for(int i=1; i<7; i++){
    nMCB += MCBs[i-1]*2/3.; 
    ndataBs += dataBs[i-1];
  }
  TString cuts[] = {"(abs(MCD)==421||abs(MCD)==411)&&MCType<13&&MCTaumode==-1",
		    "(abs(MCD)==423||abs(MCD)==413)&&MCType<13&&MCTaumode==-1",
		    "(abs(MCD)==421||abs(MCD)==411)&&MCType<13&&MCTaumode>-1",
		    "(abs(MCD)==423||abs(MCD)==413)&&MCType<13&&MCTaumode>-1",
		    "MCType>12&&MCTaumode==-1","MCType>12&&MCTaumode>-1"};

  TString texName = "babar_code/Tables/texExpectedEvents"; texName+=Dss;texName+=".txt";
  fstream tex;
  tex.open(texName,fstream::out);
  tex<<"\\begin{tabular}{|l||c|c|c|c|}"<<endl<<"\\hline"<<endl;
  tex<<texType<<" & $D^0$ & $D^{*0}$ & $D^+$ & $D^{*+}$ \\\\"<<endl;
  tex<<"\\hline \\hline"<<endl;
  TCanvas can;
  for(int i=0; i<6; i++){
    if(i==2 || i==4) tex<<"\\hline"<<endl;
    tex<<row[i];
    cout<<row[i]<<":    \t";
    int ifix = i;
    for(int j=0; j<4; j++){ 
      TString cut="candType=="; cut+=j+1; cut+="&&"; cut+=cuts[i];
      if(j<2) cut+="&&isBzero==0";
      else cut+="&&isBzero==1";
      TCut tot = selCut; tot+=cut; if(i<2||i>3) tot *= "wBF";
      sem.Draw("candM2>>pdf(100,-4,12)",tot);
      TH1F *h = (TH1F*)gDirectory->Get("pdf");
      
      double entries = h->Integral()*ndataBs/nMCB;
      tex<<" & "<<round(entries,0);
      cout<<round(entries,0)<<"\t";
    }
    tex<<" \\\\"<<endl;
    cout<<endl;
  }
  tex<<"\\hline\\hline"<<endl;
  double BFe[] = {0.0151*0.00079*0.03, 0.0151*0.0002*0.03,0.0137*0.00063*0.012,0.0137*0.00006*0.012,
		  0.0151*0.0006*0.03, 0.0151*0.00043*0.03,0.0137*0.00039*0.012,0.0137*0.00044*0.012};
  double BFepi0[] = {0.0151*0.00084*0.03, 0.0151*0.00030*0.03,0.0137*0.00133*0.012,0.0137*0.00014*0.012,
		     0.0151*0.00064*0.03, 0.0151*0.00032*0.03,0.0137*0.00088*0.012,0.0137*0.00041*0.012};
  for(int i=0; i<2; i++){
    tex<<row[i+6];
    cout<<row[i+6]<<"\t";
    for(int j=0; j<4; j++){ 
      double bf = BFe[j+(i==1)*4];
      if(Dss="dss") bf = BFepi0[j+(i==1)*4];
      tex<<" & "<<round(ndataBs/2*bf,0);
      cout<<round(ndataBs/2*bf,0)<<"\t";
    }
    tex<<" \\\\"<<endl;
    cout<<endl;
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
  if(!result.Contains(".") && e!=0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}


