TString Bmode(int theid);
TString Dmode(int theid);

void tagNames(){
  ofstream tag("tagNames.txt");
  tag << "Seed       Decay\n================================\n";

  for(int i=0; i<87; i++){
    if(i==2)  {i = 10; tag<<endl;} if(i==19) {i = 20; tag<<endl;} if(i==27) {i = 30; tag<<endl;}
    if(i==39) tag<<endl;
    if(i==49) {i = 50; tag<<endl;} if(i==59) {i = 60; tag<<endl;} if(i==62) {i = 71; tag<<endl;}
    if(i==73) {i = 80; tag<<endl;} if(i==85) {i = 86;}
    tag<<" "<<100+i<<"     "<<Dmode((100+i)*100)<<endl;
  }

}

TString Bmode(int theid){
  
  int bmode = theid%100;  

  TString base = "B->JPsi";
  if(theid>11000) base = "B->D0";
  if(theid>12000) base = "B->Dc";
  if(theid>13000) base = "B->Dc*";
  if(theid>14000) base = "B->D0*";
  if(theid>16000) base = "B->Ds";
  if(theid>17000) base = "B->Ds*";
  if(theid>18000) base = "B->Dc*";

  TString mode[]={"pi",//1
		   "K",
		   "pipi0_<1.5GeV",
		   "Kpi0_<1.5GeV",
		   "piKs",
		   "KKs",
		   "pi2pi0_<1.5GeV",
		   "K2pi0_<1.5GeV",
		   "3pi_<1.5GeV",
		   "K2pi_<1.5GeV",//10
		   "2Kpi_Ds",
		   "omegah",
		   "K2pipi0_<2.2GeV",
		   "2Kpipi0_Ds*",
		   "pipi0Ks",//15
		   "Kpi0Ks_<1.8GeV",
		   "K2pi0Ks_1.8-2.2GeV",
		   "2KsX",
		   "3pi2pi0_<2.2GeV",
		   "K2pi2pi0_<2.2GeV",//20
		   "2Kpi2pi0_Ds*",
		   "5pi_<2.3GeV",
		   "K4p_<2.7GeV",
		   "2K3pi_<2.7GeV",
		   "5pipi0_<2.2GeV",
		   "K4pipi0_<2.2GeV",
		   "2K3pipi0_<2.5GeV",
		   "3piKs_D*",
		   "3piKspi0_D*",
		   "K2piKs_D*",//30
		   "D*_Dpi0",
		   "pipi0_>1.5GeV",
		   "Kpi0_>1.5GeV",
		   "pi2pi0_1.5-2GeV",
		   "K2pi0_>1.5GeV",
		   "3pi_1.5-2GeV",
		   "K2pi_>1.5GeV",
		   "2Kpi_K*",
		   "2Kpi_other",
		   "3pipi0_<1.6GeV",//40
		   "3pipi0_1.6-2.2GeV",
		   "K2pipi0_>2.2GeV",
		   "2Kpipi0_other",
		   "Kpi0Ks_>1.8GeV",
		   "3pi2pi0_>2.2GeV",
		   "K2pi2pi0_>2.2GeV",
		   "2Kpi2pi0_other",
		   "5pi_>2.3GeV",
		   "K4p_>2.7GeV",
		   "2K3pi_>2.7GeV",//50
		   "5pipi0_>2.2GeV",
		   "3piKs_noD*",
		   "3piKspi0_noD*",
		   "Kpi_Reson",
		   "2pi_Reson",
		   "2K",
		   "pi0",
		   "Kpi_Bad",
                   "2pi_Bad"};//59


  base += ", ";
  base += mode[bmode-1];

  return base;
}

TString Dmode(int theid){
  
  int dmode = theid/100;
  TString modename;

  if(dmode == 100)   modename = "JPsi->ee";
  if(dmode == 101)   modename = "JPsi->mumu";

  if(dmode == 110)   modename = "D0->Kpi";
  if(dmode == 111)   modename = "D0->Kpipi0";
  if(dmode == 112)   modename = "D0->K3pi";
  if(dmode == 113)   modename = "D0->Kspipi";
  if(dmode == 114)   modename = "D0->Kspipipi0";
  if(dmode == 115)   modename = "D0->KK";
  if(dmode == 116)   modename = "D0->pipipi0";
  if(dmode == 117)   modename = "D0->pipi";
  if(dmode == 118)   modename = "D0->Kspi0";

  if(dmode == 130)   modename = "D*,D0->Kpi";
  if(dmode == 131)   modename = "D*,D0->Kpipi0";
  if(dmode == 132)   modename = "D*,D0->K3pi";
  if(dmode == 133)   modename = "D*,D0->Kspipi";
  if(dmode == 134)   modename = "D*,Dc->Kspipi0";
  if(dmode == 135)   modename = "D*,Dc->Kpipi";
  if(dmode == 136)   modename = "D*,Dc->Kpipipi0";
  if(dmode == 137)   modename = "D*,Dc->Kspipi0";
  if(dmode == 138)   modename = "D*,Dc->Kspipipi";
  if(dmode == 139)   modename = "D*,Dc->KKpi";

  if(dmode == 120)   modename = "Dc->Kspi";
  if(dmode == 121)   modename = "Dc->Kpipi";
  if(dmode == 122)   modename = "Dc->Kspipi0";
  if(dmode == 123)   modename = "Dc->Kpipipi0";
  if(dmode == 124)   modename = "Dc->Kspipipi";
  if(dmode == 125)   modename = "Dc->KKpi";
  if(dmode == 126)   modename = "Dc->KKpipi0";

  if(dmode == 140)   modename = "D*0->D0pi0,D0->Kpi";
  if(dmode == 141)   modename = "D*0->D0pi0,D0->Kpipi0";
  if(dmode == 142)   modename = "D*0->D0pi0,D0->K3pi";
  if(dmode == 143)   modename = "D*0->D0pi0,D0->Kspipi";
  if(dmode == 144)   modename = "D*0->D0pi0,D0->Kspipipi0";
  if(dmode == 145)   modename = "D*0->D0pi0,D0->KK";
  if(dmode == 146)   modename = "D*0->D0pi0,D0->pipipi0";
  if(dmode == 147)   modename = "D*0->D0pi0,D0->pipi";
  if(dmode == 148)   modename = "D*0->D0pi0,D0->Kspi0";

  if(dmode == 150)   modename = "D*0->D0gamma,D0->Kpi";
  if(dmode == 151)   modename = "D*0->D0gamma,D0->Kpipi0";
  if(dmode == 152)   modename = "D*0->D0gamma,D0->K3pi";
  if(dmode == 153)   modename = "D*0->D0gamma,D0->Kspipi";
  if(dmode == 154)   modename = "D*0->D0gamma,D0->Kspipipi0";
  if(dmode == 155)   modename = "D*0->D0gamma,D0->KK";
  if(dmode == 156)   modename = "D*0->D0gamma,D0->pipipi0";
  if(dmode == 157)   modename = "D*0->D0gamma,D0->pipi";
  if(dmode == 158)   modename = "D*0->D0gamma,D0->Kspi0";

  if(dmode == 160)   modename = "Ds->phipi";
  if(dmode == 161)   modename = "Ds->KsK";
  if(dmode == 171)   modename = "Ds*->Dsgamma,Ds->phipi";
  if(dmode == 172)   modename = "Ds*->Dsgamma,Ds->KsK";

  if(dmode == 180)   modename = "D*->D0pi,D0->Kspipipi0";
  if(dmode == 181)   modename = "D*->D0pi,D0->KK";
  if(dmode == 182)   modename = "D*->D0pi,D0->pipipi0";
  if(dmode == 183)   modename = "D*->D0pi,D0->pipi";
  if(dmode == 184)   modename = "D*->D0pi,D0->Kspi0";
  if(dmode == 186)   modename = "D*,Dc->KKpipi0";


  return modename;
}
