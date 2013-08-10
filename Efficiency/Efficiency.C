TString RoundNumber(double n, int e, double d=1);

void Efficiency(){
  double SigSample[3][3][4] = {{{25325.2, 44262.6, 5977.5, 23677.9},    // Ds: No weights
				{24474.7, 40175.5, 5764.6, 21498.8},    // Ds: Weights All
				{23771.9, 38907.0, 5551.6, 20802.6}},   // Ds: Weight BDT 
				{{15421.6, 1607.9, 8010.9, 93.8},       // D: No weights 
				 {14487.6, 1467.8, 7561.5, 85.1},       // D: Weights All
				 {14044.5, 1421.6, 7291.6, 82.6}},      // D: Weight BDT 
			       {{1532.4, 1646.9, 727.3, 683.0},         // Dss: No weights 
				{1460.7, 1561.2, 692.9, 656.4},         // Dss: Weights All
				{1416.5, 1482.7, 663.3, 619.4}}};       // Dss: Weight BDT 

  double Dpi0Sample[3][3][4] = {{{2463.5, 1821.1, 865.8, 1010.0},       // Ds: No weights
				 {2901.4, 2041.4, 1021.1, 1125.8},      // Ds: Weights All
				 {2906.5, 2069.9, 1022.5, 1146.4}},     // Ds: Weight BDT 
				{{832.7, 201.1, 417.5, 32.4},           // D: No weights 
				 {975.0, 229.3, 488.9, 36.1},           // D: Weights All
				 {976.9, 231.6, 489.2, 36.6}},          // D: Weight BDT 
				{{3683.8, 2021.6, 2142.0, 1044.9},      // Dss: No weights 
				 {4477.4, 2328.0, 2602.9, 1204.3},      // Dss: Weights All
				 {4475.7, 2343.8, 2602.6, 1213.3}}};    // Dss: Weight BDT 

  TString Component[3] = {"D*","D","D**"}; 
  TString NameW[3] = {"No weights","Weights All","Weights BDT"}; 
  double f[3][3][4], ef[3][3][4];

  for(int comp=0; comp<3; comp++){
    cout<<endl<<endl<<"Doing component "<<Component[comp]<<endl<<"=====================";
    for(int weight=0; weight<3; weight++){
      cout<<endl<<NameW[weight]<<":\t";
      for(int chan=0; chan<4; chan++){
	f[comp][weight][chan] = SigSample[comp][weight][chan]/Dpi0Sample[comp][weight][chan];
	ef[comp][weight][chan] = sqrt(SigSample[comp][weight][chan]/pow(Dpi0Sample[comp][weight][chan],2)+
				      pow(SigSample[comp][weight][chan],2)/pow(Dpi0Sample[comp][weight][chan],3));
	cout<<RoundNumber(f[comp][weight][chan],3)<<" +- "<<RoundNumber(ef[comp][weight][chan],3)<<", ";
      }
    }
  }

  cout<<endl<<endl<<"RATIOS";
  for(int comp=0; comp<3; comp++){
    cout<<endl<<endl<<"Doing component "<<Component[comp]<<endl<<"=====================";
    for(int weight=1; weight<3; weight++){
      cout<<endl<<NameW[weight]<<":\t";
      for(int chan=0; chan<4; chan++){
	cout<<RoundNumber(f[comp][weight][chan]*100,1,f[comp][0][chan])<<" %, ";
      }
    }
  }
  cout<<endl<<endl<<"BDT Bias";
  for(int comp=0; comp<3; comp++){
    cout<<endl<<endl<<"Doing component "<<Component[comp]<<endl<<"=====================";
    for(int weight=2; weight<3; weight++){
      cout<<endl<<NameW[weight]<<":\t";
      for(int chan=0; chan<4; chan++){
	cout<<RoundNumber(f[comp][weight][chan]*100,1,f[comp][1][chan])<<" %, ";
      }
    }
  }
  cout<<endl<<endl;
}

TString RoundNumber(double n, int e, double d){
  if(d==0) return " - ";
  double neg = 1; if(n*d<0) neg = -1;
  double b = (int)(neg*n/d*pow(10.,(double)e)+0.5);
  b /= pow(10.,(double)e)*neg;
  TString result; result+= b;
  result.ReplaceAll(" ","");
  if(!result.Contains(".") && e != 0) result += ".";
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<e-afterdot.Length(); i++)
    result += "0";
  return result;
}
