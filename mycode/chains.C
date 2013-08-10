{
  TCut MMiss = "candPMiss>.2";
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
  TCut basic;
  basic = MMiss+Q2+ee;

  TChain ch("ntp1");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run4-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run4-R22d-v06/*");
  TChain ch1("ntp1");
  ch1.Add("AWG82/ntuples/NNVeryLoose_Cleaner/SP-1235-BSemiExcl-Run2-R22d-v06/*");
  ch1.Add("AWG82/ntuples/NNVeryLoose_Cleaner/SP-1235-BSemiExcl-Run4-R22d-v06/*");
  ch1.Add("AWG82/ntuples/NNVeryLoose_Cleaner/SP-1237-BSemiExcl-Run2-R22d-v06/*");
  ch1.Add("AWG82/ntuples/NNVeryLoose_Cleaner/SP-1237-BSemiExcl-Run4-R22d-v06/*");
  ch.SetLineColor(2);
  TChain ch2("ntp1");
  ch2.Add("AWG82/ntuples/BDTVeryLoose/SP-1235-BSemiExcl-Run2-R22d-v06/*");
  ch2.Add("AWG82/ntuples/BDTVeryLoose/SP-1235-BSemiExcl-Run4-R22d-v06/*");
  ch2.Add("AWG82/ntuples/BDTVeryLoose/SP-1237-BSemiExcl-Run2-R22d-v06/*");
  ch2.Add("AWG82/ntuples/BDTVeryLoose/SP-1237-BSemiExcl-Run4-R22d-v06/*");
  ch2.SetLineColor(4);


}

