void plots(){

// stacking("candPMiss",0,4,"Missing momentum cut","p_{miss} (GeV)","PMiss",.2,"",true)
// stacking("candQ2",-4,14,"q^{2} cut","q^{2} (GeV^{2})","Q2",4,"left",true)
// stacking("candFitProb",0,.01,"Fit Probability cut","Fit Probability","FitProb",.001,"",true)
// stacking("candM2",-6,14,"m_{miss}^{2} cut","m_{miss}^{2} (GeV^{2})","Mmiss",-4,"",true)
// stacking("candPstarLep",0,2.7,"p_{lep}^{*} cut for D^{0,+}","p_{lep}^{*} (GeV)","PstarLepD",2.31,"right",true)
// stacking("candPstarLep",0,2.7,"p_{lep}^{*} cut for D*^{0,}*^{+}","p_{lep}^{*} (GeV)","PstarLepDstar",2.26,"right",true)
// stacking("candEExtra",0,3,"E_{Extra} cut (0.2 GeV) for the D^{0}","E_{Extra} (GeV)","EExtraD0",.2,"",true)
// stacking("candEExtra",0,3,"E_{Extra} cut (0.2 GeV) for the D*^{0}","E_{Extra} (GeV)","EExtraDstar0",.2,"",true)
// stacking("candEExtra",0,3,"E_{Extra} cut (0.15 GeV) for the D^{+}","E_{Extra} (GeV)","EExtraDplus",.15,"",true)
// stacking("candEExtra",0,3,"E_{Extra} cut (0.3 GeV) for the D*^{+}","E_{Extra} (GeV)","EExtraDstarplus",.3,"",true)



  // Loading files
  TChain ch("ntp1");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1235-BSemiExcl-Run4-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run2-R22d-v06/*");
  ch.Add("AWG82/ntuples/SP-1237-BSemiExcl-Run4-R22d-v06/*");

  TCanvas c("Yields","Yields",600,400);
  gStyle->SetOptStat(10);
  TLine *line= new TLine(0,0,0,0);
  line->SetLineColor(4);
  line->SetLineStyle(2);

  ch.Draw("candM2>>h(60,-6,14)");
  line->DrawLine(-4,0,-4,1.05*h->GetMaximum());
  line->DrawLine(12,0,12,1.05*h->GetMaximum());
  h->SetTitle("m_{miss}^{2} cut");
  h->SetXTitle("m_{miss}^{2} (GeV^{2})");
  c.SaveAs("mycode/Mmiss.eps");

  ch.Draw("candPstarLep>>h(60,0,2.6)","(candType==1||candType==3)");
  line->DrawLine(2.31,0,2.31,1.05*h->GetMaximum());
  h->SetTitle("p_{lep}^{*} cut for D^{0,+}");
  h->SetXTitle("p_{lep}^{*} (GeV)");
  c.SaveAs("mycode/PstarLep_D0.eps");

  ch.Draw("candPstarLep>>h(60,0,2.6)","(candType==2||candType==4)");
  line->DrawLine(2.26,0,2.26,1.05*h->GetMaximum());
  h->SetTitle("p_{lep}^{*} cut for D^{*0,*+}");
  h->SetXTitle("p_{lep}^{*} (GeV)");
  c.SaveAs("mycode/PstarLep_Dstar.eps");

  ch.Draw("candFitProb>>h(60,0,.01)");
  line->DrawLine(.001,0,.001,pow(h->GetMaximum(),1.05));
  h->SetTitle("Fit Probability cut");
  h->SetXTitle("Fit Probability");
  gPad->SetLogy(1);
  c.SaveAs("mycode/FitProb.eps");

  ch.Draw("candPMiss>>h(60,0,4)");
  line->DrawLine(.2,0,.2,1.05*h->GetMaximum());
  h->SetTitle("Missing momentum cut");
  h->SetXTitle("p_{miss} (GeV)");
  gPad->SetLogy(0);
  c.SaveAs("mycode/PMiss.eps");

  ch.Draw("candQ2>>h(60,-4,14)");
  line->DrawLine(4,0,4,1.05*h->GetMaximum());
  h->SetTitle("q^{2} cut");
  h->SetXTitle("q^{2} (GeV^{2})");
  gPad->SetLogy(0);
  c.SaveAs("mycode/Q2.eps");

  ch.Draw("candEExtra>>h(60,0,3)","candType==3");
  line->DrawLine(.15,0,.15,pow(h->GetMaximum(),1.05));
  h->SetTitle("E_{Extra} cut (0.15 GeV) for the D^{+}");
  h->SetXTitle("E_{Extra} (GeV)");
  gPad->SetLogy(1);
  c.SaveAs("mycode/EExtraDplus.eps");

  ch.Draw("candEExtra>>h(60,0,3)","candType==4");
  line->DrawLine(.3,0,.3,pow(h->GetMaximum(),1.05));
  h->SetTitle("E_{Extra} cut (0.3 GeV) for the D*^{+}");
  h->SetXTitle("E_{Extra} (GeV)");
  gPad->SetLogy(1);
  c.SaveAs("mycode/EExtraDstarplus.eps");

  ch.Draw("candEExtra>>h(60,0,3)","candType==1");
  line->DrawLine(.2,0,.2,pow(h->GetMaximum(),1.05));
  h->SetTitle("E_{Extra} cut (0.2 GeV) for the D^{0}");
  h->SetXTitle("E_{Extra} (GeV)");
  gPad->SetLogy(1);
  c.SaveAs("mycode/EExtraD0.eps");

  ch.Draw("candEExtra>>h(60,0,3)","candType==2");
  line->DrawLine(.2,0,.2,pow(h->GetMaximum(),1.05));
  h->SetTitle("E_{Extra} cut (0.2 GeV) for the D*^{0}");
  h->SetXTitle("E_{Extra} (GeV)");
  gPad->SetLogy(1);
  c.SaveAs("mycode/EExtraDstar0.eps");


}

