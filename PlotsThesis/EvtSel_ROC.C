void EvtSel_ROC(){
   TCanvas *c = new TCanvas("c", "TMVA ROC",1000,520);
   gStyle->SetOptStat(0);
   c->Range(-0.128266,0.0769231,1.05938,1.10256);
   c->SetFillColor(170);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetGridx();
   c->SetGridy();
   c->SetTickx();
   c->SetTicky();
   c->SetRightMargin(0.05);
   c->SetBottomMargin(0.12);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#fffffd");
   c->SetFrameFillColor(ci);
   c->SetFrameBorderMode(0);

   ci = TColor::GetColor("#fffffd");
   c->SetFrameFillColor(ci);
   c->SetFrameBorderMode(0);
   
   TH1 *frame = new TH2F("frame","Background rejection versus Signal efficiency",500,0,1,500,0.2,1);
   frame->SetStats(0);
   frame->SetLineWidth(2);
   frame->SetMarkerStyle(21);
   frame->SetMarkerSize(0.3);
   frame->GetXaxis()->SetTitle("Signal efficiency");
   frame->GetXaxis()->SetLabelOffset(0.012);
   frame->GetXaxis()->SetTitleSize(0.045);
   frame->GetXaxis()->SetTitleOffset(1.25);
   frame->GetYaxis()->SetTitle("Background rejection");
   frame->GetYaxis()->SetLabelOffset(0.012);
   frame->GetYaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleOffset(1.22);
   frame->Draw("");
   
   TH1 *MVA_BDT_rejBvsS = new TH1F("MVA_BDT_rejBvsS","MVA_BDT",100,0,1);
   MVA_BDT_rejBvsS->SetBinContent(1,0.999757);
   MVA_BDT_rejBvsS->SetBinContent(2,0.998992);
   MVA_BDT_rejBvsS->SetBinContent(3,0.99796);
   MVA_BDT_rejBvsS->SetBinContent(4,0.9967);
   MVA_BDT_rejBvsS->SetBinContent(5,0.996205);
   MVA_BDT_rejBvsS->SetBinContent(6,0.995513);
   MVA_BDT_rejBvsS->SetBinContent(7,0.994214);
   MVA_BDT_rejBvsS->SetBinContent(8,0.993197);
   MVA_BDT_rejBvsS->SetBinContent(9,0.992519);
   MVA_BDT_rejBvsS->SetBinContent(10,0.9919);
   MVA_BDT_rejBvsS->SetBinContent(11,0.991039);
   MVA_BDT_rejBvsS->SetBinContent(12,0.990365);
   MVA_BDT_rejBvsS->SetBinContent(13,0.98932);
   MVA_BDT_rejBvsS->SetBinContent(14,0.987882);
   MVA_BDT_rejBvsS->SetBinContent(15,0.987512);
   MVA_BDT_rejBvsS->SetBinContent(16,0.986445);
   MVA_BDT_rejBvsS->SetBinContent(17,0.984871);
   MVA_BDT_rejBvsS->SetBinContent(18,0.983454);
   MVA_BDT_rejBvsS->SetBinContent(19,0.982267);
   MVA_BDT_rejBvsS->SetBinContent(20,0.979838);
   MVA_BDT_rejBvsS->SetBinContent(21,0.978194);
   MVA_BDT_rejBvsS->SetBinContent(22,0.975945);
   MVA_BDT_rejBvsS->SetBinContent(23,0.97426);
   MVA_BDT_rejBvsS->SetBinContent(24,0.972234);
   MVA_BDT_rejBvsS->SetBinContent(25,0.969192);
   MVA_BDT_rejBvsS->SetBinContent(26,0.967101);
   MVA_BDT_rejBvsS->SetBinContent(27,0.96483);
   MVA_BDT_rejBvsS->SetBinContent(28,0.961545);
   MVA_BDT_rejBvsS->SetBinContent(29,0.959971);
   MVA_BDT_rejBvsS->SetBinContent(30,0.958524);
   MVA_BDT_rejBvsS->SetBinContent(31,0.955725);
   MVA_BDT_rejBvsS->SetBinContent(32,0.953177);
   MVA_BDT_rejBvsS->SetBinContent(33,0.951323);
   MVA_BDT_rejBvsS->SetBinContent(34,0.948241);
   MVA_BDT_rejBvsS->SetBinContent(35,0.945286);
   MVA_BDT_rejBvsS->SetBinContent(36,0.942486);
   MVA_BDT_rejBvsS->SetBinContent(37,0.938514);
   MVA_BDT_rejBvsS->SetBinContent(38,0.935591);
   MVA_BDT_rejBvsS->SetBinContent(39,0.931696);
   MVA_BDT_rejBvsS->SetBinContent(40,0.928399);
   MVA_BDT_rejBvsS->SetBinContent(41,0.923816);
   MVA_BDT_rejBvsS->SetBinContent(42,0.92019);
   MVA_BDT_rejBvsS->SetBinContent(43,0.917445);
   MVA_BDT_rejBvsS->SetBinContent(44,0.914074);
   MVA_BDT_rejBvsS->SetBinContent(45,0.909883);
   MVA_BDT_rejBvsS->SetBinContent(46,0.90615);
   MVA_BDT_rejBvsS->SetBinContent(47,0.902021);
   MVA_BDT_rejBvsS->SetBinContent(48,0.897983);
   MVA_BDT_rejBvsS->SetBinContent(49,0.893959);
   MVA_BDT_rejBvsS->SetBinContent(50,0.88748);
   MVA_BDT_rejBvsS->SetBinContent(51,0.883318);
   MVA_BDT_rejBvsS->SetBinContent(52,0.878563);
   MVA_BDT_rejBvsS->SetBinContent(53,0.873076);
   MVA_BDT_rejBvsS->SetBinContent(54,0.867647);
   MVA_BDT_rejBvsS->SetBinContent(55,0.860999);
   MVA_BDT_rejBvsS->SetBinContent(56,0.856797);
   MVA_BDT_rejBvsS->SetBinContent(57,0.850621);
   MVA_BDT_rejBvsS->SetBinContent(58,0.844718);
   MVA_BDT_rejBvsS->SetBinContent(59,0.83838);
   MVA_BDT_rejBvsS->SetBinContent(60,0.832226);
   MVA_BDT_rejBvsS->SetBinContent(61,0.824744);
   MVA_BDT_rejBvsS->SetBinContent(62,0.820169);
   MVA_BDT_rejBvsS->SetBinContent(63,0.813739);
   MVA_BDT_rejBvsS->SetBinContent(64,0.806751);
   MVA_BDT_rejBvsS->SetBinContent(65,0.801329);
   MVA_BDT_rejBvsS->SetBinContent(66,0.792065);
   MVA_BDT_rejBvsS->SetBinContent(67,0.781413);
   MVA_BDT_rejBvsS->SetBinContent(68,0.773035);
   MVA_BDT_rejBvsS->SetBinContent(69,0.763219);
   MVA_BDT_rejBvsS->SetBinContent(70,0.751574);
   MVA_BDT_rejBvsS->SetBinContent(71,0.741596);
   MVA_BDT_rejBvsS->SetBinContent(72,0.7271);
   MVA_BDT_rejBvsS->SetBinContent(73,0.717824);
   MVA_BDT_rejBvsS->SetBinContent(74,0.704425);
   MVA_BDT_rejBvsS->SetBinContent(75,0.696524);
   MVA_BDT_rejBvsS->SetBinContent(76,0.682949);
   MVA_BDT_rejBvsS->SetBinContent(77,0.671874);
   MVA_BDT_rejBvsS->SetBinContent(78,0.658909);
   MVA_BDT_rejBvsS->SetBinContent(79,0.644928);
   MVA_BDT_rejBvsS->SetBinContent(80,0.627895);
   MVA_BDT_rejBvsS->SetBinContent(81,0.61611);
   MVA_BDT_rejBvsS->SetBinContent(82,0.597258);
   MVA_BDT_rejBvsS->SetBinContent(83,0.583402);
   MVA_BDT_rejBvsS->SetBinContent(84,0.569894);
   MVA_BDT_rejBvsS->SetBinContent(85,0.554513);
   MVA_BDT_rejBvsS->SetBinContent(86,0.537413);
   MVA_BDT_rejBvsS->SetBinContent(87,0.512538);
   MVA_BDT_rejBvsS->SetBinContent(88,0.490974);
   MVA_BDT_rejBvsS->SetBinContent(89,0.468371);
   MVA_BDT_rejBvsS->SetBinContent(90,0.442445);
   MVA_BDT_rejBvsS->SetBinContent(91,0.406539);
   MVA_BDT_rejBvsS->SetBinContent(92,0.374898);
   MVA_BDT_rejBvsS->SetBinContent(93,0.352486);
   MVA_BDT_rejBvsS->SetBinContent(94,0.309766);
   MVA_BDT_rejBvsS->SetBinContent(95,0.271402);
   MVA_BDT_rejBvsS->SetBinContent(96,0.237676);
   MVA_BDT_rejBvsS->SetBinContent(97,0.195688);
   MVA_BDT_rejBvsS->SetBinContent(98,0.140793);
   MVA_BDT_rejBvsS->SetBinContent(99,0.0906321);
   MVA_BDT_rejBvsS->SetBinContent(100,0.0457354);
   MVA_BDT_rejBvsS->SetEntries(100);
   MVA_BDT_rejBvsS->SetStats(0);
   MVA_BDT_rejBvsS->SetLineWidth(3);
   MVA_BDT_rejBvsS->SetMarkerStyle(21);
   MVA_BDT_rejBvsS->SetMarkerSize(0.3);
   MVA_BDT_rejBvsS->GetXaxis()->SetTitle("signal eff");
   MVA_BDT_rejBvsS->GetYaxis()->SetTitle("backgr rejection (1-eff)");
   MVA_BDT_rejBvsS->Draw("csame");
   
   TH1 *MVA_Fisher_rejBvsS = new TH1F("MVA_Fisher_rejBvsS","MVA_Fisher",100,0,1);
   MVA_Fisher_rejBvsS->SetBinContent(1,0.999517);
   MVA_Fisher_rejBvsS->SetBinContent(2,0.99876);
   MVA_Fisher_rejBvsS->SetBinContent(3,0.997573);
   MVA_Fisher_rejBvsS->SetBinContent(4,0.995966);
   MVA_Fisher_rejBvsS->SetBinContent(5,0.99457);
   MVA_Fisher_rejBvsS->SetBinContent(6,0.992677);
   MVA_Fisher_rejBvsS->SetBinContent(7,0.992072);
   MVA_Fisher_rejBvsS->SetBinContent(8,0.990324);
   MVA_Fisher_rejBvsS->SetBinContent(9,0.98917);
   MVA_Fisher_rejBvsS->SetBinContent(10,0.987145);
   MVA_Fisher_rejBvsS->SetBinContent(11,0.986491);
   MVA_Fisher_rejBvsS->SetBinContent(12,0.985541);
   MVA_Fisher_rejBvsS->SetBinContent(13,0.984646);
   MVA_Fisher_rejBvsS->SetBinContent(14,0.982344);
   MVA_Fisher_rejBvsS->SetBinContent(15,0.979214);
   MVA_Fisher_rejBvsS->SetBinContent(16,0.977092);
   MVA_Fisher_rejBvsS->SetBinContent(17,0.974054);
   MVA_Fisher_rejBvsS->SetBinContent(18,0.972071);
   MVA_Fisher_rejBvsS->SetBinContent(19,0.968897);
   MVA_Fisher_rejBvsS->SetBinContent(20,0.968153);
   MVA_Fisher_rejBvsS->SetBinContent(21,0.965785);
   MVA_Fisher_rejBvsS->SetBinContent(22,0.963124);
   MVA_Fisher_rejBvsS->SetBinContent(23,0.960515);
   MVA_Fisher_rejBvsS->SetBinContent(24,0.957593);
   MVA_Fisher_rejBvsS->SetBinContent(25,0.955433);
   MVA_Fisher_rejBvsS->SetBinContent(26,0.952896);
   MVA_Fisher_rejBvsS->SetBinContent(27,0.950163);
   MVA_Fisher_rejBvsS->SetBinContent(28,0.947994);
   MVA_Fisher_rejBvsS->SetBinContent(29,0.946201);
   MVA_Fisher_rejBvsS->SetBinContent(30,0.94376);
   MVA_Fisher_rejBvsS->SetBinContent(31,0.940612);
   MVA_Fisher_rejBvsS->SetBinContent(32,0.938321);
   MVA_Fisher_rejBvsS->SetBinContent(33,0.933889);
   MVA_Fisher_rejBvsS->SetBinContent(34,0.930632);
   MVA_Fisher_rejBvsS->SetBinContent(35,0.926064);
   MVA_Fisher_rejBvsS->SetBinContent(36,0.920407);
   MVA_Fisher_rejBvsS->SetBinContent(37,0.916375);
   MVA_Fisher_rejBvsS->SetBinContent(38,0.912577);
   MVA_Fisher_rejBvsS->SetBinContent(39,0.908476);
   MVA_Fisher_rejBvsS->SetBinContent(40,0.903299);
   MVA_Fisher_rejBvsS->SetBinContent(41,0.899642);
   MVA_Fisher_rejBvsS->SetBinContent(42,0.89658);
   MVA_Fisher_rejBvsS->SetBinContent(43,0.892669);
   MVA_Fisher_rejBvsS->SetBinContent(44,0.88809);
   MVA_Fisher_rejBvsS->SetBinContent(45,0.884205);
   MVA_Fisher_rejBvsS->SetBinContent(46,0.880766);
   MVA_Fisher_rejBvsS->SetBinContent(47,0.876381);
   MVA_Fisher_rejBvsS->SetBinContent(48,0.869708);
   MVA_Fisher_rejBvsS->SetBinContent(49,0.865738);
   MVA_Fisher_rejBvsS->SetBinContent(50,0.860085);
   MVA_Fisher_rejBvsS->SetBinContent(51,0.854962);
   MVA_Fisher_rejBvsS->SetBinContent(52,0.850326);
   MVA_Fisher_rejBvsS->SetBinContent(53,0.844698);
   MVA_Fisher_rejBvsS->SetBinContent(54,0.84022);
   MVA_Fisher_rejBvsS->SetBinContent(55,0.835128);
   MVA_Fisher_rejBvsS->SetBinContent(56,0.829885);
   MVA_Fisher_rejBvsS->SetBinContent(57,0.823556);
   MVA_Fisher_rejBvsS->SetBinContent(58,0.818264);
   MVA_Fisher_rejBvsS->SetBinContent(59,0.812738);
   MVA_Fisher_rejBvsS->SetBinContent(60,0.807467);
   MVA_Fisher_rejBvsS->SetBinContent(61,0.801566);
   MVA_Fisher_rejBvsS->SetBinContent(62,0.792718);
   MVA_Fisher_rejBvsS->SetBinContent(63,0.787314);
   MVA_Fisher_rejBvsS->SetBinContent(64,0.779798);
   MVA_Fisher_rejBvsS->SetBinContent(65,0.771163);
   MVA_Fisher_rejBvsS->SetBinContent(66,0.762074);
   MVA_Fisher_rejBvsS->SetBinContent(67,0.7551);
   MVA_Fisher_rejBvsS->SetBinContent(68,0.745752);
   MVA_Fisher_rejBvsS->SetBinContent(69,0.735424);
   MVA_Fisher_rejBvsS->SetBinContent(70,0.725324);
   MVA_Fisher_rejBvsS->SetBinContent(71,0.71299);
   MVA_Fisher_rejBvsS->SetBinContent(72,0.705165);
   MVA_Fisher_rejBvsS->SetBinContent(73,0.694113);
   MVA_Fisher_rejBvsS->SetBinContent(74,0.683108);
   MVA_Fisher_rejBvsS->SetBinContent(75,0.670822);
   MVA_Fisher_rejBvsS->SetBinContent(76,0.660654);
   MVA_Fisher_rejBvsS->SetBinContent(77,0.651417);
   MVA_Fisher_rejBvsS->SetBinContent(78,0.633455);
   MVA_Fisher_rejBvsS->SetBinContent(79,0.623424);
   MVA_Fisher_rejBvsS->SetBinContent(80,0.608087);
   MVA_Fisher_rejBvsS->SetBinContent(81,0.592989);
   MVA_Fisher_rejBvsS->SetBinContent(82,0.580502);
   MVA_Fisher_rejBvsS->SetBinContent(83,0.559074);
   MVA_Fisher_rejBvsS->SetBinContent(84,0.5437);
   MVA_Fisher_rejBvsS->SetBinContent(85,0.525862);
   MVA_Fisher_rejBvsS->SetBinContent(86,0.510847);
   MVA_Fisher_rejBvsS->SetBinContent(87,0.489605);
   MVA_Fisher_rejBvsS->SetBinContent(88,0.470206);
   MVA_Fisher_rejBvsS->SetBinContent(89,0.448106);
   MVA_Fisher_rejBvsS->SetBinContent(90,0.420391);
   MVA_Fisher_rejBvsS->SetBinContent(91,0.393161);
   MVA_Fisher_rejBvsS->SetBinContent(92,0.370534);
   MVA_Fisher_rejBvsS->SetBinContent(93,0.337573);
   MVA_Fisher_rejBvsS->SetBinContent(94,0.303176);
   MVA_Fisher_rejBvsS->SetBinContent(95,0.267573);
   MVA_Fisher_rejBvsS->SetBinContent(96,0.232942);
   MVA_Fisher_rejBvsS->SetBinContent(97,0.194267);
   MVA_Fisher_rejBvsS->SetBinContent(98,0.152642);
   MVA_Fisher_rejBvsS->SetBinContent(99,0.102358);
   MVA_Fisher_rejBvsS->SetBinContent(100,0.0372351);
   MVA_Fisher_rejBvsS->SetEntries(100);
   MVA_Fisher_rejBvsS->SetStats(0);
   MVA_Fisher_rejBvsS->SetLineColor(2);
   MVA_Fisher_rejBvsS->SetLineWidth(3);
   MVA_Fisher_rejBvsS->SetMarkerStyle(21);
   MVA_Fisher_rejBvsS->SetMarkerSize(0.3);
   MVA_Fisher_rejBvsS->GetXaxis()->SetTitle("signal eff");
   MVA_Fisher_rejBvsS->GetYaxis()->SetTitle("backgr rejection (1-eff)");
   MVA_Fisher_rejBvsS->Draw("csame");
   
   TH1 *MVA_Likelihood_rejBvsS = new TH1F("MVA_Likelihood_rejBvsS","MVA_Likelihood",100,0,1);
   MVA_Likelihood_rejBvsS->SetBinContent(1,0.999641);
   MVA_Likelihood_rejBvsS->SetBinContent(2,0.999178);
   MVA_Likelihood_rejBvsS->SetBinContent(3,0.998604);
   MVA_Likelihood_rejBvsS->SetBinContent(4,0.997648);
   MVA_Likelihood_rejBvsS->SetBinContent(5,0.996876);
   MVA_Likelihood_rejBvsS->SetBinContent(6,0.996069);
   MVA_Likelihood_rejBvsS->SetBinContent(7,0.99475);
   MVA_Likelihood_rejBvsS->SetBinContent(8,0.99373);
   MVA_Likelihood_rejBvsS->SetBinContent(9,0.993112);
   MVA_Likelihood_rejBvsS->SetBinContent(10,0.991892);
   MVA_Likelihood_rejBvsS->SetBinContent(11,0.991218);
   MVA_Likelihood_rejBvsS->SetBinContent(12,0.990278);
   MVA_Likelihood_rejBvsS->SetBinContent(13,0.988945);
   MVA_Likelihood_rejBvsS->SetBinContent(14,0.988119);
   MVA_Likelihood_rejBvsS->SetBinContent(15,0.98732);
   MVA_Likelihood_rejBvsS->SetBinContent(16,0.985731);
   MVA_Likelihood_rejBvsS->SetBinContent(17,0.984675);
   MVA_Likelihood_rejBvsS->SetBinContent(18,0.982762);
   MVA_Likelihood_rejBvsS->SetBinContent(19,0.981429);
   MVA_Likelihood_rejBvsS->SetBinContent(20,0.979483);
   MVA_Likelihood_rejBvsS->SetBinContent(21,0.978237);
   MVA_Likelihood_rejBvsS->SetBinContent(22,0.976168);
   MVA_Likelihood_rejBvsS->SetBinContent(23,0.973868);
   MVA_Likelihood_rejBvsS->SetBinContent(24,0.971381);
   MVA_Likelihood_rejBvsS->SetBinContent(25,0.968858);
   MVA_Likelihood_rejBvsS->SetBinContent(26,0.967278);
   MVA_Likelihood_rejBvsS->SetBinContent(27,0.965332);
   MVA_Likelihood_rejBvsS->SetBinContent(28,0.96208);
   MVA_Likelihood_rejBvsS->SetBinContent(29,0.958921);
   MVA_Likelihood_rejBvsS->SetBinContent(30,0.956396);
   MVA_Likelihood_rejBvsS->SetBinContent(31,0.954661);
   MVA_Likelihood_rejBvsS->SetBinContent(32,0.950913);
   MVA_Likelihood_rejBvsS->SetBinContent(33,0.949357);
   MVA_Likelihood_rejBvsS->SetBinContent(34,0.945956);
   MVA_Likelihood_rejBvsS->SetBinContent(35,0.943117);
   MVA_Likelihood_rejBvsS->SetBinContent(36,0.939745);
   MVA_Likelihood_rejBvsS->SetBinContent(37,0.937627);
   MVA_Likelihood_rejBvsS->SetBinContent(38,0.933636);
   MVA_Likelihood_rejBvsS->SetBinContent(39,0.929935);
   MVA_Likelihood_rejBvsS->SetBinContent(40,0.926329);
   MVA_Likelihood_rejBvsS->SetBinContent(41,0.920709);
   MVA_Likelihood_rejBvsS->SetBinContent(42,0.916536);
   MVA_Likelihood_rejBvsS->SetBinContent(43,0.912573);
   MVA_Likelihood_rejBvsS->SetBinContent(44,0.909263);
   MVA_Likelihood_rejBvsS->SetBinContent(45,0.905056);
   MVA_Likelihood_rejBvsS->SetBinContent(46,0.902064);
   MVA_Likelihood_rejBvsS->SetBinContent(47,0.898747);
   MVA_Likelihood_rejBvsS->SetBinContent(48,0.893631);
   MVA_Likelihood_rejBvsS->SetBinContent(49,0.889311);
   MVA_Likelihood_rejBvsS->SetBinContent(50,0.883735);
   MVA_Likelihood_rejBvsS->SetBinContent(51,0.878458);
   MVA_Likelihood_rejBvsS->SetBinContent(52,0.873827);
   MVA_Likelihood_rejBvsS->SetBinContent(53,0.867891);
   MVA_Likelihood_rejBvsS->SetBinContent(54,0.862876);
   MVA_Likelihood_rejBvsS->SetBinContent(55,0.857492);
   MVA_Likelihood_rejBvsS->SetBinContent(56,0.852259);
   MVA_Likelihood_rejBvsS->SetBinContent(57,0.846485);
   MVA_Likelihood_rejBvsS->SetBinContent(58,0.83956);
   MVA_Likelihood_rejBvsS->SetBinContent(59,0.832632);
   MVA_Likelihood_rejBvsS->SetBinContent(60,0.826106);
   MVA_Likelihood_rejBvsS->SetBinContent(61,0.819874);
   MVA_Likelihood_rejBvsS->SetBinContent(62,0.811465);
   MVA_Likelihood_rejBvsS->SetBinContent(63,0.80226);
   MVA_Likelihood_rejBvsS->SetBinContent(64,0.795147);
   MVA_Likelihood_rejBvsS->SetBinContent(65,0.788993);
   MVA_Likelihood_rejBvsS->SetBinContent(66,0.779508);
   MVA_Likelihood_rejBvsS->SetBinContent(67,0.77212);
   MVA_Likelihood_rejBvsS->SetBinContent(68,0.765169);
   MVA_Likelihood_rejBvsS->SetBinContent(69,0.756221);
   MVA_Likelihood_rejBvsS->SetBinContent(70,0.746177);
   MVA_Likelihood_rejBvsS->SetBinContent(71,0.735721);
   MVA_Likelihood_rejBvsS->SetBinContent(72,0.72492);
   MVA_Likelihood_rejBvsS->SetBinContent(73,0.717486);
   MVA_Likelihood_rejBvsS->SetBinContent(74,0.707751);
   MVA_Likelihood_rejBvsS->SetBinContent(75,0.694377);
   MVA_Likelihood_rejBvsS->SetBinContent(76,0.68487);
   MVA_Likelihood_rejBvsS->SetBinContent(77,0.671831);
   MVA_Likelihood_rejBvsS->SetBinContent(78,0.658245);
   MVA_Likelihood_rejBvsS->SetBinContent(79,0.647587);
   MVA_Likelihood_rejBvsS->SetBinContent(80,0.628551);
   MVA_Likelihood_rejBvsS->SetBinContent(81,0.610466);
   MVA_Likelihood_rejBvsS->SetBinContent(82,0.593383);
   MVA_Likelihood_rejBvsS->SetBinContent(83,0.576184);
   MVA_Likelihood_rejBvsS->SetBinContent(84,0.558991);
   MVA_Likelihood_rejBvsS->SetBinContent(85,0.539532);
   MVA_Likelihood_rejBvsS->SetBinContent(86,0.52456);
   MVA_Likelihood_rejBvsS->SetBinContent(87,0.505848);
   MVA_Likelihood_rejBvsS->SetBinContent(88,0.484198);
   MVA_Likelihood_rejBvsS->SetBinContent(89,0.459002);
   MVA_Likelihood_rejBvsS->SetBinContent(90,0.436671);
   MVA_Likelihood_rejBvsS->SetBinContent(91,0.412721);
   MVA_Likelihood_rejBvsS->SetBinContent(92,0.391933);
   MVA_Likelihood_rejBvsS->SetBinContent(93,0.362264);
   MVA_Likelihood_rejBvsS->SetBinContent(94,0.328588);
   MVA_Likelihood_rejBvsS->SetBinContent(95,0.288286);
   MVA_Likelihood_rejBvsS->SetBinContent(96,0.251312);
   MVA_Likelihood_rejBvsS->SetBinContent(97,0.212345);
   MVA_Likelihood_rejBvsS->SetBinContent(98,0.159233);
   MVA_Likelihood_rejBvsS->SetBinContent(99,0.101972);
   MVA_Likelihood_rejBvsS->SetBinContent(100,0.043166);
   MVA_Likelihood_rejBvsS->SetEntries(100);
   MVA_Likelihood_rejBvsS->SetStats(0);
   MVA_Likelihood_rejBvsS->SetLineColor(3);
   MVA_Likelihood_rejBvsS->SetLineWidth(3);
   MVA_Likelihood_rejBvsS->SetMarkerStyle(21);
   MVA_Likelihood_rejBvsS->SetMarkerSize(0.3);
   MVA_Likelihood_rejBvsS->GetXaxis()->SetTitle("signal eff");
   MVA_Likelihood_rejBvsS->GetYaxis()->SetTitle("backgr rejection (1-eff)");
   MVA_Likelihood_rejBvsS->Draw("csame");
   
   TH1 *MVA_MLP_rejBvsS = new TH1F("MVA_MLP_rejBvsS","MVA_MLP",100,0,1);
   MVA_MLP_rejBvsS->SetBinContent(1,0.999068);
   MVA_MLP_rejBvsS->SetBinContent(2,0.99856);
   MVA_MLP_rejBvsS->SetBinContent(3,0.997362);
   MVA_MLP_rejBvsS->SetBinContent(4,0.996121);
   MVA_MLP_rejBvsS->SetBinContent(5,0.995426);
   MVA_MLP_rejBvsS->SetBinContent(6,0.99379);
   MVA_MLP_rejBvsS->SetBinContent(7,0.993316);
   MVA_MLP_rejBvsS->SetBinContent(8,0.992019);
   MVA_MLP_rejBvsS->SetBinContent(9,0.991035);
   MVA_MLP_rejBvsS->SetBinContent(10,0.990064);
   MVA_MLP_rejBvsS->SetBinContent(11,0.988578);
   MVA_MLP_rejBvsS->SetBinContent(12,0.986876);
   MVA_MLP_rejBvsS->SetBinContent(13,0.985599);
   MVA_MLP_rejBvsS->SetBinContent(14,0.983714);
   MVA_MLP_rejBvsS->SetBinContent(15,0.982259);
   MVA_MLP_rejBvsS->SetBinContent(16,0.980299);
   MVA_MLP_rejBvsS->SetBinContent(17,0.978724);
   MVA_MLP_rejBvsS->SetBinContent(18,0.977748);
   MVA_MLP_rejBvsS->SetBinContent(19,0.976287);
   MVA_MLP_rejBvsS->SetBinContent(20,0.974392);
   MVA_MLP_rejBvsS->SetBinContent(21,0.971946);
   MVA_MLP_rejBvsS->SetBinContent(22,0.970244);
   MVA_MLP_rejBvsS->SetBinContent(23,0.96787);
   MVA_MLP_rejBvsS->SetBinContent(24,0.966512);
   MVA_MLP_rejBvsS->SetBinContent(25,0.964559);
   MVA_MLP_rejBvsS->SetBinContent(26,0.961885);
   MVA_MLP_rejBvsS->SetBinContent(27,0.959453);
   MVA_MLP_rejBvsS->SetBinContent(28,0.957505);
   MVA_MLP_rejBvsS->SetBinContent(29,0.954738);
   MVA_MLP_rejBvsS->SetBinContent(30,0.951504);
   MVA_MLP_rejBvsS->SetBinContent(31,0.949791);
   MVA_MLP_rejBvsS->SetBinContent(32,0.946907);
   MVA_MLP_rejBvsS->SetBinContent(33,0.9434);
   MVA_MLP_rejBvsS->SetBinContent(34,0.941844);
   MVA_MLP_rejBvsS->SetBinContent(35,0.938901);
   MVA_MLP_rejBvsS->SetBinContent(36,0.93726);
   MVA_MLP_rejBvsS->SetBinContent(37,0.934413);
   MVA_MLP_rejBvsS->SetBinContent(38,0.930402);
   MVA_MLP_rejBvsS->SetBinContent(39,0.925513);
   MVA_MLP_rejBvsS->SetBinContent(40,0.922843);
   MVA_MLP_rejBvsS->SetBinContent(41,0.920417);
   MVA_MLP_rejBvsS->SetBinContent(42,0.91779);
   MVA_MLP_rejBvsS->SetBinContent(43,0.912466);
   MVA_MLP_rejBvsS->SetBinContent(44,0.908053);
   MVA_MLP_rejBvsS->SetBinContent(45,0.904962);
   MVA_MLP_rejBvsS->SetBinContent(46,0.901399);
   MVA_MLP_rejBvsS->SetBinContent(47,0.896701);
   MVA_MLP_rejBvsS->SetBinContent(48,0.892511);
   MVA_MLP_rejBvsS->SetBinContent(49,0.888742);
   MVA_MLP_rejBvsS->SetBinContent(50,0.882579);
   MVA_MLP_rejBvsS->SetBinContent(51,0.880225);
   MVA_MLP_rejBvsS->SetBinContent(52,0.875466);
   MVA_MLP_rejBvsS->SetBinContent(53,0.87132);
   MVA_MLP_rejBvsS->SetBinContent(54,0.86573);
   MVA_MLP_rejBvsS->SetBinContent(55,0.859795);
   MVA_MLP_rejBvsS->SetBinContent(56,0.856475);
   MVA_MLP_rejBvsS->SetBinContent(57,0.850415);
   MVA_MLP_rejBvsS->SetBinContent(58,0.843168);
   MVA_MLP_rejBvsS->SetBinContent(59,0.83778);
   MVA_MLP_rejBvsS->SetBinContent(60,0.829608);
   MVA_MLP_rejBvsS->SetBinContent(61,0.823767);
   MVA_MLP_rejBvsS->SetBinContent(62,0.816876);
   MVA_MLP_rejBvsS->SetBinContent(63,0.808367);
   MVA_MLP_rejBvsS->SetBinContent(64,0.802159);
   MVA_MLP_rejBvsS->SetBinContent(65,0.793444);
   MVA_MLP_rejBvsS->SetBinContent(66,0.784533);
   MVA_MLP_rejBvsS->SetBinContent(67,0.777028);
   MVA_MLP_rejBvsS->SetBinContent(68,0.769886);
   MVA_MLP_rejBvsS->SetBinContent(69,0.761204);
   MVA_MLP_rejBvsS->SetBinContent(70,0.748396);
   MVA_MLP_rejBvsS->SetBinContent(71,0.740148);
   MVA_MLP_rejBvsS->SetBinContent(72,0.730887);
   MVA_MLP_rejBvsS->SetBinContent(73,0.720931);
   MVA_MLP_rejBvsS->SetBinContent(74,0.71042);
   MVA_MLP_rejBvsS->SetBinContent(75,0.701005);
   MVA_MLP_rejBvsS->SetBinContent(76,0.690861);
   MVA_MLP_rejBvsS->SetBinContent(77,0.677795);
   MVA_MLP_rejBvsS->SetBinContent(78,0.663915);
   MVA_MLP_rejBvsS->SetBinContent(79,0.650732);
   MVA_MLP_rejBvsS->SetBinContent(80,0.632826);
   MVA_MLP_rejBvsS->SetBinContent(81,0.619294);
   MVA_MLP_rejBvsS->SetBinContent(82,0.607352);
   MVA_MLP_rejBvsS->SetBinContent(83,0.588325);
   MVA_MLP_rejBvsS->SetBinContent(84,0.573247);
   MVA_MLP_rejBvsS->SetBinContent(85,0.556829);
   MVA_MLP_rejBvsS->SetBinContent(86,0.539586);
   MVA_MLP_rejBvsS->SetBinContent(87,0.518965);
   MVA_MLP_rejBvsS->SetBinContent(88,0.49396);
   MVA_MLP_rejBvsS->SetBinContent(89,0.464501);
   MVA_MLP_rejBvsS->SetBinContent(90,0.435213);
   MVA_MLP_rejBvsS->SetBinContent(91,0.404139);
   MVA_MLP_rejBvsS->SetBinContent(92,0.375079);
   MVA_MLP_rejBvsS->SetBinContent(93,0.345841);
   MVA_MLP_rejBvsS->SetBinContent(94,0.316698);
   MVA_MLP_rejBvsS->SetBinContent(95,0.283725);
   MVA_MLP_rejBvsS->SetBinContent(96,0.245958);
   MVA_MLP_rejBvsS->SetBinContent(97,0.198987);
   MVA_MLP_rejBvsS->SetBinContent(98,0.15001);
   MVA_MLP_rejBvsS->SetBinContent(99,0.0875898);
   MVA_MLP_rejBvsS->SetBinContent(100,0.0314129);
   MVA_MLP_rejBvsS->SetEntries(100);
   MVA_MLP_rejBvsS->SetStats(0);
   MVA_MLP_rejBvsS->SetLineColor(4);
   MVA_MLP_rejBvsS->SetLineWidth(3);
   MVA_MLP_rejBvsS->SetMarkerStyle(21);
   MVA_MLP_rejBvsS->SetMarkerSize(0.3);
   MVA_MLP_rejBvsS->GetXaxis()->SetTitle("signal eff");
   MVA_MLP_rejBvsS->GetYaxis()->SetTitle("backgr rejection (1-eff)");
   MVA_MLP_rejBvsS->Draw("csame");
   
   TH1 *MVA_SVM_Gauss_rejBvsS = new TH1F("MVA_SVM_Gauss_rejBvsS","MVA_SVM_Gauss",100,0,1);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(1,0.999185);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(2,0.995499);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(3,0.990925);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(4,0.985588);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(5,0.981931);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(6,0.976917);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(7,0.97338);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(8,0.967957);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(9,0.962884);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(10,0.95567);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(11,0.950721);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(12,0.944938);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(13,0.93885);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(14,0.932912);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(15,0.928356);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(16,0.924095);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(17,0.919248);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(18,0.91516);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(19,0.909624);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(20,0.903913);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(21,0.897304);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(22,0.891347);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(23,0.883604);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(24,0.87742);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(25,0.871895);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(26,0.866421);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(27,0.860418);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(28,0.852779);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(29,0.845251);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(30,0.838116);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(31,0.830981);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(32,0.823847);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(33,0.816712);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(34,0.809577);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(35,0.802443);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(36,0.795308);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(37,0.788173);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(38,0.781039);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(39,0.77092);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(40,0.758563);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(41,0.746207);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(42,0.73385);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(43,0.721493);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(44,0.709136);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(45,0.69678);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(46,0.684423);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(47,0.672066);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(48,0.659709);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(49,0.647353);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(50,0.634996);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(51,0.622639);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(52,0.610282);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(53,0.597926);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(54,0.585569);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(55,0.573212);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(56,0.560856);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(57,0.548499);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(58,0.536142);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(59,0.523785);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(60,0.511428);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(61,0.499072);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(62,0.486715);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(63,0.474358);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(64,0.462001);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(65,0.449645);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(66,0.437288);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(67,0.424931);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(68,0.412574);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(69,0.400218);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(70,0.387861);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(71,0.375504);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(72,0.363147);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(73,0.350791);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(74,0.338434);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(75,0.326077);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(76,0.31372);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(77,0.301364);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(78,0.289007);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(79,0.27665);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(80,0.264293);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(81,0.251937);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(82,0.23958);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(83,0.227223);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(84,0.214866);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(85,0.20251);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(86,0.190153);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(87,0.177796);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(88,0.165439);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(89,0.153083);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(90,0.140726);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(91,0.128369);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(92,0.116012);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(93,0.103656);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(94,0.091299);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(95,0.0775716);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(96,0.0597602);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(97,0.0445914);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(98,0.0297268);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(99,0.0148898);
   MVA_SVM_Gauss_rejBvsS->SetBinContent(100,0.00420487);
   MVA_SVM_Gauss_rejBvsS->SetEntries(100);
   MVA_SVM_Gauss_rejBvsS->SetStats(0);
   MVA_SVM_Gauss_rejBvsS->SetLineColor(6);
   MVA_SVM_Gauss_rejBvsS->SetLineWidth(3);
   MVA_SVM_Gauss_rejBvsS->SetMarkerStyle(21);
   MVA_SVM_Gauss_rejBvsS->SetMarkerSize(0.3);
   MVA_SVM_Gauss_rejBvsS->GetXaxis()->SetTitle("signal eff");
   MVA_SVM_Gauss_rejBvsS->GetYaxis()->SetTitle("backgr rejection (1-eff)");
   MVA_SVM_Gauss_rejBvsS->Draw("csame");
   
   TH1 *frame = new TH2F("frame","Background rejection versus Signal efficiency",500,0,1,500,0.2,1);
   frame->SetStats(0);
   frame->SetLineWidth(2);
   frame->SetMarkerStyle(21);
   frame->SetMarkerSize(0.3);
   frame->GetXaxis()->SetTitle("Signal efficiency");
   frame->GetXaxis()->SetLabelOffset(0.012);
   frame->GetXaxis()->SetTitleSize(0.045);
   frame->GetXaxis()->SetTitleOffset(1.25);
   frame->GetYaxis()->SetTitle("Background rejection");
   frame->GetYaxis()->SetLabelOffset(0.012);
   frame->GetYaxis()->SetTitleSize(0.045);
   frame->GetYaxis()->SetTitleOffset(1.22);
   frame->Draw("sameaxis");
   
   TLegend *leg = new TLegend(0.15,0.171,0.5,0.501,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("NULL","MVA Method:","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("MVA_BDT_rejBvsS","BDT","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextAlign(12);
   entry->SetTextColor(1);
   entry=leg->AddEntry("MVA_Likelihood_rejBvsS","Likelihood","l");
   entry->SetLineColor(3);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextAlign(12);
   entry->SetTextColor(1);
   entry=leg->AddEntry("MVA_MLP_rejBvsS","MLP","l");
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextAlign(12);
   entry->SetTextColor(1);
   entry=leg->AddEntry("MVA_Fisher_rejBvsS","Fisher","l");
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextAlign(12);
   entry->SetTextColor(1);
   entry=leg->AddEntry("MVA_SVM_Gauss_rejBvsS","SVM_Gauss","l");
   entry->SetLineColor(6);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextAlign(12);
   entry->SetTextColor(1);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.01,0.939068,0.71,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(1);
   pt->SetFillColor(183);
   TText *text = pt->AddText("Background rejection versus Signal efficiency");
   pt->Draw();
  
// ------------>Primitives in pad: imgpad
   TPad *imgpad = new TPad("imgpad", "imgpad",0.789257,0.91,0.95,0.965);
   imgpad->Draw();
   imgpad->cd();
   imgpad->Range(0,0,1,1);
   imgpad->SetFillColor(0);
   imgpad->SetBorderMode(0);
   imgpad->SetBorderSize(2);
   imgpad->SetTickx();
   imgpad->SetTicky();
   imgpad->SetLeftMargin(0);
   imgpad->SetRightMargin(0);
   imgpad->SetTopMargin(0);
   imgpad->SetBottomMargin(0);

   ci = TColor::GetColor("#fffffd");
   imgpad->SetFrameFillColor(ci);
   imgpad->SetFrameBorderMode(0);
   imgpad->SetFrameLineColor(0);
   imgpad->SetFrameBorderMode(0);

/* XPM */
 char *xpm_tmva_logo_gif_[] = {
/* columns rows colors chars-per-pixel */
"455 107 6 1",
"  c #FF0707",
". c #383899",
"X c #3A3AA0",
"o c #FFAC28",
"O c #FFFFFF",
"+ c None",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                    OX...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X.OOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                    OX...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X.OOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                   OOX...X...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X..OOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                   OOX...X...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X..OOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                   OOX...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...OOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                    OOX...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...OOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                    OOX...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...XOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                   OOOX...X...X...X...X...XOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...XOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                   OOOX...X...X...X...X...XOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X.OOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                    OOOX...X...X...X...X...X.OOOOOOOOOOOOOOOOOOX...X...X...X...X...X.OOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                                                    OOOX...X...X...X...X...X.OOOOOOOOOOOOOOOOOOX...X...X...X...X...X.OOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X..OOOOOOOOOOOOOOOOX...X...X...X...X...X..OOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X..OOOOOOOOOOOOOOOOX...X...X...X...X...X..OOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOOX...X...X...X...X...X...OOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...XOOOOOOOOOOOOX...X...X...X...X...X...XOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...XOOOOOOOOOOOOX...X...X...X...X...X...XOOOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X.OOOOOOOOOOX...X...X...X...X...X...X.OOOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X.OOOOOOOOOOX...X...X...X...X...X...X.OOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...XOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X..OOOOOOOOX...X...X...X...X...X...X..OOOOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OX...X...X...OOOOOOOOX...X...X..OX...X...X...X..OOOOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OX...X...X...OOOOOOOX...X...X...OX...X...X...X..OOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOX...X...X...OOOOOOX...X...X..OOX...X...X...X..OOOOOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOX...X...X...OOOOOX...X...X...OOX...X...X...X..OOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOooooooooooooooooooooooooooooOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOX...X...X...XOOOOX...X...X...OOX...X...X...X..OOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOoooooooooooooooooooooooooOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOX...X...X...OOOX...X...X...OOOX...X...X...X..OOOOOOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOOOOoooooooooooooooooooooOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOX...X...X...XOOX...X...X...OOOX...X...X...X..OOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOOOOOOOOoooooooooooooooOOOOOOOOOOOOooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOX...X...X...OOX...X...X..OOOOX...X...X...X..OOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooOOOOOOOOOOOOOOOOooooOOOOOOOOOOOOOOOOOOooooOOOOOOOOOOOOOOOOOOooooOOOOOOOOOOOOOOOooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOX...X...X...X...X...X...XOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOX...X...X...X...OOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooOOOOOOOOOOOoooOOOOOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOOOOOOOOoooOOOOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOX...X...X...X...X...X..OOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooOOOOOOOOOOOOOOOOOoooOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOOOOOooOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOX...X...X...X...X...X..OOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOOooOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOX...X...X...X...X...XOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOX...X...X...X...OOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOooOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOX...X...X...X...X...XOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOX...X...X...X...X...XOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOX...X...X...X...X..OOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOX...X...X...X...X..OOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOX...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOX...X...X...X...XOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOX...X...X...X...OOOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOX...X...X...X...XOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOX...X...X...X..OOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOX...X...X...X..OOX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOX...X...X...X..OOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOX...X...X...X...OX...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...X.OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO                OOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOX...X...X...XOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...XOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOX...X...X...XOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...XOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO               OOOOOOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X.OOOOOOOOOOX...X...X...XOOOOOOOOOOX...X...X...X..OOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...XOOOOOOOOOOOOOOOOOOOOOOOX...X...X...X...X...X...X...X...X...X...X...X...OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
"OoooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooo                oooooooooooooooooooooooooooo...X...X...X..ooooooooooo..X...X...Xooooooooooo...X...X...X...oooooooooooooooooooooooo..X...X...X...X...X...X...Xooooooooooooooooooooooo...X...X...X...X...X...X...X...X...X...X...X...X..ooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooo                oooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xoooooooooooooooooooooooo...X...X...X...X...X...X...ooooooooooooooooooooooo..X...X...X...X...X...X...X...X...X...X...X...X...ooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooo               ooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooo...X...X...X...X...X...X.oooooooooooooooooooooooo..X...X...X...X...X...X...X...X...X...X...X...X...ooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooo               ooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooo...X...X...X...X...X...X.ooooooooooooooooooooooo..X...X...X...X...X...X...X...X...X...X...X...X...X.oooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooo               ooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooo...X...X...X...X...X...X.ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooo...X...X...X...XoooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooo                ooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xoooooooooooooooooooooooooo...X...X...X...X...X...oooooooooooooooooooooooo..X...X...X...ooooooooooooooooooooooo..X...X...X...XoooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooo               oooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xoooooooooooooooooooooooooo...X...X...X...X...X...ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooo...X...X...X...XooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooo               oooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xoooooooooooooooooooooooooo...X...X...X...X...X...ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooo...X...X...X...XooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooo               oooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooo...X...X...X...X...X.ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooo...X...X...X...XoooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooo                oooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooo...X...X...X...X...X.ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooo...X...X...X...XoooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooo               ooooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooo...X...X...X...X...X.ooooooooooooooooooooooo..X...X...X...oooooooooooooooooooooooooo..X...X...X...X.oooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooo               ooooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xoooooooooooooooooooooooooooo...X...X...X...X...ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooo...X...X...X...XooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooo               ooooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xoooooooooooooooooooooooooooo...X...X...X...X...ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooo...X...X...X...XooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooo                ooooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooooo...X...X...X...X.oooooooooooooooooooooooo..X...X...X...oooooooooooooooooooooooooooo..X...X...X...X.ooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooo                ooooooooooooooooooooooooooooooo...X...X...X..ooooooooooooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooooo...X...X...X...X.ooooooooooooooooooooooo..X...X...X...Xooooooooooooooooooooooooooooo...X...X...X...XoooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOOOOOOOOOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"OoooooooooooooooooooooooooooooooooooooooooooooooOOOoooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooOOOOooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooO",
"+OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOoooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"
};


   TImage *tmva_logo_gif_ = TImage::Create();
   tmva_logo_gif_->SetImageBuffer((char**)&xpm_tmva_logo_gif_, TImage::kXpm);
   tmva_logo_gif_->Draw();
   imgpad->Modified();
   c->cd();
   c->Modified();
   c->cd();
   c->SetSelected(c);
   TString pName = "public_html/EvtSel_ROC.eps"; 
   c->SaveAs(pName);
}
