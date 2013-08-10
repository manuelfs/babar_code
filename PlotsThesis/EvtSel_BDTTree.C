void EvtSel_BDTTree(){
   TCanvas *c1 = new TCanvas("c1", "Tree 151",220,136,1044,522);
   gStyle->SetOptStat(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(170);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetTickx();
   c1->SetTicky();
   c1->SetLeftMargin(0.12);
   c1->SetRightMargin(0.05);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#fffffd");
   c1->SetFrameFillColor(ci);
   c1->SetFrameBorderMode(0);
   TLine *line = new TLine(0.375,0.9,0.25,0.8);
   line->SetLineWidth(2);
   line->Draw();
   line = new TLine(0.19,0.7,0.16,0.59);
   line->SetLineWidth(2);
   line->Draw();
   line = new TLine(0.1,0.49,0.07,0.38);
   line->SetLineWidth(2);
   line->Draw();
   
   TPaveText *pt = new TPaveText(0.02,0.24,0.16,0.38,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#dd0032");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   TText *text = pt->AddText("N=656");
   text = pt->AddText("S/(S+B)=0.460");
   pt->Draw();
   line = new TLine(0.24,0.49,0.27,0.38);
   line->SetLineWidth(2);
   line->Draw();
   line = new TLine(0.18,0.27,0.15,0.17);
   line->SetLineWidth(2);
   line->Draw();
   
   pt = new TPaveText(0.10,0.03,0.24,0.17,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#dd0032");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=7540");
   text = pt->AddText("S/(S+B)=0.492");
   pt->Draw();
   line = new TLine(0.32,0.27,0.35,0.17);
   line->SetLineWidth(2);
   line->Draw();
   
   pt = new TPaveText(0.26,0.03,0.4,0.17,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#2144a5");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=454");
   text = pt->AddText("S/(S+B)=0.555");
   pt->Draw();
   
   pt = new TPaveText(0.18,0.18,0.32,0.38,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#32aa76");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=7994");
   text = pt->AddText("S/(S+B)=0.495");
   text = pt->AddText("E_{Extra} < 0.195");
   pt->Draw();
   
   pt = new TPaveText(0.1,0.39,0.24,0.59,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#32aa76");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=8651");
   text = pt->AddText("S/(S+B)=0.493");
   text = pt->AddText("cos(#Delta#theta_{T}) > -0.75");
   pt->Draw();
   line = new TLine(0.33,0.7,0.36,0.59);
   line->SetLineWidth(2);
   line->Draw();
   
   pt = new TPaveText(0.28,0.45,0.42,0.59,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#2144a5");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=73");
   text = pt->AddText("S/(S+B)=0.643");
   pt->Draw();
   
   pt = new TPaveText(0.19,0.6,0.33,0.8,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#32aa76");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=8725");
   text = pt->AddText("S/(S+B)=0.494");
   text = pt->AddText("#DeltaE > 0.0683");
   pt->Draw();
   line = new TLine(0.625,0.9,0.75,0.8);
   line->SetLineWidth(2);
   line->Draw();
   line = new TLine(0.69,0.7,0.66,0.59);
   line->SetLineWidth(2);
   line->Draw();
   
   pt = new TPaveText(0.6,0.45,0.74,0.59,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#dd0032");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=798");
   text = pt->AddText("S/(S+B)=0.466");
   pt->Draw();
   line = new TLine(0.83,0.7,0.86,0.59);
   line->SetLineWidth(2);
   line->Draw();
   
   pt = new TPaveText(0.78,0.45,0.92,0.59,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#2144a5");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=6481");
   text = pt->AddText("S/(S+B)=0.522");
   pt->Draw();
   
   pt = new TPaveText(0.69,0.6,0.83,0.8,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#32aa76");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=7275");
   text = pt->AddText("S/(S+B)=0.516");
   text = pt->AddText("#DeltaE < 0.0395");
   pt->Draw();
   
   pt = new TPaveText(0.375,0.78,0.625,0.985,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#32aa76");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("N=16000");
   text = pt->AddText("S/(S+B)=0.504");
   text = pt->AddText("E_{Extra} < 0.18");
   pt->Draw();
   
   pt = new TPaveText(0.82,0.92,0.98,0.987,"brNDC");
   pt->SetBorderSize(1);
   pt->SetFillColor(155);
   pt->SetFillStyle(1001);
   text = pt->AddText("Decision Tree no.: 151");
   pt->Draw();
   
   pt = new TPaveText(0.015,0.91,0.17,0.985,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#32aa76");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("Intermediate Nodes");
   pt->Draw();
   
   pt = new TPaveText(0.015,0.825,0.17,0.9,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#2144a5");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("Signal Leaf Nodes");
   pt->Draw();
   
   pt = new TPaveText(0.015,0.735,0.17,0.815,"brNDC");
   pt->SetBorderSize(1);

   ci = TColor::GetColor("#dd0032");
   pt->SetFillColor(ci);
   pt->SetFillStyle(1001);
   pt->SetTextColor(10);
   text = pt->AddText("Backgr. Leaf Nodes");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   TString pName = "public_html/EvtSel_BDTTree.eps"; 
   c1->SaveAs(pName);
}
