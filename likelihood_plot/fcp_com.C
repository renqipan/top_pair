void fcp_com(){
	gStyle->SetOptStat(0);
	TCanvas *c=new TCanvas("c","",8, 30, 600, 600);
	TCanvas *c2=new TCanvas("c2","",8,30,600,600);
	c2->SetFillColor(0);
	c2->SetBorderMode(0);
	c2->SetBorderSize(2);
	c2->SetTickx(1);
	c2->SetTicky(1);
	c2->SetLeftMargin(0.17);
	c2->SetRightMargin(0.05);
	c2->SetTopMargin(0.07);
	c2->SetBottomMargin(0.13);
	c2->SetFrameFillStyle(0);
	c2->SetFrameBorderMode(0);
	c2->SetFrameFillStyle(0);
	c2->SetFrameBorderMode(0);
	TLegend *leg = new TLegend(0.3,0.67,0.7,0.90);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	int color[14]={1,kViolet,4,2,2,kCyan, kMagenta,kBlue-10, 3,kYellow,93};
	 TH1 *frame = c2->DrawFrame(-1,0,1,28); //Draw an empty pad frame 
	 //DrawFrame (Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, const char *title="")
	 TString fileName[]={"higgsCombinettbar_fcp.MultiDimFit.mH125.root",
	 					"tHq_fCP.root","ttH_scan.root",
	 					"higgsCombinethw_fcp_nosys22.MultiDimFit.mH125.root",
	 					"higgsCombinethw_fcp_nosys21.MultiDimFit.mH125.root",
						};
	TString text[]={"t#bar{t}","tHq","ttH","tHW","tHW(no ttH)",};

	 

	for(int i=0;i<5;i++){
		TGraph *gr0;
		TTree* limit;
		if(i==1||i==2){
			TFile fthq(fileName[i],"read");
			TCanvas *cthq = (TCanvas*)fthq.Get("can1");
			gr0 = (TGraph*)cthq->GetListOfPrimitives()->FindObject("Exp1D");
			gr0->GetXaxis()->SetTitle("");
			//c2->cd();
			//gr0->Draw();
		}
		else{
			TFile *file=TFile::Open(fileName[i]);
			limit=(TTree*)file->Get("limit");
			 c->cd();
			 limit->Draw("2*deltaNLL:x");
			 gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();		 
			 }
		gr0->Sort();
		gr0->SetLineColor(color[i]);
		gr0->SetLineWidth(2);
		c2->cd();
		gr0->Draw("csame");
		if(i==4) gr0->SetLineStyle(7);
		TLegendEntry* leg_entry0=leg->AddEntry(gr0,text[i],"l");
	    leg_entry0->SetTextColor(color[i]);
	    leg->Draw("same");
	    if(fileName[i].Contains("ttbar")){
	    	c->cd();
		 	limit->Draw("2*deltaNLL:-x");
		 	TGraph *gr1 = (TGraph*) gROOT->FindObject("Graph")->Clone();
		 	gr1->Sort();
			gr1->SetLineColor(color[i]);
			gr1->SetLineWidth(2);
			c2->cd();
			gr1->Draw("csame");
		
		 }

	}


	frame->GetXaxis()->SetNdivisions(505);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelOffset(0.007);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetTitleSize(0.06);
	frame->GetXaxis()->SetTitleOffset(0.9);
	frame->GetXaxis()->SetTitleFont(42);

	frame->GetYaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelOffset(0.007);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.06);
	frame->GetYaxis()->SetTitleOffset(1.1);
	frame->GetYaxis()->SetTitleFont(42);

	frame->GetXaxis()->CenterTitle();
	frame->GetYaxis()->CenterTitle();
	frame->GetXaxis()->SetTitle("f_{CP}");
	frame->GetYaxis()->SetTitle("-2#Deltaln L");


  TPaveText *oneSig = new TPaveText(0.50,0.16,0.65,0.21,"NDC");
  //TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
  oneSig->SetFillColor(0);
  oneSig->SetTextFont(42);
  oneSig->SetTextSize(0.04);
  oneSig->SetTextColor(1);
  oneSig->SetBorderSize(0);
  oneSig->AddText("68% CL"); 
  oneSig->Draw();
  TLine *l1=new TLine();
  l1->SetLineStyle(9);
  l1->SetLineWidth(2);
  l1->SetLineColor(1);
  l1->DrawLine(-1.0,1.0,1.0,1.0);
  l1->Draw("same");

  TPaveText *twoSig = new TPaveText(0.50,0.24,0.65,0.29,"NDC");
  //TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
  twoSig->SetFillColor(0);
  twoSig->SetTextFont(42);
  twoSig->SetTextSize(0.04);
  twoSig->SetTextColor(9);
  twoSig->SetBorderSize(0);
  twoSig->AddText("95% CL"); 
  twoSig->Draw();
  TLine *l2=new TLine();
  l2->SetLineStyle(9);
  l2->SetLineWidth(2);
  l2->SetLineColor(9);
  l2->DrawLine(-1.0,3.84,1.0,3.84);
  l2->Draw("same");
 // TLatex lumi;
 // lumi.DrawLatex(-0.35,15,"\\mathscr{L}=300 \\mathrm{fb^{-1}}");
  TPaveText *lumi = new TPaveText(0.68,0.80,0.78,0.90,"NDC");
  //TPaveText (Double_t x1, Double_t y1, Double_t x2, Double_t y2, Option_t *option="br")
  lumi->SetFillColor(0);
  lumi->SetTextFont(42);
  lumi->SetTextSize(0.055);
  lumi->SetTextColor(9);
  lumi->SetBorderSize(0);
  lumi->AddText("\\mathscr{L}=300 \\mathrm{fb^{-1}}"); 
  lumi->Draw();

  c2->Print("com_fcp.png");
  c2->Print("com_fcp.pdf");
  c2->Print("com_fcp.eps");


}
