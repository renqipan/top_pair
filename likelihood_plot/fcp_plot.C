void fcp_plot(){
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
	TLegend *leg = new TLegend(0.3,0.7,0.7,0.93);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	 TH1 *frame = c2->DrawFrame(-1,0,1,8); //Draw an empty pad frame 
	 //DrawFrame (Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, const char *title="")
	 TString input1="./higgsCombinethw_fcp_nosys22.MultiDimFit.mH125.root";
	 TFile *file1=TFile::Open(input1);
	 TTree* limit1=(TTree*)file1->Get("limit");
	 
	 c->cd();
	 limit1->Draw("2*deltaNLL:x");
	 TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
	 
	 gr0->Sort();
	 gr0->SetLineColor(4);
	 gr0->SetLineWidth(2);
	 c2->cd();
	 gr0->Draw("c");

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
	//TLegendEntry* leg_entry0=leg->AddEntry(gr0,"Fix #kappa=1","l");
	//leg_entry0->SetTextColor(4);
	//leg->Draw();

  TPaveText *oneSig = new TPaveText(0.65,0.25,0.80,0.30,"NDC");
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
  c2->Print("thw_deltaNLL.png");
}
