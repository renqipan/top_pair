void format_canvas(TCanvas *c) {
  	c->SetFillColor(0);
 	c->SetBorderMode(0);
  	c->SetBorderSize(2);
  	c->SetTickx(1);
  	c->SetTicky(1);
  	c->SetLeftMargin(0.18);
  	c->SetRightMargin(0.05);
  	c->SetTopMargin(0.03);
  	c->SetBottomMargin(0.15);
  	c->SetFrameFillStyle(0);
  	c->SetFrameBorderMode(0);
  	c->SetFrameFillStyle(0);
  	c->SetFrameBorderMode(0);
}
void plot_signal() {
  	TString fileNames[] = {"new_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root"};
  	TString legend[] = {"wrong reco", "non reco", "background", "right reco"};
  	TString xvars[] = {"rapidity_tt", "mass_tt","likelihood","rectop_pt","neutrino_pz"};
  	TString xtitle[] = {"#Deltay_{t#bar{t}}", "M_{t#bar{t}} [GeV]","likelihood","P_{T}(t) [GeV]","P_{z}^{#nu} [GeV]"};
  	TString title[] = {"detlay", "mtt","likelihood","pt","neutrino_pz"};
  	float xlow[] = {-5.0, 0.0,10.0,0.0,-1100.0};
  	float xup[] = {5.0, 2000.0,30.0,1000.0,1100.0};
  	int color[] = {2, 1, 4, 226, 6, kOrange + 2, 3, kYellow, 93};
  	int xbins[] = {40, 50, 30, 50,110};
  	TString cuts[] = {"tt_efficiency==2||tt_efficiency==4||tt_efficiency==5",
                    "tt_efficiency==1", "tt_efficiency==0", "tt_efficiency==3"};
  //	TString seletName[] = {"wrong_reco", "non_reco", "background", "right_reco"};
  	for (int i = 0; i < 5; i++) {                         // loop over variables
   	 	auto c1 = new TCanvas("c1", "c1", 8, 30, 600, 600); // temporary canvas
    	auto c2 = new TCanvas("c2", "c2", 8, 30, 600, 600); // draw on this canvas
   		format_canvas(c2);
    	c2->cd();
    	auto leg = new TLegend(.65, .65, .85, .95);
    	leg->SetFillColor(0);
   		leg->SetLineColor(0);
    	leg->SetBorderSize(0);
    	leg->SetFillStyle(0);
    	leg->SetTextFont(42);
    	leg->SetTextSize(0.052);
    /*	double range=0;
    	for (int m = 0; m < 4; m++) {   // select: all, right, wrong
        	TFile *file = TFile::Open(fileNames[0]);
        	TH1F *h2 = new TH1F("h2", "", xbins[i], xlow[i], xup[i]);
          	TTree *tree = (TTree *)file->Get("mytree"); // draw from tree
          	c1->cd();
          	tree->Draw(xvars[i] + ">>h2", cuts[m]);
       		if(i==0){
       			range=h2->GetMaximum();
       		}
       		else if(h2->GetMaximum()>range){
       			range=h2->GetMaximum();
       		}
       	}*/
    	for (int m = 0; m < 4; m++) {   // select: all, right, wrong
      		for (int k = 0; k < 1; k++) { // loop over files
        		TFile *file = TFile::Open(fileNames[k]);
        		TH1F *h1 = new TH1F("h1", "", xbins[i], xlow[i], xup[i]);
        		
          			TTree *tree = (TTree *)file->Get("mytree"); // draw from tree
          			c1->cd();
          			tree->Draw(xvars[i] + ">>h1", cuts[m]);
       			
       			h1->SetLineColor(color[m]);
        		h1->SetLineWidth(2);
        		h1->SetStats(kFALSE);
        		h1->GetXaxis()->SetTitle(xtitle[i]);
        		h1->GetYaxis()->SetTitle("Normalized");
        		h1->SetTitle("");
        		h1->GetXaxis()->CenterTitle();
        		h1->GetYaxis()->CenterTitle();
        		h1->GetXaxis()->SetTitleSize(0.06);
       			h1->GetYaxis()->SetTitleSize(0.06);
        		h1->GetXaxis()->SetTitleOffset(0.9);
	        	h1->GetYaxis()->SetTitleOffset(1.3);
	        	h1->GetXaxis()->SetLabelSize(0.05);
	        	h1->GetYaxis()->SetLabelSize(0.05);
	        	if(i==2){
	        		h1->GetYaxis()->SetRangeUser(0., h1->GetMaximum() * 2.5);
	        	}
	        	else{
	        		h1->GetYaxis()->SetRangeUser(0., h1->GetMaximum() * 1.3);
	        	}
	        	h1->GetXaxis()->SetRangeUser(xlow[i], xup[i]);
				c2->cd();
	        	if (m < 1) {
	          		h1->DrawNormalized("hist");
				} else
	          		h1->DrawNormalized("histsame");
	        	TLegendEntry *leg_entry0 = leg->AddEntry(h1, legend[m]);
	        	leg_entry0->SetTextColor(color[m]);
      		}
    	}
	    c2->cd();
	    leg->Draw("same");
	    c2->Print(title[i] + ".png");
	//    c2->Print(title[i] + ".pdf");
  	}
}
