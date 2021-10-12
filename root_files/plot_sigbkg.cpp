//plot some distributions for signal and background in normalized graphs 
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

void plot_sigbkg() {
  	TString fileNames[] = {	"new_DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM.root",
  							"new_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM.root",
  							"new_DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM.root",
  							"new_DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM.root",
  							"new_DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM.root",
  							"new_DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM.root",
							"new_DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM.root",														
							"new_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM.root",
							"new_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root",
							//"new_ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo.root",
							//"new_ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin.root",
							//"new_ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin.root",
							"new_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg.root",
							"new_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg.root",																																																						
							"new_QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							//"new_QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							//"new_TTTo2L2Nu_TuneCP5_13TeV-powheg.root",
							//"new_TTToHadronic_TuneCP5_13TeV-powheg.root",							
							//"new_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM.root",							
							//"new_WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM.root",							
							"new_WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM.root",
							"new_WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM.root",
							"new_WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM.root",
							"new_WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM.root",
							"new_WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM.root",
							"new_WW_TuneCP5_13TeV.root",
							"new_WZ_TuneCP5_13TeV.root",
							"new_ZZ_TuneCP5_13TeV.root"};
  	//TString legend[] = {"wrong reco", "non reco", "background", "right reco"};
	float weight[]={1.37,0.771,0.305,0.096,0.027,0.014,0.0,0.001,
			96.738,
			0.604,0.509,
		199211.838,7118.803,76.752,19.688,10.556,1.549,3.578,
		0.982,0.223,0.103,0.026,0.001,
		1.096,0.579,0.932};
  	TString xvars[] = {"rapidity_tt", "mass_tt","likelihood","rectop_pt","neutrino_pz","jet_num","mass_thad","mass_tlep","jet_pt","jet_pt[0]"};
  	TString xtitle[] = {"#Deltay_{t#bar{t}}", "M_{t#bar{t}} [GeV]","likelihood","P_{T}(t) [GeV]","P_{z}^{#nu} [GeV]","njet","M_{thad}","M_{tlep}","jet_P_{T}","leading_jet_P_{T}"};
  	TString title[] = {"detlay", "mtt","likelihood","pt","neutrino_pz","jet_num","mass_thad","mass_tlep","jet_pt","leading_jet_pt"};
  	TString legend[]={"DY","tt_semi","single_top","QCD","V+jets","VV","V+jets","VV"};
  	int edge[]={0,8,9,11,18,23,26,31,34};
  	float xlow[] = {-5.0, 0.0, 17.0, 0.0, -1100.0,4,50,100,0,0};
  	float xup[] = {5.0, 1500.0, 20.0, 1000.0, 1100.0,15,300,300,600,600};
  	int color[] = {2, 1, 4, 226, 6, kOrange + 2, 3, kYellow, 93};
  	int xbins[] = {40, 50, 60, 50,110,11,25,20,50,50};
  	float max;
  //	TString seletName[] = {"wrong_reco", "non_reco", "background", "right_reco"};
  	for (int i = 0; i < 10; i++) {                         // loop over variables
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
    	max=0;
    	//find max
    	for (int k = 0; k < 6; k++) { 
      			TH1F *h2 = new TH1F("h2", "", xbins[i], xlow[i], xup[i]);
      			for(int j=edge[k];j<edge[k+1];j++){// loop over files
        		TFile *file = TFile::Open(fileNames[j]);
        				TH1F* hist=new TH1F("hist", "", xbins[i], xlow[i], xup[i]);
          			TTree *tree = (TTree *)file->Get("mytree"); // draw from tree
          			c1->cd();
          			tree->Draw(xvars[i] + ">>hist");
       			if(j==edge[k]){
       					h2=hist;
       					h2->Scale(weight[j]);
       			}
       			else{
       					h2->Add(hist,weight[j]);
       			}
       		}
       		if(max<h2->GetMaximum()/h2->Integral())
       			max=h2->GetMaximum()/h2->Integral();
       	}
      	for (int k = 0; k < 6; k++) { 
      		TH1F *h1 = new TH1F("h1", "", xbins[i], xlow[i], xup[i]);
      		for(int j=edge[k];j<edge[k+1];j++){// loop over files
        		TFile *file = TFile::Open(fileNames[j]);
        		TH1F* hist=new TH1F("hist", "", xbins[i], xlow[i], xup[i]);
          		TTree *tree = (TTree *)file->Get("mytree"); // draw from tree
          		c1->cd();
          		tree->Draw(xvars[i] + ">>hist");
       			if(j==edge[k]){
       				h1=hist;
       				h1->Scale(weight[j]);
       			}
       			else{
       					h1->Add(hist,weight[j]);
       			}
       		}
       		h1->Scale(1/h1->Integral());
       		
       		h1->SetLineColor(color[k]);
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
	       		//h1->GetYaxis()->SetRangeUser(0., max*1.4);
	        h1->GetXaxis()->SetRangeUser(xlow[i], xup[i]);
			c2->cd();
	        if (k < 1) {
	          	h1->Draw("hist");
			}else{
	          	h1->Draw("histsame");
	        }
	        h1->GetYaxis()->SetRangeUser(0., max*1.4);
	        TLegendEntry *leg_entry0 = leg->AddEntry(h1, legend[k]);
	        leg_entry0->SetTextColor(color[k]);

      	}	
	    c2->cd();
	    leg->Draw("same");
	    c2->Print(title[i] + ".png");
	//  c2->Print(title[i] + ".pdf");
  	}
}