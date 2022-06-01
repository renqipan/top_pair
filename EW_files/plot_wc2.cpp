void format_canvas(TCanvas* c){
		c->SetFillColor(0);
		c->SetBorderMode(0);
		c->SetBorderSize(2);
		c->SetTickx(1);
		c->SetTicky(1);
		c->SetLeftMargin(0.17);
		c->SetRightMargin(0.12);
		c->SetTopMargin(0.07);
		c->SetBottomMargin(0.13);
		c->SetFrameFillStyle(0);
		c->SetFrameBorderMode(0);
		c->SetFrameFillStyle(0);
		c->SetFrameBorderMode(0);

	}

void plot_wc2(){
	TString fileName[]={"ci1.00.00.00.0.root","ci2.00.00.00.0.root","ci0.01.00.00.0.root",
						"ci0.02.00.00.0.root","ci0.00.01.00.0.root","ci0.00.02.00.0.root",
					"ci0.00.00.01.0.root","ci0.00.00.02.0.root"};
    TString legName[]={"(1.0, 0.0, 0.0, 0.0)",
    					"(2.0, 0.0, 0.0, 0.0)",
   						"(0.0, 1.0, 0.0, 0.0)",
   						"(0.0, 2.0, 0.0, 0.0)",
   						"(0.0, 0.0, 1.0, 0.0)",
   					    "(0.0, 0.0, 2.0, 0.0)",
   					    "(0.0, 0.0, 0.0, 1.0)",
   						"(0.0, 0.0, 0.0, 2.0)"};
    TString var_ew[]={"id7","id10","id13","id16"};//mtt,pt,ytt, yt
    float xlow[]={350,0,-4.0,-4.0};
	float xup[]={1200,1200,4.0,4.0};
	float ylow[]={-0.6,-0.8,-0.6,-0.6};
	float yup[]={0.3,0.3,0.3,0.3};
	TString xtitle[]={"M_{t#bar{t}} [GeV]","p_{T}^{t} [GeV]","#Deltay_{t#bar{t}}","y_{t}"};
	TString ytitle[]={"(d#delta#sigma_{weak}/dM_{t#bar{t}})/(d#sigma_{LO}/dM_{t#bar{t}})",
	          "(d#delta#sigma_{weak}/dp_{T}^{t})/(d#sigma_{LO}/dp_{T}^{t})",
				"(d#delta#sigma_{weak}/d#Deltay_{t#bar{t}})/(d#sigma_{LO}/d#Deltay_{t#bar{t}})",
				"(d#delta#sigma_{weak}/dy_{t})/(d#sigma_{LO}/dy_{t})"};
	TString title[]={"dist_mtt","dist_pt","dist_ytt","dis_yt"};
	int color[]={2,1,4,226,6,kOrange+2,kViolet+1,kAzure+10,93};
	int bins=50;

	for(int i=0;i<4;i++){ //loop over variables
		auto c2=new TCanvas("c2","",30, 30, 800, 800); //draw on this canvas
		format_canvas(c2);
		c2->cd();
		auto leg=new TLegend(.28,.72,.85,.9,"#frac{v^{2}}{#Lambda^{2}}(C^{#varphiq3}, C^{#varphiu}, Re[C^{u#varphi}], Im[C^{u#varphi}])");	
		leg->SetFillColor(0);
		leg->SetLineColor(0);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.021);
		leg->SetNColumns(2);
			
		for(int k=0;k<8;k++){ //loop over files 
		//TH1F *ew_hist = (TH1F*) file->Get(varName[i]);
		TFile* file=TFile::Open(fileName[k]);
		TH1F *ew_hist = (TH1F*) file->Get(var_ew[i]);
		ew_hist->Smooth();
		int nbins=ew_hist->GetNbinsX();
		double x[nbins],y[nbins];
		for(int ii=0;ii<nbins;ii++){
			x[ii]=ew_hist->GetXaxis()->GetBinCenter(ii+1);
			y[ii]=ew_hist->GetBinContent(ii+1);
		//	cout<<x[ii]<<" "<<y[ii]<<endl;
		}
		TGraph* h1=new TGraph(nbins,x,y);
	
		//TH1F* h1=(TH1F*)ew_hist->Clone();
		h1->SetLineColor(color[k]);
		h1->SetLineWidth(3);
		//h1->SetStats(kFALSE);
		h1->GetXaxis()->SetTitle(xtitle[i]);
		h1->GetYaxis()->SetTitle(ytitle[i]);
		h1->SetTitle("");
		h1->GetXaxis()->CenterTitle();
		h1->GetYaxis()->CenterTitle();
		h1->GetXaxis()->SetTitleSize(0.05);
		h1->GetYaxis()->SetTitleSize(0.05);
		h1->GetYaxis()->SetTitleOffset(1.4);
		h1->GetYaxis()->SetRangeUser(ylow[i],yup[i]);
		h1->GetXaxis()->SetRangeUser(xlow[i],xup[i]);

		c2->cd();
		if(k<1)
			h1->Draw("AC");
		else
			h1->Draw("Csame");
		TLegendEntry* leg_entry0=leg->AddEntry(h1,legName[k]);
		leg_entry0->SetTextColor(color[k]);
		}
		c2->cd();
		leg->Draw("same");
	    //c2->Print(title[i]+"_wc.png");
	    c2->Print(title[i]+"_wc2.pdf");

	}
	
        
}
