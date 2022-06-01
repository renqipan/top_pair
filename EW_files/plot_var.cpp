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

void plot_var(){
	TString fileName[]={"ci0.00.01.00.0.root","ci0.00.00.00.0.root","ci0.00.01.01.0.root"};
    TString legName[]={"#frac{v^{2}}{#Lambda^{2}}Re[C^{u#varphi}]=1(#kappa=0)","#frac{v^{2}}{#Lambda^{2}}Re[C^{u#varphi}]=0(#kappa=1)",
         "#frac{v^{2}}{#Lambda^{2}}Im[C^{u#varphi}]=1(#tilde{#kappa}=1)"};
    //TString var_lo[]={"id4","id7","id10"};//mtt,pt,ytt
  //  TString var_ew[]={"id5","id8","id11"};//mtt,pt,ytt
    TString varName[]={"id7","id10","id13"};//mtt,pt,ytt
    float xlow[]={350,0,-5.0};
	float xup[]={1200,1000,5.0};
	float ylow[]={-0.15,-0.15,-0.08};
	float yup[]={0.15,0.15,0.10};
	TString xtitle[]={"M_{t#bar{t}} [GeV]","p_{T}^{t} [GeV]","#Deltay_{t#bar{t}}"};
	TString ytitle[]={"(d#delta#sigma_{weak}/dM_{t#bar{t}})/(d#sigma_{LO}/dM_{t#bar{t}})",
	          "(d#delta#sigma_{weak}/dp_{T}^{t})/(d#sigma_{LO}/dp_{T}^{t})",
				"(d#delta#sigma_{weak}/d#Deltay_{t#bar{t}})/(d#sigma_{LO}/d#Deltay_{t#bar{t}})"};
	TString title[]={"dist_mtt","dist_pt","dist_ytt"};
	int color[]={2,1,4,226,6,kOrange+2,3,kYellow,93};
	//int color[]={2,1,4,226,93,kOrange+2,3,kYellow,93};
	int bins=50;

	for(int i=0;i<3;i++){ //loop over variables
		//auto c1=new TCanvas("c1","",8, 30, 700, 700);//temporary canvas

		auto c2=new TCanvas("c2","",30, 30, 800, 800); //draw on this canvas
		format_canvas(c2);
		c2->cd();
		auto leg=new TLegend(.50,.62,.75,.9);	
		leg->SetFillColor(0);
		leg->SetLineColor(0);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.04);
			
		for(int k=0;k<3;k++){ //loop over files 
		TFile* file=TFile::Open(fileName[k]);
		TH1F *h1 = (TH1F*) file->Get(varName[i]);

		h1->SetLineColor(color[k]);
		h1->SetLineWidth(2);
		h1->SetStats(kFALSE);
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
			h1->Draw("hist");
		else
			h1->Draw("samehist");
		TLegendEntry* leg_entry0=leg->AddEntry(h1,legName[k]);
		leg_entry0->SetTextColor(color[k]);
		}
		c2->cd();
		leg->Draw("same");
	    //c2->Print(title[i]+"_veri.png");
	    c2->Print(title[i]+"_veri.pdf");

	}
	
        
}
