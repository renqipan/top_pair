void format_canvas(TCanvas* c){
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

void plot_2Dtt(){
	TString fileNames[]={"new_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root"};
    TString legend[]={"SM: t#bar{t}","c_{G}=0,c_{tG}=0","c_{G}=1,c_{tG}=0","c_{G}=0,c_{tG}=1","SM: NLO","Mad: NLO","SM: NLO"};
	TString xvars[]={"mass_thad","mass_tlep"};
	TString yvars[]={"mass_whad","mass_wlep"};
	TString xtitle[]={"M_{t_{had}} [GeV]","M_{t_{lep}} [GeV]"};
	TString ytitle[]={"M_{w_{had}} [GeV]","M_{w_{lep}} [GeV]"};
	TString title[]={"mtwhad","mtwlep"};
	float xlow[]={0.0,0.0};
	float xup[]={600.0,600.0,};
	float ylow[]={0.0,0.0};
	float yup[]={600,600};
	int color[]={2,1,4,226,6,kOrange+2,3,kYellow,93};
	int xbins=50,ybins=40;
	for(int i=0;i<2;i++){ //loop over variables
		auto c1=new TCanvas("c1","",8, 30, 700, 700);//temporary canvas
		auto c2=new TCanvas("c2","",8, 30, 700, 700); //draw on this canvas
		format_canvas(c2);
		c2->cd();
		auto leg=new TLegend(.65,.65,.85,.95);	
		leg->SetFillColor(0);
		leg->SetLineColor(0);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.052);
			
		for(int k=0;k<1;k++){ //loop over files 
			TFile* file=TFile::Open(fileNames[k]);
			TH2F *h1=new TH2F("h1","",xbins,xlow[i],xup[i],ybins,ylow[i],yup[i]);
			
			if(k<1){	
				TTree* tree=(TTree*)file->Get("mytree");// draw from tree
				c1->cd();
				tree->Draw(yvars[i]+":"+xvars[i]+">>h1");
			}
			
			h1->SetLineColor(kBlue);
			h1->SetLineWidth(2);
			h1->SetStats(kFALSE);
			h1->GetXaxis()->SetTitle(xtitle[i]);
			h1->GetYaxis()->SetTitle(ytitle[i]);
			h1->SetTitle("");
			h1->GetXaxis()->CenterTitle();
			h1->GetYaxis()->CenterTitle();
			h1->GetXaxis()->SetTitleSize(0.06);
			h1->GetYaxis()->SetTitleSize(0.06);
			h1->GetXaxis()->SetTitleOffset(1.5);
			h1->GetYaxis()->SetTitleOffset(1.6);
			h1->GetXaxis()->SetLabelSize(0.05);
			h1->GetYaxis()->SetLabelSize(0.05);
			//h1->GetYaxis()->SetRangeUser(0.,h1->GetMaximum()*1.3);
			//h1->GetXaxis()->SetRangeUser(xlow[i],xup[i]);

			c2->cd();
			if(k<1)
				{h1->DrawNormalized("surf2");
				}
			else
				{h1->DrawNormalized("samehist");
			}
			//TLegendEntry* leg_entry0=leg->AddEntry(h1,legend[k]);
			//leg_entry0->SetTextColor(color[k]);
		}
		c2->cd();
		leg->Draw("same");
	   c2->Print(title[i]+"_surf2.png");
	   c2->Print(title[i]+"_surf2.pdf");

	}
	
        
}
