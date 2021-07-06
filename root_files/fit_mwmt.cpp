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
void format_hist(TH1F* h1){
	h1->SetLineWidth(2);
	h1->SetTitle("");
	h1->GetXaxis()->CenterTitle();
	h1->GetYaxis()->CenterTitle();
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetXaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetLabelSize(0.05);
	h1->GetYaxis()->SetRangeUser(0.,h1->GetMaximum()*1.3);	
}

void fit_mwmt(){
	gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetOptFit(0111);
    gStyle->SetStatX(.92); gStyle->SetStatY(.93);
    auto c1= new TCanvas("c1","c1");
	auto c2= new TCanvas("c2","c2");
	auto c3= new TCanvas("c3","c3");
	auto c4= new TCanvas("c4","c4");
	format_canvas(c1);
	format_canvas(c2);
	format_canvas(c3);
	format_canvas(c4);

	float xi_thad=18.0,x0_thad=179,xi_wlep=2.0,x0_wlep=80,xi_tlep=8.5,x0_tlep=169,xi_whad=14.0,x0_whad=84;
    float mw_lep=85,mt_had=193, mt_lep=178,mw_had=94.0,sigmaw_lep=5.5,sigmat_lep=21,sigmaw_had=28.0,sigmat_had=32,rho=0.3;
	TString inputFile="new_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root";
	TFile *file=TFile::Open(inputFile);
	TTree *mytree=(TTree*) file->Get("mytree");
	TH1F *h1=new TH1F("h1","mt_had",80,0,1200);//for mt_had
	TH1F *h2=new TH1F("h2","mw_lep",30,0,300);//for mw_lep
	TH1F *h3=new TH1F("h3","mt_lep",60,0,800);//for mt_lep
	TH1F *h4=new TH1F("h4","mw_had",40,0,600);//for mw_had
	mytree->Draw("mass_thad>>h1");
	mytree->Draw("mass_wlep>>h2");
	mytree->Draw("mass_tlep>>h3");
	mytree->Draw("mass_whad>>h4");
	float scale1=1.0/h1->Integral("width");
	h1->Scale(scale1);
	float scale2=1.0/h2->Integral("width");
	h2->Scale(scale2);
	float scale3=1.0/h3->Integral("width");
	h3->Scale(scale3);
	float scale4=1.0/h4->Integral("width");
	h4->Scale(scale4);
	//ROOT::Math::landau_pdf (double x, double xi=1, double x0=0)
	bool landau_flag=false;
	bool gaus_flag=true;
	TF1 *f1, *f2, *f3, *f4;
	if(landau_flag){
		f1=new TF1("f1","ROOT::Math::landau_pdf(x,[0],[1])",0,1200);
		f2=new TF1("f2","ROOT::Math::landau_pdf(x,[0],[1])",0,300);
		f3=new TF1("f3","ROOT::Math::landau_pdf(x,[0],[1])",0,800);
	    f4=new TF1("f4","ROOT::Math::landau_pdf(x,[0],[1])",0,600);
		f1->SetParNames( "#xi","x0");
		f1->SetParameters(xi_thad,x0_thad);//for t_had
		f2->SetParNames("#xi","x0");
		f2->SetParameters(xi_wlep,x0_wlep);//for w_lep
		f3->SetParNames("#xi","x0");
		f3->SetParameters(xi_tlep,x0_tlep); //for t_lep
		f4->SetParNames("#xi","x0");
		f4->SetParameters(xi_whad,x0_whad); //for w_had
		
		
		h1->Fit("f1");
		h2->Fit("f2");
		h3->Fit("f3");
		h4->Fit("f4");

		xi_thad=f1->GetParameter(0);
		x0_thad=f1->GetParameter(1);
		xi_wlep=f2->GetParameter(0);
		x0_wlep=f2->GetParameter(1);
		xi_tlep=f3->GetParameter(0);
		x0_tlep=f3->GetParameter(1);
		xi_whad=f4->GetParameter(0);
		x0_whad=f4->GetParameter(1);
		cout<<"-------------------------------------------------------"<<endl;
		cout<<" xi_thad: "<<xi_thad<<" x0_thad: "<<x0_thad<<endl;
		cout<<" xi_wlep: "<<xi_wlep<<" x0_wlep: "<<x0_wlep<<endl;
		cout<<" xi_tlep: "<<xi_tlep<<" x0_tlep: "<<x0_tlep<<endl;
		cout<<" xi_whad: "<<xi_whad<<" x0_whad: "<<x0_whad<<endl;
		cout<<"-------------------------------------------------------"<<endl;

	}
////////////////////////////////////////////////////////////////////////
	if(gaus_flag){
		//ROOT::Math::gaussian_pdf (double x, double sigma=1, double x0=0)
		f1=new TF1("f1","ROOT::Math::gaussian_pdf (x,[0],[1])",0,1200);
		f2=new TF1("f2","ROOT::Math::gaussian_pdf (x,[0],[1])",0,300);
		f3=new TF1("f3","ROOT::Math::gaussian_pdf (x,[0],[1])",0,800);
	    f4=new TF1("f4","ROOT::Math::gaussian_pdf (x,[0],[1])",0,600);
		f1->SetParNames( "#sigma_{thad}","#mu_{thad}");
		f1->SetParameters(sigmat_had,mt_had);//for t_had
		f2->SetParNames("#sigma_{wlep}","#mu_{wlep}");
		f2->SetParameters(sigmaw_lep,mw_lep);//for w_lep
		f3->SetParNames("#sigma_{tlep}","#mu_{tlep}");
		f3->SetParameters(sigmat_lep,mt_lep); //for t_lep
		f4->SetParNames("#sigma_{whad}","#mu_{whad}");
		f4->SetParameters(sigmaw_had,mw_had); //for w_had
		
		
		h1->Fit("f1");
		h2->Fit("f2");
		h3->Fit("f3");
		h4->Fit("f4");

		sigmat_had=f1->GetParameter(0);
		mt_had=f1->GetParameter(1);
		sigmaw_lep=f2->GetParameter(0);
		mw_lep=f2->GetParameter(1);
		sigmat_lep=f3->GetParameter(0);
		mt_lep=f3->GetParameter(1);
		sigmaw_had=f4->GetParameter(0);
		mw_had=f4->GetParameter(1);
		cout<<"-------------------------------------------------------"<<endl;
		cout<<" sigmat_had: "<<sigmat_had<<" mt_had: "<<mt_had<<endl;
		cout<<" sigmaw_lep: "<<sigmaw_lep<<" mw_lep: "<<mw_lep<<endl;
		cout<<" sigmat_lep: "<<sigmat_lep<<" mt_lep: "<<mt_lep<<endl;
		cout<<" sigmaw_had: "<<sigmaw_had<<" mw_had: "<<mw_had<<endl;
		cout<<"-------------------------------------------------------"<<endl;


	}

		h1->GetXaxis()->SetTitle("m_{t_{had}}");
		h2->GetXaxis()->SetTitle("m_{W_{lep}}");
		h3->GetXaxis()->SetTitle("m_{t_{lep}}");
		h4->GetXaxis()->SetTitle("m_{W_{had}}");
		format_hist(h1);
		format_hist(h2);
		format_hist(h3);
		format_hist(h4);
		c1->cd();	
		h1->Draw("hist");
		f1->Draw("same");
		c2->cd();
		h2->Draw("hist");
		f2->Draw("same");
		c3->cd();
		h3->Draw("hist");
		f3->Draw("same");
		c4->cd();
		h4->Draw("hist");
		f4->Draw("same");
		
		c1->Print("mass_thad.pdf");
		c2->Print("mass_wlep.pdf");
		c3->Print("mass_tlep.pdf");
		c4->Print("mass_whad.pdf");
	

}

	/*TPaveStats *st = (TPaveStats*)h2->FindObject("stats");
		st->SetX1NDC(0.65); //new x start position
		st->SetX2NDC(0.93); //new x end position
		st->SetY1NDC(0.6);
		st->SetY2NDC(0.95);
	*/