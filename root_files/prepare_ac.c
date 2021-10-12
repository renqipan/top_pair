
using namespace std;
using namespace RooFit;
const int extraSysn_all = 3;
const int extraCat= 2;
const int extraSysn_zjet = 3;
TString extraSysName[extraCat]={"jec","btag"};
TString extraSys_handName_zjet [extraSysn_zjet]={"zjet_4mu","zjet_4e","zjet_2e2mu"};
TString extraSys_handName_all[extraSysn_all]={"lumi_13TeV","CMS_eff_m","CMS_eff_e"};
float extraSys_zjet_up [3][extraSysn_zjet]= {
	1.104, 1.314,    1.152,
	1.32,	1.38,	1.33,  
	1.30,	1.37,	1.24 
};

float extraSys_zjet_dn [3][extraSysn_zjet]= {
	0.899,   0.728,      0.868,
	0.68, 	0.64,	0.67,
	0.70,   0.63,      0.76
};

float extraSys_all_up [3][extraSysn_all][3]= {
	{
		1.026,  1.026,  1.026,
		1.046,  1.,     1.025,
		1.,     1.082,  1.039
	},
	{
		1.023,	1.023,	1.023,
		1.056,	1.,	1.03,
		1.,	1.125,	1.058
	},
	{
		1.025,	1.025,	1.025,
		1.016,	1.,	1.011,
		1.,	1.161,	1.074
	}
};
float extraSys_all_dn [3][extraSysn_all][3]= {
	{
		0.974,  0.974,  0.974,
		0.953,  1.,     0.975,
		1.,     0.914,  0.96
	},
	{
		0.977,	0.977,	0.977,
		0.937,	1.,	0.968,
		1.,	0.862,	0.939
	},
	{
		0.975,	0.975,	0.975,
		0.978,	1.,	0.992,
		1.,	0.850,	0.928
	}
};

void setlabel(vector<TString> xlabel,vector<TString> ylabel,TH2F *h2yield){
	for(int i=0;i<xlabel.size();i++){
		h2yield->GetXaxis()->SetBinLabel(i+1,xlabel[i]);
		for (int j=0;j<h2yield->GetNbinsY();j++)
			h2yield->GetYaxis()->SetBinLabel(j+1,ylabel[j]);
	}
}
void setlabel(vector<TString> xlabel,TH2F *h2yield){
	for(int i=0;i<xlabel.size();i++){
		h2yield->GetXaxis()->SetBinLabel(i+1,xlabel[i]);
	}
}

void Floor(TH3F* histo){
	//TH3F *histo_new = (TH3F*)histo->Clone(Form("%s_new",histo->GetName())); 
	for (int i=0;i<histo->GetNbinsX();i++){
		for (int j=0;j<histo->GetNbinsY();j++){
			for (int k=0;j<histo->GetNbinsZ();k++){
				if(histo->GetBinContent(i+1,j+1,k+1)==0)
					histo->SetBinContent(i+1,j+1,k+1,1.E-10);
			}
		}
	}
}

void Floor(TH2F* histo){
	for (int i=0;i<histo->GetNbinsX();i++){
		for (int j=0;j<histo->GetNbinsY();j++){
			if(histo->GetBinContent(i+1,j+1)==0)
				histo->SetBinContent(i+1,j+1,1.E-10);
		}
	}
}

void Floor(TH1F* histo){
	//	histo->Smooth();
	for (int i=0;i<histo->GetNbinsX();i++){
		if(histo->GetBinContent(i+1)==0)
			histo->SetBinContent(i+1,1.E-10);
	}
}
void writeline(vector<TString> arr , ofstream &card){
	for (int i=0;i<arr.size();i++)
		card<< arr[i]<<"\t";
	card<<endl;
}

void writeline(vector<int> arr , ofstream &card){
	for (int i=0;i<arr.size();i++)
		card<< arr[i]<<"\t";
	card<<endl;
}
void writeline(vector<float> arr , ofstream &card){
	//	gStyle->SetPaintTextFormat("2.2f");
	card<< std::fixed;
	card<< std::setprecision(3);
	for (int i=0;i<arr.size();i++)
		card<< arr[i]<<"\t";
	card<<endl;
}
void prepare_ac(TString category="htxs_stage1_reco_cat", TString categoryName="htxs_stage1_reco_catName",TString outputDir="default",TString file_app = "_newcate",int doSys=1, int readPara=0,TString stageName = "htxs_stage1_red_cat", TString year="2018", int flav=0){
	gStyle->SetOptStat(0);
	TString chanName[3] = {"4mu" , "4e", "2e2mu"}; 
	float lumi = 41.5;
	int yearn = 1;
	if (year=="2016"){
		lumi = 35.9;
		yearn=0;
	}
	else if (year=="2018"){
		lumi = 59.7;
		yearn=2;
	}
	TString inputDir = "/eos/home-x/xiaomeng/ac_SIP"; 
	TChain *tsys= new TChain ("SelectedTree");
	TChain *tsig = new TChain ("SelectedTree");
	TChain *tggh = new TChain ("SelectedTree");
	TChain *tzx = new TChain ("SelectedTree");
	TChain *tqqzz = new TChain ("SelectedTree");
	TChain *tggzz = new TChain ("SelectedTree");
	TChain *tdata = new TChain("SelectedTree");
	tdata->Add(inputDir+"/AllData"+file_app+"_"+year+"_bdt.root");
	tzx->Add(inputDir+"/AllData_ZX"+file_app+"*"+year+"_bdt.root");
	tzx->Add(inputDir+"/ZZTo4lext"+file_app+"*"+year+"_bdt.root");
	tzx->Add(inputDir+"/ggTo*"+file_app+"*"+year+"_bdt.root");
	tsys->Add(inputDir+"/ggH125_minloHJJ"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/VBFH125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ggH0PM_M125_mcanlo"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ggH0M_M125_mcanlo"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ggH0Mf05_M125_mcanlo"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/WplusH125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/WminusH125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ZH125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ggZH125"+file_app+"*"+year+"_bdt.root");
	//	tsig->Add(inputDir+"/ttH0PM_M125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ttH125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/ttH0M_M125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/tqH125"+file_app+"*"+year+"_bdt.root");
	tsig->Add(inputDir+"/bbH125"+file_app+"*"+year+"_bdt.root");
	tggh->Add(inputDir+"/ggH125"+file_app+"*"+year+"_bdt.root");

	RooRealVar* dbkg=new RooRealVar("dbkg","",0,1);
	RooRealVar* dbkg_2j=new RooRealVar("dbkg_2j","",0,1);
	RooRealVar* BDTG=new RooRealVar("BDTG","",-1,1);
	RooRealVar* DVHDEC=new RooRealVar("DVHDEC","",0,1);
	RooRealVar* d_2j_qq=new RooRealVar("d_2j_qq","",0,1);
	//	RooRealVar* d_2j_int=new RooRealVar("d_2j_int","",-1,1);
	RooRealVar* D2jet=new RooRealVar("D2jet","",0,1);

	dbkg->setBins(20);
	DVHDEC->setBins(20);
	BDTG->setBins(4);
	d_2j_qq->setBins(4);
	//	d_2j_int->setBins(2);
	D2jet->setBins(4);

	//tsig->Draw(Form("%s:%s>>h_reco",category.Data(),categoryName.Data()),"chan==3");
	tsig->Draw(Form("%s>>h_reco",categoryName.Data()),"chan==3");
	TH1F *h=(TH1F*)gROOT->FindObject("h_reco")->Clone();
	int reco_cat_all = h->GetNbinsX();
	vector<int> reco_cat_num;
	vector<TString> ylabel;
	cout<< "reconstructed categories "<<reco_cat_all <<endl;
	for(int i=0;i<reco_cat_all;i++){
		cout<< h->GetXaxis()->GetBinLabel(i+1)<<endl;
		TString reco_name = h->GetXaxis()->GetBinLabel(i+1);
		ylabel.push_back(reco_name);
		tsig->Draw(Form("%s",category.Data()),Form("chan==3&&%s==\"%s\"",categoryName.Data(),reco_name.Data()));
		int reco_cat = *(tsig->GetV1());
		reco_cat_num.push_back(reco_cat);
	}

	tsig->Draw(Form("%sName>>h_stage1",stageName.Data()),Form("%s>0 && %s%100!=0",stageName.Data(),stageName.Data()));
	TH1F *h_stage1=(TH1F*)gROOT->FindObject("h_stage1")->Clone();
	const int stage1_cat = h_stage1->GetNbinsX();
	TTree *tsig_prod [stage1_cat+3];
	cout<< "stage1 categories "<< stage1_cat <<endl;
	vector<TString> xlabel; 
	TH2F *h2bkg = new TH2F("h2bkg","",3,0,3,20,0,1);
	tzx->Draw("dbkg:htxs_stage1_red_catName>>h2bkg",Form("weight*%f/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,lumi));
	for(int i=0;i<stage1_cat+3;i++){
		if(i<stage1_cat){
			cout<< h_stage1->GetXaxis()->GetBinLabel(i+1)<<endl;
			xlabel.push_back( h_stage1->GetXaxis()->GetBinLabel(i+1));
			tsig_prod[i] = (TTree*)tsig->CopyTree(Form("%sName==\"%s\"",stageName.Data(), h_stage1->GetXaxis()->GetBinLabel(i+1)));
		}
		else{
			cout<< h2bkg->GetXaxis()->GetBinLabel(i-stage1_cat+1)<<endl;
			xlabel.push_back( h2bkg->GetXaxis()->GetBinLabel(i-stage1_cat+1));
			tsig_prod[i] = (TTree*)tzx->CopyTree(Form("htxs_stage1_red_catName==\"%s\"", h2bkg->GetXaxis()->GetBinLabel(i-stage1_cat+1)));
		}
	}
	gStyle->SetPaintTextFormat(".3f");
	//setColZGradient_OneColor(5);


	TCanvas *c=new TCanvas ("c","",1200,600);
	TH2F *h2 = new TH2F("h2","",xlabel.size(),0,xlabel.size(),20,0,1);
	TH2F *h2yield = new TH2F("h2yield","",xlabel.size(),0,xlabel.size(),reco_cat_all,0,reco_cat_all);

	setlabel(xlabel,h2);
	setlabel(xlabel,ylabel,h2yield);
	vector<TH2F> *h2yield_sys_up=new vector<TH2F>;
	vector<TH2F> *h2yield_sys_dn=new vector<TH2F>;


	tsig->Draw(Form("dbkg:%sName>>h2",stageName.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
	//	tsig->Draw(Form("%s:%sName>>h2yield",categoryName.Data(),stageName.Data()),Form("weight*%f*(chan==%d && (abs(JetEta1)<2.5 ||abs(JetEta1)>3)&&(abs(JetEta2)<2.5 ||abs(JetEta2)>3) )/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
	tsig->Draw(Form("%s:%sName>>h2yield",categoryName.Data(),stageName.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
	tzx->Draw("dbkg:htxs_stage1_red_catName>>+h2",Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
	tzx->Draw(Form("%s:htxs_stage1_red_catName>>+h2yield",categoryName.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
	h2->Draw("colztext");
	gPad->Print("stageVsDbkg.png");
	h2yield->Draw("colztext");
	gPad->Print("yield.png");
	vector <TString> sysnames ;
	if(doSys){

		TH2F *h_sys=new TH2F("h_sys","",xlabel.size(),0,xlabel.size(),100,0,100);
		setlabel(xlabel,h_sys);
		//		tsig->Draw("weight_name:"+stageName+"Name>>h_sys","!weight_name.Contains(\"qcd_ggH\")");
		tsig->Draw("weight_name:"+stageName+"Name>>h_sys","!weight_name.Contains(\"QCDscale_muR_ggH\")&&!weight_name.Contains(\"QCDscale_muF_ggH\")");
		tzx->Draw("weight_name:htxs_stage1_red_catName>>+h_sys","!weight_name.Contains(\"_ggZZ\")");
		for(int i=0;i<h_sys->GetNbinsY();i++){
			TString obl = h_sys->GetYaxis()->GetBinLabel(i+1);
			sysnames.push_back(obl);
			cout<<"sys "<<obl<<endl;
		}	       
		h_sys->Draw("colz");
		gPad->Print("systematic.png");
		cout<<"Systematic uncertainties "<<sysnames.size()<<endl;
		//		tsig->Draw("weight_name>>+h_sys");
		for(int i=0;i<sysnames.size()+extraCat;i++){
			TString sysname_this;
			if(i<sysnames.size())
				sysname_this=h_sys->GetYaxis()->GetBinLabel(i+1);
			else
				sysname_this = extraSysName[i-sysnames.size()];
			TH2F *h2yield_sys_uptmp =new TH2F("h2yield_"+  sysname_this+"_up","",xlabel.size(),0,xlabel.size(),reco_cat_all,0,reco_cat_all);
			TH2F *h2yield_sys_dntmp =new TH2F("h2yield_"+  sysname_this+"_dn","",xlabel.size(),0,xlabel.size(),reco_cat_all,0,reco_cat_all);
			setlabel(xlabel,ylabel,h2yield_sys_uptmp);
			setlabel(xlabel,ylabel,h2yield_sys_dntmp);

			cout<<"sys "<<sysname_this<<endl;
			if(i<sysnames.size()){
				for(int istage=0;istage<xlabel.size();istage++){
					if(h_sys->GetBinContent(istage+1,i+1)!=0){
						if(xlabel[istage]=="ggH"){
							tggh->Draw(Form("%s:%sName>>+h2yield_%s_up",categoryName.Data(),stageName.Data(),sysname_this.Data()),Form("weight_up*%f*(chan==%d&&weight_name==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,sysname_this.Data(),lumi));
							tggh->Draw(Form("%s:%sName>>+h2yield_%s_dn",categoryName.Data(),stageName.Data(),sysname_this.Data()),Form("weight_dn*%f*(chan==%d&&weight_name==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,sysname_this.Data(),lumi));
						}
						else{
							tsig_prod[istage]->Draw(Form("%s:%sName>>+h2yield_%s_up",categoryName.Data(),stageName.Data(),sysname_this.Data()),Form("weight_up*%f*(chan==%d&&weight_name==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,sysname_this.Data(),lumi));
							tsig_prod[istage]->Draw(Form("%s:%sName>>+h2yield_%s_dn",categoryName.Data(),stageName.Data(),sysname_this.Data()),Form("weight_dn*%f*(chan==%d&&weight_name==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,sysname_this.Data(),lumi));
						}
						cout<<"per"<<endl;
					}
				}
			}
			else{
				tsig->Draw(Form("%s_%s_up:%sName>>+h2yield_%s_up",categoryName.Data(),sysname_this.Data(),stageName.Data(),sysname_this.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
				tsig->Draw(Form("%s_%s_dn:%sName>>+h2yield_%s_dn",categoryName.Data(),sysname_this.Data(),stageName.Data(),sysname_this.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
				tzx->Draw(Form("%s_%s_up:%sName>>+h2yield_%s_up",categoryName.Data(),sysname_this.Data(),stageName.Data(),sysname_this.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
				tzx->Draw(Form("%s_%s_dn:%sName>>+h2yield_%s_dn",categoryName.Data(),sysname_this.Data(),stageName.Data(),sysname_this.Data()),Form("weight*%f*(chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,flav,lumi));
			}

			h2yield_sys_uptmp->Draw("colztext");
			gPad->Print(sysname_this+".png");
			h2yield_sys_up->push_back(*h2yield_sys_uptmp);
			h2yield_sys_dn->push_back(*h2yield_sys_dntmp);
		}
		for(int i=0;i<extraCat;i++)
			sysnames.push_back(extraSysName[i]);
	}

	const int nsys = sysnames.size(); 
	cout<< "Systematic uncertainties "<<nsys<<endl;

	for (int j =0;j< reco_cat_all;j++){

		TString reco_name = ylabel[j]; 
//		if(!reco_name.Contains("VBF_2j")){
//			continue;
//		}
		ofstream card;
		TString cat_name = reco_name + "_"+chanName[flav-1];
		card.open (outputDir+"/"+cat_name+year+".txt");
		card <<"Datacard for event category: "<< cat_name<<endl;
		card<< "imax 1 number of bins"<<endl;
		card<< "jmax * number of processes minus 1"<<endl;
		card<< "kmax * number of nuisance parameters"<<endl;
		card<<"---------------------------------"<<endl;
		card<<endl;
		card<< "shapes * "<< cat_name << " "<< cat_name+year+".root w:$PROCESS w:$PROCESS_$SYSTEMATIC" <<endl;
		card<<"---------------------------------"<<endl;
		card<< "bin "<< cat_name <<endl;
		card<< "observation -1 "<<endl;
		card<<"---------------------------------"<<endl;
		card<<endl;

		vector<TString> bin_arr;
		vector<TString> processName_arr;
		vector<int> process_arr;
		vector<float> yield_arr;
		vector<vector<float>> yield_arr_up;
		vector<vector<float>> yield_arr_dn;
		vector<TString> sysnames_final;
		for (int isys=0;isys<nsys;isys++)
			sysnames_final.push_back(sysnames.at(isys));
		for (int nextra = 0; nextra < extraSysn_all; nextra++)
			sysnames_final.push_back(extraSys_handName_all[nextra]);
		for (int nextra = flav-1; nextra < flav; nextra++)
			sysnames_final.push_back(extraSys_handName_zjet[nextra]);

		sysnames_final.push_back("hzz_br");
		sysnames_final.push_back("QCDscale_ggZZ");

		TFile *f=new TFile(outputDir+"/"+cat_name+year+".root","recreate");
		RooWorkspace w("w");

		for (int i=0;i<xlabel.size();i++){
			float yield = h2yield->GetBinContent(i+1,j+1);
			TString stage1_name = xlabel[i];
			if(yield<0.005)
				continue;
			bin_arr.push_back(cat_name);
			processName_arr.push_back(stage1_name);
			if(i<stage1_cat)
				process_arr.push_back ( i-stage1_cat);
			else
				process_arr.push_back ( i-stage1_cat+1);
			//				if(stage1_name.Contains("Mix"))
			//					yield=yield*2;
			yield_arr.push_back(yield);
			if(doSys){
				vector<float> yield_arr_up_tmp;
				vector<float> yield_arr_dn_tmp;
				for (int isys=0;isys<nsys;isys++){
					float yield_up = h2yield_sys_up->at(isys).GetBinContent(i+1,j+1);
					float yield_dn = h2yield_sys_dn->at(isys).GetBinContent(i+1,j+1);
					//if(stage1_name.Contains("ggH")&&reco_name.Contains("VBF_2j")){
					//	yield_up=yield_up*0.083;
					//	yield_dn=yield_dn*0.083;
					//}
					yield_arr_up_tmp.push_back(yield_up/yield);
					yield_arr_dn_tmp.push_back(yield_dn/yield);

				}
				for (int nextra = 0; nextra < extraSysn_all; nextra++){
					float up=0; 
					float dn=0; 
					if (stage1_name!="ZX"){
						up = extraSys_all_up[yearn][nextra][flav-1];
						dn = extraSys_all_dn[yearn][nextra][flav-1];
					}
					yield_arr_up_tmp.push_back(up);	
					yield_arr_dn_tmp.push_back(dn);	
				}
				for (int nextra = flav-1; nextra < flav; nextra++){
					float up=0; 
					float dn=0; 
					if (stage1_name=="ZX"){
						up = extraSys_zjet_up[yearn][flav-1];
						dn = extraSys_zjet_dn[yearn][flav-1];
					}
					yield_arr_up_tmp.push_back(up);	
					yield_arr_dn_tmp.push_back(dn);	
				}
				if(i<stage1_cat){
					yield_arr_up_tmp.push_back(1.02);	
					yield_arr_dn_tmp.push_back(0.98);	
				}
				else{
					yield_arr_up_tmp.push_back(0);	
					yield_arr_dn_tmp.push_back(0);	
				}
				if(stage1_name=="ggZZ"){
					yield_arr_up_tmp.push_back(1.1);	
					yield_arr_dn_tmp.push_back(0.9);	
				}
				else{
					yield_arr_up_tmp.push_back(0);	
					yield_arr_dn_tmp.push_back(0);	
				}

				cout<<"overall sys "<<sysnames_final.size()<<"\t"<< yield_arr_up_tmp.size()<<"\t"<<yield_arr_dn_tmp.size()<<endl;
				yield_arr_up.push_back(yield_arr_up_tmp);
				yield_arr_dn.push_back(yield_arr_dn_tmp);
			}

			TH1F *hstage_dbkg = (TH1F*)h2->ProjectionY(stage1_name+"_h1",i+1,i+1);
			RooDataHist *datahist_stage_dbkg = new RooDataHist(stage1_name+cat_name+"_datahistdbkg"+year,stage1_name+cat_name+"_datahistdbkg"+year,*dbkg,hstage_dbkg);
			RooHistPdf *histpdf_stage_dbkg = new RooHistPdf(stage1_name+"_dbkg"+year,stage1_name+"_dbkg"+year,*dbkg,*datahist_stage_dbkg);

			if(reco_name.Contains("VBF_2j")){
				//TH3F *hstage = new TH3F (stage1_name+"h3",stage1_name+"h3",5,0,1,4,0.5,1,4,0,1); 
				TH2F *hstage = new TH2F (stage1_name+"h3",stage1_name+"h3",4,0.,1,4,0,1); 
				TH2F *hstage_sys = new TH2F (stage1_name+"h3_sys",stage1_name+"h3_sys",4,0.,1,4,0,1); 
				int entries=0;
				if(i<stage1_cat)
					entries=tsig->Draw("d_2j_qq:DVBFDEC>>"+stage1_name+"h3",Form("weight*%f*(%s==%d&&chan==%d&&%sName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stageName.Data(),stage1_name.Data(),lumi)); 
				else
					entries=tzx->Draw("d_2j_qq:DVBFDEC>>"+stage1_name+"h3",Form("weight*%f*(%s==%d&&chan==%d&&htxs_stage1_red_catName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stage1_name.Data(),lumi)); 

				if (stage1_name=="ggH")
					tsys->Draw("d_2j_qq:DVBFDEC>>"+stage1_name+"h3_sys",Form("weight*%f*(%s==%d&&chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,lumi));
				else if (stage1_name=="ggH_ALT")
					tsys->Draw("d_2j_qq:DVBFDEC>>"+stage1_name+"h3_sys",Form("weight*(1./d_2jgen-1)*%f*(%s==%d&&chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,lumi));
				else if (stage1_name=="ggH_Mix")
					tsys->Draw("d_2j_qq:DVBFDEC>>"+stage1_name+"h3_sys",Form("weight*(2-d_2j_intgen)/d_2jgen*%f*(%s==%d&&chan==%d)/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,lumi));
				cout<<Form("weight*%f*(%s==%d&&chan==%d&&%sName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stageName.Data(),stage1_name.Data(),lumi)<<endl;
				cout<<entries<<endl;
				Floor(hstage);
				hstage->Draw("colz");
				gPad->Print(stage1_name+"h3_vbf2j.png");
				RooDataHist *datahist_stage = new RooDataHist(stage1_name+cat_name+"_datahist"+year,stage1_name+cat_name+"_datahist"+year,RooArgSet(*D2jet,*d_2j_qq),hstage);
				RooHistPdf *histpdf_stage = new RooHistPdf(stage1_name+cat_name+"_dis"+year,stage1_name+cat_name+"_dis"+year,RooArgSet(*D2jet,*d_2j_qq),*datahist_stage);
				if(stage1_name.Contains("ggH")){
					Floor(hstage_sys);
					hstage->Scale(1./hstage->Integral());
					hstage_sys->Scale(1./hstage_sys->Integral());
					TH2F *hstage_sys_mirr = (TH2F*)hstage->Clone();
					hstage_sys_mirr->Add(hstage,1);
					hstage_sys_mirr->Add(hstage_sys,-1);
					RooDataHist *datahist_stage_sys = new RooDataHist(stage1_name+cat_name+"_datahist_minloUp"+year,stage1_name+cat_name+"_datahist_minloUp"+year,RooArgSet(*D2jet,*d_2j_qq),hstage_sys);
					RooHistPdf *histpdf_stage_sys = new RooHistPdf(stage1_name+cat_name+"_dis_minloUp"+year,stage1_name+cat_name+"_dis_minloUp"+year,RooArgSet(*D2jet,*d_2j_qq),*datahist_stage_sys);

					RooDataHist *datahist_stage_sys_mirr = new RooDataHist(stage1_name+cat_name+"_datahist_minloDown"+year,stage1_name+cat_name+"_datahist_minloDown"+year,RooArgSet(*D2jet,*d_2j_qq),hstage_sys_mirr);
					RooHistPdf *histpdf_stage_sys_mirr = new RooHistPdf(stage1_name+cat_name+"_dis_minloDown"+year,stage1_name+cat_name+"_dis_minloDown"+year,RooArgSet(*D2jet,*d_2j_qq),*datahist_stage_sys_mirr);
					RooProdPdf *prodpdf_stage_sys = new RooProdPdf(stage1_name+"_minlo_"+year+"Up",stage1_name+"_minloUp",RooArgSet(*histpdf_stage_dbkg,*histpdf_stage_sys));
					RooProdPdf *prodpdf_stage_sys_mirr = new RooProdPdf(stage1_name+"_minlo_"+year+"Down",stage1_name+"_minloDown",RooArgSet(*histpdf_stage_dbkg,*histpdf_stage_sys_mirr));
					w.import(*prodpdf_stage_sys,RecycleConflictNodes());
					w.import(*prodpdf_stage_sys_mirr,RecycleConflictNodes());
					//  RooArgList pdflist("pdflist");
					//    RooArgList coeffs("coeffs");
					//    pdflist.add(*histpdf_stage);
					//    pdflist.add(*histpdf_stage_sys);
					//    pdflist.add(*histpdf_stage_sys_mirr);
					//     RooRealVar *sysrr=new RooRealVar("minlo","minlo",-3,3);
					//                    coeffs.add(*sysrr);

					//VerticalInterpPdf *histpdf_vert= new VerticalInterpPdf(stage1_name+"_vert","",RooArgSet(*D2jet,*d_2j_qq),pdflist,coeffs);
					//RooProdPdf *prodpdf_stage_vertical = new RooProdPdf(stage1_name,stage1_name,RooArgSet(*histpdf_stage_dbkg,*histpdf_vert));
					// w.import(*prodpdf_stage_vertical,RecycleConflictNodes());
				}
				RooProdPdf *prodpdf_stage = new RooProdPdf(stage1_name,stage1_name,RooArgSet(*histpdf_stage_dbkg,*histpdf_stage));
				w.import(*prodpdf_stage,RecycleConflictNodes());
			}
			else if(reco_name.Contains("TTH")){
				TH1F *hstage = new TH1F (stage1_name+"h2",stage1_name+"h2",4,-1,1); 
				if(i<stage1_cat)
					tsig->Draw("BDTG>>"+stage1_name+"h2",Form("weight*%f*(%s==%d&&chan==%d&&%sName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stageName.Data(),stage1_name.Data(),lumi)); 
				else
					tzx->Draw("BDTG>>"+stage1_name+"h2",Form("weight*%f*(%s==%d&&chan==%d&&htxs_stage1_red_catName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stage1_name.Data(),lumi)); 
				Floor(hstage);
				hstage->Draw("colz");
				gPad->Print(stage1_name+"h2_tth.png");
				RooDataHist *datahist_stage = new RooDataHist(stage1_name+cat_name+"_datahist"+year,stage1_name+cat_name+"_datahist"+year,RooArgSet(*BDTG),hstage);
				RooHistPdf *histpdf_stage = new RooHistPdf(stage1_name+cat_name+"_dis"+year,stage1_name+cat_name+"_dis"+year,RooArgSet(*BDTG),*datahist_stage);
				RooProdPdf *prodpdf_stage = new RooProdPdf(stage1_name,stage1_name,RooArgSet(*histpdf_stage_dbkg,*histpdf_stage));
				w.import(*prodpdf_stage,RecycleConflictNodes());
			}
			else if(reco_name.Contains("VH_Had")){
				TH1F *hstage = new TH1F (stage1_name+"h2",stage1_name+"h2",20,0,1); 
				if(i<stage1_cat)
					tsig->Draw("DVHDEC>>"+stage1_name+"h2",Form("weight*%f*(%s==%d&&chan==%d&&%sName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stageName.Data(),stage1_name.Data(),lumi)); 
				else
					tzx->Draw("DVHDEC>>"+stage1_name+"h2",Form("weight*%f*(%s==%d&&chan==%d&&htxs_stage1_red_catName==\"%s\")/(htxs_stage1_red_cat==-2? %f: 35.9)",lumi,category.Data(),reco_cat_num[j],flav,stage1_name.Data(),lumi)); 
				Floor(hstage);
				hstage->Draw("colz");
				gPad->Print(stage1_name+"h2_vh.png");
				RooDataHist *datahist_stage = new RooDataHist(stage1_name+cat_name+"_datahist"+year,stage1_name+cat_name+"_datahist"+year,RooArgSet(*DVHDEC),hstage);
				RooHistPdf *histpdf_stage = new RooHistPdf(stage1_name,stage1_name,RooArgSet(*DVHDEC),*datahist_stage);
				w.import(*histpdf_stage,RecycleConflictNodes());
			}
			else{
				histpdf_stage_dbkg->SetName(stage1_name);
				histpdf_stage_dbkg->SetTitle(stage1_name);
				w.import(*histpdf_stage_dbkg,RecycleConflictNodes());
			}
		}

		RooDataSet *data;
			TTree* reducetree= tdata->CopyTree(Form("chan==%d&&%s==%d",flav,category.Data(),reco_cat_num[j]));
		if(reco_name.Contains("VBF_2j"))
			data=new RooDataSet("data_obs","",reducetree,RooArgList(*dbkg,*D2jet,*d_2j_qq));
		else if(reco_name.Contains("TTH"))
			data=new RooDataSet("data_obs","",reducetree,RooArgList(*dbkg,*BDTG));
		else if(reco_name.Contains("VH_Had"))
			data=new RooDataSet("data_obs","",reducetree,RooArgList(*DVHDEC));
		else
			data=new RooDataSet("data_obs","",reducetree,RooArgList(*dbkg));

		w.import(*data);


		card<<"bin \t";
		writeline(bin_arr,card);
		card<<"process \t";
		writeline(processName_arr,card);
		card<<"process \t";
		writeline(process_arr,card);
		card<<"rate \t";
		writeline(yield_arr,card);
		if(doSys){
			for(int isys=0;isys<sysnames_final.size();isys++){
				float dn =0; float up=0;
				card<<sysnames_final[isys]<<"\t lnN \t";
				for (int prc=0;prc<yield_arr_dn.size();prc++){
					if (yield_arr_dn.at(prc).at(isys)==0)	
						card<<"- \t";
					else if(processName_arr.at(prc)=="TTH_ALT" )
						if( processName_arr.at(prc-1)=="TTH")
							card<<yield_arr_dn.at(prc-1).at(isys)<<"/"<<yield_arr_up.at(prc-1).at(isys)<<"\t";
						else
							card<<"1/1"<<"\t";
					else if (processName_arr.at(prc)=="ggH"){
						dn = yield_arr_dn.at(prc).at(isys); up = yield_arr_up.at(prc).at(isys);
						card<<yield_arr_dn.at(prc).at(isys)<<"/"<<yield_arr_up.at(prc).at(isys)<<"\t";
					}
					else if (processName_arr.at(prc).Contains("ggH")){
						card<<dn<<"/"<<up<<"\t";
					}
					else
						card<<yield_arr_dn.at(prc).at(isys)<<"/"<<yield_arr_up.at(prc).at(isys)<<"\t";
				}
				card<<endl;
			}
			if(reco_name.Contains("VBF_2j")){
				//card<<"minlo param 0 1 [-3,3]"<<endl;
				card<< "minlo"<<"\t shape \t";
				for (int prc=0;prc<yield_arr_dn.size();prc++){
					if (processName_arr.at(prc).Contains("ggH"))	
						card<<"1 \t";
					else
						card<<"- \t";
				}
				card<<endl;
			}
		}
		card.close();
		f->cd();
		w.Write();
		f->Close();
	}


}
