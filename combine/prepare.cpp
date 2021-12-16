//build 2D RooHistPdf and datacard file
//for ttbar and its background
using namespace std;
using namespace RooFit;

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
void Floor(TH2D* histo){
	for (int i=0;i<histo->GetNbinsX();i++){
		for (int j=0;j<histo->GetNbinsY();j++){
			if(!(histo->GetBinContent(i+1,j+1)>1.E-8)){
				histo->SetBinContent(i+1,j+1,1.E-8);
				float xx=histo->GetXaxis()->GetBinCenter(i+1);
				float yy=histo->GetYaxis()->GetBinCenter(j+1);
				//cout<<"warning!!!!! No events in x: "<<xx<<" y: "<<yy<<endl;
			}
		}
	}
}
void prepare(){
	TString dosys_th_th=true;
	const int nsample=26;
	const int nsignal=6;

	TString fileNames[nsample]={"new_TTToSemiLeptonic_TuneCP5_13TeV-powheg.root",
                            "new_TTTo2L2Nu_TuneCP5_13TeV-powheg.root",
                            "new_TTToHadronic_TuneCP5_13TeV-powheg.root",

                            "new_DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                            "new_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                            "new_DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                            "new_DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                            "new_DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                            "new_DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                            "new_DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM.root",                           
                            "new_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM.root",
                        
                            "new_ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo.root",
                            "new_ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin.root",
                            "new_ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin.root",
                            "new_ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg.root",
                            "new_ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg.root",                                                                               
                            
                            "new_WWTo1L1Nu2Q_TuneCP5_13TeV-amcatnloFXFX.root",
                            "new_ZZTo2L2Nu_TuneCP5_13TeV_powheg.root",
                            "new_WZTo1L1Nu2Q_TuneCP5_13TeV-amcatnloFXFX.root",
                            
                            "new_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM.root",

                         /*   "new_QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                          */  
                   
};
							
	Float_t cross_sections[nsample]={366.91, 89.05, 377.96,
									169.9, 147.4, 41.0, 5.7, 1.4, 0.63, 0.15, 0.0036,
									3.36, 136.02, 80.95, 35.6, 35.6,
									45.68, 0.9738, 11.66,
									1345.7, 359.7, 48.9, 12.1, 5.5, 1.3, 0.032,
								    //27990000, 1712000, 347700, 32100, 6831, 1207, 119.9, 25.2,
								    };
	Float_t K_Factor[nsample]={1.0, 1.0, 1.0,
								1.23,1.23,1.23,1.23,1.23,1.23,1.23,1.23,
								1.0,1.0,1.0,1.0,1.0,
								1.0,1.0,1.0,
								1.21,1.21,1.21,1.21,1.21,1.21,1.21,
								//1.0, 1.0, 1.0,1.0, 1.0, 1.0,1.0, 1.0,
							};		
	TString dir="./output2/";
	TString process[]={"ttbar","DYJets","STop","VV","WJets","QCD"};
	Int_t sample_id[]={2, 10, 15, 18, 25, 33};
	Int_t Cpq3[6]={ 0, 0, 0, 0, 0, 0};
	Int_t Cpu[6]={  0, 1, 0, 0, 2, 0};
	Int_t ReCup[6]={0, 0, 1, 0, 0, 2};
	Int_t ImCup[6]={0, 0, 0, 1, 0, 0};
	float lumi=137.1;
	TString outputDir="datacard2";

	Double_t mtt_edges[9]={0,370,420,500,600,700,800,950,2000};
	Double_t ytt_edges[10]={-5.0,-1.4,-0.9,-0.5,-0.15,0.15,0.5,0.9,1.4,5.0};
	RooRealVar* mtt=new RooRealVar("mass_tt","mass_tt",0,2000);
	RooRealVar* ytt=new RooRealVar("rapidity_tt","rapidity_tt",-5,5);
	const int xbin=8, ybin=9;
	mtt->setBins(xbin);
	ytt->setBins(ybin);
	TString cuts[]={"(jet_num == 3 && likelihood<20.0)","(jet_num >= 4 && likelihood<20.0 )"};
	TString cutsName[]={"3jets","4jets"};
	Float_t entries[2][nsample];// number of events in 3jets and 4jets final states
	for(int s=0; s<2; s++){ //loop over final states
		int nprocess=0; //count process was dealed with
		TString category="ttbar_"+cutsName[s];
		TFile *file=new TFile(outputDir+"/"+category+".root","recreate");
		RooWorkspace w("w");

		std::vector<TString> process_names; //names of sigal and bkg
		std::vector<int> process_id;  //process ID; minus and zero for sigal; positive for bkg
		std::vector<TString> bin_arr;   //category name
		std::vector<float> yield_array;  //rate(event yeild)
		std::vector<TString> STop_norm;  //single top  normlization uncertainty
		std::vector<TString> VV_norm;  //diboson  normlization uncertainty
		std::vector<TString> DYJets_norm;  //Drell-Yan+jets  normlization uncertainty
		std::vector<TString> WJets_norm;  // WJets normlization uncertainty
		std::vector<TString> sig_norm;   //signal norlization uncertainty
		std::vector<TString> cms_lumi;
		std::vector<TString> ew_weight;
    	TH2D* h2dist[nsignal+4];//6 signal + 4 background
    	std::vector<TH2D> h2sys_up[nsignal+4];//2D array [process][sys]
		std::vector<TH2D> h2sys_dn[nsignal+4];//2D array [process][sys]
    	//////////////////////////////////////////////
    	//for systematic uncentaitny
		  std::vector<TString> sysNames;

		for(int i=0;i<nsample;i++) { //loop over samples
			TChain* chain=new TChain("mytree");
			TChain* chain2=new TChain("rawtree");
			if(s==0)
				fileNames[i]=fileNames[i].ReplaceAll(".root","_*.root");
			chain->Add(dir+fileNames[i]);
			chain2->Add(dir+fileNames[i]);
			Int_t nMC, ncut;
			nMC=chain2->GetEntries();
			ncut=chain->GetEntries();
			cout<<nMC<<" events simulated and "<<ncut<<" events selected in "<<fileNames[i]<<endl;
			float global_weight=cross_sections[i]*1000*lumi/nMC*K_Factor[i];
			TString gen_weight="Generator_weight/abs(Generator_weight)";
			TString condition="(mass_tt<=2000.0)&&(abs(rapidity_tt)<=5.0)&&(likelihood <20.0)";
			chain->Draw("mass_tt",Form("%s*%s*%s",cuts[s],condition,gen_weight));
			Float_t entry_cut=chain->GetSumOfWeights();
			entries[s][i]=entry_cut*global_weight; //number of events in each channel
			TString sample_name=fileNames[i];
			sample_name=sample_name.ReplaceAll("_*.root","_hist");
			sample_name=sample_name.ReplaceAll("new_","");
			//for sysmetic uncentainties with weights in Nanoaod		
			if(i==0){ //get the systematic weight names at first sample
					TH1D* hname=new TH1D("hname","hname",20,0,20);
					chain->Draw("weight_name>>hname");
					for(int k=0;k<hname->GetNbinsX();k++){
					TString sysname=hname->GetBinLabel(k+1);
					sysNames.push_back(sysname);
					cout<<"systematic: "<<sysname<<endl;
				}

			}
			std::vector<TH2D> h2sample_sysup;
			std::vector<TH2D> h2sample_sysdn;
				
			if(i <= sample_id[0]){
				for(int k=0;k<nsignal;k++){ //loop over EW weights
			        TString weight_EW=Form("weight_ci%d%d%d%d",Cpq3[k],Cpu[k],ReCup[k],ImCup[k]);
			        TString weight=Form("%f*%s*%s",global_weight,weight_EW.Data(),gen_weight);
			        TString sample_weighted=sample_name+"_"+weight_EW;
			        TH2D* h2sample=new TH2D(sample_weighted,sample_weighted,xbin,mtt_edges, ybin, ytt_edges);
			        h2sample->Sumw2();
					chain->Draw("rapidity_tt:mass_tt>>"+sample_weighted, weight+"*"+cuts[s] );
					
					if(dosys_th){	
						for(int n=0;n < sysNames.size();n++){
							TString h2sysName=sample_weighted+"_"+sysNames[n];
							TH2D* h2sys_up=new TH2D(h2sysName+"_up",h2sysName+"_up",xbin,mtt_edges,ybin, ytt_edges);
							TH2D* h2sys_dn=new TH2D(h2sysName+"_dn",h2sysName+"_dn",xbin,mtt_edges,ybin, ytt_edges);
			        		h2sys_up->Sumw2();
			        		h2sys_dn->Sumw2();
							chain->Draw("rapidity_tt:mass_tt>>"h2sysName+"_up", Form("%s*%s*weight_up*(weight_name==%s)",weight,cuts[s],sysNames[n]));
							chain->Draw("rapidity_tt:mass_tt>>"h2sysName+"_dn", Form("%s*%s*weight_down*(weight_name==%s)",weight,cuts[s],sysNames[n]));
							h2sample_sysup.push_back(h2sys_up);
							h2sample_sysdn.push_back(h2sys_dn);
							delete h2sys_dn, h2sys_up;
						}
					}

					if(i==0){
						h2dist[k]=(TH2D*)h2sample->Clone();
						h2dist[k]->SetName("ttbar_"+weight_EW);
						h2dist[k]->SetTitle("ttbar_"+weight_EW);
						if(dosys_th){
							for(int n=0;n<sysNames.size();n++){
								h2sys_up[k][n]=(TH2D*)h2sample_sysup->Clone();
								h2sys_dn[k][n]=(TH2D*)h2sample_sysdn->Clone();
								h2sys_up[k][n]->SetName("ttbar_"+weight_EW+"_"+sysNames[n]+"Up");
								h2sys_dn[k][n]->SetName("ttbar_"+weight_EW+"_"+sysNames[n]+"Down");
							}
						}
					}
					else
						h2dist[k]->Add(h2sample);
					    if(dosys_th){
							for(int n=0;n<sysNames.size();n++){
								h2sys_up[k][n]->Add(h2sample_sysup[n]);
								h2sys_dn[k][n]->Add(h2sample_sysdn[n]);
							}
						}
					if (i==sample_id[0]){
						process_id.push_back(-k);
						TString process_name=h2dist[k]->GetName();
						process_names.push_back(process_name);
						Floor(h2dist[k]);
						h2dist[k]->Draw("colz text");
						gPad->Print(outputDir+"/"+process_name+"_"+cutsName[s]+"_hist.png");
						float yield=h2dist[k]->GetSumOfWeights();
						yield_array.push_back(yield);
						bin_arr.push_back(category);
						sig_norm.push_back("1.05");
						STop_norm.push_back("-");
						VV_norm.push_back("-");
						WJets_norm.push_back("-");
						DYJets_norm.push_back("-");
						cms_lumi.push_back("1.025");
						ew_weight.push_back("1.01");
						cout<<"after reweight there are "<<yield<<" events in "<<process_name<<" in "<<cutsName[s]<<endl;
						
						RooDataHist* datahist=new RooDataHist(process_name+"_datahist",process_name+"_datahist",RooArgSet(*mtt,*ytt),h2dist[k]);
						RooHistPdf* hist_pdf=new RooHistPdf(process_name,process_name,RooArgSet(*mtt,*ytt),*datahist);
						w.import(*hist_pdf,RecycleConflictNodes());
						if(dosys_th){
							for(int n=0;n<sysNames.size();n++){
								TString pro_sysup_name=h2sys_up[k][n]->GetName();
								TString pro_sysdn_name=h2sys_dn[k][n]->GetName();
								RooDataHist* hist_up=new RooDataHist(pro_sysup_name+"_hist",pro_sysup_name+"_hist",RooArgSet(*mtt,*ytt),h2sys_up[k][n]);
								RooHistPdf* pdf_up=new RooHistPdf(pro_sysup_name,pro_sysup_name,RooArgSet(*mtt,*ytt),*hist_up);
								RooDataHist* hist_dn=new RooDataHist(pro_sysdn_name+"_hist",pro_sysdn_name+"_hist",RooArgSet(*mtt,*ytt),h2sys_dn[k][n]);
								RooHistPdf* pdf_dn=new RooHistPdf(pro_sysdn_name,pro_sysdn_name,RooArgSet(*mtt,*ytt),*hist_dn);
								w.import(*pdf_up,RecycleConflictNodes());
								w.import(*pdf_dn,RecycleConflictNodes());
								delete hist_up,hist_dn,pdf_up,pdf_dn;
							}
						}

					}
                    delete h2sample, h2sample_sysup, h2sample_sysdn;

				}  //end of loop over EW weights
                if (i==sample_id[0]) nprocess++;
			}
			else{
				TH2D* h2sample=new TH2D(sample_name,sample_name,xbin,mtt_edges, ybin, ytt_edges);
				h2sample->Sumw2();
				chain->Draw("rapidity_tt:mass_tt>>"+sample_name, Form("%f*%s*%s",global_weight,cuts[s].Data(),gen_weight));
				if(dosys_th){	
						for(int n=0;n< sysNames.size();n++){
							TString h2sysName=sample_name+"_"+sysNames[n];
							TH2D* h2sys_up=new TH2D(h2sysName+"_up",h2sysName+"_up",xbin,mtt_edges,ybin, ytt_edges);
							TH2D* h2sys_dn=new TH2D(h2sysName+"_dn",h2sysName+"_dn",xbin,mtt_edges,ybin, ytt_edges);
			        		h2sys_up->Sumw2();
			        		h2sys_dn->Sumw2();
							chain->Draw("rapidity_tt:mass_tt>>"h2sysName+"_up", Form("%f*%s*%s*weight_up*(weight_name==%s)",global_weight,gen_weight,cuts[s],sysNames[n]));
							chain->Draw("rapidity_tt:mass_tt>>"h2sysName+"_dn", Form("%f*%s*%s*weight_down*(weight_name==%s)",global_weight,gen_weight,cuts[s],sysNames[n]));
							h2sample_sysup.push_back(h2sys_up);
							h2sample_sysdn.push_back(h2sys_dn);
							delete h2sys_dn, h2sys_up;
						}
					}

				if(i>sample_id[nprocess-1] && i <= sample_id[nprocess]){
					if(i==sample_id[nprocess-1]+1){
						h2dist[nsignal-1+nprocess]=(TH2D*)h2sample->Clone();
						h2dist[nsignal-1+nprocess]->SetName(process[nprocess]);
						h2dist[nsignal-1+nprocess]->SetTitle(process[nprocess]);

						if(dosys_th){
							for(int n=0;n<sysNames.size();n++){
								h2sys_up[nsignal-1+nprocess][n]=(TH2D*)h2sample_sysup->Clone();
								h2sys_dn[nsignal-1+nprocess][n]=(TH2D*)h2sample_sysdn->Clone();
								h2sys_up[nsignal-1+nprocess][n]->SetName(process[nprocess]+"_"+sysNames[n]+"Up");
								h2sys_dn[nsignal-1+nprocess][n]->SetName(process[nprocess]+"_"+sysNames[n]+"Down");
							}
						}
					}
					else
						{
							h2dist[nsignal-1+nprocess]->Add(h2sample);
	                        if(dosys_th){
								for(int n=0;n<sysNames.size();n++){
									h2sys_up[nsignal-1+nprocess][n]->Add(h2sample_sysup[n]);
									h2sys_dn[nsignal-1+nprocess][n]->Add(h2sample_sysdn[n]);
								}
							}
                        }
					if(i==sample_id[nprocess]){
						process_id.push_back(nprocess);
						TString process_name=h2dist[nsignal-1+nprocess]->GetName();
						process_names.push_back(process_name);
						Floor(h2dist[nsignal-1+nprocess]);
						h2dist[nsignal-1+nprocess]->Draw("colz text");
						gPad->Print(outputDir+"/"+process_name+"_"+cutsName[s]+"_hist.png");
						float yield=h2dist[nsignal-1+nprocess]->GetSumOfWeights();
						yield_array.push_back(yield);
						bin_arr.push_back(category);
						sig_norm.push_back("-");
						cms_lumi.push_back("1.025");
						ew_weight.push_back("-");
						if (process_name.Contains("STop")){
							STop_norm.push_back("1.15");
							VV_norm.push_back("-");
							WJets_norm.push_back("-");
							DYJets_norm.push_back("-");
						}
						else if(process_name.Contains("WJets")){
							STop_norm.push_back("-");
							VV_norm.push_back("-");
							WJets_norm.push_back("1.30");
							DYJets_norm.push_back("-");
						}
						else if(process_name.Contains("DYJets")){
							STop_norm.push_back("-");
							VV_norm.push_back("-");
							WJets_norm.push_back("-");
							DYJets_norm.push_back("1.30");
						}
						else if(process_name.Contains("VV")){
							STop_norm.push_back("-");
							VV_norm.push_back("1.30");
							WJets_norm.push_back("-");
							DYJets_norm.push_back("-");
						}
						cout<<"after reweight there are "<<yield<<" events in "<<process_name<<" in "<<cutsName[s]<<endl;
						
						RooDataHist* datahist=new RooDataHist(process_name+"_datahist",process_name+"_datahist",RooArgSet(*mtt,*ytt),h2dist[nsignal-1+nprocess]);
						RooHistPdf* hist_pdf=new RooHistPdf(process_name,process_name,RooArgSet(*mtt,*ytt),*datahist);
						w.import(*hist_pdf,RecycleConflictNodes());
						if(dosys_th){
							for(int n=0;n<sysNames.size();n++){
								TString pro_sysup_name=h2sys_up[nsignal-1+nprocess][n]->GetName();
								TString pro_sysdn_name=h2sys_dn[nsignal-1+nprocess][n]->GetName();
								RooDataHist* hist_up=new RooDataHist(pro_sysup_name+"_hist",pro_sysup_name+"_hist",RooArgSet(*mtt,*ytt),h2sys_up[nsignal-1+nprocess][n]);
								RooHistPdf* pdf_up=new RooHistPdf(pro_sysup_name,pro_sysup_name,RooArgSet(*mtt,*ytt),*hist_up);
								RooDataHist* hist_dn=new RooDataHist(pro_sysdn_name+"_hist",pro_sysdn_name+"_hist",RooArgSet(*mtt,*ytt),h2sys_dn[nsignal-1+nprocess][n]);
								RooHistPdf* pdf_dn=new RooHistPdf(pro_sysdn_name,pro_sysdn_name,RooArgSet(*mtt,*ytt),*hist_dn);
								w.import(*pdf_up,RecycleConflictNodes());
								w.import(*pdf_dn,RecycleConflictNodes());
								delete hist_up,hist_dn,pdf_up,pdf_dn;
							}
						}
						nprocess++;
					}
                    delete h2sample, h2sample_sysdn, h2sample_sysup;

				}
			}
		}
			
		ofstream card;
		card.open (outputDir+"/"+category+".txt");
		card <<"Datacard for event category: "<< category<<endl;
		card<< "imax 1 number of channels"<<endl;
		card<< "jmax 9 number of processes minus 1"<<endl;
		card<< "kmax * number of nuisance parameters"<<endl;
		card<<"---------------------------------"<<endl;
		card<<endl;
		card<< "shapes * "<< category << " "<< category+".root w:$PROCESS w:$PROCESS_$SYSTEMATIC" <<endl;
		card<<"---------------------------------"<<endl;
		card<< "bin           "<< category <<endl;
		card<< "observation   "<< "-1"<<endl;
		card<<"---------------------------------"<<endl;
		card<<endl;
		card<<"bin \t";
		writeline(bin_arr,card);
		card<<"process \t";
		writeline(process_names,card);
		card<<"process \t";
		writeline(process_id,card);
		card<<"rate \t";
		writeline(yield_array,card);
		card<<"sig_norm"<<"\t lnN \t";
		writeline(sig_norm,card);
		card<<"DYJets_norm"<<"\t lnN \t";
		writeline(DYJets_norm,card);
		card<<"STop_norm"<<"\t lnN \t";
		writeline(STop_norm,card);
		card<<"VV_normli"<<"\t lnN \t";
		writeline(VV_norm,card);
		card<<"WJets_norm"<<"\t lnN \t";
		writeline(WJets_norm,card);
		card<<"cms_lumi"<<"\t lnN \t";
		writeline(cms_lumi,card);
		card<<"EW_weight"<<"\t lnN \t";
		writeline(ew_weight,card);
		if(dosys_th){
			for(int n=0;n<sysNames.size();n++){
				card<<sysNames[n]<<"\t shape \t"
				for(int p=0;p<nsignal+4;p++){
					card<<"1"<<"\t"
				}
				card<<endl;

			}
		}
		//build dataset, but for expected results dataset is not needed.
		TChain *chain_data=new TChain("mytree");
		chain_data->Add(dir+fileNames[0]); //fake data now
	    RooDataHist *data;
	    TString hist_data_name="hist_data";
        TH2D* hist_data=new TH2D(hist_data_name,hist_data_name,8,mtt_edges, 9, ytt_edges);
	    chain_data->Draw("abs(rapidity_tt):mass_tt>>"+hist_data_name);
	    data=new RooDataHist("data_obs","",RooArgSet(*mtt,*ytt),Import(*hist_data));
	    w.import(*data);
	    w.Print();

	    card.close();
		file->cd();
		w.Write();
		file->Close();
		delete hist_data;
		delete h2dist, h2sys_up, h2sys_dn;
	
	}// end of loop over final states
	cout.setf(ios::fixed, ios::floatfield); // set fixed floating format
	cout.precision(2); // for fixed format, two decimal places
	//cout << fixed << setprecision(2); // same effects, but using manipulators	
	cout<<entries[0][0]<<" events in 3jets from ttbar semiletonic channel"<<endl;
	cout<<entries[1][0]<<" events in 4jets from ttbar semiletonic channel"<<endl;
	Float_t nttbar_3jet=entries[0][0]+entries[0][1]+entries[0][2];
	Float_t nttbar_4jet=entries[1][0]+entries[1][1]+entries[1][2];
	cout<<nttbar_3jet<<" events in ttbar 3jets "<<endl;
	cout<<nttbar_4jet<<" events in ttbar 4jets "<<endl;
	cout<<"expect "<<nttbar_4jet+nttbar_3jet<<" events in ttbar"<<endl;

}
