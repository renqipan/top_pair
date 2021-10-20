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
void Floor(TH2F* histo){
	for (int i=0;i<histo->GetNbinsX();i++){
		for (int j=0;j<histo->GetNbinsY();j++){
			if(histo->GetBinContent(i+1,j+1)==0)
				histo->SetBinContent(i+1,j+1,1.E-8);
		}
	}
}
void prepare(){
	Int_t nsample=34;
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
                            
                            "new_WW_TuneCP5_13TeV.root",
                            "new_WZ_TuneCP5_13TeV.root",
                            "new_ZZ_TuneCP5_13TeV.root",
                            
                            "new_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM.root",
                            "new_WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM.root",

                            "new_QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
							"new_QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            "new_QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM.root",
                            
                   
};
							
	Float_t cross_sections[nsample]={366.91, 89.05, 377.96,
									169.9, 147.4, 41.0, 5.7, 1.4, 0.63, 0.15, 0.0036,
									3.36, 136.02, 80.95, 35.6, 35.6,
									118.7, 16.5, 47.1,
									1345.7, 359.7, 48.9, 12.1, 5.5, 1.3, 0.032,
								    27990000, 1712000, 347700, 32100, 6831, 1207, 119.9, 25.2,
								    };
	Float_t K_Factor[nsample]={1.0, 1.0, 1.0,
								1.23,1.23,1.23,1.23,1.23,1.23,1.23,1.23,
								1.0,1.0,1.0,1.0,1.0,
								1.0,1.0,1.0,
								1.21,1.21,1.21,1.21,1.21,1.21,1.21,
								1.0, 1.0, 1.0,1.0, 1.0, 1.0,1.0, 1.0,};		
	TString dir="/afs/cern.ch/user/r/repan/work/top_pair/condor/output/";
	TString process[]={"ttbar","DYJets","STop","VV","WJets","QCD"};
	Int_t sample_id[]={2, 10, 15, 18, 25, 33};
	Int_t nsignal=9;
	Int_t Cpq3[9]={ 0, 1, 0, 0, 0, 2, 0, 0, 1};
	Int_t Cpu[9]={  0, 0, 1, 0, 0, 0, 2, 0, 1};
	Int_t ReCup[9]={0, 0, 0, 1, 0, 0, 0, 2, 0};
	Int_t ImCup[9]={0, 0, 0, 0, 1, 0, 0, 0, 0};
	float lumi=137.1;
	Double_t mtt_edges[10]={0,350,400,500,600,700,800,950,1200,2000};
	Double_t ytt_edges[12]={-5.0,-3.0,-1.6,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.6,3.0,5.0};
	TH2F* h2dist[nsignal+5];//9 signal + 5 background
	TString outputDir="datacard";

	RooRealVar* mtt=new RooRealVar("mass_tt","mass_tt",0,2000);
	RooRealVar* ytt=new RooRealVar("rapidity_tt","rapidity_tt",-5,5);
	mtt->setBins(9);
	ytt->setBins(11);
	TString cuts[]={"(jet_num == 3)","(jet_num >= 4)"};
	TString cutsName[]={"3jets","4jets"};
	Float_t entries[2][nsample];// number of events in 3jets and 4jets final states
	for(int s=0; s<2; s++){
		int nprocess=0; //count process was dealed with
		TString category="ttbar_"+cutsName[s];
		TFile *file=new TFile(outputDir+"/"+category+".root","recreate");
		RooWorkspace w("w");

		std::vector<TString> process_names; //names of sigal and bkg
		std::vector<int> process_id;  //process ID; minus and zero for sigal; positive for bkg
		std::vector<TString> bin_arr;   //category name
		std::vector<float> yield_array;  //rate(event yeild)
		std::vector<TString> bkg_norm;  //background  normlization uncertainty
		std::vector<TString> sig_norm;   //signal norlization uncertainty

		for(int i=0;i<nsample;i++) { //loop over samples
			TChain* chain=new TChain("mytree");
			TChain* chain2=new TChain("rawtree");
			chain->Add(dir+fileNames[i]);
			chain2->Add(dir+fileNames[i]);
			Int_t nMC, ncut;
			nMC=chain2->GetEntries();
			ncut=chain->GetEntries();
			cout<<nMC<<"events simulated and "<<ncut<<" events selected in "<<fileNames[i]<<endl;
			float global_weight=cross_sections[i]*1000*lumi/nMC*K_Factor[i];
			Int_t entry_cut=chain->Draw("mass_tt",cuts[s]);
			entries[s][nsample]=entry_cut*global_weight; //number of events in each channel
			TString sample_name=fileNames[i].ReplaceAll(".root","_hist");
			sample_name=sample_name.ReplaceAll("new_","");
			if(i <= sample_id[0]){
				for(int k=0;k<nsignal;i++){
			        TString weight_EW=Form("ci%d%d%d%d",Cpq3[k],Cpu[k],ReCup[k],ImCup[k]);
			        TString weight=Form("%f*%s",global_weight,weight_EW);
			        TString sample_weighted=sample_name+"_"+weight_EW;
			        TH2F* h2sample=new TH2F(sample_weighted,sample_weighted,9,mtt_edges, 11, ytt_edges);
			        h2sample->SumW2();
					chain->Draw("abs(rapidity_tt):mass_tt>>"+sample_weighted, weight+"*"+cuts[s] );
					chain->Draw("-abs(rapidity_tt):mass_tt>>+"+sample_weighted,weight+"*"+cuts[s]);
					h2sample->Scale(1/2.0);
					if(i==0){
						h2dist[k]=h2sample->Clone();
						h2dist[k]->SumW2();
						h2dist[k]->SetName("ttbar_"+weight_EW);
						h2dist[k]->SetTitle("ttbar_"+weight_EW);
					}
					else
						h2dist[k]->Add(h2sample);
					if (i==sample_id[0]){
						process_id.push_back(-k);
						process_name=h2dist[k]->GetName();
						process_names.push_back(process_name);
						Floor(h2dist[k]);
						h2dist[k]->Draw("colz text");
						gPad->Print(outputDir+"/"+process_name+"_"+cutsName[s]+"_hist.png");
						float yield=h2dist[k]->GetSumOfWeights();
						yield_array.push_back(yield);
						bin_arr.push_back(category);
						sig_norm.push_back("1.05");
						bkg_norm.push_back("-");
						cout<<"after reweight there are "<<yield<<" events in "<<process_name<<" in "<<cutsName[s]<<endl;
						
						RooDataHist* datahist=new RooDataHist(process_name+"_datahist",process_name+"_datahist",RooArgSet(*mtt,*ytt),h2dist[k]);
						RooHistPdf* hist_pdf=new RooHistPdf(process_name,process_name,RooArgSet(*mtt,*ytt),*datahist);
						w.import(*hist_pdf,RecycleConflictNodes());
						nprocess++;

					}

				}  

			}
			else{
				TH2F* h2sample=new TH2F(sample_name,sample_name,9,mtt_edges, 11, ytt_edges);
				h2sample->SumW2();
				chain->Draw("abs(rapidity_tt):mass_tt>>"+sample_name, Form("%f*%s",global_weight,cuts[s]));
				chain->Draw("-abs(rapidity_tt):mass_tt>>+"+sample_name,Form("%f*%s",global_weight,cuts[s]));
				h2sample->Scale(1/2.0);
				if(i>sample_id[nprocess-1] && i <= sample_id[nprocess]){
					if(i==sample_id[nprocess-1]+1){
						h2dist[nsignal-1+nprocess]=h2sample->Clone();
						h2dist[nsignal-1+nprocess]->SumW2();
						h2dist[nsignal-1+nprocess]->SetName(process[nprocess]);
						h2dist[nsignal-1+nprocess]->SetTitle(process[nprocess]);
					}
					else
						h2dist[nsignal-1+nprocess]->Add(h2sample);
					if(i==sample_id[nprocess]){
						process_id.push_back(nprocess);
						process_name=h2dist[nsignal-1+nprocess]->GetName();
						process_names.push_back(process_name);
						Floor(h2dist[nsignal-1+nprocess]);
						h2dist[nsignal-1+nprocess]->Draw("colz text");
						gPad->Print(outputDir+"/"+process_name+"_"+cutsName[s]+"_hist.png");
						float yield=h2dist[nsignal-1+nprocess]->GetSumOfWeights();
						yield_array.push_back(yield);
						bin_arr.push_back(category);
						sig_norm.push_back("-");
						if (process_name.Contains("STop")){
							bkg_norm.push_back("1.15");
						}
						else if(process_name.Contains("WJets")){
							bkg_norm.push_back("1.30");
						}
						else if(process_name.Contains("QCD")){
							bkg_norm.push_back("1.30");
						}
						cout<<"after reweight there are "<<yield<<" events in "<<process_name<<" in "<<cutsName[s]<<endl;
						
						RooDataHist* datahist=new RooDataHist(process_name+"_datahist",process_name+"_datahist",RooArgSet(*mtt,*ytt),h2dist[nsignal-1+nprocess]);
						RooHistPdf* hist_pdf=new RooHistPdf(process_name,process_name,RooArgSet(*mtt,*ytt),*datahist);
						w.import(*hist_pdf,RecycleConflictNodes());
						nprocess++;
					}

				}
			}
		}
			
		ofstream card;
		card.open (outputDir+"/"+category+".txt");
		card <<"Datacard for event category: "<< category<<endl;
		card<< "imax 1 number of channels"<<endl;
		card<< "jmax 13 number of processes minus 1"<<endl;
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
		card<<"bkg_norm"<<"\t lnN \t";
		writeline(bkg_norm,card);
		
		//build dataset, but for expected results dataset is not needed.
		TChain *chain_data=new TChain("mytree");
		chain_data->Add(fileNames[0]); 
	    RooDataHist *data;
	    TH2F* hist_data=new TH2F("hist_data","hist_data",9,mtt_edges, 11, ytt_edges)
	    chain_data->Draw("abs(rapidity_tt):mass_tt>>"+"hist_data");
		chain_data->Draw("-abs(rapidity_tt):mass_tt>>+"+"hist_data");
		hist_data->Scale(1/2.0)
	    data=new RooDataHist("data_obs","",RooArgSet(*mtt,*ytt),Import(*hist_data));
	    w.import(*data);
	    w.Print();

	    card.close();
		file->cd();
		w.Write();
		file->Close();
	
	}
	cout.setf(ios::fixed, ios::floatfield); // set fixed floating format
	cout.precision(2); // for fixed format, two decimal places
	//cout << fixed << setprecision(2); // same effects, but using manipulators	
	cout<<entries[0][0]<<" events in 3jets from ttbar semiletonic channel"<<endl;
	cout<<entries[1][0]<<" events in 4jets from ttbar semiletonic channel"<<endl;
	Float_t nttbar_3jet=entries[0][0]+entries[0][1]+entries[0][2];
	Float_t nttbar_4jet=entries[1][0]+entries[1][1]+entries[1][2];
	cout<<nttbar_3jet<<"events in ttbar 3jets "<endl;
	cout<<nttbar_4jet<<"events in ttbar 4jets "<endl;
	cout<<"expect "<<nttbar_4jet+nttbar_3jet<<" events in ttbar"<<endl;

}