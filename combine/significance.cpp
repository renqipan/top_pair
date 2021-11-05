//calculate significance of ttbar by s/sqrt(s+b)
//written by Ren-Qi Pan in Nov 5, 2021.
#include <cmath>
using namespace std;

void significance(){
	const int nsample=26;
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
									118.7, 16.5, 47.1,
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
	TString dir="./output/";
	TString process[]={"ttbar","DYJets","STop","VV","WJets","QCD"};
	Int_t sample_id[]={2, 10, 15, 18, 25, 33};
	const int nsignal=9;
	Int_t Cpq3[9]={ 0, 1, 0, 0, 0, 2, 0, 0, 1};
	Int_t Cpu[9]={  0, 0, 1, 0, 0, 0, 2, 0, 1};
	Int_t ReCup[9]={0, 0, 0, 1, 0, 0, 0, 2, 0};
	Int_t ImCup[9]={0, 0, 0, 0, 1, 0, 0, 0, 0};
	float lumi=137.1;
	Double_t mtt_edges[9]={0,370,420,500,600,700,800,950,2000};
	Double_t ytt_edges[10]={-5.0,-1.4,-0.9,-0.5,-0.15,0.15,0.5,0.9,1.4,5.0};
//	TH2F* h2dist[nsignal+5];//9 signal + 5 background
	TString outputDir="datacard";

	TString cuts[]={"(jet_num == 3)","(jet_num >= 4)"};
	TString cutsName[]={"3jets","4jets"};
	Float_t entries[2][nsample];// number of events in 3jets and 4jets final states 

	ofstream card;
	card.open ("sig.txt");
	card<<std::fixed<<std::setprecision(2); // same effects, but using manipulators	
	card<<"likelihood"<<", "<<cutsName[0]<<", "<<cutsName[1]<<", "<<"semiletonic"<<endl;
	float max_sig,max_likeL;//maximum significance and the corresponds likelihood criteria
	int initial=0;
	for(float likeL=16;likeL<40;likeL++){


		for(int s=0; s<2; s++){ //loop over final states
			int nprocess=0; //count process was dealed with
			TString category="ttbar_"+cutsName[s];
			std::vector<TString> process_names; //names of sigal and bkg
			std::vector<int> process_id;  //process ID; minus and zero for sigal; positive for bkg
			std::vector<TString> bin_arr;   //category name
			std::vector<float> yield_array;  //rate(event yeild)
			std::vector<TString> bkg_norm;  //background  normlization uncertainty
			std::vector<TString> sig_norm;   //signal norlization uncertainty
	    	TH2F* h2dist[nsignal+5];//9 signal + 5 background
			for(int i=0;i<nsample;i++) { //loop over samples
				TChain* chain=new TChain("mytree");
				TChain* chain2=new TChain("rawtree");
				if(s==0 && initial==0)
					fileNames[i]=fileNames[i].ReplaceAll(".root","_*.root");
				chain->Add(dir+fileNames[i]);
				chain2->Add(dir+fileNames[i]);
				Int_t nMC, ncut;
				nMC=chain2->GetEntries();
				ncut=chain->GetEntries();
				//cout<<nMC<<" events simulated and "<<ncut<<" events selected in "<<fileNames[i]<<endl;
				float global_weight=cross_sections[i]*1000*lumi/nMC*K_Factor[i];
				TString condition="(mass_tt<=2000)&&(abs(rapidity_tt)<=5)";
				TString likelihood=Form("(likelihood < %f)",likeL);
				Int_t entry_cut=chain->GetEntries(cuts[s]+"&&"+condition+"&&"+likelihood);
				entries[s][i]=entry_cut*global_weight; //number of events in each channel
				TString sample_name=fileNames[i];
				sample_name=sample_name.ReplaceAll("_*.root","_hist");
				sample_name=sample_name.ReplaceAll("new_","");
				//cout<<Form("entries[%d][%d]: ",s,i)<<entries[s][i]<<endl;

				delete chain;
				delete chain2;
				
			}//end of samples loop
			
			
			
		}//end of cut loop

		Float_t nttbar_3jet=entries[0][0]+entries[0][1]+entries[0][2];
		Float_t nttbar_4jet=entries[1][0]+entries[1][1]+entries[1][2];
		
		float total[2]={0.0,0.0};
		for(int i=0; i<2; i++){
			for(int j=0; j<nsample; j++){
				total[i]=total[i]+entries[i][j];

			}
		}
	    float sig_3jet=nttbar_3jet/sqrt(total[0]);
	    float sig_4jet=nttbar_4jet/sqrt(total[1]);
	    float sig_semi=(nttbar_4jet+nttbar_3jet)/sqrt((total[0]+total[1]));
		
		card<<likeL<<", "<<sig_3jet<<", "<<sig_4jet<<", "<<sig_semi<<endl;
		cout<<std::fixed<<std::setprecision(2); 
		cout<<likeL<<", "<<sig_3jet<<", "<<sig_4jet<<", "<<sig_semi<<endl;

		if(initial==0){
			max_sig=sig_semi;
			max_likeL=likeL;
			initial++;
		}
		else{
			if(max_sig < sig_semi)
				max_sig=sig_semi;
				max_likeL=likeL;
		}

	}	//end of likelihood scan
	card.close();
	cout<<"significance is summarized in sig.txt"<<endl;
	cout<<"at likelihood creteria "<<max_likeL<<", the significance of semi has maximum"<<endl;


}