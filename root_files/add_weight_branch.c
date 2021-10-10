#include <iostream>
#include <cstdio>
using namespace std;
void add_weight_branch(){
	//add weight(kappa=2) to generator level tree
	
    TH2F* hist[17];
    TString dir="/Users/renqi/Documents/top_pairs/EW_2Dcorrection/";
    Int_t Cpq3[17]={ 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1 };
    Int_t Cpu[17]={  0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 1, 0, 2, 1 };
    Int_t ReCup[17]={0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 2, 1 };
    Int_t ImCup[17]={0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 1, 1, 2, 1 };
                  
    for(Int_t i=0;i<17;i++){
      TString file;
      file.Form("EWci%d%d%d%d.root",Cpq3[i],Cpu[i],ReCup[i],ImCup[i]);
      TFile* fhist=TFile::Open(dir+file);
      hist[i]=(TH2F*)fhist->Get("h2");
    }
    std::vector <TString> files={"new_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root",
                                "semi_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root"};
  for(int k=0;k<files.size();k++){
      TString fileName=files[k];
      if (fileName.Contains("TT")){
          TFile *file=new TFile(fileName,"update");
          TTree *mytree=(TTree*) file->Get("mytree");
          Float_t M_tt,delta_rapidity;
          mytree->SetBranchAddress("M_tt_gen",&M_tt);
          mytree->SetBranchAddress("delta_rapidity_gen",&delta_rapidity);
          
          Float_t weight[17];TBranch* branch[17];
          for(int i=0;i<17;i++){
              TString weight_name=Form("weight_ci%d%d%d%d",Cpq3[i],Cpu[i],ReCup[i],ImCup[i]);
              branch[i]=mytree->Branch(weight_name,&weight[i],weight_name+"/F");
              cout<<weight_name<<endl;
           }  
          Int_t entries=mytree->GetEntries();
          cout<<"total number of events: "<<entries<<endl;
          for(Int_t i=0;i<entries;i++){
            mytree->GetEntry(i);
            for(Int_t i=0;i<17;i++){
                Int_t nbin=hist[i]->FindBin(M_tt,delta_rapidity);
                weight[i]=1.0+hist[i]->GetBinContent(nbin);
                branch[i]->Fill();
                //cout<<"weight[i]: "<<weight[i]<<endl;
              }

          }
          mytree->Write("",TObject::kOverwrite);

            cout<<"number of entries in mytree: "<<mytree->GetEntries()<<endl;
            cout<<"EW weights has been added."<<endl;
            cout<<fileName<<" has been updated"<<endl;
            file->Close();
            
      }
      else
          perror("the sample is not a TT process");       

  }
      
}