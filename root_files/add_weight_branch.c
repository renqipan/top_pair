#include <iostream>
#include <cstdio>
using namespace std;
void add_weight_branch(TString fileName){
	//add weight(kappa=2) to generator level tree
	
    TH2F* hist[9];
    //TString dir="/afs/cern.ch/user/r/repan/work/top_pair/correction_roots/";
    TString dir="/Users/renqi/Documents/top_pairs/EW_2Dcorrection/";
    Int_t Cpq3[9]={ 0, 1, 0, 0, 0, 2, 0, 0, 0 };
    Int_t Cpu[9]={  0, 0, 1, 0, 0, 0, 2, 0, 0 };
    Int_t ReCup[9]={0, 0, 0, 1, 0, 0, 0, 2, 0 };
    Int_t ImCup[9]={0, 0, 0, 0, 1, 0, 0, 0, 2 };
                  
    for(Int_t i=0;i<9;i++){
      TString file;
      file.Form("EWci%d%d%d%d.root",Cpq3[i],Cpu[i],ReCup[i],ImCup[i]);
      TFile* fhist=TFile::Open(dir+file);
      hist[i]=(TH2F*)fhist->Get("h2");
    }
      TString treeName[]={"mytree","jesUp","jesDown","jerUp","jerDown","unclusup","unclusdown"};
      if (fileName.Contains("TT")){
          TFile *file=new TFile(fileName,"update");
          for(int t=0;t<7;t++){ //loop over trees in a same file
              TTree *mytree=(TTree*) file->Get(treeName[t]);
              Float_t M_tt,delta_rapidity;
              mytree->SetBranchAddress("M_tt_gen",&M_tt);
              mytree->SetBranchAddress("delta_rapidity_gen",&delta_rapidity);
              
              Float_t weight[9];TBranch* branch[9];
              for(int i=0;i<9;i++){
                  TString weight_name=Form("weight_ci%d%d%d%d",Cpq3[i],Cpu[i],ReCup[i],ImCup[i]);
                  branch[i]=mytree->Branch(weight_name,&weight[i],weight_name+"/F");
                  //cout<<weight_name<<endl;
               }  
              Int_t entries=mytree->GetEntries();
              cout<<"total number of events: "<<entries<<endl;
              for(Int_t i=0;i<entries;i++){
                mytree->GetEntry(i);
                for(Int_t i=0;i<9;i++){
                    Int_t nbin=hist[i]->FindBin(M_tt,delta_rapidity);
                    weight[i]=1.0+hist[i]->GetBinContent(nbin);
                    branch[i]->Fill();
                    //cout<<"weight[i]: "<<weight[i]<<endl;
                  }

              }
              mytree->Write("",TObject::kOverwrite);

              cout<<"the tree "<<mytree->GetName()<<" has "<<mytree->GetEntries()<<" entries"<<endl;
              delete mytree;
          }
            cout<<"EW weights has been added."<<endl;
            cout<<fileName<<" has been updated"<<endl;
            file->Close();
            
      }
      else
          perror("the sample is not a TT process");       

 //}
      
}
