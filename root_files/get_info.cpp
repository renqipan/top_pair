//get useful information from NanoAOD files for ttbar analysis
// each event contains the information for top quark pairs and their decay products
// written by Renqi Pan in 16th June, 2021.

void get_info(){
    TChain chain("Events");
	TString inputFile="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root";
	chain.Add(inputFile);	
	TString output="new_"+inputFile;
  	TFile *file=new TFile(output,"RECREATE");
  	TTree *mytree=new TTree("mytree"," tree with branches");
  	TTree *rawtree=new TTree("rawtree","tree without selection");
  	Int_t nevents=0, nevents2=0; //count the number of events written in tree

    cout<<inputFile<<" is reading and processing"<<endl;
    cout<<"total number of events: "<<chain.GetEntries()<<endl;
    Float_t LHEPart_eta[9],LHEPart_mass[9],LHEPart_phi[9],LHEPart_pt[9];
	Int_t LHEPart_pdgId[9],LHEPart_status[9];
	UInt_t nLHEPart;
	Float_t M_tt_gen,delta_rapidity_gen,lep_charge;
	Float_t top_pt,top_eta,top_mass,top_phi,antitop_pt,antitop_eta,antitop_phi,antitop_mass;
	Float_t b_pt,b_eta,b_mass,b_phi,antib_pt,antib_eta,antib_phi,antib_mass;
	Float_t lep_pt,lep_eta,lep_mass,lep_phi,nu_pt,nu_eta,nu_phi,nu_mass;
	Float_t up_pt,up_eta,up_mass,up_phi,down_pt,down_eta,down_phi,down_mass;
	//difine branch for ttbar process information at parton level
  	if(inputFile.Contains("ttbar_semi")||inputFile.Contains("TTToSemiLeptonic")){ 
	    chain.SetBranchAddress("LHEPart_eta", LHEPart_eta);
	    chain.SetBranchAddress("LHEPart_mass", LHEPart_mass);
	    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
	    chain.SetBranchAddress("LHEPart_pt", LHEPart_pt);
	    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
	    chain.SetBranchAddress("LHEPart_pdgId",LHEPart_pdgId);
	    chain.SetBranchAddress("nLHEPart",&nLHEPart);
	   // chain.SetBranchAddress("LHEPart_status",LHEPart_status);
	    mytree->Branch("M_tt_gen",&M_tt_gen,"M_tt_gen/F");
	    mytree->Branch("delta_rapidity_gen",&delta_rapidity_gen,"delta_rapidity_gen/F");
	    mytree->Branch("lep_charge",&lep_charge,"lep_charge/F");
	    mytree->Branch("top_pt",&top_pt,"top_pt/F");
	    mytree->Branch("top_eta",&top_eta,"top_eta/F");
	    mytree->Branch("top_phi",&top_phi,"top_phi/F");
	    mytree->Branch("top_mass",&top_mass,"top_mass/F");
	    mytree->Branch("antitop_pt",&antitop_pt,"antitop_pt/F");
	    mytree->Branch("antitop_eta",&antitop_eta,"antitop_eta/F");
	    mytree->Branch("antitop_phi",&antitop_phi,"antitop_phi/F");
	    mytree->Branch("antitop_mass",&antitop_mass,"antitop_mass/F");
	    mytree->Branch("up_pt",&up_pt,"up_pt/F");
	    mytree->Branch("up_eta",&up_eta,"up_eta/F");
	    mytree->Branch("up_phi",&up_phi,"up_phi/F");
	    mytree->Branch("up_mass",&up_mass,"up_mass/F");    
	    mytree->Branch("down_pt",&down_pt,"down_pt/F");
	    mytree->Branch("down_eta",&down_eta,"down_eta/F");
	    mytree->Branch("down_phi",&down_phi,"down_phi/F");
	    mytree->Branch("down_mass",&down_mass,"down_mass/F");
	    mytree->Branch("b_pt",&b_pt,"b_pt/F");
	    mytree->Branch("b_eta",&b_eta,"b_eta/F");
	    mytree->Branch("b_phi",&b_phi,"b_phi/F");
	    mytree->Branch("b_mass",&b_mass,"b_mass/F");
	    mytree->Branch("antib_pt",&antib_pt,"antib_pt/F");
	    mytree->Branch("antib_eta",&antib_eta,"antib_eta/F");
	    mytree->Branch("antib_phi",&antib_phi,"antib_phi/F");
	    mytree->Branch("antib_mass",&antib_mass,"antib_mass/F");
	    mytree->Branch("lep_pt",&lep_pt,"lep_pt/F");
	    mytree->Branch("lep_eta",&lep_eta,"lep_eta/F");
	    mytree->Branch("lep_phi",&lep_phi,"lep_phi/F");
	    mytree->Branch("lep_mass",&lep_mass,"lep_mass/F");
	    mytree->Branch("nu_pt",&nu_pt,"nu_pt/F");
	    mytree->Branch("nu_eta",&nu_eta,"nu_eta/F");
	    mytree->Branch("nu_phi",&nu_phi,"nu_phi/F");
	    mytree->Branch("nu_mass",&nu_mass,"nu_mass/F");
	}
/////////////////////////////////////////////////////////
  //difine branch for final state at detector level
	Float_t MET_pt, MET_phi;
	Float_t Electron_eta[8], Electron_mass[8],Electron_pt[8],Electron_phi[8];
	Float_t Muon_mass[9],Muon_phi[9],Muon_pt[9],Muon_eta[9];
	Float_t lepton_mass[17], lepton_phi[17],lepton_eta[17],lepton_pt[17];
	UInt_t nMuon, nElectron, nJet,nlepton,Jet_btaged[35],nBtag;
	Int_t Jet_partonFlavour[35],Muon_charge[9],Electron_charge[8],lepton_charge[17];
	Float_t Jet_btagCSVV2[35], Jet_btagDeepB[35],Jet_eta[35],Jet_mass[35],Jet_phi[35],Jet_pt[35];
	chain.SetBranchAddress("MET_pt",&MET_pt);
	chain.SetBranchAddress("MET_phi",&MET_phi);
	chain.SetBranchAddress("Electron_phi",Electron_phi);
	chain.SetBranchAddress("Electron_pt",Electron_pt);
	chain.SetBranchAddress("Electron_mass",Electron_mass);
	chain.SetBranchAddress("Electron_eta",Electron_eta);
	chain.SetBranchAddress("nElectron",&nElectron);
	chain.SetBranchAddress("Electron_charge",Electron_charge);
	chain.SetBranchAddress("nMuon",&nMuon);
	chain.SetBranchAddress("nJet",&nJet);
	chain.SetBranchAddress("Muon_eta",Muon_eta);
	chain.SetBranchAddress("Muon_pt",Muon_pt);
	chain.SetBranchAddress("Muon_phi",Muon_phi);
	chain.SetBranchAddress("Muon_mass",Muon_mass);
	chain.SetBranchAddress("Muon_charge",Muon_charge);
	chain.SetBranchAddress("Jet_partonFlavour",Jet_partonFlavour);
	chain.SetBranchAddress("Jet_btagCSVV2",Jet_btagCSVV2);
	chain.SetBranchAddress("Jet_btagDeepB",Jet_btagDeepB);
	chain.SetBranchAddress("Jet_eta",Jet_eta);
	chain.SetBranchAddress("Jet_pt",Jet_pt);
	chain.SetBranchAddress("Jet_phi",Jet_phi);
	chain.SetBranchAddress("Jet_mass",Jet_mass);
	mytree->Branch("MET_phi",&MET_phi,"MET_phi/F");
	mytree->Branch("MET_pt",&MET_pt,"MET_pt/F");
	mytree->Branch("nlepton",&nlepton,"nlepton/I");
	mytree->Branch("lepton_eta",lepton_eta,"lepton_eta[nlepton]/F");
	mytree->Branch("lepton_pt",lepton_pt,"lepton_pt[nlepton]/F");
	mytree->Branch("lepton_mass",lepton_mass,"lepton_mass[nlepton]/F");
	mytree->Branch("lepton_phi",lepton_phi,"lepton_phi[nlepton]/F");
	mytree->Branch("lepton_charge",lepton_charge,"lepton_charge[nlepton]/I");
	//mytree->Branch("nMuon",&nMuon,"nMuon/F");
	//mytree->Branch("nElectron",&nElectron,"nElectron/F");
	mytree->Branch("nJet",&nJet,"nJet/I");
	mytree->Branch("nBtag",&nBtag,"nBtag/I");
	mytree->Branch("Jet_btaged",Jet_btaged,"Jet_btaged[nJet]/I");
	mytree->Branch("Jet_btagCSVV2",Jet_btagCSVV2,"Jet_btagCSVV2[nJet]/F");
	mytree->Branch("Jet_btagDeepB",Jet_btagDeepB,"Jet_btagDeepB[nJet]/F");
	mytree->Branch("Jet_partonFlavour",Jet_partonFlavour,"Jet_partonFlavour[nJet]/I");
	mytree->Branch("Jet_eta",Jet_eta,"Jet_eta[nJet]/F");
	mytree->Branch("Jet_mass",Jet_mass,"Jet_mass[nJet]/F");
	mytree->Branch("Jet_phi",Jet_phi,"Jet_phi[nJet]/F");
	mytree->Branch("Jet_pt",Jet_pt,"Jet_pt[nJet]/F");
	rawtree->Branch("nJet",&nJet,"nJet/I");
	rawtree->Branch("nlepton",&nlepton,"nlepton/I");
	rawtree->Branch("Jet_pt",Jet_pt,"Jet_pt[nJet]/F");
		   //loop over entry
	    cout<<"infomation is writing. Please wait for a while"<<endl;
	    cout<<"infomation is writing. Please wait for a while"<<endl;
        Int_t njet_need=4;// the at least number of jet of semileptonic final satate

	    for(Int_t entry = 0; entry < chain.GetEntries(); entry++){
	    	chain.GetEntry(entry);
	    	if(inputFile.Contains("ttbar_semi")||inputFile.Contains("TTToSemiLeptonic")){ 
	    	//get information of ttbar process at parton level from LHEPart	
		    	int index_b,index_antib,index_up,index_down,index_lep,index_nu;
		    	for(int i=0;i<nLHEPart;i++){
		    //		cout<<Form("LHEPart_pdgId[%d]: ",i)<<LHEPart_pdgId[i]<<endl;
		    		if(LHEPart_pdgId[i]==5) index_b=i;
		    		else if(LHEPart_pdgId[i]==-5) index_antib=i;
		    		else if(LHEPart_pdgId[i]==2||LHEPart_pdgId[i]==4||LHEPart_pdgId[i]==-2||LHEPart_pdgId[i]==-4) 
		    				index_up=i;
		    		else if(LHEPart_pdgId[i]==1||LHEPart_pdgId[i]==3||LHEPart_pdgId[i]==-1||LHEPart_pdgId[i]==-3) 
		    				index_down=i;
		    		else if(LHEPart_pdgId[i]==11||LHEPart_pdgId[i]==13||LHEPart_pdgId[i]==15||LHEPart_pdgId[i]==-11||LHEPart_pdgId[i]==-13||LHEPart_pdgId[i]==-15)
		    				index_lep=i;
		    		else if(LHEPart_pdgId[i]==12||LHEPart_pdgId[i]==14||LHEPart_pdgId[i]==16||LHEPart_pdgId[i]==-12||LHEPart_pdgId[i]==-14||LHEPart_pdgId[i]==-16)
		    				index_nu=i;
	    	   }
	    	   TLorentzVector p4_b, p4_antib,p4_up,p4_down,p4_lep,p4_nu,p4_top,p4_antitop;
	    	   p4_b.SetPtEtaPhiM(LHEPart_pt[index_b],LHEPart_eta[index_b],LHEPart_phi[index_b],LHEPart_mass[index_b]);
	    	   p4_antib.SetPtEtaPhiM(LHEPart_pt[index_antib],LHEPart_eta[index_antib],LHEPart_phi[index_antib],LHEPart_mass[index_antib]);
	    	   p4_up.SetPtEtaPhiM(LHEPart_pt[index_up],LHEPart_eta[index_up],LHEPart_phi[index_up],LHEPart_mass[index_up]);
	    	   p4_down.SetPtEtaPhiM(LHEPart_pt[index_down],LHEPart_eta[index_down],LHEPart_phi[index_down],LHEPart_mass[index_down]);
	    	   p4_lep.SetPtEtaPhiM(LHEPart_pt[index_lep],LHEPart_eta[index_lep],LHEPart_phi[index_lep],LHEPart_mass[index_lep]);
	    	   p4_nu.SetPtEtaPhiM(LHEPart_pt[index_nu],LHEPart_eta[index_nu],LHEPart_phi[index_nu],LHEPart_mass[index_nu]);
	    	   if(LHEPart_pdgId[index_lep] > 0){
		    	   	p4_antitop=p4_antib+p4_lep+p4_nu;
		    	   	p4_top=p4_b+p4_up+p4_down;
		    	   	lep_charge=-1.0;

	    	   }
	    	   else{p4_top=p4_b+p4_lep+p4_nu;
	    	   	    p4_antitop=p4_antib+p4_up+p4_down;
	    	   	    lep_charge=1.0;
				}
				top_pt=p4_top.Pt();
				top_mass=p4_top.M();
				top_phi=p4_top.Phi();
				top_eta=p4_top.Eta();
				antitop_pt=p4_antitop.Pt();
				antitop_mass=p4_antitop.M();
				antitop_phi=p4_antitop.Phi();
				antitop_eta=p4_antitop.Eta();
				up_pt=p4_up.Pt();
				up_mass=p4_up.M();
				up_phi=p4_up.Phi();
				up_eta=p4_up.Eta();
				down_pt=p4_down.Pt();
				down_mass=p4_down.M();
				down_phi=p4_down.Phi();
				down_eta=p4_down.Eta();
				lep_pt=p4_lep.Pt();
				lep_mass=p4_lep.M();
				lep_phi=p4_lep.Phi();
				lep_eta=p4_lep.Eta();
				nu_pt=p4_nu.Pt();
				nu_mass=p4_nu.M();
				nu_phi=p4_nu.Phi();
				nu_eta=p4_nu.Eta();
				b_pt=p4_b.Pt();
				b_mass=p4_b.M();
				b_phi=p4_b.Phi();
				b_eta=p4_b.Eta();
				antib_pt=p4_antib.Pt();
				antib_mass=p4_antib.M();
				antib_phi=p4_antib.Phi();
				antib_eta=p4_antib.Eta();
				M_tt_gen=(p4_top+p4_antitop).M();
				delta_rapidity_gen=p4_top.Rapidity() - p4_antitop.Rapidity();
			}
			/////////////////////////////////////////////////////

            //get information for final state at detector level
			nlepton=nMuon+nElectron;
			TLorentzVector p4_lepton[17];
			for(int i=0;i<nlepton;i++){
				if(i<nElectron){
				p4_lepton[i].SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
				lepton_charge[i]=Electron_charge[i];
				} 
				else {
					p4_lepton[i].SetPtEtaPhiM(Muon_pt[i-nElectron],Muon_eta[i-nElectron],Muon_phi[i-nElectron],Muon_mass[i-nElectron]);
					lepton_charge[i]=Muon_charge[i-nElectron];			
				}
			}
			//sort leptons according to their pt and make the first jet has maximum pt
			for(int k=0;k<nlepton;k++){
				for(int j=k+1;j<nlepton;j++)
				{
					if(p4_lepton[k].Pt() < p4_lepton[j].Pt())
					{   TLorentzVector temp; Int_t temp_charge;
						temp=p4_lepton[k];
						p4_lepton[k]=p4_lepton[j];
						p4_lepton[j]=temp;
						temp_charge=lepton_charge[k];
						lepton_charge[k]=lepton_charge[j];
						lepton_charge[j]=temp_charge;


					}
				}
			}
			for(int i=0;i<nlepton;i++){
				lepton_pt[i]=p4_lepton[i].Pt();
				lepton_eta[i]=p4_lepton[i].Eta();
				lepton_phi[i]=p4_lepton[i].Phi();
				lepton_mass[i]=p4_lepton[i].M();
			}
			nBtag=0;//count number of bjet 
			for(int i=0; i<nJet;i++){
				if(Jet_btagCSVV2[i] > 0.8)
				{
					Jet_btaged[i]=1;
					nBtag++;
				}
				else Jet_btaged[i]=0;
			}

	        nevents++;
	        rawtree->Fill();
	       //////////////////////////////////////////////////////////////////
	        //select ttbar semiletopnic final state
	        //selec jets
	        bool jet_flag=false;//if true pass the selection
	        int jet_satisfy=0,Nbtag3=0;
	        if(nJet>=njet_need){ //njet_need=4 by default
	        	for(int i=0; i<njet_need;i++){
	        		if(Jet_eta[i] < 2.4 && Jet_pt[i] >30)
	        			{jet_satisfy++;
	        			 Nbtag3=Nbtag3+Jet_btaged[i];
		        		}
	        		else
	        			break;
	        	}
	        	if(jet_satisfy==njet_need && Nbtag3 >=2)
	        		jet_flag=true;
	        }
	        //select lepton
	        bool lepton_flag=false; //if true pass the selction
	        int nlepton_satisfy=1;
	        if(nlepton >=1 && lepton_pt[0] > 30 && lepton_eta[0] < 2.4){
	        	for (int i=1; i<nlepton; i++){ //start from the second lepton
	        		if(lepton_pt[i] > 15 &&lepton_eta[i] > 2.4)
	        			break;
	        		else
	        			nlepton_satisfy++;
	        	}
	        	if(nlepton_satisfy==nlepton)
	        		lepton_flag=true;
	        }

	        //select satisfied events
	        if(jet_flag==true && lepton_flag==true){
	        	mytree->Fill();
	        	nevents2++;
	        }
	    	 
	    	//cout<<"nlepton: "<<nlepton<<" lepton_pt[0]: "<<lepton_pt[0]<<" lepton_pt[1]: "<<lepton_pt[1]<<" lepton_pt[2]: "<<lepton_pt[2]<<endl;
	    	//cout<<"nJet: "<<nJet<<" Jet_pt[0]: "<<Jet_pt[0]<<" Jet_pt[1]: "<<Jet_pt[1]<<" Jet_pt[2]: "<<Jet_pt[2]<<endl;
	    }

	mytree->Write();
	rawtree->Write();
    cout<<inputFile<<" has "<<chain.GetEntries()<<" events"<<endl;
    cout<<output<<" is created"<<endl;
    cout<<nevents<<" events are written into "<<"rawtree."<<endl;
    cout<<nevents2<<" events are written into "<<"mytree."<<endl;

}



