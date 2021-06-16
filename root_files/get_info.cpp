//get useful information from NanoAOD files for ttbar analysis
// each event contains the information for top quark pairs and their decay products
// written by Renqi Pan in 16th June, 2021.

void get_info(){
	TChain chain("Events");
	TString inputFile="ttbar_semi1.root";
	chain.Add("ttbar_semi1.root");
	TString output="new_"+inputFile;
  	TFile *file=new TFile(output,"RECREATE");
  	TTree *mytree=new TTree("mytree"," tree with branches");
  	Int_t nevents=0; //count the number of events written in tree

    cout<<inputFile<<" is reading and processing"<<endl;
    cout<<"total number of events: "<<chain.GetEntries()<<endl;
    cout<<"get information of ttbar process at parton level"<<endl;
    //get parton level infomation from LHEPart for ttbar process
  	if(inputFile.Contains("ttbar_semi")||inputFile.Contains("TTToSemiLeptonic")){

	    Float_t LHEPart_eta[9],LHEPart_mass[9],LHEPart_phi[9],LHEPart_pt[9];
	    Int_t LHEPart_pdgId[9],LHEPart_status[9];
	    chain.SetBranchAddress("LHEPart_eta", LHEPart_eta);
	    chain.SetBranchAddress("LHEPart_mass", LHEPart_mass);
	    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
	    chain.SetBranchAddress("LHEPart_pt", LHEPart_pt);
	    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
	    chain.SetBranchAddress("LHEPart_pdgId",LHEPart_pdgId);
	    chain.SetBranchAddress("LHEPart_status",LHEPart_status);

	    Float_t M_tt_gen,delta_rapidity_gen,lep_charge;
	    Float_t top_pt,top_eta,top_mass,top_phi,antitop_pt,antitop_eta,antitop_phi,antitop_mass;
	    Float_t b_pt,b_eta,b_mass,b_phi,antib_pt,antib_eta,antib_phi,antib_mass;
	    Float_t lep_pt,lep_eta,lep_mass,lep_phi,nu_pt,nu_eta,nu_phi,nu_mass;
	    Float_t up_pt,up_eta,up_mass,up_phi,down_pt,down_eta,down_phi,down_mass;

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
	   //loop over entry
	    cout<<"parton level infomation is writing. Please wait for a while"<<endl;
	    for(Int_t entry = 0; entry < chain.GetEntries(); entry++){
	    	    chain.GetEntry(entry);
		    	int index_b,index_antib,index_up,index_down,index_lep,index_nu;
		    	for(int i=0;i<9;i++){
		    		//cout<<Form("LHEPart_pdgId[%d]: ",i)<<LHEPart_pdgId[i]<<endl;
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

	            nevents++;
	    	    mytree->Fill();
	    	}
	}  

	//get final state information at detector level
	mytree->Write();
    cout<<inputFile<<" has "<<chain.GetEntries()<<" events"<<endl;
    cout<<output<<" is created"<<endl;
    cout<<nevents<<" events are written into "<<output<<endl;
}



