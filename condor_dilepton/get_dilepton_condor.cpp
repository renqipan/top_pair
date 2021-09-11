#include "dineu.cpp"
#include "intersect.cpp"
void getbleppair(TLorentzVector mom_b[],TLorentzVector mom_l[],UInt_t bpair[],double metx,double mety,bool& flag){
	bool pair11,pair12,pair21,pair22,pair1,pair2;//1:0to0,1to1;2:0to1,1to0
	TMatrixD N11(3,3),N12(3,3),N21(3,3),N22(3,3);
	TMatrixD gamma(3,3);
	gamma(0,0)=-1;
	gamma(1,1)=-1;
	gamma(2,2)=1;
	gamma(0,2)=metx;
	gamma(1,2)=mety;
	gamma(1,0)=0;
	gamma(0,1)=0;
	gamma(2,0)=0;
	gamma(2,1)=0;
	pair11=NeutrinoSolver(&mom_l[0],&mom_b[0],80.0,172.5,N11);
	pair12=NeutrinoSolver(&mom_l[1],&mom_b[1],80.0,172.5,N12);
	pair21=NeutrinoSolver(&mom_l[0],&mom_b[1],80.0,172.5,N21);
	pair22=NeutrinoSolver(&mom_l[1],&mom_b[0],80.0,172.5,N22);
	pair1=pair11&&pair12;
	pair2=pair21&&pair22;
	if(pair1&&pair2){
		TMatrixD n12(3,3),n22(3,3),tmp(gamma);
		tmp.T();
		n12=tmp*N12*gamma;
		n22=tmp*N22*gamma;
		pair1=intersect1(N11,n12);
		pair2=intersect1(N21,n22);
		if((pair1&&pair2)||(!pair1&&!pair2)){
			if(mom_l[0].DeltaR(mom_b[0])+mom_l[1].DeltaR(mom_b[1])<mom_l[0].DeltaR(mom_b[1])+mom_l[1].DeltaR(mom_b[0])){
				bpair[0]=0;
				bpair[1]=1;
			}
			else{
				bpair[0]=1;
				bpair[1]=0;
			}
		}
		else if(pair1){
			bpair[0]=0;
			bpair[1]=1;
		}
		else{
			bpair[0]=1;
			bpair[1]=0;
		}
	}
	else if(pair1){
		bpair[0]=0;
		bpair[1]=1;
	}
	else if(pair2){
		bpair[0]=1;
		bpair[1]=0;
	}
	else{
		flag=false;
	}
	return;
}
void get_dilepton_condor(TString inputFile){
	TString dir="root://cms-xrd-global.cern.ch/";
	TChain chain("Events");
	ifstream infile;
	infile.open(inputFile);
	if(!infile.is_open()){
	    cout<<"error: reading input file failed!!! "<<endl;
	    exit(1);
	}
	else{
	       std::string line;
	       while (getline(infile,line)){
	       TString file_dir(line);
	        if(file_dir.BeginsWith("/store"))
	          chain.Add(dir+file_dir);
	        else throw std::invalid_argument("error dasgoclient root file lists");
	        }
	}
	TString process=inputFile.ReplaceAll(".txt",".root");
	TString output = "dilepton_"+process;
  	TFile *file = new TFile(output, "RECREATE");
  	TTree *mytree = new TTree("mytree", " tree with branches");
  	TTree *rawtree = new TTree("rawtree", "tree without selection");
	Int_t nevents = 0, nevents2 = 0; // count the number of events written in tree

	cout << inputFile << " is reading and processing" << endl;
	cout << "total number of events: " << chain.GetEntries() << endl;
	TH2F* hist[8];// hists for  weights of EW corrections 
	if(inputFile.Contains("TTToSemiLeptonic")||inputFile.Contains("TTTo2L2Nu")||inputFile.Contains("TTToHadronic")){
		TString dir="/afs/cern.ch/user/r/repan/work/top_pair/correction_roots/";
		TString files[8]={"correction_kappa10","correction_kappa20","correction_kappa30",
	              "correction_kappa01","correction_kappa02","correction_kappa11",
	              "correction_kappa22","correction_kappa00"};
		for(Int_t i=0;i<8;i++){
	    	TFile* fhist=TFile::Open(dir+files[i]+".root");
	   		hist[i]=(TH2F*)fhist->Get("h2");
		}             
	}
  ////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////
	Float_t LHEPart_eta[9], LHEPart_mass[9], LHEPart_phi[9], LHEPart_pt[9];
	Int_t LHEPart_pdgId[9], LHEPart_status[9];
	UInt_t nLHEPart;

	UInt_t LHE_nlep = 0, LHE_nhad = 0, LHE_tao = 0;
	Float_t M_tt_gen, delta_rapidity_gen;
	Float_t deltaybb_gen,deltayll_gen,deltaybl_gen,deltayblbar_gen,deltaybbll_gen,Mbbll_gen;
	Float_t top_pt, top_eta, top_mass, top_phi, antitop_pt, antitop_eta,
	  antitop_phi, antitop_mass;
	Float_t b_pt, b_eta, b_mass, b_phi, antib_pt, antib_eta, antib_phi,
	  antib_mass;
	Float_t lep_pt[2], lep_eta[2], lep_mass[2], lep_phi[2], nu_pt[2], nu_eta[2], nu_phi[2], nu_mass[2];//0:lepn,nup;1:lepp,nun
	// Int_t LHE_had[6];
	// difine branch for ttbar process information at parton level
	if(inputFile.Contains("TTToSemiLeptonic")||inputFile.Contains("TTTo2L2Nu")||inputFile.Contains("TTToHadronic")){
    chain.SetBranchAddress("LHEPart_eta", LHEPart_eta);
    chain.SetBranchAddress("LHEPart_mass", LHEPart_mass);
    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
    chain.SetBranchAddress("LHEPart_pt", LHEPart_pt);
    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
    chain.SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId);
    chain.SetBranchAddress("nLHEPart", &nLHEPart);
    mytree->Branch("M_tt_gen", &M_tt_gen, "M_tt_gen/F");
    mytree->Branch("delta_rapidity_gen", &delta_rapidity_gen,"delta_rapidity_gen/F");
  }
	if (inputFile.Contains("TTTo2L2Nu")) {
		mytree->Branch("top_pt", &top_pt, "top_pt/F");
		mytree->Branch("top_eta", &top_eta, "top_eta/F");
		mytree->Branch("top_phi", &top_phi, "top_phi/F");
		mytree->Branch("top_mass", &top_mass, "top_mass/F");
		mytree->Branch("antitop_pt", &antitop_pt, "antitop_pt/F");
		mytree->Branch("antitop_eta", &antitop_eta, "antitop_eta/F");
		mytree->Branch("antitop_phi", &antitop_phi, "antitop_phi/F");
		mytree->Branch("antitop_mass", &antitop_mass, "antitop_mass/F");
		mytree->Branch("b_pt", &b_pt, "b_pt/F");
		mytree->Branch("b_eta", &b_eta, "b_eta/F");
		mytree->Branch("b_phi", &b_phi, "b_phi/F");
		mytree->Branch("b_mass", &b_mass, "b_mass/F");
		mytree->Branch("antib_pt", &antib_pt, "antib_pt/F");
		mytree->Branch("antib_eta", &antib_eta, "antib_eta/F");
		mytree->Branch("antib_phi", &antib_phi, "antib_phi/F");
		mytree->Branch("antib_mass", &antib_mass, "antib_mass/F");
		mytree->Branch("lep_pt", lep_pt, "lep_pt[2]/F");
		mytree->Branch("lep_eta", lep_eta, "lep_eta[2]/F");
		mytree->Branch("lep_phi", lep_phi, "lep_phi[2]/F");
		mytree->Branch("lep_mass", lep_mass, "lep_mass[2]/F");
		mytree->Branch("nu_pt", nu_pt, "nu_pt[2]/F");
		mytree->Branch("nu_eta", nu_eta, "nu_eta[2]/F");
		mytree->Branch("nu_phi", nu_phi, "nu_phi[2]/F");
		mytree->Branch("nu_mass", nu_mass, "nu_mass[2]/F");
		mytree->Branch("deltaybb_gen",&deltaybb_gen,"deltabb_gen/F");
		mytree->Branch("deltayll_gen",&deltayll_gen,"deltall_gen/F");
		mytree->Branch("deltaybl_gen",&deltaybl_gen,"deltabl_gen/F");
		mytree->Branch("deltayblbar_gen",&deltayblbar_gen,"deltayblbar_gen/F");
		mytree->Branch("deltaybbll_gen",&deltaybbll_gen,"deltaybbll_gen/F");
		mytree->Branch("Mbbll_gen",&Mbbll_gen,"Mbbll_gen/F");
	}
	Float_t weight[8];
	if(inputFile.Contains("TTToSemiLeptonic")||inputFile.Contains("TTTo2L2Nu")||inputFile.Contains("TTToHadronic")){	
		mytree->Branch("weight_kappa10",&weight[0],"weight_kappa10/F");
		mytree->Branch("weight_kappa20",&weight[1],"weight_kappa20/F");
		mytree->Branch("weight_kappa30",&weight[2],"weight_kappa30/F");
		mytree->Branch("weight_kappa01",&weight[3],"weight_kappa01/F");
		mytree->Branch("weight_kappa02",&weight[4],"weight_kappa02/F");
		mytree->Branch("weight_kappa11",&weight[5],"weight_kappa11/F");
		mytree->Branch("weight_kappa22",&weight[6],"weight_kappa22/F");
		mytree->Branch("weight_kappa00",&weight[7],"weight_kappa00/F");

	}
	/////////////////////////////////////////////////////////
	// difine branch for final state at detector level
	Float_t MET_pt, MET_phi;
	Float_t Electron_eta[9], Electron_mass[9], Electron_pt[9], Electron_phi[9];
	Float_t Muon_mass[11], Muon_phi[11], Muon_pt[11], Muon_eta[11];
	Float_t lepton_mass[20], lepton_phi[20], lepton_eta[20], lepton_pt[20];
	Float_t slepton_mass[2], slepton_phi[2], slepton_eta[2], slepton_pt[2];
	UInt_t nMuon, nElectron, nJet, nlepton, Jet_btaged[45], nBtag;
	Int_t Jet_partonFlavour[45], Muon_charge[11], Electron_charge[9],lepton_charge[20],slepton_charge[2];
	Float_t Jet_btagCSVV2[45], Jet_eta[45], Jet_mass[45], Jet_phi[45], Jet_pt[45],Jet_btagDeepB[45];
	Float_t jet_btagCSVV2[2], jet_eta[2], jet_mass[2], jet_phi[2], jet_pt[2],jet_btagDeepB[2];
	Float_t MtW;
	Int_t jet_index[45];
	Float_t deltaybb,deltayll,deltaybbll,Mbbll;
	Float_t Electron_deltaEtaSC[9],Electron_dxy[9], Electron_dz[9],Muon_pfRelIso04_all[11];
  	Int_t Electron_cutBased[9], Jet_jetId[45];
  	Int_t PV_npvsGood;
  	Bool_t Muon_tightId[11], Muon_looseId[11];
	chain.SetBranchAddress("Electron_deltaEtaSC",Electron_deltaEtaSC);
	chain.SetBranchAddress("Electron_dz",Electron_dz);
	chain.SetBranchAddress("Electron_dxy",Electron_dxy);
	chain.SetBranchAddress("PV_npvsGood",&PV_npvsGood);
	chain.SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all);
	chain.SetBranchAddress("Muon_tightId",Muon_tightId);
	chain.SetBranchAddress("Muon_looseId",Muon_looseId);
	chain.SetBranchAddress("Electron_cutBased",Electron_cutBased);
	chain.SetBranchAddress("Jet_jetId",Jet_jetId);
	chain.SetBranchAddress("MET_pt", &MET_pt);
	chain.SetBranchAddress("MET_phi", &MET_phi);
	chain.SetBranchAddress("Electron_phi", Electron_phi);
	chain.SetBranchAddress("Electron_pt", Electron_pt);
	chain.SetBranchAddress("Electron_mass", Electron_mass);
	chain.SetBranchAddress("Electron_eta", Electron_eta);
	chain.SetBranchAddress("nElectron", &nElectron);
	chain.SetBranchAddress("Electron_charge", Electron_charge);
	chain.SetBranchAddress("nMuon", &nMuon);
	chain.SetBranchAddress("nJet", &nJet);
	chain.SetBranchAddress("Muon_eta", Muon_eta);
	chain.SetBranchAddress("Muon_pt", Muon_pt);
	chain.SetBranchAddress("Muon_phi", Muon_phi);
	chain.SetBranchAddress("Muon_mass", Muon_mass);
	chain.SetBranchAddress("Muon_charge", Muon_charge);
	chain.SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour);
	chain.SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2);
	chain.SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB);
	chain.SetBranchAddress("Jet_eta", Jet_eta);
	chain.SetBranchAddress("Jet_pt", Jet_pt);
	chain.SetBranchAddress("Jet_phi", Jet_phi);
	chain.SetBranchAddress("Jet_mass", Jet_mass);
	mytree->Branch("MET_phi", &MET_phi, "MET_phi/F");
	mytree->Branch("MET_pt", &MET_pt, "MET_pt/F");
	mytree->Branch("slepton_eta", slepton_eta, "slepton_eta[2]/F");
	mytree->Branch("slepton_pt", slepton_pt, "slepton_pt[2]/F");
	mytree->Branch("slepton_mass", slepton_mass, "slepton_mass[2]/F");
	mytree->Branch("slepton_phi", slepton_phi, "slepton_phi[2]/F");
	mytree->Branch("slepton_charge", slepton_charge, "slepton_charge[2]/I");
	mytree->Branch("jet_btagCSVV2", jet_btagCSVV2, "jet_btagCSVV2[2]/F");
	mytree->Branch("jet_btagDeepB", jet_btagDeepB, "jet_btagDeepB[2]/F");
	mytree->Branch("jet_eta", jet_eta, "jet_eta[2]/F");
	mytree->Branch("jet_mass", jet_mass, "jet_mass[2]/F");
	mytree->Branch("jet_phi", jet_phi, "jet_phi[2]/F");
	mytree->Branch("jet_pt", jet_pt, "jet_pt[2]/F");
	mytree->Branch("deltaybb",&deltaybb,"deltabb/F");
	mytree->Branch("deltayll",&deltayll,"deltall/F");
	mytree->Branch("deltaybbll",&deltaybbll,"deltabbll/F");
	//mytree->Branch("deltayblbar",&deltayblbar,"deltablbar/F");
	mytree->Branch("Mbbll",&Mbbll,"Mbbll/F");
	rawtree->Branch("nJet", &nJet, "nJet/I");
	rawtree->Branch("nlepton", &nlepton, "nlepton/I");
	rawtree->Branch("Jet_pt", Jet_pt, "Jet_pt[nJet]/F");

	////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////
	// loop over entry
	cout << "infomation is writing. Please wait for a while" << endl;
	cout << "infomation is writing. Please wait for a while" << endl;
	/////////////////////////////////////////////////////////////////////
	int total_entry=chain.GetEntries();
  	if(inputFile.Contains("TTTo2L2Nu")||inputFile.Contains("TT_Tune")){
    	if(total_entry > 8.E+7)
       		total_entry=8.E+7;
   }
    else if (inputFile.Contains("QCD_HT100to200")||inputFile.Contains("QCD_HT200to300")||
    	inputFile.Contains("QCD_HT300to500")||inputFile.Contains("QCD_Pt-15to7000")){
        if(total_entry > 1.E+9)
          total_entry=1.E+9;
  	}
    else{
    	if(total_entry > 1.E+8)
      		total_entry=1.E+8;
  	}
	for (Int_t entry = 0; entry < total_entry; entry++){
	    chain.GetEntry(entry);
	    int index_b, index_antib, index_lepn, index_nun, index_lepp, index_nup;
	    TLorentzVector p4_b, p4_antib, p4_lepn, p4_nun, p4_lepp, p4_nup, p4_top,p4_antitop;
	    if(inputFile.Contains("TTTo2L2Nu")){
	      // get information of ttbar process at parton level from LHEPart
	      //  int index_b,index_antib,index_up,index_down,index_lep,index_nu;
	      	for(int i = 0; i < nLHEPart; i++) {
	        //    cout<<Form("LHEPart_pdgId[%d]:
	        //",i)<<LHEPart_pdgId[i]<<endl;
				if (LHEPart_pdgId[i] == 5) index_b = i;
				else if (LHEPart_pdgId[i] == -5) index_antib = i;
				else if (LHEPart_pdgId[i] == 11|| LHEPart_pdgId[i] == 13||LHEPart_pdgId[i] == 15)
					index_lepn = i;
				else if (LHEPart_pdgId[i] == 12 || LHEPart_pdgId[i] == 14 || LHEPart_pdgId[i] == 16)
				  	index_nun = i;
				else if (LHEPart_pdgId[i] == -11 || LHEPart_pdgId[i] == -13 || LHEPart_pdgId[i] == -15)
				  	index_lepp = i;
				else if (LHEPart_pdgId[i] == -12 ||LHEPart_pdgId[i] == -14 || LHEPart_pdgId[i] == -16)
				  	index_nup = i;
	        }	
	      //  TLorentzVector p4_b,
	      //  p4_antib,p4_up,p4_down,p4_lep,p4_nu,p4_top,p4_antitop;
			p4_b.SetPtEtaPhiM(LHEPart_pt[index_b], LHEPart_eta[index_b],LHEPart_phi[index_b], LHEPart_mass[index_b]);
			p4_antib.SetPtEtaPhiM(LHEPart_pt[index_antib], LHEPart_eta[index_antib],LHEPart_phi[index_antib],LHEPart_mass[index_antib]);
			p4_lepp.SetPtEtaPhiM(LHEPart_pt[index_lepp], LHEPart_eta[index_lepp],LHEPart_phi[index_lepp], LHEPart_mass[index_lepp]);
			p4_lepn.SetPtEtaPhiM(LHEPart_pt[index_lepn], LHEPart_eta[index_lepn],LHEPart_phi[index_lepn], LHEPart_mass[index_lepn]);
			p4_nup.SetPtEtaPhiM(LHEPart_pt[index_nup], LHEPart_eta[index_nup],LHEPart_phi[index_nup], LHEPart_mass[index_nup]);
			p4_nun.SetPtEtaPhiM(LHEPart_pt[index_nun], LHEPart_eta[index_nun],LHEPart_phi[index_nun], LHEPart_mass[index_nun]);
			p4_top = p4_b + p4_lepp + p4_nun;
			p4_antitop = p4_antib + p4_lepn + p4_nup;
			top_pt = p4_top.Pt();
			top_mass = p4_top.M();
			top_phi = p4_top.Phi();
			top_eta = p4_top.Eta();
			antitop_pt = p4_antitop.Pt();
			antitop_mass = p4_antitop.M();
			antitop_phi = p4_antitop.Phi();
			antitop_eta = p4_antitop.Eta();
			lep_pt[0] = p4_lepn.Pt();
			lep_mass[0] = p4_lepn.M();
			lep_phi[0] = p4_lepn.Phi();
			lep_eta[0] = p4_lepn.Eta();
			nu_pt[0] = p4_nup.Pt();
			nu_mass[0] = p4_nup.M();
			nu_phi[0] = p4_nup.Phi();
			nu_eta[0] = p4_nup.Eta();
			lep_pt[1] = p4_lepp.Pt();
			lep_mass[1] = p4_lepp.M();
			lep_phi[1] = p4_lepp.Phi();
			lep_eta[1] = p4_lepp.Eta();
			nu_pt[1] = p4_nun.Pt();
			nu_mass[1] = p4_nun.M();
			nu_phi[1] = p4_nun.Phi();
			nu_eta[1] = p4_nun.Eta();
			b_pt = p4_b.Pt();
			b_mass = p4_b.M();
			b_phi = p4_b.Phi();
			b_eta = p4_b.Eta();
			antib_pt = p4_antib.Pt();
			antib_mass = p4_antib.M();
			antib_phi = p4_antib.Phi();
			antib_eta = p4_antib.Eta();
			M_tt_gen = (p4_top + p4_antitop).M();
			delta_rapidity_gen = p4_top.Rapidity() - p4_antitop.Rapidity();
			deltaybb_gen=p4_b.Rapidity()-p4_antib.Rapidity();
			deltaybl_gen=p4_b.Rapidity()-p4_lepp.Rapidity();
			deltaybbll_gen=(p4_b+p4_lepp).Rapidity()-(p4_antib+p4_lepn).Rapidity();
			deltayll_gen=p4_lepn.Rapidity()-p4_lepp.Rapidity();
			deltayblbar_gen=p4_antib.Rapidity()-p4_lepn.Rapidity();
			Mbbll_gen=(p4_b+p4_antib+p4_lepn+p4_lepp).M();
	    }
	 	if (inputFile.Contains("TTToSemiLeptonic")) {
	    	int index_b, index_antib, index_up, index_down, index_lep, index_nu;
	    	TLorentzVector p4_b, p4_antib, p4_up, p4_down, p4_lep, p4_nu, p4_top,p4_antitop;
	      // get information of ttbar process at parton level from LHEPart
	      //  int index_b,index_antib,index_up,index_down,index_lep,index_nu;
	      for (int i = 0; i < nLHEPart; i++) {
	        //    cout<<Form("LHEPart_pdgId[%d]: ",i)<<LHEPart_pdgId[i]<<endl;
	        if (LHEPart_pdgId[i] == 5)
	          index_b = i;
	        else if (LHEPart_pdgId[i] == -5)
	          index_antib = i;
	        else if (LHEPart_pdgId[i] == 2 || LHEPart_pdgId[i] == 4 ||
	                 LHEPart_pdgId[i] == -2 || LHEPart_pdgId[i] == -4)
	          index_up = i;
	        else if (LHEPart_pdgId[i] == 1 || LHEPart_pdgId[i] == 3 ||
	                 LHEPart_pdgId[i] == -1 || LHEPart_pdgId[i] == -3)
	          index_down = i;
	        else if (LHEPart_pdgId[i] == 11 || LHEPart_pdgId[i] == 13 ||
	                 LHEPart_pdgId[i] == 15 || LHEPart_pdgId[i] == -11 ||
	                 LHEPart_pdgId[i] == -13 || LHEPart_pdgId[i] == -15)
	          index_lep = i;
	        else if (LHEPart_pdgId[i] == 12 || LHEPart_pdgId[i] == 14 ||
	                 LHEPart_pdgId[i] == 16 || LHEPart_pdgId[i] == -12 ||
	                 LHEPart_pdgId[i] == -14 || LHEPart_pdgId[i] == -16)
	          index_nu = i;
      		}

	      p4_b.SetPtEtaPhiM(LHEPart_pt[index_b], LHEPart_eta[index_b],LHEPart_phi[index_b], LHEPart_mass[index_b]);
	      p4_antib.SetPtEtaPhiM(LHEPart_pt[index_antib], LHEPart_eta[index_antib],LHEPart_phi[index_antib],LHEPart_mass[index_antib]);
	      p4_up.SetPtEtaPhiM(LHEPart_pt[index_up], LHEPart_eta[index_up],LHEPart_phi[index_up], LHEPart_mass[index_up]);
	      p4_down.SetPtEtaPhiM(LHEPart_pt[index_down], LHEPart_eta[index_down],LHEPart_phi[index_down], LHEPart_mass[index_down]);
	      p4_lep.SetPtEtaPhiM(LHEPart_pt[index_lep], LHEPart_eta[index_lep],LHEPart_phi[index_lep], LHEPart_mass[index_lep]);
	      p4_nu.SetPtEtaPhiM(LHEPart_pt[index_nu], LHEPart_eta[index_nu],LHEPart_phi[index_nu], LHEPart_mass[index_nu]);
	      if (LHEPart_pdgId[index_lep] > 0) {
	        p4_antitop = p4_antib + p4_lep + p4_nu;
	        p4_top = p4_b + p4_up + p4_down;

	      } 
	      else {
	        p4_top = p4_b + p4_lep + p4_nu;
	        p4_antitop = p4_antib + p4_up + p4_down;
	      }
	      
	      M_tt_gen = (p4_top + p4_antitop).M();
	      delta_rapidity_gen = p4_top.Rapidity() - p4_antitop.Rapidity();
    }

    if(inputFile.Contains("TTToHadronic")){
        int index_b, index_antib, index_upbar, index_up, index_downbar, index_down;
        TLorentzVector p4_b, p4_antib, p4_upbar, p4_up, p4_downbar, p4_down, p4_top,p4_antitop;
        for(int i = 0; i < nLHEPart; i++) {
            if (LHEPart_pdgId[i] == 5) index_b = i;
            else if (LHEPart_pdgId[i] == -5) index_antib = i;
            else if (LHEPart_pdgId[i] == 2|| LHEPart_pdgId[i] == 4)
                index_up = i;
            else if (LHEPart_pdgId[i] == -2 || LHEPart_pdgId[i] == -4)
                index_upbar = i;
            else if (LHEPart_pdgId[i] == 1 || LHEPart_pdgId[i] == 3)
                index_down = i;
            else if (LHEPart_pdgId[i] == -1 ||LHEPart_pdgId[i] == -3)
                index_downbar = i;
          } 
      p4_b.SetPtEtaPhiM(LHEPart_pt[index_b], LHEPart_eta[index_b],LHEPart_phi[index_b], LHEPart_mass[index_b]);
      p4_antib.SetPtEtaPhiM(LHEPart_pt[index_antib], LHEPart_eta[index_antib],LHEPart_phi[index_antib],LHEPart_mass[index_antib]);
      p4_up.SetPtEtaPhiM(LHEPart_pt[index_up], LHEPart_eta[index_up],LHEPart_phi[index_up], LHEPart_mass[index_up]);
      p4_upbar.SetPtEtaPhiM(LHEPart_pt[index_upbar], LHEPart_eta[index_upbar],LHEPart_phi[index_upbar], LHEPart_mass[index_upbar]);
      p4_down.SetPtEtaPhiM(LHEPart_pt[index_down], LHEPart_eta[index_down],LHEPart_phi[index_down], LHEPart_mass[index_down]);
      p4_downbar.SetPtEtaPhiM(LHEPart_pt[index_downbar], LHEPart_eta[index_downbar],LHEPart_phi[index_downbar], LHEPart_mass[index_downbar]);
      p4_top = p4_b + p4_up + p4_downbar;
      p4_antitop = p4_antib + p4_upbar + p4_down;      
      M_tt_gen = (p4_top + p4_antitop).M();
      delta_rapidity_gen = p4_top.Rapidity() - p4_antitop.Rapidity();
    }
	    /////////////////////////////////////////////////////
	    // get information for final state at detector level
	    nlepton = nMuon + nElectron;
	    TLorentzVector p4_lepton[20];
	    for (int i = 0; i < nlepton; i++) {
	      	if (i < nElectron) {
	        	p4_lepton[i].SetPtEtaPhiM(Electron_pt[i], Electron_eta[i],Electron_phi[i], Electron_mass[i]);
	        	lepton_charge[i] = Electron_charge[i];
	      	}
	      	else{
	        	p4_lepton[i].SetPtEtaPhiM(Muon_pt[i-nElectron], Muon_eta[i-nElectron],Muon_phi[i-nElectron], Muon_mass[i-nElectron]);
	        	lepton_charge[i] = Muon_charge[i-nElectron];
	      	}
	    }
	    for (int i = 0; i < nlepton; i++) {
			lepton_pt[i] = p4_lepton[i].Pt();
			lepton_eta[i] = p4_lepton[i].Eta();
			lepton_phi[i] = p4_lepton[i].Phi();
			lepton_mass[i] = p4_lepton[i].M();
	    }
	    nevents++;
	    rawtree->Fill();
	    ////////////////////////////////////////////////////////////////////
	    bool lepton_flag = false; // if true pass the selction
    	UInt_t lepton_index[20];
    	TLorentzVector mom_lep[2];
    	int num_select1=0, num_select2=0,num_select3=0;
    	for (int i = 0; i < nlepton; i++) {
      		if(i<nElectron){
          		if(Electron_cutBased[i]==4&& abs(Electron_eta[i]) <2.4 && (abs(Electron_eta[i])<1.4442||abs(Electron_eta[i])>1.5660)&&Electron_pt[i]>20){
            		if((abs(Electron_deltaEtaSC[i]+Electron_eta[i])<1.479&&abs(Electron_dxy[i])<0.05&&abs(Electron_dz[i])<0.1)
						||(abs(Electron_deltaEtaSC[i]+Electron_eta[i])>=1.479&&abs(Electron_dxy[i])<0.1&&abs(Electron_dz[i])<0.2)){	
            				lepton_index[num_select2]=i;
            				num_select2++;
            				if(Electron_pt[i]>30)
            					num_select3++;
            		}
            	}
            }       
     		else {
             	if(Muon_tightId[i-nElectron]==1&&Muon_pfRelIso04_all[i-nElectron]<0.15&&Muon_pt[i-nElectron]>20&&abs(Muon_eta[i-nElectron])<2.4){ 
              		lepton_index[num_select2]=i;
              		num_select2++;
              		if(Muon_pt[i-nElectron]>30)
              			num_select3++;
              	}
            }
        }
        if(num_select2>=2&&num_select3>=1) {
      		lepton_flag=true;
  		}
        for(int i = 0; i < nlepton; i++){
      		if(i<nElectron){
        		if(Electron_cutBased[i]>1&& abs(Electron_eta[i]) <2.4 && (abs(Electron_eta[i])<1.4442||abs(Electron_eta[i])>1.5660)&&Electron_pt[i]>15){
          			if(lepton_flag==true&&i!=lepton_index[0]&&i!=lepton_index[1]){
          				if((abs(Electron_deltaEtaSC[i]+Electron_eta[i])<1.479&&abs(Electron_dxy[i])<0.05&&abs(Electron_dz[i])<0.1)
						||(abs(Electron_deltaEtaSC[i]+Electron_eta[i])>=1.479&&abs(Electron_dxy[i])<0.1&&abs(Electron_dz[i])<0.2))
          					num_select1++;
          			}
          		}
            }       
     		else{
        		if(Muon_looseId[i-nElectron]==1&&Muon_pfRelIso04_all[i-nElectron]<0.25&&Muon_pt[i-nElectron]>15&&abs(Muon_eta[i-nElectron])<2.4){    
          			if(lepton_flag==true&&i!=lepton_index[0]&&i!=lepton_index[1])
          				num_select1++;
          		}
            }
        }
        if(lepton_flag==true&&num_select1>0)
        	lepton_flag=false;
        //cout<<num_select1<<" "<<num_select2<<" "<<num_select3<<endl;
	    ////////////////////////////////////////////////////////////////////
	    bool jet_flag=false;
	    TLorentzVector mom_jets[45],mom_jet[2];
	    UInt_t index[45];
	    if(lepton_flag==true){
	    	int jet_num = 0; // count number fot bjets satisfy the selection criteria
		    for (int i=0; i<nJet; i++) {
		      	mom_jets[i].SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
		      	if(abs(Jet_eta[i]) < 2.4 && Jet_pt[i] > 30&&Jet_btagDeepB[i]>0.14&&Jet_jetId[i]==6&&mom_jets[i].DeltaR(p4_lepton[lepton_index[0]])>0.4&&mom_jets[i].DeltaR(p4_lepton[lepton_index[1]])>0.4) {
						jet_index[jet_num]=i;
						jet_num=jet_num+1;			
		    	}
		    }
		    if(jet_num>=2){
		    	jet_flag=true;
		    }
		    ///sort b-tagged jets with DeepCSV
		   	for(int i=0;i<jet_num;i++)
	      		index[i]=jet_index[i];
		   	for (int kk=0;kk<2;kk++) {
	      		int max=kk;
	      		for(int tt=kk+1;tt<jet_num;tt++){
	        		if(Jet_btagDeepB[index[tt]]>Jet_btagDeepB[index[max]]) {
	          			max=tt;
	        		}		
	      		}
	      		int tmp=index[max];
	      		index[max]=index[kk];
	      		index[kk]=tmp;
	    	}
	    }
	    /////////////////////////////////////////////////////////////////////////
      	double METx=MET_pt*cos(MET_phi);
      	double METy=MET_pt*sin(MET_phi);
      	UInt_t bpairl[2];
      	//////////////////////////////////////////////////////////////
      	if (jet_flag==true&&MET_pt>=30&&PV_npvsGood>=1){
      		for(int i=0;i<2;i++){
				jet_eta[i]=Jet_eta[index[i]];
				jet_pt[i]=Jet_pt[index[i]];
				jet_mass[i]=Jet_mass[index[i]];
				jet_phi[i]=Jet_phi[index[i]];
				jet_btagCSVV2[i]=Jet_btagCSVV2[index[i]];
				jet_btagDeepB[i]=Jet_btagDeepB[index[i]];
				mom_jet[i]=mom_jets[index[i]];
	    	}
	    	for(int i=0;i<2;i++){
	      		slepton_pt[i]=lepton_pt[lepton_index[i]];
	      		slepton_eta[i]=lepton_eta[lepton_index[i]];
	      		slepton_phi[i]=lepton_phi[lepton_index[i]];
	      		slepton_mass[i]=lepton_mass[lepton_index[i]];
	      		slepton_charge[i]=lepton_charge[lepton_index[i]];
	      		mom_lep[i]=p4_lepton[lepton_index[i]];
      		}
      		bool flag=true;
      		getbleppair(mom_jet,mom_lep,bpairl,METx,METy,flag);
      		if(flag){
	      		deltaybb=mom_jet[0].Rapidity()-mom_jet[1].Rapidity();	
				deltayll=mom_lep[0].Rapidity()-mom_lep[1].Rapidity();
				deltaybbll=(mom_jet[0]+mom_lep[bpairl[0]]).Rapidity()-(mom_jet[1]+mom_lep[bpairl[1]]).Rapidity();
				Mbbll=(mom_jet[0]+mom_jet[1]+mom_lep[0]+mom_lep[1]).M();
				if(inputFile.Contains("TTToSemiLeptonic")||inputFile.Contains("TTTo2L2Nu")||inputFile.Contains("TTToHadronic")){
					for(Int_t i=0;i<8;i++){
					   Int_t nbin=hist[i]->FindBin(M_tt_gen,delta_rapidity_gen);
					   weight[i]=1.0+hist[i]->GetBinContent(nbin);
					 //  cout<<"weight[i]: "<<weight[i]<<endl;
					}
	            }
	            mytree->Fill();
				nevents2++;
			}
      	}
    }
    file->cd();
    mytree->Write();
  	rawtree->Write();
  	cout << inputFile << " has " << total_entry << " events" << endl;
  	cout << output << " is created" << endl;
  	cout << nevents << " events are written into "<< "rawtree." << endl;
  	cout << nevents2 << " events are written into "<< "mytree." << endl;
}
