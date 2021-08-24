// get useful information from NanoAOD files for ttbar analysis
// each event contains the information for top quark pairs and their decay
// products written by Renqi Pan in 16th June, 2021.
//////////////////////////////////////////////////////////////////////////////
// declare global vraibles
Float_t nu_px, nu_py, neutrino_pz;
TLorentzVector mom_nu, mom_lep, mom_jets[45];
Int_t btag_num; // btag_num: number of bjets among the leading four jets
Float_t LepCharge, Jet_btagDeepB[45],jet_btagDeepB[45]; // charge of the satisfied leading lepton
UInt_t jet_num;     // count the number of jets satisfy pt and eta criteria
int bjets_index[2]; // store the indexes of two bjets with the hightest btaged
                  // score.
int jets_index[45]; // store the indexes of light jets
int jet_index[45];
int bjet_lep, bjet_had, min_j1, min_j2; // denote the minimum likelihood case
Float_t mass_wlep, mass_whad, mass_tlep, mass_thad;
Double_t pro_wlep, pro_tlep, pro_thad, pro_whad, pro_twlep;
int bindex; // for loop over b-jet index

float xi_thad = 18.0, x0_thad = 179, xi_wlep = 2.0, x0_wlep = 80, xi_tlep = 8.5,
      x0_tlep = 169, xi_whad = 14.0, x0_whad = 84;
float mw_lep = 80, mt_had = 173, mt_lep = 173, mw_had = 80.0,
	  sigmat_lep = 40, sigmaw_had = 35.0, sigmat_had = 50, sigmaw_lep = 20.0;
TLorentzVector mom_top, mom_antitop;
float rectop_mass, recantitop_mass, rectop_pt, mass_tt, rapidity_tt,
      rectop_eta,rectop_rapidity,rectop_costheta;
Double_t nupz_min = -1000.0, nupz_max = 1000.0;
Double_t minimum;//minimum likelihood

////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//define likelihood function with two b-jets
Double_t likelihood(Double_t *pz,Double_t *pars){
    Double_t nu_pz=pz[0];
    int j1=pars[0];
    int j2=pars[1];
    Float_t nu_E=sqrt(nu_px*nu_px+nu_py*nu_py+nu_pz*nu_pz);
    mom_nu= TLorentzVector(nu_px,nu_py,nu_pz,nu_E);
    mass_wlep=(mom_nu+mom_lep).M();
    mass_whad=(mom_jets[jets_index[j1]]+mom_jets[jets_index[j2]]).M();
    mass_tlep=(mom_nu+mom_lep+mom_jets[bjets_index[bindex]]).M();
    if(bindex==0){
          mass_thad=(mom_jets[jets_index[j1]]+mom_jets[jets_index[j2]]+mom_jets[bjets_index[1]]).M();
          }
    else
          mass_thad=(mom_jets[jets_index[j1]]+mom_jets[jets_index[j2]]+mom_jets[bjets_index[0]]).M();
    
    pro_wlep=ROOT::Math::gaussian_pdf(mass_wlep,sigmaw_lep,mw_lep);
    pro_tlep=ROOT::Math::gaussian_pdf(mass_tlep,sigmat_lep,mt_lep);
    pro_whad=ROOT::Math::gaussian_pdf(mass_whad,sigmaw_had,mw_had);
    pro_thad=ROOT::Math::gaussian_pdf(mass_thad,sigmat_had,mt_had);
    Double_t log_nupz;
    log_nupz=-TMath::Log(pro_tlep)-TMath::Log(pro_wlep)-TMath::Log(pro_thad)-TMath::Log(pro_whad); 
     
    return log_nupz;
}
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

// reconstruct ttbar pair using likelihood function
// accorind to the leading four jets and two bjets among them.
void recons_tt() {
  if (btag_num >= 2) {
    //reorder bjet scores acrroding to their scores
    int index[jet_num];
    for (int i = 0; i < jet_num; i++)
      index[i] = jet_index[i];

    for (int kk = 0; kk < 2; kk++) {
      int max = kk;
      for (int tt = kk + 1; tt < jet_num; tt++) {
        if (Jet_btagDeepB[index[tt]] > Jet_btagDeepB[index[max]]) {
          max = tt;
        }
      }
      int tmp = index[max];
      index[max] = index[kk];
      index[kk] = tmp;
    }
    bjets_index[0] = index[0];
    bjets_index[1] = index[1];
    int light_jets = 0;
    for (int k = 0; k < jet_num; k++) {
      if (jet_index[k] != bjets_index[0] && jet_index[k] != bjets_index[1]) {
        jets_index[light_jets] = jet_index[k];
        light_jets++; // light flavour jets
      }
    }

    Double_t minimum_likelihood, nupz;
    bjet_lep = 0, bjet_had = 1, min_j1 = 0, min_j2 = 1;
    for (bindex = 0; bindex < 2; bindex++) {
      for (int j1 = 0; j1 < jet_num - 2; j1++) {
        for (int j2 = j1 + 1; j2 < jet_num - 2; j2++) {
          TF1 *likelihood_fun =
              new TF1("likelihood_fun", likelihood, nupz_min, nupz_max, 2);
          likelihood_fun->SetParameters(j1, j2); // pass the index of j1 and j2 as parameters to a
                                                 // function with type TF1
          minimum_likelihood = likelihood_fun->GetMinimum(nupz_min, nupz_max);
          nupz = likelihood_fun->GetMinimumX(nupz_min, nupz_max);

          if (bindex == 0 && j1 == 0 && j2 == 1) {
            minimum = minimum_likelihood;
            neutrino_pz = nupz;
          } else if (minimum_likelihood < minimum) {
            minimum = minimum_likelihood;
            neutrino_pz = nupz;
            bjet_lep = bindex;
            bjet_had = bindex == 0 ? 1 : 0;
            min_j1 = j1;
            min_j2 = j2;
          }
          
        }
      }
    }

    Float_t nu_E =
        sqrt(nu_px * nu_px + nu_py * nu_py + neutrino_pz * neutrino_pz);
    mom_nu = TLorentzVector(nu_px, nu_py, neutrino_pz, nu_E);
    mass_wlep = (mom_nu + mom_lep).M();
    mass_whad =
        (mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]]).M();
    mass_tlep = (mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]]).M();
    mass_thad = (mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]] +
                 mom_jets[bjets_index[bjet_had]])
                    .M();
    if (LepCharge > 0) {
      mom_top = mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]];
      mom_antitop = mom_jets[jets_index[min_j1]] +
                    mom_jets[jets_index[min_j2]] +
                    mom_jets[bjets_index[bjet_had]];
    } else {
      mom_top = mom_jets[jets_index[min_j1]] + mom_jets[jets_index[min_j2]] +
                mom_jets[bjets_index[bjet_had]];
      mom_antitop = mom_nu + mom_lep + mom_jets[bjets_index[bjet_lep]];
    }
    rectop_mass = mom_top.M();
    recantitop_mass = mom_antitop.M();
    rectop_pt = mom_top.Pt();
    rectop_rapidity=mom_top.Rapidity();
    rectop_eta=mom_top.Eta();
    rectop_costheta=mom_top.CosTheta();
    rapidity_tt = mom_top.Rapidity() - mom_antitop.Rapidity();
    mass_tt = (mom_antitop + mom_top).M();
  }
}
///////////////////////////////////////////////////////////////////////
// select the semileptonic final states and reconstruct top quark pairs.
void get_info() {
  TChain chain("Events");
  TString inputFile ="TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv6p1_2018.root";
  chain.Add(inputFile);
  chain.Add("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1TopNanoAODv5p1_2018.root");
  TString output = "new_new_"+inputFile;
  TFile *file = new TFile(output, "RECREATE");
  TTree *mytree = new TTree("mytree", " tree with branches");
  TTree *rawtree = new TTree("rawtree", "tree without selection");
  Int_t nevents = 0, nevents2 = 0; // count the number of events written in tree

  cout << inputFile << " is reading and processing" << endl;
  cout << "total number of events: " << chain.GetEntries() << endl;
  TH2F* hist[8];// hists for  weights of EW corrections 
   if(inputFile.Contains("TTToSemiLeptonic")){
     TString dir="/Users/renqi/Documents/top_pairs/correction_roots/";
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
  UInt_t tt_efficiency;
  UInt_t LHE_nlep = 0, LHE_nhad = 0, LHE_tao = 0;
  Float_t M_tt_gen, delta_rapidity_gen, lep_charge;
  Float_t top_pt, top_eta, top_mass, top_phi, antitop_pt, antitop_eta,
      antitop_phi, antitop_mass;
  Float_t b_pt, b_eta, b_mass, b_phi, antib_pt, antib_eta, antib_phi,
      antib_mass;
  Float_t lep_pt, lep_eta, lep_mass, lep_phi, nu_pt, nu_eta, nu_phi, nu_mass;
  Float_t up_pt, up_eta, up_mass, up_phi, down_pt, down_eta, down_phi, down_mass;
  // Int_t LHE_had[6];
  // difine branch for ttbar process information at parton level
  if (inputFile.Contains("ttbar_semi") ||
      inputFile.Contains("TTToSemiLeptonic")) {
    chain.SetBranchAddress("LHEPart_eta", LHEPart_eta);
    chain.SetBranchAddress("LHEPart_mass", LHEPart_mass);
    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
    chain.SetBranchAddress("LHEPart_pt", LHEPart_pt);
    chain.SetBranchAddress("LHEPart_phi", LHEPart_phi);
    chain.SetBranchAddress("LHEPart_pdgId", LHEPart_pdgId);
    chain.SetBranchAddress("nLHEPart", &nLHEPart);
    // chain.SetBranchAddress("LHEPart_status",LHEPart_status);
    mytree->Branch("M_tt_gen", &M_tt_gen, "M_tt_gen/F");
    mytree->Branch("delta_rapidity_gen", &delta_rapidity_gen,
                   "delta_rapidity_gen/F");
    mytree->Branch("lep_charge", &lep_charge, "lep_charge/F");
    mytree->Branch("top_pt", &top_pt, "top_pt/F");
    mytree->Branch("top_eta", &top_eta, "top_eta/F");
    mytree->Branch("top_phi", &top_phi, "top_phi/F");
    mytree->Branch("top_mass", &top_mass, "top_mass/F");
    mytree->Branch("antitop_pt", &antitop_pt, "antitop_pt/F");
    mytree->Branch("antitop_eta", &antitop_eta, "antitop_eta/F");
    mytree->Branch("antitop_phi", &antitop_phi, "antitop_phi/F");
    mytree->Branch("antitop_mass", &antitop_mass, "antitop_mass/F");
    mytree->Branch("up_pt", &up_pt, "up_pt/F");
    mytree->Branch("up_eta", &up_eta, "up_eta/F");
    mytree->Branch("up_phi", &up_phi, "up_phi/F");
    mytree->Branch("up_mass", &up_mass, "up_mass/F");
    mytree->Branch("down_pt", &down_pt, "down_pt/F");
    mytree->Branch("down_eta", &down_eta, "down_eta/F");
    mytree->Branch("down_phi", &down_phi, "down_phi/F");
    mytree->Branch("down_mass", &down_mass, "down_mass/F");
    mytree->Branch("b_pt", &b_pt, "b_pt/F");
    mytree->Branch("b_eta", &b_eta, "b_eta/F");
    mytree->Branch("b_phi", &b_phi, "b_phi/F");
    mytree->Branch("b_mass", &b_mass, "b_mass/F");
    mytree->Branch("antib_pt", &antib_pt, "antib_pt/F");
    mytree->Branch("antib_eta", &antib_eta, "antib_eta/F");
    mytree->Branch("antib_phi", &antib_phi, "antib_phi/F");
    mytree->Branch("antib_mass", &antib_mass, "antib_mass/F");
    mytree->Branch("lep_pt", &lep_pt, "lep_pt/F");
    mytree->Branch("lep_eta", &lep_eta, "lep_eta/F");
    mytree->Branch("lep_phi", &lep_phi, "lep_phi/F");
    mytree->Branch("lep_mass", &lep_mass, "lep_mass/F");
    mytree->Branch("nu_pt", &nu_pt, "nu_pt/F");
    mytree->Branch("nu_eta", &nu_eta, "nu_eta/F");
    mytree->Branch("nu_phi", &nu_phi, "nu_phi/F");
    mytree->Branch("nu_mass", &nu_mass, "nu_mass/F");
    mytree->Branch("tt_efficiency", &tt_efficiency, "tt_efficiency/I");
    
  }

  Float_t weight[8];
  if(inputFile.Contains("TTToSemiLeptonic")){
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
  Float_t Muon_mass[9], Muon_phi[9], Muon_pt[9], Muon_eta[9];
  Float_t lepton_mass[18], lepton_phi[18], lepton_eta[18], lepton_pt[18];
  UInt_t nMuon, nElectron, nJet, nlepton, Jet_btaged[45], nBtag;
  Int_t Jet_partonFlavour[45], Muon_charge[9], Electron_charge[9],
      lepton_charge[18];
  Float_t Jet_btagCSVV2[45], Jet_eta[45], Jet_mass[45], Jet_phi[45], Jet_pt[45];
  Float_t jet_btagCSVV2[45], jet_eta[45], jet_mass[45], jet_phi[45], jet_pt[45],jet_partonFlavour[45],jet_btaged[45];
  Float_t MtW;
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
  mytree->Branch("nlepton", &nlepton, "nlepton/I");
  mytree->Branch("lepton_eta", lepton_eta, "lepton_eta[nlepton]/F");
  mytree->Branch("lepton_pt", lepton_pt, "lepton_pt[nlepton]/F");
  mytree->Branch("lepton_mass", lepton_mass, "lepton_mass[nlepton]/F");
  mytree->Branch("lepton_phi", lepton_phi, "lepton_phi[nlepton]/F");
  mytree->Branch("lepton_charge", lepton_charge, "lepton_charge[nlepton]/I");
  mytree->Branch("jet_num", &jet_num,"jet_num/I"); // number of jets satisfy the  seletion criteria
  mytree->Branch("nBtag", &nBtag, "nBtag/I");
  mytree->Branch("jet_btaged", jet_btaged, "jet_btaged[jet_num]/I");
  mytree->Branch("jet_btagCSVV2", jet_btagCSVV2, "jet_btagCSVV2[jet_num]/F");
  mytree->Branch("jet_btagDeepB", jet_btagDeepB, "jet_btagDeepB[jet_num]/F");
  mytree->Branch("jet_partonFlavour", jet_partonFlavour,
                 "jet_partonFlavour[jet_num]/I");
  mytree->Branch("jet_eta", jet_eta, "jet_eta[jet_num]/F");
  mytree->Branch("jet_mass", jet_mass, "jet_mass[jet_num]/F");
  mytree->Branch("jet_phi", jet_phi, "jet_phi[jet_num]/F");
  mytree->Branch("Jet_pt", jet_pt, "jet_pt[jet_num]/F");
  rawtree->Branch("nJet", &nJet, "nJet/I");
  rawtree->Branch("nlepton", &nlepton, "nlepton/I");
  rawtree->Branch("Jet_pt", Jet_pt, "Jet_pt[nJet]/F");
  ////////////////////////////////////////////////////////////////
  // add information at reconstruction level.
  mytree->Branch("rectop_pt", &rectop_pt, "rectop_pt/F");
  mytree->Branch("rectop_eta",&rectop_eta,"rectop_eta/F");
  mytree->Branch("rectop_rapidity",&rectop_rapidity,"rectop_rapidity/F");
  mytree->Branch("rectop_costheta",&rectop_costheta,"rectop_costheta/F");
  mytree->Branch("rectop_mass", &rectop_mass, "rectop_mass/F");
  mytree->Branch("recantitop_mass", &recantitop_mass, "recantitop_mass/F");
  mytree->Branch("rapidity_tt", &rapidity_tt, "rapidity_tt/F");
  mytree->Branch("mass_tt", &mass_tt, "mass_tt/F");
  mytree->Branch("mass_wlep", &mass_wlep, "mass_wlep/F");
  mytree->Branch("mass_whad", &mass_whad, "mass_whad/F");
  mytree->Branch("neutrino_pz", &neutrino_pz, "neutrino_pz/F");
  mytree->Branch("mass_thad", &mass_thad, "mass_thad/F");
  mytree->Branch("mass_tlep", &mass_tlep, "mass_tlep/F");
  mytree->Branch("MtW",&MtW,"MtW/F");
  mytree->Branch("likelihood",&minimum,"minimum/D" );
  //////////////////////////////////////////////////////////////////////
  // loop over entry
  cout << "infomation is writing. Please wait for a while" << endl;
  cout << "infomation is writing. Please wait for a while" << endl;
  Int_t njet_need =4; // the at least number of jet of semileptonic final satate

  for (Int_t entry = 0; entry < chain.GetEntries(); entry++) {
    chain.GetEntry(entry);
    int index_b, index_antib, index_up, index_down, index_lep, index_nu;
    TLorentzVector p4_b, p4_antib, p4_up, p4_down, p4_lep, p4_nu, p4_top,
        p4_antitop;
    if (inputFile.Contains("ttbar_semi") ||
        inputFile.Contains("TTToSemiLeptonic")) {
      // get information of ttbar process at parton level from LHEPart
      //  int index_b,index_antib,index_up,index_down,index_lep,index_nu;
      for (int i = 0; i < nLHEPart; i++) {
        //    cout<<Form("LHEPart_pdgId[%d]:
        //",i)<<LHEPart_pdgId[i]<<endl;
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

      //  TLorentzVector p4_b,
      //  p4_antib,p4_up,p4_down,p4_lep,p4_nu,p4_top,p4_antitop;
      p4_b.SetPtEtaPhiM(LHEPart_pt[index_b], LHEPart_eta[index_b],
                        LHEPart_phi[index_b], LHEPart_mass[index_b]);
      p4_antib.SetPtEtaPhiM(LHEPart_pt[index_antib], LHEPart_eta[index_antib],
                            LHEPart_phi[index_antib],
                            LHEPart_mass[index_antib]);
      p4_up.SetPtEtaPhiM(LHEPart_pt[index_up], LHEPart_eta[index_up],
                         LHEPart_phi[index_up], LHEPart_mass[index_up]);
      p4_down.SetPtEtaPhiM(LHEPart_pt[index_down], LHEPart_eta[index_down],
                           LHEPart_phi[index_down], LHEPart_mass[index_down]);
      p4_lep.SetPtEtaPhiM(LHEPart_pt[index_lep], LHEPart_eta[index_lep],
                          LHEPart_phi[index_lep], LHEPart_mass[index_lep]);
      p4_nu.SetPtEtaPhiM(LHEPart_pt[index_nu], LHEPart_eta[index_nu],
                         LHEPart_phi[index_nu], LHEPart_mass[index_nu]);
      if (LHEPart_pdgId[index_lep] > 0) {
        p4_antitop = p4_antib + p4_lep + p4_nu;
        p4_top = p4_b + p4_up + p4_down;
        lep_charge = -1.0;

      } else {
        p4_top = p4_b + p4_lep + p4_nu;
        p4_antitop = p4_antib + p4_up + p4_down;
        lep_charge = 1.0;
      }
      top_pt = p4_top.Pt();
      top_mass = p4_top.M();
      top_phi = p4_top.Phi();
      top_eta = p4_top.Eta();
      antitop_pt = p4_antitop.Pt();
      antitop_mass = p4_antitop.M();
      antitop_phi = p4_antitop.Phi();
      antitop_eta = p4_antitop.Eta();
      up_pt = p4_up.Pt();
      up_mass = p4_up.M();
      up_phi = p4_up.Phi();
      up_eta = p4_up.Eta();
      down_pt = p4_down.Pt();
      down_mass = p4_down.M();
      down_phi = p4_down.Phi();
      down_eta = p4_down.Eta();
      lep_pt = p4_lep.Pt();
      lep_mass = p4_lep.M();
      lep_phi = p4_lep.Phi();
      lep_eta = p4_lep.Eta();
      nu_pt = p4_nu.Pt();
      nu_mass = p4_nu.M();
      nu_phi = p4_nu.Phi();
      nu_eta = p4_nu.Eta();
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
    }
    /////////////////////////////////////////////////////
    // get information for final state at detector level
    nlepton = nMuon + nElectron;
    TLorentzVector p4_lepton[18];
    for (int i = 0; i < nlepton; i++) {
      if (i < nElectron) {
        p4_lepton[i].SetPtEtaPhiM(Electron_pt[i], Electron_eta[i],
                                  Electron_phi[i], Electron_mass[i]);
        lepton_charge[i] = Electron_charge[i];
      } else {
        p4_lepton[i].SetPtEtaPhiM(
            Muon_pt[i - nElectron], Muon_eta[i - nElectron],
            Muon_phi[i - nElectron], Muon_mass[i - nElectron]);
        lepton_charge[i] = Muon_charge[i - nElectron];
      }
    }
    // sort leptons according to their pt and make the first jet has maximum pt
    for (int k = 0; k < nlepton; k++) {
      for (int j = k + 1; j < nlepton; j++) {
        if (p4_lepton[k].Pt() < p4_lepton[j].Pt()) {
          TLorentzVector temp;
          Int_t temp_charge;
          temp = p4_lepton[k];
          p4_lepton[k] = p4_lepton[j];
          p4_lepton[j] = temp;
          temp_charge = lepton_charge[k];
          lepton_charge[k] = lepton_charge[j];
          lepton_charge[j] = temp_charge;
        }
      }
    }
    for (int i = 0; i < nlepton; i++) {
      lepton_pt[i] = p4_lepton[i].Pt();
      lepton_eta[i] = p4_lepton[i].Eta();
      lepton_phi[i] = p4_lepton[i].Phi();
      lepton_mass[i] = p4_lepton[i].M();
    }
    nBtag = 0;   // count number of bjet among all the jets
    jet_num = 0; // count number fot jets satisfy the selection criteria
    for (int i = 0; i < nJet; i++) {
      mom_jets[i].SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i],
                                       Jet_mass[i]);
      if (abs(Jet_eta[i]) < 2.4 && Jet_pt[i] > 30) {
   //     mom_jets[jet_num].SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i],
     //                                  Jet_mass[i]);
        // jets_btag[i]=Jet_btaged[i]; //for reconstruct
        jet_index[jet_num] = i;
        jet_num = jet_num + 1;
        if (Jet_btagDeepB[i] > 0.14) {
          Jet_btaged[i] = 1;
          nBtag++;
        } else
          Jet_btaged[i] = 0;
      }
    }
    for(int i=0;i<jet_num;i++){
      jet_btaged[i]=Jet_btaged[jet_index[i]];
      jet_eta[i]=Jet_eta[jet_index[i]];
      jet_pt[i]=Jet_pt[jet_index[i]];
      jet_mass[i]=Jet_mass[jet_index[i]];
      jet_phi[i]=Jet_phi[jet_index[i]];
      jet_btagCSVV2[i]=Jet_btagCSVV2[jet_index[i]];
      jet_btagDeepB[i]=Jet_btagDeepB[jet_index[i]];
      jet_partonFlavour[i]=Jet_partonFlavour[jet_index[i]];

    }
    nevents++;
    rawtree->Fill();
    //////////////////////////////////////////////////////////////////
    // select ttbar semiletopnic final state
    // selec jets
    bool jet_flag = false; // if true pass the selection
    int jet_satisfy = 0;
    btag_num = 0;            // count bjet number in the leading four jets.
    if (nJet >= njet_need) { // njet_need=4 by default
      for (int i = 0; i < nJet; i++) {
        if (abs(Jet_eta[i]) < 2.4 && Jet_pt[i] > 30) {
          jet_satisfy++;
          btag_num = btag_num + Jet_btaged[i];

        } 
        else
            break;
      }
      if (jet_satisfy >= njet_need && btag_num >= 2)
        jet_flag = true;
    }
    // select lepton
    bool lepton_flag = false; // if true pass the selction
    int nlepton_satisfy = 1;
    if (nlepton >= 1 && lepton_pt[0] > 30 && abs(lepton_eta[0]) < 2.4) {
      for (int i = 1; i < nlepton; i++) { // start from the second lepton
        if (lepton_pt[i] > 15 && abs(lepton_eta[i]) < 2.4)
          break;
        else
          nlepton_satisfy++;
      }
      if (nlepton_satisfy == nlepton) {
        lepton_flag = true;
        mom_lep = p4_lepton[0];       // the lepton momenton for reconstrut
        LepCharge = lepton_charge[0]; // for reconstruct
      }
    }

    // select satisfied events
    if (jet_flag == true && lepton_flag == true ) {

      nu_px = MET_pt * cos(MET_phi);
      nu_py = MET_pt * sin(MET_phi);
      MtW=sqrt(2*(mom_lep.Pt()*MET_pt-mom_lep.Px()*nu_px-mom_lep.Py()*nu_py));
      
        recons_tt();

        if( minimum < 20.0  ){ 
        /////////////////////////////////////////
        //add weights according to invariant mass and rapidity difference at generator level.
          if(inputFile.Contains("TTToSemiLeptonic")){
                  for(Int_t i=0;i<8;i++){
                       Int_t nbin=hist[i]->FindBin(M_tt_gen,delta_rapidity_gen);
                       weight[i]=1.0+hist[i]->GetBinContent(nbin);
                     //  cout<<"weight[i]: "<<weight[i]<<endl;
                     }
               }
        //////////////////////////////////////////////  
        // calculate tt efficiency
        if (inputFile.Contains("ttbar_semi") ||inputFile.Contains("TTToSemiLeptonic")) {
	        int b_flag = 0;
	        int l_flag = 0;
	         int flag = 0;
	        for (int iso = 0; iso < nJet; iso++) {
	          for (int jso = 0; jso < nJet; jso++) {
	            for (int bso = 0; bso < nJet; bso++) {
	              for (int aso = 0; aso < nJet; aso++) {
	                if (p4_b.DeltaR(mom_jets[bso]) < 0.4 &&
	                    p4_antib.DeltaR(mom_jets[aso]) < 0.4 &&
	                    p4_up.DeltaR(mom_jets[iso]) < 0.4 &&
	                    p4_down.DeltaR(mom_jets[jso]) < 0.4) {
	                  if (bso != aso && iso != jso && iso != bso && iso != aso &&
	                      bso != jso && aso != jso) {
	                      flag = 1;
                        break;
	                  }
	                }
	              }
                if(flag==1) break;
	            }
              if(flag==1) break;
	          }
            if(flag==1) break;
	        }
        
          LHE_nhad = 0;
          LHE_nlep = 0;
          LHE_tao = 0;
          for (int i = 1; i < nLHEPart; i++) {
            if (LHEPart_pdgId[i] == 2 || LHEPart_pdgId[i] == 4 ||
                LHEPart_pdgId[i] == -2 || LHEPart_pdgId[i] == -4) {
              LHE_nhad++;
            }
            if (LHEPart_pdgId[i] == 1 || LHEPart_pdgId[i] == 3 ||
                LHEPart_pdgId[i] == -1 || LHEPart_pdgId[i] == -3) {
              LHE_nhad++;
            }
            if (LHEPart_pdgId[i] == 11 || LHEPart_pdgId[i] == 13 ||
                LHEPart_pdgId[i] == -11 || LHEPart_pdgId[i] == -13)
              LHE_nlep++;
            if (LHEPart_pdgId[i] == 15 || LHEPart_pdgId[i] == -15)
              LHE_tao++;
          }
          if (!(LHE_nhad == 2 && LHE_tao ==0 &&LHE_nlep==1)) {
                      tt_efficiency = 0;//tt background
                } 
                else if (flag == 0) {
                    tt_efficiency = 1;//non reco
                }
                else if (!((p4_b.DeltaR(mom_jets[bjets_index[0]]) < 0.4 &&
                            p4_antib.DeltaR(mom_jets[bjets_index[1]]) < 0.4) ||
                            (p4_b.DeltaR(mom_jets[bjets_index[1]]) < 0.4 &&
                            p4_antib.DeltaR(mom_jets[bjets_index[0]]) < 0.4))) {
                    tt_efficiency = 2; //wrong reco
                }
                else if (!((p4_up.DeltaR(mom_jets[jets_index[min_j1]]) < 0.4 &&
                            p4_down.DeltaR(mom_jets[jets_index[min_j2]]) < 0.4) ||
                            (p4_down.DeltaR(mom_jets[jets_index[min_j1]]) < 0.4 &&
                            p4_up.DeltaR(mom_jets[jets_index[min_j2]]) < 0.4))) {
                    tt_efficiency = 2; //wrong reco
                }
                else if (!((p4_b.DeltaR(mom_jets[bjets_index[bjet_lep]]) < 0.4 &&
                            p4_antib.DeltaR(mom_jets[bjets_index[bjet_had]]) < 0.4 && LepCharge>0) ||
                            (p4_b.DeltaR(mom_jets[bjets_index[bjet_had]]) < 0.4 &&
                            p4_antib.DeltaR(mom_jets[bjets_index[bjet_lep]]) < 0.4 && LepCharge<0))) {
                    tt_efficiency = 2; //wrong reco
                }
              
              else {
                    tt_efficiency = 3; //right
                 }
        }
        mytree->Fill();
        nevents2++;
      } // end of cut of minimum
    } //end of satisfied lepton and jets criteria  

  } // end loop over entries
  file->cd();
  mytree->Write();
  rawtree->Write();
  cout << inputFile << " has " << chain.GetEntries() << " events" << endl;
  cout << output << " is created" << endl;
  cout << nevents << " events are written into "
       << "rawtree." << endl;
  cout << nevents2 << " events are written into "
       << "mytree." << endl;
}
