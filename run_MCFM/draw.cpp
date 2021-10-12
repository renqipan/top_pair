void draw(){
gStyle->SetOptStat(0);
const int nbinsx=13;
const int nbinsy=40;
const double ylow=-5.0;
const double yhigh=5.0;
double xbins[]={345,360,380,400,450,500,550,600,650,700,750,800,950,2000};
TString inputFileNames[nbinsx]={
"kappa02_345_360.root",
"kappa02_360_380.root",
"kappa02_380_400.root",
"kappa02_400_450.root",
"kappa02_450_500.root",
"kappa02_500_550.root",
"kappa02_550_600.root",
"kappa02_600_650.root",
"kappa02_650_700.root",
"kappa02_700_750.root",
"kappa02_750_800.root",
"kappa02_800_950.root",
"kappa02_950_2000.root"
};

TH2F *h2 = new TH2F("h2","", nbinsx, xbins, nbinsy, ylow, yhigh);
for (int i =0;i<nbinsx; i++){
	for (int j =0;j<nbinsy; j++){
		TFile *f = new TFile(inputFileNames[i]);
		TH1F *htemp = (TH1F*) f->Get("id6");
		h2 ->SetBinContent ( i+1, j+1, htemp->GetBinContent(j+1) );
	}
}
h2->GetXaxis()->SetTitle("M_ttbar");
h2->GetYaxis()->SetTitle("y_ttbar");
h2->GetYaxis()->CenterTitle();
h2->GetXaxis()->CenterTitle();
h2->Draw("colz");
TFile *fnew = new TFile("correction_tilde2.root","recreate");
fnew->cd();
h2->Write();
fnew->Close();
}