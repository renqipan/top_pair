void draw(TString couplings){
gStyle->SetOptStat(0);
const int nbinsx=14;
const int nbinsy=40;
const double ylow=-5.0;
const double yhigh=5.0;
double xbins[]={345,360,380,400,450,500,550,600,650,700,750,800,900,1200,2000};
TString mtt[nbinsx]={
"345_360.root",
"360_380.root",
"380_400.root",
"400_450.root",
"450_500.root",
"500_550.root",
"550_600.root",
"600_650.root",
"650_700.root",
"700_750.root",
"750_800.root",
"800_900.root",
"900_1200.root",
"1200_2000.root"
};

TH2F *h2 = new TH2F("h2","", nbinsx, xbins, nbinsy, ylow, yhigh);
for (int i =0;i<nbinsx; i++){
	for (int j =0;j<nbinsy; j++){
		TFile *f = new TFile(couplings+"_"+mtt[i]);
		TH1F *htemp = (TH1F*) f->Get("id13");
		h2 ->SetBinContent ( i+1, j+1, htemp->GetBinContent(j+1) );
	}
}
h2->GetXaxis()->SetTitle("Mtt");
h2->GetYaxis()->SetTitle("ytt");
h2->GetYaxis()->CenterTitle();
h2->GetXaxis()->CenterTitle();
h2->Draw("colz");
TFile *fnew = new TFile("EW"+couplings+".root","recreate");
fnew->cd();
h2->Write();
fnew->Close();
cout<<"EW"+couplings+".root is created"<<endl;
}
