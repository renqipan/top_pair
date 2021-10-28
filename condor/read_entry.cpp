using namespace std;
void read_entry(TString fileName){
//	TString dir="/Users/renqi/Documents/top_pairs/root_files/";
	 dir="root://cms-xrd-global.cern.ch/";
	TFile* file=TFile::Open(dir+fileName,"READ");
	TTree* tree=(TTree*) file->Get("Events");
	int entries=tree->GetEntries();
	cout<<"has "<<setprecision(9)<<entries<<" entries"<<endl;
	file->Close();
	delete file;
}
