using namespace std;
void get_entries(TString fileName){
//	TString dir="/Users/renqi/Documents/top_pairs/root_files/";
	//TString dir="root://cms-xrd-global.cern.ch/";
	TFile* file=TFile::Open(fileName,"READ");
	TTree* tree=(TTree*) file->Get("mytree");
	int entries=tree->GetEntries();
	cout<<"has "<<setprecision(9)<<entries<<" entries"<<endl;
	file->Close();
	delete file;
}
