using namespace std;
void read_entry(TString fileName){
//	TString dir="/Users/renqi/Documents/top_pairs/root_files/";
	 dir="root://cms-xrd-global.cern.ch/";
	TFile* file=TFile::Open(dir+fileName);
	TTree* tree=(TTree*) file->Get("Events");
	float entries=tree->GetEntries();
	cout<<"has "<<entries<<" entries"<<endl;
	file->Close();
	delete file;
}
