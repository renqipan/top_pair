void read_pdfhist(){
	bool do_pdfsys=true;
	TString cutsName[]={"3jets","4jets"};
	for(int s=0;s<2; s++){
		if(do_pdfsys){
			TString pdf_dir="/eos/user/y/yuekai/output3";
			TString pdf_file="pdfweight_"+cutsName[s]+".root";
			cout<<"pdf weight file: "<<pdf_file<<endl;
			TFile *fpdf = TFile::Open(pdf_dir+"/"+pdf_file);
			TFile *file = new TFile("pdfhist_"+cutsName[s]+".root","recreate");
		    TIter keyList(fpdf->GetListOfKeys());
		    TKey *key;
		    while ((key = (TKey*)keyList())) {
		        TClass *cl = gROOT->GetClass(key->GetClassName());
		       if (!cl->InheritsFrom("TH1D")) 
		      		cout<<"the format is not TH1D"<<endl;
		       else{
		      		TH1D *hpdf = (TH1D*)key->ReadObj();
		      		cout<<"hist name: "<<hpdf->GetName()<<endl;
				    file->cd();
				    hpdf->Write();
				    delete hpdf;
				}

		    }
		}
	}
	
}