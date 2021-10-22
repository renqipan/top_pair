//#include "CMS_lumi.C"
#include "plotstyle.h"
TGraph* bestFit(TTree *t, TString x, TString y, TCut cut, double minv) {
	int nfind = t->Draw(y+":"+x, cut + Form("abs(deltaNLL -%f)<0.001",minv));
	//   cout<<"find "<<nfind<<endl;
	//int nfind = t->Draw(y+":"+x, cut + Form("deltaNLL == %f",0));
	if (nfind == 0) {
		TGraph *gr0 = new TGraph(1);
		gr0->SetPoint(0,-999,-999);
		gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
		return gr0;
	} else {
		TGraph *gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone();
		gr0->SetMarkerStyle(34); gr0->SetMarkerSize(2.0);
		cout<< gr0->GetX()[0]<<endl; 
		if (gr0->GetN() > 1) gr0->Set(1);
		return gr0;
	}
}

TH2 *treeToHist2D(TTree *t, TString x, TString y, TString name, TCut cut, double xmin, double xmax, double ymin, double ymax, int xbins, int ybins, double minv) {
//	t->Draw(Form("2*(deltaNLL-%f):%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)", minv,y.Data(), x.Data(), name.Data(), xbins, xmin, xmax, ybins, ymin, ymax), cut+Form( "deltaNLL-%f != 0",minv), "PROF");
		t->Draw(Form("2*(deltaNLL):%s:%s>>%s_prof(%d,%10g,%10g,%d,%10g,%10g)", y.Data(), x.Data(), name.Data(), xbins, xmin, xmax, ybins, ymin, ymax), cut + "deltaNLL !=0", "PROF");
	TH2 *prof = (TH2*) gROOT->FindObject(name+"_prof");
	TH2D *h2d = new TH2D(name, name, xbins, xmin, xmax, ybins, ymin, ymax);
	for (int ix = 1; ix <= xbins; ++ix) {
		for (int iy = 1; iy <= ybins; ++iy) {
			double z = prof->GetBinContent(ix,iy);
			if (z != z) z = (name.Contains("bayes") ? 0 : 999); // protect agains NANs
			if (z==0 &&  prof->GetBinContent(ix+1,iy)!=0)
				z = (prof->GetBinContent(ix+1,iy)+prof->GetBinContent(ix-1,iy))/2;
			h2d->SetBinContent(ix, iy, z);
		}
	}
	h2d->GetXaxis()->SetTitle(x);
	h2d->GetYaxis()->SetTitle(y);
	h2d->GetZaxis()->SetTitle("-2#Deltaln L");
	h2d->SetDirectory(0);
	return h2d;
}
TH2D* frameTH2D(TH2D *in, double threshold){
	// NEW LOGIC:
	//   - pretend that the center of the last bin is on the border if the frame
	//   - add one tiny frame with huge values
	double frameValue = 1000;
	if (TString(in->GetName()).Contains("bayes")) frameValue = -1000;

	Double_t xw = in->GetXaxis()->GetBinWidth(1);
	Double_t yw = in->GetYaxis()->GetBinWidth(1);

	Int_t nx = in->GetNbinsX();
	Int_t ny = in->GetNbinsY();

	Double_t x0 = in->GetXaxis()->GetXmin();
	Double_t x1 = in->GetXaxis()->GetXmax();

	Double_t y0 = in->GetYaxis()->GetXmin();
	Double_t y1 = in->GetYaxis()->GetXmax();
	Double_t xbins[999], ybins[999]; 
	double eps = 0.1;

	xbins[0] = x0 - eps*xw - xw; xbins[1] = x0 + eps*xw - xw;
	for (int ix = 2; ix <= nx; ++ix) xbins[ix] = x0 + (ix-1)*xw;
	xbins[nx+1] = x1 - eps*xw + 0.5*xw; xbins[nx+2] = x1 + eps*xw + xw;

	ybins[0] = y0 - eps*yw - yw; ybins[1] = y0 + eps*yw - yw;
	for (int iy = 2; iy <= ny; ++iy) ybins[iy] = y0 + (iy-1)*yw;
	ybins[ny+1] = y1 - eps*yw + yw; ybins[ny+2] = y1 + eps*yw + yw;

	TH2D *framed = new TH2D(
			Form("%s framed",in->GetName()),
			Form("%s framed",in->GetTitle()),
			nx + 2, xbins,
			ny + 2, ybins 
			);

	//Copy over the contents
	for(int ix = 1; ix <= nx ; ix++){
		for(int iy = 1; iy <= ny ; iy++){
			framed->SetBinContent(1+ix, 1+iy, in->GetBinContent(ix,iy));
		}
	}
	//Frame with huge values
	nx = framed->GetNbinsX();
	ny = framed->GetNbinsY();
	for(int ix = 1; ix <= nx ; ix++){
		framed->SetBinContent(ix,  1, frameValue);
		framed->SetBinContent(ix, ny, frameValue);
	}
	for(int iy = 2; iy <= ny-1 ; iy++){
		framed->SetBinContent( 1, iy, frameValue);
		framed->SetBinContent(nx, iy, frameValue);
	}

	return framed;
}

TList* contourFromTH2(TH2 *h2in, double threshold, int minPoints=20) {
	std::cout << "Getting contour at threshold " << threshold << " from " << h2in->GetName() << std::endl;
	//http://root.cern.ch/root/html/tutorials/hist/ContourList.C.html
	Double_t contours[1];
	contours[0] = threshold;
	if (h2in->GetNbinsX() * h2in->GetNbinsY() > 10000) minPoints = 50;
	if (h2in->GetNbinsX() * h2in->GetNbinsY() <= 100) minPoints = 10;

	TH2D *h2 = frameTH2D((TH2D*)h2in,threshold);

	h2->SetContour(1, contours);

	// Draw contours as filled regions, and Save points
	h2->Draw("CONT Z LIST");
	gPad->Update(); // Needed to force the plotting and retrieve the contours in TGraphs


	// Get Contours
	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	TList* contLevel = NULL;

	if (conts == NULL || conts->GetSize() == 0){
		printf("*** No Contours Were Extracted!\n");
		return 0;
	}

	TList *ret = new TList();
	for(int i = 0; i < conts->GetSize(); i++){
		contLevel = (TList*)conts->At(i);
		//printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
		for (int j = 0, n = contLevel->GetSize(); j < n; ++j) {
			TGraph *gr1 = (TGraph*) contLevel->At(j);
			//printf("\t Graph %d has %d points\n", j, gr1->GetN());
			if (gr1->GetN() > minPoints) ret->Add(gr1->Clone());
			//break;
		}
	}
	return ret;
}

void styleMultiGraph(TList *tmg, int lineColor, int lineWidth, int lineStyle) {
	for (int i = 0; i < tmg->GetSize(); ++i) {
		TGraph *g = (TGraph*) tmg->At(i);
		g->SetLineColor(lineColor); g->SetLineWidth(lineWidth); g->SetLineStyle(lineStyle);
	}
}
void styleMultiGraphMarker(TList *tmg, int markerColor, int markerSize, int markerStyle) {
	for (int i = 0; i < tmg->GetSize(); ++i) {
		TGraph *g = (TGraph*) tmg->At(i);
		g->SetMarkerColor(markerColor); g->SetMarkerSize(markerSize); g->SetMarkerStyle(markerStyle);
	}
}


/** Make a 2D contour plot from the output of MultiDimFit
 *
 * Inputs:
 *  - gFile should point to the TFile containing the 'limit' TTree
 *  - xvar should be the variable to use on the X axis, with xbins bins in the [xmin, xmax] range
 *  - yvar should be the variable to use on the Y axis, with ybins bins in the [ymin, ymax] range
 *  - (smx, smy) are the coordinates at which to put a diamond representing the SM expectation
 *  - if fOut is not null, then the output objects will be saved to fOut:
 *     - the 2D histogram will be saved as a TH2 with name name+"_h2d"
 *     - the 68% CL contour will be saved as a TList of TGraphs with name name+"_c68"
 *     - the 95% CL contour will be saved as a TList of TGraphs with name name+"_c95"
 *     - the 99.7% CL contour will be saved as a TList of TGraphs with name name+"_c997"
 *     - the best fit point will be saved as a TGraph with name name+"_best"
 *
 * Notes:
 *     - it's up to you to make sure that the binning you use for this plot matches the one used
 *       when running MultiDimFit (but you can just plot a subset of the points; e.g. if you had
 *       100x100 points in [-1,1]x[-1,1] you can make a 50x50 plot for [0,1]x[0,1])
 *     - the 99.7 contour is not plotted by default
 *     - the SM marker is not saved
 */
void contour2D(TString datafile,TString xvar, int xbins, float xmin, float xmax, TString yvar, int ybins, float ymin, float ymax, float smx=1.0, float smy=1.0, TString name="contour2D",TString xtitle="#mu_{ggH,bbH,ttH}", TString ytitle="#mu_{VBF,VH}",TFile *fOut=0 ) {
	//    TTree *tree = (TTree*) gFile->Get("limit") ;
	gStyle->SetOptStat(0);
	TChain *tree = new TChain("limit","") ;
	tree->Add(datafile);
	//    cout<< tree->GetEntries()<<endl;
	//	int ent =tree->Draw("deltaNLL:rv:rf","abs(deltaNLL)<20");
	int ent =tree->Draw("deltaNLL:"+yvar+":"+xvar,"");
	double *dnllV = tree->GetV1();
	double *rv = tree->GetV2();
	double *rf = tree->GetV3();
	//    for(int i =0;i<2000;i++)
	//	    cout<< *(dnllV+i)<<"\t"<< *(rv+i)<<"\t"<<*(rf+i)<<endl;
	double minv = *std::min_element(dnllV, dnllV+ent);

	cout<< minv<<endl;
	//tree->Scan("2*deltaNLL:rv:rf","abs(2*deltaNLL-1)<0.03");
	TH2 *hist2d = treeToHist2D(tree, xvar, yvar, "", "", xmin, xmax, ymin, ymax, xbins, ybins,minv);
	//hist2d->SetContour(200);
	//hist2d->GetZaxis()->SetRangeUser(-0,100);
	for (int i=0;i<hist2d->GetNbinsX();i++){
	for (int j=0;j<hist2d->GetNbinsY();j++){
		float cc=hist2d->GetBinContent(i+1,j+1);
		if (cc<20&&cc>0)
		cout<<hist2d->GetXaxis()->GetBinCenter(i+1)<<"\t"<< hist2d->GetYaxis()->GetBinCenter(j+1)<< "\t"<<cc<<endl;
	}}

	hist2d->Draw("colz");
	hist2d->GetXaxis()->CenterTitle();
	hist2d->GetYaxis()->CenterTitle();
	gPad->Print(name+"_2d.png");
	TGraph *fit = bestFit(tree, xvar, yvar, "",minv);
	TCanvas *c = new TCanvas ("c","",600,600);
	TList *c68 = contourFromTH2(hist2d, 2.30);
	TList *c95 = contourFromTH2(hist2d, 5.99);
	//	c68->Print("v");
	//	c95->Print("v");
	//	TList *c997 = contourFromTH2(hist2d, 11.83);
	styleMultiGraph(c68,  /*color=*/1, /*width=*/3, /*style=*/2);
	styleMultiGraph(c95,  /*color=*/1, /*width=*/3, /*style=*/1);
	//	styleMultiGraph(c997, /*color=*/1, /*width=*/3, /*style=*/2);
	
	//hist2d->GetXaxis()->SetTitle("#mu_{#lower[-0.2]{ggH,b#bar{b}H,t#bar{t}H,tH}}"); 
	//hist2d->GetYaxis()->SetTitle("#mu_{#lower[-0.1]{VBF,VH}}"); 
	hist2d->GetXaxis()->SetTitle(xtitle); 
	hist2d->GetYaxis()->SetTitle(ytitle); 
	hist2d->GetXaxis()->SetTitleOffset(0.4); 
	hist2d->GetYaxis()->SetTitleOffset(0.4);
	c->SetRightMargin(0.01);
	c->SetTopMargin(0.07);
	c->SetBottomMargin(0.13);
	c->SetLeftMargin(0.14);
	hist2d->Draw("col"); 
	gPad->Update();
	//TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
	//palette->SetX1NDC(palette->GetX1NDC()-0.045);
	//palette->SetX2NDC(palette->GetX2NDC()-0.05);
	//palette->SetX2NDC(0.60);
	//palette->SetX2NDC(0.90);
	//palette->SetY1NDC(0.13);
	//palette->SetY2NDC(0.93);
	//	cout<< c68->size()<<endl;
	//	cout<< c95->size()<<endl;
	//if(c68->At(0)){
	  c68->Draw("C SAME");
	  c95->Draw("C SAME");

		
	//}
	//if(c95->At(0)){
	//     c95->At(0)->Draw("C SAME");
	//}

	TMarker m;
	m.SetMarkerSize(3.0); m.SetMarkerColor(97); m.SetMarkerStyle(33); 
	m.DrawMarker(smx,smy);
	//fit->Draw("P SAME");

	//    m.SetMarkerSize(1.8); m.SetMarkerColor(89); m.SetMarkerStyle(33); 
	//    m.DrawMarker(smx,smy);
	//
	//CMS_lumi( c, 4, 0, true);
	TLegend *leg = new TLegend(0.76,0.70,0.96,0.88);

	leg->SetFillColor(0);
	leg->SetLineColor(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.05);
	gStyle->SetEndErrorSize(10);
	leg->AddEntry(&m, "SM","p");
	//leg->AddEntry(fit, "best fit","p");
	leg->AddEntry(c68->At(0), "68% CL","l");
	leg->AddEntry(c95->At(0), "95% CL","l");
	leg->Draw();
	TPaveText *ptc= new TPaveText(0.65,0.8,0.8,0.9,"brNDC");
	ptc->SetBorderSize(0);
	//        ptc->SetTextAlign(12);
	ptc->SetFillStyle(0);
	ptc->SetTextFont(42);
	ptc->SetTextSize(0.04);
	ptc->AddText("H#rightarrowZZ#rightarrow4l");
//	ptc->AddText("m_{H} profiled");
	//ptc->Draw();
	gPad->SetGrid(1,1);
	//applycanvasstyle((TCanvas*)gPad);
	hist2d->GetZaxis()->SetRangeUser(0,20); 

	applylegendstyle(leg);
	applyaxesstyle(hist2d);
	setColZGradient_TwoColors();//set color style
	gPad->Print(name+".png");
	gPad->Print(name+".pdf");
	gPad->Print(name+".eps");

	if (fOut != 0) {
		/*hist2d->SetName(name+"_h2d");*/  fOut->WriteTObject(hist2d,0);
		fit->SetName(name+"_best");    fOut->WriteTObject(fit,0);
		c68->SetName(name+"_c68");     fOut->WriteTObject(c68,0,"SingleKey");
		c95->SetName(name+"_c95");     fOut->WriteTObject(c95,0,"SingleKey");
		//c997->SetName(name+"_c997");   fOut->WriteTObject(c997,0,"SingleKey");
	}
}


