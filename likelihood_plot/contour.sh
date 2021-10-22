 #contour2D(TString datafile,TString xvar, int xbins, float xmin, float xmax, TString yvar, int ybins, float ymin, 
 #float ymax, float smx=1.0, float smy=1.0, TString name="contour2D",TString xtitle="#mu_{ggH,bbH,ttH}", 
 #TString ytitle="#mu_{VBF,VH}",TFile *fOut=0 )

root -l -q -b contour2D.cxx"(\"higgsCombinethw_expect.MultiDimFit.mH125.root\",\"kappa\",100,-2.0,2.0,\"kappa_t\",50,-3.5,3.5,1.,0.,\"scan2D\",\"#kappa\",\"#tilde{#kappa}\")"
