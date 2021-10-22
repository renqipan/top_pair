void applycanvasstyle(TCanvas *c){
    c->SetFillColor(0);
    c->SetBorderMode(0);
    c->SetBorderSize(2);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.07);
    c->SetBottomMargin(0.13);
    c->SetFrameFillStyle(0);
    c->SetFrameBorderMode(0);
    c->SetFrameFillStyle(0);
    c->SetFrameBorderMode(0);
}
void applylegendstyle(TLegend *l2D){
    l2D->SetBorderSize(0);
    l2D->SetTextFont(42);
    l2D->SetTextSize(0.04);
    l2D->SetLineColor(1);
    l2D->SetLineStyle(1);
    l2D->SetLineWidth(1);
    l2D->SetFillColor(0);
    l2D->SetFillStyle(0);
}
void applyaxesstyle(TH1 *h){
    h->GetXaxis()->SetNdivisions(505);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelOffset(0.007);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleSize(0.08);
    h->GetXaxis()->SetTitleOffset(0.7);
    h->GetXaxis()->SetTitleFont(42);
    h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelOffset(0.007);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetYaxis()->SetTitleOffset(0.7);
    h->GetYaxis()->SetTitleFont(42);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetZaxis()->CenterTitle();

    h->GetZaxis()->SetTitleSize(0.07);
    h->GetZaxis()->SetTitleOffset(0.6);
    h->GetZaxis()->SetTitleFont(42);
}

void setColZGradient_TwoColors() {
    const Int_t NRGBs = 4;
    const Int_t NCont = 255;
//    Double_t stops[NRGBs] = { 0.,0.023,0.0599,1.0};
   Double_t stops[NRGBs] = { 0.0,0.230,0.599,1.0};
   // Double_t red[NRGBs]   = { 60./256., 1.00 };
   // Double_t green[NRGBs] = { 140./256., 1.00 };
   // Double_t blue[NRGBs]  = { 1.00, 218./256. };
    Double_t red[NRGBs]   = { 25./256.,246./256.,1,0};
    Double_t green[NRGBs] = { 121./256.,198./256.,1,153./256.};
    Double_t blue[NRGBs]  = { 218./256.,108./256.,1,150./256.};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
