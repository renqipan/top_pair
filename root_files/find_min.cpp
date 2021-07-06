Double_t myfunction(Double_t *x, Double_t *par)
		{
		   Float_t xx =x[0];
		   Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
		   return f;
		}
void find_min(){
	
    TF1 *fun= new TF1("fun",myfunction,0,10,2);
    fun->SetParameters(1,2);
    fun->Draw();
    double min=fun->GetMinimum(0.001,2);
    double minx=fun->GetMinimumX(0.001,2);
    cout<<"min: "<<min<<" minx: "<<minx<<endl;

}