#include <TMatrixD.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TMath.h>
using namespace std;
using namespace TMath;
TMatrixD RotationX(double a)
{
	double ca = cos(a);
	double sa = sin(a);
	TMatrixD D(3, 3);
	D(0, 0) = 1.; 
	D(1, 0) = 0.; 
	D(2, 0) = 0.; 
	D(0, 1) = 0.; 
	D(1, 1) = ca; 
	D(2, 1) = sa; 
	D(0, 2) = 0.; 
	D(1, 2) = -sa; 
	D(2, 2) = ca;
	return(D);
}

TMatrixD RotationY(double a)
{
	double ca = cos(a);
	double sa = sin(a);
	TMatrixD D(3, 3);
	D(0, 0) = ca; 
	D(1, 0) = 0.; 
	D(2, 0) = -sa; 
	D(0, 1) = 0.; 
	D(1, 1) = 1.; 
	D(2, 1) = 0.; 
	D(0, 2) = sa; 
	D(1, 2) = 0.; 
	D(2, 2) = ca;
	return(D);
}

TMatrixD RotationZ(double a)
{
	double ca = cos(a);
	double sa = sin(a);
	TMatrixD D(3, 3);
	D(0, 0) = ca; 
	D(1, 0) = sa; 
	D(2, 0) = 0.; 
	D(0, 1) = -sa; 
	D(1, 1) = ca; 
	D(2, 1) = 0.; 
	D(0, 2) = 0; 
	D(1, 2) = 0.; 
	D(2, 2) = 1.;
	return(D);
}
bool NeutrinoSolver(TLorentzVector* lep, TLorentzVector* bjet, double MW, double MT,TMatrixD &N) 
{	
	double Mt = MT;
	double Mw = MW;
	//Ml = 0.105;
	//Mb = 4.2;
	//Mn = 0.;
	double Ml = lep->M();
	double Mb = bjet->M();
	double Mn = 0.;
	TVector3 l(lep->Vect());
	TVector3 b(bjet->Vect());

	double pl = l.Mag();
	double pb = b.Mag();
	double El = sqrt(Ml*Ml + pl*pl);
	double Eb = sqrt(Mb*Mb + pb*pb);

	double cosbl = (l*b)/(pl*pb);;
	double sinbl = sqrt(1.-cosbl*cosbl);

	double betal = pl/El;
	double betab = pb/Eb;
	double gammali = Ml/El;
	//double gammab = Eb/Mb;

	double x0 = -0.5/El*(Mw*Mw - Ml*Ml - Mn*Mn);
	double x0p = -0.5/Eb*(Mt*Mt - Mw*Mw - Mb*Mb);
	double epsilon = (Mw*Mw - Mn*Mn)*gammali*gammali;
	double Sx = (x0*betal - pl*gammali*gammali)/(betal*betal);
	double Sy = 1./sinbl*(x0p/betab - cosbl*Sx);


	double omega = 1./sinbl*(betal/betab - cosbl);
	//double omega = -1./sinbl*(betal/betab + cosbl);

	double OmegaS = omega*omega - gammali*gammali;
	double Omega = sqrt(OmegaS);

	double x1 = Sx - (Sx + omega*Sy)/OmegaS;
	double y1 = Sy - (Sx + omega*Sy)*omega/OmegaS;
	double ZS = x1*x1*OmegaS - (Sy - omega*Sx)*(Sy - omega*Sx) - Mw*Mw + x0*x0 + epsilon*epsilon;

	if(ZS < 0) return false;

	double Z = sqrt(ZS);

	TMatrixD Ht(3,3), H(3,3),h(3,3);
	Ht(0, 0) = Z/Omega; 
	Ht(1, 0) = Z*omega/Omega; 
	Ht(2, 0) = 0.; 
	Ht(0, 1) = 0.; 
	Ht(1, 1) = 0.; 
	Ht(2, 1) = Z; 
	Ht(0, 2) = x1-pl; 
	Ht(1, 2) = y1; 
	Ht(2, 2) = 0.;

	TVector3 bn(b);
	TVector3 ln(l);
	double w1 = TMath::ATan2(ln.Y(), ln.X());

	//cout << "A" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;

	bn.RotateZ(-1.*w1);
	ln.RotateZ(-1.*w1);

	//cout << "B" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;

	double w2 = TMath::ATan2(ln.Z(), ln.X());

	bn.RotateY(w2);
	ln.RotateY(w2);
	//cout << "C" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;
	double w3 = TMath::ATan2(bn.Z(), bn.Y());

	bn.RotateX(-1.*w3);
	ln.RotateX(-1.*w3);

	//cout << "D" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;

	TMatrixD R1(RotationZ(w1));
	TMatrixD R2(RotationY(-1.*w2));
	TMatrixD R3(RotationX(w3));

	H=R1*R2*R3*Ht;
	h(0,0)=H(0,0);
	h(0,1)=H(0,1);
	h(0,2)=H(0,2);
	h(1,0)=H(1,0);
	h(1,1)=H(1,1);
	h(1,2)=H(1,2);
	h(2,0)=0;
	h(2,1)=0;
	h(2,2)=1;
	double adet=h(0,0)*h(1,1)*h(2,2)+h(0,1)*h(1,2)*h(2,0)+h(1,0)*h(2,1)*h(0,2)-h(0,2)*h(2,0)*h(1,1)-h(2,1)*h(1,2)*h(0,0)-h(1,0)*h(0,1)*h(2,2);
	if(!(abs(adet)>1e-6)) return false;
	TMatrixD U(3,3);
	U(0,0)=1;
	U(0,1)=0;
	U(0,2)=0;
	U(1,0)=0;
	U(1,1)=1;
	U(1,2)=0;
	U(2,0)=0;
	U(2,1)=0;
	U(2,2)=-1;
	TMatrixD tmp(h);
	tmp.Invert();
	tmp.T();
	h.Invert();
	N=tmp*U*h;
	return true;

	
}
