#include "NeutrinoSolver.h"


TMatrixD NeutrinoSolver::RotationX(double a)
{
	double ca = Cos(a);
	double sa = Sin(a);
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

TMatrixD NeutrinoSolver::RotationY(double a)
{
	double ca = Cos(a);
	double sa = Sin(a);
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

TMatrixD NeutrinoSolver::RotationZ(double a)
{
	double ca = Cos(a);
	double sa = Sin(a);
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

NeutrinoSolver::NeutrinoSolver(TLorentzVector* lep, TLorentzVector* bjet, double MW, double MT) : ERROR(false), H(3, 3), T(3,1), MET(2, 1), VM(2, 2)
{
	Mt = MT;
	Mw = MW;
	//Ml = 0.105;
	//Mb = 4.2;
	//Mn = 0.;
	Ml = lep->M();
	Mb = bjet->M();
	Mn = 0.;
	TVector3 l(lep->Vect());
	TVector3 b(bjet->Vect());

	double pl = l.Mag();
	double pb = b.Mag();
	double El = Sqrt(Ml*Ml + pl*pl);
	double Eb = Sqrt(Mb*Mb + pb*pb);

	double cosbl = (l*b)/(pl*pb);
	double sinbl = Sqrt(1.-cosbl*cosbl);

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
	double Omega = Sqrt(OmegaS);

	double x1 = Sx - (Sx + omega*Sy)/OmegaS;
	double y1 = Sy - (Sx + omega*Sy)*omega/OmegaS;
	double ZS = x1*x1*OmegaS - (Sy - omega*Sx)*(Sy - omega*Sx) - Mw*Mw + x0*x0 + epsilon*epsilon;


	if(ZS < 0)
	{
		ERROR = true;
		//omega = -1./sinbl*(betal/betab + cosbl);

		//OmegaS = omega*omega - gammali*gammali;
		//Omega = Sqrt(OmegaS);
		//x1 = Sx - (Sx + omega*Sy)/OmegaS;
		//y1 = Sy - (Sx + omega*Sy)*omega/OmegaS;
		//ZS = x1*x1*OmegaS - (Sy - omega*Sx)*(Sy - omega*Sx) - Mw*Mw + x0*x0 + epsilon*epsilon;
	}

	double Z = Sqrt(ZS);

	TMatrixD Ht(3, 3);
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
	double w1 = ATan2(ln.Y(), ln.X());

	//cout << "A" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;

	bn.RotateZ(-1.*w1);
	ln.RotateZ(-1.*w1);

	//cout << "B" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;

	double w2 = ATan2(ln.Z(), ln.X());

	bn.RotateY(w2);
	ln.RotateY(w2);
	//cout << "C" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;
	double w3 = ATan2(bn.Z(), bn.Y());

	bn.RotateX(-1.*w3);
	ln.RotateX(-1.*w3);

	//cout << "D" << endl;
	//cout << bn.X() << " " << bn.Y() << " " << bn.Z() << endl;
	//cout << ln.X() << " " << ln.Y() << " " << ln.Z() << endl;

	TMatrixD R1(RotationZ(w1));
	TMatrixD R2(RotationY(-1.*w2));
	TMatrixD R3(RotationX(w3));

	H=R1*R2*R3*Ht;
}

void NeutrinoSolver::Solve(double t)
{
	T(0, 0) = Cos(t);
	T(1, 0) = Sin(t);
	T(2, 0) = 1.;

	T = H*T;
}

TLorentzVector NeutrinoSolver::GetSolution(double t)
{
	Solve(t);
	return TLorentzVector(T(0,0), T(1,0), T(2,0), Sqrt(T(0,0)*T(0,0) + T(1,0)*T(1,0) + T(2,0)*T(2,0) + Mn*Mn));
}

TMatrixD NeutrinoSolver::GetPtSolution(double t)
{
	Solve(t);
	TMatrixD res(2,1);
	res(0,0) = T(0,0);
	res(1,0) = T(1,0);

	return(res);
}

double NeutrinoSolver::Chi2(double t)
{
	TMatrixD SOL(GetPtSolution(t));
	return(((MET-SOL).T()*VM*(MET-SOL))(0,0));
}

pair<double, double> NeutrinoSolver::Extrem(double t, bool MIN)
{
	double sign = -1.;
	if(MIN){sign = 1.;}
	double step = 0.05;
	double old = sign*Chi2(t);
	bool right = true;
	while(Abs(step) > 0.00001)
	{
		double n = sign*Chi2(t+step);
		//cout << old << " " << n << " " << t << " " << step << endl;
		if(n < old)
		{
			t = t + step;
			old = n;
			right = false;
		}	
		else
		{
			if(right)
			{
				step = (-0.5*step);
			}
			else
			{
				step = (-1.*step);
				right = true;
			}
		}
	}
	return pair<double, double>(t, old);
}

TLorentzVector NeutrinoSolver::GetBest(double metx, double mety, double metxerr, double metyerr, double metxyrho, double& test)
{
	if(ERROR){ test = -1; return(TLorentzVector(0.,0.,0.,0.));}

	MET(0,0) = metx;
	MET(1,0) = mety;

	VM(0,0) = metxerr*metxerr;
	VM(1,1) = metyerr*metyerr;
	VM(1,0) = metxerr*metyerr*metxyrho;
	VM(0,1) = metxerr*metyerr*metxyrho;
		
	VM.Invert();

	pair<double, double> maximum = Extrem(0., false);
	pair<double, double> minimuma = Extrem(maximum.first+0.1, true);
	pair<double, double> minimumb = Extrem(maximum.first-0.1, true);

	if(minimuma.second > minimumb.second)
	{
		test = minimumb.second;
		return(GetSolution(minimumb.first));	
	}
	else
	{
		test = minimuma.second;
		return(GetSolution(minimuma.first));	
	}


//	double testmin = 1E100;
//	double t = 0.;
//	for(double tb = 0 ; tb < 6.4 ; tb+=0.1)
//	{
//		double test = Chi2(tb);
//		if(test < testmin)
//		{
//			testmin = test;
//			t = tb;
//		}
//	}
//	double step = 0.05;
//	double old = Chi2(t);
//	bool right = true;
//	while(Abs(step) > 0.00001)
//	{
//		double n = Chi2(t+step);
//		//cout << old << " " << n << " " << t << " " << step << endl;
//		if(n < old)
//		{
//			t = t + step;
//			old = n;
//			right = true;
//		}	
//		else
//		{
//			if(right)
//			{
//				step = (-0.5*step);
//			}
//			else
//			{
//				step = (-1.*step);
//				right = false;
//			}
//		}
//	}
//
//
////	if(old > 10)
////	{
////		cout << "Solution: " << t << " " << old << endl;
////		for(double tb = 0 ; tb < 7 ; tb+=0.01)
////		{
////			cout << tb << " " << Chi2(tb) << endl;
////		}
////	}
//
//	test = old;
//
//	return(GetSolution(t));	
}
