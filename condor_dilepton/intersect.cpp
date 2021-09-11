//determine whether two quadric have intersections or not
void ShengJin(double a,double b,double c,double d,vector<double> &X123)
{
	double A=b*b-3*a*c;
	double B=b*c-9*a*d;
	double C=c*c-3*b*d;
	double f=B*B-4*A*C;
	double i_value;
	double Y1,Y2;
	if (abs(A)<1e-6 && abs(B)<1e-6)//公式1
	{
		X123.push_back(-b/(3*a));
		X123.push_back(-b/(3*a));
		X123.push_back(-b/(3*a));
	}
	else if (abs(f)<1e-6)   //公式3
	{
		double K=B/A;
		X123.push_back(-b/a+K);
		X123.push_back(-K/2);
		X123.push_back(-K/2);
	}
	else if (f>1e-6)      //公式2
	{
		Y1=A*b+3*a*(-B+sqrt(f))/2;
		Y2=A*b+3*a*(-B-sqrt(f))/2;
		double Y1_value=(Y1/fabs(Y1))*pow((double)fabs(Y1),1.0/3);
		double Y2_value=(Y2/fabs(Y2))*pow((double)fabs(Y2),1.0/3);
		X123.push_back((-b-Y1_value-Y2_value)/(3*a));//虚根我不要
		//i_value=sqrt(3.0)/2*(Y1_value-Y2_value)/(3*a);
		//if (fabs(i_value)<1e-1)
		//{
		//	X123.push_back((-b+0.5*(Y1_value+Y2_value))/(3*a));
		//}
	}
	else if (f<-1e-6)   //公式4
	{
		double T=(2*A*b-3*a*B)/(2*A*sqrt(A));
		double S=acos(T);
		X123.push_back((-b-2*sqrt(A)*cos(S/3))/(3*a));
		X123.push_back((-b+sqrt(A)*(cos(S/3)+sqrt(3.0)*sin(S/3)))/(3*a));
		X123.push_back((-b+sqrt(A)*(cos(S/3)-sqrt(3.0)*sin(S/3)))/(3*a));
	}
}
bool intersect1(TMatrixD A, TMatrixD B){
	double tmp1,tmp2,tmp3,tmp4;
	TMatrixD G(3,3);
	double a,b;
	a=A(0,0)*A(1,1)*A(2,2)+A(0,1)*A(1,2)*A(2,0)+A(1,0)*A(2,1)*A(0,2)-A(0,2)*A(2,0)*A(1,1)-A(2,1)*A(1,2)*A(0,0)-A(1,0)*A(0,1)*A(2,2);
	b=B(0,0)*B(1,1)*B(2,2)+B(0,1)*B(1,2)*B(2,0)+B(1,0)*B(2,1)*B(0,2)-B(0,2)*B(2,0)*B(1,1)-B(2,1)*B(1,2)*B(0,0)-B(1,0)*B(0,1)*B(2,2);
	if(abs(a)<1e-6){
		G=B;
	}
	else if(abs(b)<1e-6){
		G=A;
		A=B;
	}
	else{
		TMatrixD E=A.Invert()*B;
		double a,b,c,d;
		a=-1;
		b=E(0,0)+E(1,1)+E(2,2);
		c=-E(0,0)*E(1,1)-E(0,0)*E(2,2)-E(1,1)*E(2,2)+E(0,2)*E(2,0)+E(2,1)*E(1,2)+E(1,0)*E(0,1);
		d=E(0,0)*E(1,1)*E(2,2)+E(0,1)*E(1,2)*E(2,0)+E(1,0)*E(2,1)*E(0,2)-E(0,2)*E(2,0)*E(1,1)-E(2,1)*E(1,2)*E(0,0)-E(1,0)*E(0,1)*E(2,2);
		vector<double> values;
		//cout<<a<<endl;
		//cout<<b<<endl;
		//cout<<c<<endl;
		//cout<<d<<endl;
		ShengJin(a,b,c,d,values);
		//cout<<endl;
	
		/*TMatrixDEigen C(E);
		TMatrixD D=C.GetEigenValues();*/
	//	cout<<values[0]<<endl;
		G=B-values[0]*A.Invert();
		
	
	}
	
	/*cout<<"init G:"<<endl;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++)
				cout<<G(i,j)<<" ";
			cout<<endl;
		}
		cout<<"init A:"<<endl;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++)
				cout<<A(i,j)<<" ";
			cout<<endl;
		}*/
	//cout<<values[0]<<endl;
	//cout<<G(0,0)<<endl;
	//cout<<G(1,0)<<endl;
	//cout<<G(1,1)<<endl;
	if(G(0,0)==0&&G(1,1)==0&&G(1,0)==0){
		/*cout<<"G:"<<endl;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++)
				cout<<G(i,j)<<" ";
			cout<<endl;
		}
		cout<<"A:"<<endl;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++)
				cout<<A(i,j)<<" ";
			cout<<endl;
		}*/
		if(G(0,2)==0&&G(1,2)==0){
			if(G(2,2)!=0)
				return false;
			else
				return true;
		}
		else if(G(1,2)){
		//	cout<<"function are:"<<endl;
		//	cout<<A(0,0)+A(1,1)*pow(G(0,2)/G(1,2),2)<<"x^2+"<<2*A(0,2)+G(0,2)*G(2,2)/pow(G(1,2),2)-2*A(1,2)*G(0,2)/G(1,2)<<"x+"<<A(2,2)+A(1,1)/4*pow(G(2,2)/G(1,2),2)-A(1,2)*G(2,2)/G(1,2)<<endl;
			double delta=pow(2*A(0,2)+G(0,2)*G(2,2)/pow(G(1,2),2)-2*A(1,2)*G(0,2)/G(1,2),2)-4*(A(0,0)+A(1,1)*pow(G(0,2)/G(1,2),2))*(A(2,2)+A(1,1)/4*pow(G(2,2)/G(1,2),2)-A(1,2)*G(2,2)/G(1,2));
			if(delta<0) return false;
			else if(A(0,0)+A(1,1)*pow(G(0,2)/G(1,2),2)==0&&delta==0) return false;
			else return true;		
		}
		else{
		//	cout<<"function are:"<<endl;
		//	cout<<A(1,1)+A(0,0)*pow(G(1,2)/G(0,2),2)<<"x^2+"<<2*A(1,2)+G(1,2)*G(2,2)/pow(G(0,2),2)-2*A(0,2)*G(1,2)/G(0,2)<<"x+"<<A(2,2)+A(0,0)/4*pow(G(2,2)/G(0,2),2)-A(0,2)*G(2,2)/G(0,2)<<endl;
			double delta=pow(2*A(1,2)+G(1,2)*G(2,2)/pow(G(0,2),2)-2*A(0,2)*G(1,2)/G(0,2),2)-4*(A(1,1)+A(0,0)*pow(G(1,2)/G(0,2),2))*(A(2,2)+A(0,0)/4*pow(G(2,2)/G(0,2),2)-A(0,2)*G(2,2)/G(0,2));
			if(delta<0) return false;
			else if(A(1,1)+A(0,0)*pow(G(1,2)/G(0,2),2)==0&&delta==0) return false;
			else return true;
		}
	}
	else{
		//rotation
		double theta;
		if(G(0,0)-G(1,1)){
			theta=atan(2*G(0,1)/(G(0,0)-G(1,1)))/2;
		}
		else theta=3.1415926/4;
		//cout<<theta<<endl;
		tmp1=G(0,2);
		tmp2=G(1,2);
		G(0,0)=G(0,0)+G(0,1)*tan(theta);
		G(1,1)=G(1,1)-G(0,1)*tan(theta);
		G(0,1)=0;
		G(1,0)=0;
		G(0,2)=tmp1*cos(theta)+tmp2*sin(theta);
		G(1,2)=tmp2*cos(theta)-tmp1*sin(theta);
		G(2,0)=G(0,2);
		G(2,1)=G(1,2);
		G(2,2)=G(2,2);
		tmp1=A(0,0);
		tmp2=A(0,1);
		tmp3=A(1,1);
		A(0,0)=tmp1*pow(cos(theta),2)+tmp2*sin(2*theta)+tmp3*pow(sin(theta),2);
		A(1,1)=tmp1*pow(sin(theta),2)-tmp2*sin(2*theta)+tmp3*pow(cos(theta),2);
		A(0,1)=1/2*(tmp3-tmp1)*sin(2*theta)+tmp2*cos(2*theta);
		A(1,0)=A(0,1);
		A(0,2)=tmp1*cos(theta)+tmp3*sin(theta);
		A(1,2)=tmp3*cos(theta)-tmp1*sin(theta);
		A(2,0)=A(0,2);
		A(2,1)=A(1,2);
		A(2,2)=A(2,2);
		//translation
		if(G(0,0)&&G(1,1)){
			tmp1=A(0,2);
			tmp2=A(1,2);
			double i=G(0,2)/G(0,0);
			double j=G(1,2)/G(1,1);
			G(0,2)=0;
			G(1,2)=0;
			G(2,0)=G(0,2);
			G(2,1)=G(1,2);
			G(2,2)=0;
			A(0,2)=tmp1-i*A(0,0)-j*A(0,1);
			A(1,2)=tmp2-j*A(1,1)-i*A(0,1);
			A(2,0)=A(0,2);
			A(2,1)=A(1,2);
			A(2,2)=A(2,2)+A(0,0)*pow(i,2)+A(1,1)*pow(j,2)-2*tmp1*i-2*tmp2*j+2*A(0,1)*i*j;
		/*	cout<<"G:"<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++)
					cout<<G(i,j)<<" ";
				cout<<endl;
			}
			cout<<"A:"<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++)
					cout<<A(i,j)<<" ";
				cout<<endl;
			}*/
			if(G(0,0)*G(1,1)>0){
				if(A(2,2)) return false;
				else return true;
			}
			else{
			//	cout<<"function are:"<<endl;
			//	cout<<A(0,0)-A(1,1)*G(0,0)/G(1,1)+2*A(0,1)*sqrt(-G(0,0)/G(1,1))<<"x^2+"<<2*A(0,2)+2*A(1,2)*sqrt(-G(0,0)/G(1,1))<<"x+"<<A(2,2)<<endl;
			//	cout<<A(0,0)-A(1,1)*G(0,0)/G(1,1)-2*A(0,1)*sqrt(-G(0,0)/G(1,1))<<"x^2+"<<2*A(0,2)-2*A(1,2)*sqrt(-G(0,0)/G(1,1))<<"x+"<<A(2,2)<<endl;
				double delta1=pow(2*A(0,2)+2*A(1,2)*sqrt(-G(0,0)/G(1,1)),2)-4*(A(0,0)-A(1,1)*G(0,0)/G(1,1)+2*A(0,1)*sqrt(-G(0,0)/G(1,1)))*A(2,2);
				double delta2=pow(2*A(0,2)-2*A(1,2)*sqrt(-G(0,0)/G(1,1)),2)-4*(A(0,0)-A(1,1)*G(0,0)/G(1,1)-2*A(0,1)*sqrt(-G(0,0)/G(1,1)))*A(2,2);
				if(delta1<0&&delta2<0) return false;
				else if((A(0,0)-A(1,1)*G(0,0)/G(1,1)+2*A(0,1)*sqrt(-G(0,0)/G(1,1))==0)&&(A(0,0)-A(1,1)*G(0,0)/G(1,1)-2*A(0,1)*sqrt(-G(0,0)/G(1,1))==0)&&delta1==0&&delta2==0) return false;
				else return true;
			}
		}
		else if(G(0,0)){
			tmp1=G(0,2);
			tmp2=G(1,2);
			tmp3=A(0,2);
			tmp4=A(1,2);
			double i=G(0,2)/G(0,0);
			G(0,2)=0;
			G(1,2)=0;
			G(2,0)=G(0,2);
			G(2,1)=G(1,2);
			G(2,2)=G(2,2)-pow(tmp1,2)/G(0,0);
			A(0,2)=tmp3-i*A(0,0);
			A(1,2)=tmp4-i*A(0,1);
			A(2,0)=A(0,2);
			A(2,1)=A(1,2);
			A(2,2)=A(2,2)+A(0,0)*pow(i,2)-2*tmp3*i;
			/*cout<<"G:"<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++)
					cout<<G(i,j)<<" ";
				cout<<endl;
			}
			cout<<"A:"<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++)
					cout<<A(i,j)<<" ";
				cout<<endl;
			}*/
			if(G(0,0)*G(2,2)>0) return false;
			else{
				//cout<<"function are:"<<endl;
				//cout<<A(1,1)<<"x^2+"<<2*A(1,2)-2*A(0,1)*sqrt(-G(2,2)/G(0,0))<<"x+"<<A(2,2)-A(0,0)*G(2,2)/G(0,0)-2*A(0,2)*sqrt(-G(2,2)/G(0,0))<<endl;
				//cout<<A(1,1)<<"x^2+"<<2*A(1,2)+2*A(0,1)*sqrt(-G(2,2)/G(0,0))<<"x+"<<A(2,2)-A(0,0)*G(2,2)/G(0,0)+2*A(0,2)*sqrt(-G(2,2)/G(0,0))<<endl;
				double delta1=pow(2*A(1,2)+2*A(0,1)*sqrt(-G(2,2)/G(0,0)),2)-4*A(1,1)*(A(2,2)-A(0,0)*G(2,2)/G(0,0)+2*A(0,2)*sqrt(-G(2,2)/G(0,0)));
				double delta2=pow(2*A(1,2)-2*A(0,1)*sqrt(-G(2,2)/G(0,0)),2)-4*A(1,1)*(A(2,2)-A(0,0)*G(2,2)/G(0,0)-2*A(0,2)*sqrt(-G(2,2)/G(0,0)));
				if(delta1<0&&delta2<0) return false;
				else if(A(1,1)&&delta1==0&&delta2==0) return false;
				else return true;
			}
		}
		else{
			tmp1=G(0,2);
			tmp2=G(1,2);
			tmp3=A(0,2);
			tmp4=A(1,2);
			double j=G(1,2)/G(1,1);
			G(0,2)=0;
			G(1,2)=0;
			G(2,0)=G(0,2);
			G(2,1)=G(1,2);
			G(2,2)=G(2,2)-pow(tmp2,2)/G(1,1);
			A(0,2)=tmp3-j*A(0,1);
			A(1,2)=tmp4-j*A(1,1);
			A(2,0)=A(0,2);
			A(2,1)=A(1,2);
			A(2,2)=A(2,2)+A(1,1)*pow(j,2)-2*tmp4*j;
		/*	cout<<"G:"<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++)
					cout<<G(i,j)<<" ";
				cout<<endl;
			}
			cout<<"A:"<<endl;
			for(int i=0;i<3;i++){
				for(int j=0;j<3;j++)
					cout<<A(i,j)<<" ";
				cout<<endl;
			}*/
			if(G(1,1)*G(2,2)>0) return false;
			else{
			//	cout<<"function are:"<<endl;
			//	cout<<A(0,0)<<"x^2+"<<2*A(0,2)+2*A(0,1)*sqrt(-G(2,2)/G(1,1))<<"x+"<<A(2,2)-A(1,1)*G(2,2)/G(1,1)+2*A(1,2)*sqrt(-G(2,2)/G(1,1))<<endl;
			//	cout<<A(0,0)<<"x^2+"<<2*A(0,2)-2*A(0,1)*sqrt(-G(2,2)/G(1,1))<<"x+"<<A(2,2)-A(1,1)*G(2,2)/G(1,1)-2*A(1,2)*sqrt(-G(2,2)/G(1,1))<<endl;
				double delta1=pow(2*A(0,2)+2*A(0,1)*sqrt(-G(2,2)/G(1,1)),2)-4*A(0,0)*(A(2,2)-A(1,1)*G(2,2)/G(1,1)+2*A(1,2)*sqrt(-G(2,2)/G(1,1)));
				double delta2=pow(2*A(0,2)-2*A(0,1)*sqrt(-G(2,2)/G(1,1)),2)-4*A(0,0)*(A(2,2)-A(1,1)*G(2,2)/G(1,1)-2*A(1,2)*sqrt(-G(2,2)/G(1,1)));
				if(delta1<0&&delta2<0) return false;
				else if(A(0,0)&&delta1==0&&delta2==0) return false;
				else return true;
			}
		}
	}
}

/*void intersect(){
	TMatrixD A(3,3);
	TMatrixD B(3,3);
	A(0, 0) = 0; 
	A(1, 0) = 0.5; 
	A(2, 0) = 0; 
	A(0, 1) = 0.5; 
	A(1, 1) = 0; 
	A(2, 1) = 0; 
	A(0, 2) = 0; 
	A(1, 2) = 0; 
	A(2, 2) = -1;
	B(0, 0) = 1; 
	B(1, 0) = 0.; 
	B(2, 0) = -1; 
	B(0, 1) = 0; 
	B(1, 1) = 1; 
	B(2, 1) = 0; 
	B(0, 2) = -1; 
	B(1, 2) = 0; 
	B(2, 2) = 0;
	//A.Determinant();
	int i=intersect1(A,B);
	cout<<i<<endl;
}*/