#ifndef _DSP_H_
#define _DSP_H_

#include "Complex.h"

#include <iostream>
#include <math.h>

#define BLACKMAN 0
#define HANNING 1
#define HAMMING 2

using namespace std;

const double PIE  = M_PI;
const double PI2 = PIE*2.0;
const double SQRT2 = M_SQRT2;

Vector<Complex> rdft(const Vector<double> x) {
  const int    NT = x.Dim();
  const int    NF = NT/2+1;
  const double N  = (double)NT;
  const double w  = PI2/N;
  
  Vector<Complex> X(NF);
  
  for (int ki=0; ki<NF; ki++) {
    double k = (double)ki;
    for (int ti=0; ti<NT; ti++) {
      double t = (double)ti;
      Complex z = Complex::Polar(x(ti),-k*w*t);
      X(ki)+=z;
    }
  } 
  return X;
}

Vector<double> irdft(const Vector<Complex> X) {
  const int    NF = X.Dim();
  const int    NT = (NF-1)*2;
  const double N  = (double)NT;
  const double w  = PI2/N;

  Vector<Complex> F(X);  
  Vector<double> x(NT);
  
  for (int f=0; f<NF; f++) {
    F(f).Real(F(f).Real()/(N/2.0));
    F(f).Imag(F(f).Imag()/(N/-2.0));
  }
  
  F(0).Real(F(0).Real()/2.0);
  F(NF-1).Real(F(0).Real()/2.0);
  
  for (int f=0; f<NF; f++) {
    double k = (double)f;
    for (int t=0; t<NT; t++) { 
      double i = (double)t;
      x(t) += F(f).Real() * cos(k*w*t) + F(f).Imag() * sin(k*w*t);
    }
  }
  return x;
}

Vector<Complex> cdft(const Vector<Complex>& z, bool inverse) {
	const int N = z.Dim();
	const double norm = 1.0/(double)N;
	const double w = (2*PIE)*norm;
	double sign = -1.0;
	if (inverse) sign=1.0;
	
	Vector<Complex> Z(N);
		
	for (int k=0; k<N; k++) {
		const double kw = (double)k*w*sign;
		for (int t=0; t<N; t++) {
			const double kwt = kw*(double)t;
			Z(k)+=Complex::Polar(1.0,kwt)*z(t);
		}	
		if (inverse) Z(k)*=norm;
	}
	return Z;
}

void fft(Vector<Complex>& z, bool inverse) {
	const int N = z.Dim();
	const double sign = inverse ? 1.0 : -1.0;
	
    int i, j, k, l, m, M, M2, ip;
    Complex tmp,u,w;   

	m = (int)log2((double)N);	    							//number of 
	j = 0;
	
    for (i=0; i<(N-1); i++) {									//bit reversal
        if (i<j) {
            tmp=z(j); 
            z(j)=z(i); 
            z(i)=tmp;
        }
        k = N/2;
        while (k<=j) {
            j-=k;
            k/=2;
        }
        j+=k;
    }
    
    M=1;
    for (l=1; l<=m; l++)  {										//loop each stage
        M2=M;
        M*=2;
        u = Complex::Polar(1.0,0.0);							//complex sinusoid
        w = Complex::Polar(1.0,(sign*(M_PI/(double)M2)));		//base freq for sub matrix
        
        for (j=0; j<M2; j++) {									//loop each sub dft
            i = j;
            while(i < N) {										//loop each butterfly
                ip=i+M2;
                tmp=(z(ip)*u);									//calculate butterfly
                z(ip)=(z(i)-tmp);								//positive frequency
                z(i)+=tmp;										//negative frequency
                i+=M;
            }
            u*=w;
        }
    }
}


//Blackman. N specifies the number of samples, must be odd.
Vector<double> blackmanWindow(const int N) {
	const double M = N-1;
	Vector<double> h(N);
	for (int i=0; i<N; i++) h(i) = 0.42 - (0.5 * cos((2*M_PI)*(double)i/M)) + (0.08*cos((4*M_PI)*(double)i/M));
	return h;
}

//Hanning. N specifies the number of samples, must be odd.
Vector<double> hanningWindow(const int N) {
	const double M = N-1;
	Vector<double> h(N);
	for (int i=0; i<N; i++) h(i) = 0.5 * (1.0 - cos((2*M_PI)*(double)i/M));
	return h;
}

//Hamming. N specifies the number of samples, must be odd.
Vector<double> hammingWindow(const int N) {
	const double M = (double)(N-1);
	Vector<double> h(N);
	for (int i=0; i<N; i++) h(i) = 0.54 - (0.46*cos((2*M_PI)*(double)i/M));
	return h;
}

//sinc fuction. N specifies the number of samples, must be odd.
Vector<double> genSinc(const int N, const double cutf) {
	const double M = N-1;
	Vector<double> s(N);

	for (int i=0; i<N; i++) {
		double t = (double)i-M/2.0; 
		if (t == 0.0) s(i)=2.0*cutf;
		else s(i)=sin((2*M_PI)*cutf*t)/(M_PI*t);
	}        
	return s;
}

//low-pass FIR sinc filter. N: wraplen(even), M: filterlen(odd) and 0 < lcutf < .5. 
Vector<double> lpassSinc(const int N, const int M, const double cutf, const int window) {
	Vector<double> r = blackmanWindow(M); 
	Vector<double> s = genSinc(M,cutf);
	Vector<double> h(N);
	for (int i=0; i<M; i++) h(i)=s(i)*r(i);
	return h;
}

//high-pass FIR sinc filter. N: wraplen(even), M: filterlen(odd) and 0 < hcutf < .5. 
Vector<double> hpassSinc(const int N, const int M, const double cutf, const int window) {
	Vector<double> h = lpassSinc(N,M,cutf,window);
	for (int i=0; i<M; i++) h(i)*=-1.0;
	h(M/2)+=1.0;
	return h;	
}

//band-stop FIR sinc filter. N: wraplen(even), M: filterlen(odd) and 0 < l/hcutf < .5. 
Vector<double> bstopSinc(const int N, const int M, const double lcutf, const double hcutf, const int window) {
  Vector<double> hl = lpassSinc(N,M,lcutf,window);
  Vector<double> hh = hpassSinc(N,M,hcutf,window);
  Vector<double> h(N);  
  for (int i=0; i<M; i++) h(i)=hl(i)+hh(i);
  return h;
}

//band-pass FIR sinc filter. N: wraplen(even), M: filterlen(odd) and 0 < l/hcutf < .5.  
Vector<double> bpassSinc(const int N, const int M, const double lcutf, const double hcutf, const int window) {
  Vector<double> h = bstopSinc(N,M,lcutf,hcutf,window);
  for (int i=0; i<M; i++) h(i)*=-1.0;
  h(M/2)+=1.0;
  return h;
}

Vector<double> ramp(const int N, const double a, const double b) {
	Vector<double> rmp(N);
	for (int i=0; i<N; i++) rmp(i)=(double)i*a+b;
	return rmp;
}

Matrix<double> drt(const Matrix<double>& src) {
	const double M = (double)src.Dimx();
	const double T = PIE*(M-1);
	const double R = 2*M-1;
	const double dt = PIE/T;
	const double dp = 1.0/sqrt(2.0);//.5*sqrt(1^2+1^2)
	const double xymin = -(M-1)/2;
	const double pmin = -(R-1)/2;
		
	Matrix<double> tgt((int)ceil(T),(int)R,src.Dimz());
	
	for (int t=0; t<(int)ceil(T); t++) {
		const double theta = t*dt;
		const double cost = cos(theta);
		const double sint = sin(theta);
		const double poff = xymin*(cost+sint);
		if (sint > 1.0/SQRT2) {
			const double a = -cost/sint;//x/y
			for (int r=0; r<(int)R; r++) {
				const double p = (pmin+(double)r)*dp;
				const double b = (p-poff)/sint;//scale rho to function with period
				Vector<double> sum(src.Dimz());
				for (int m=0; m<(int)M; m++) {
					unsigned int n = (unsigned int)round((double)m*a+b);
					if (n>=0 && n<(int)M) sum+=src(m,n);
				}
				tgt(t,r)=sum/sint;
			}
		}
		else {
			const double a = -sint/cost;//y/x
			for (int r=0; r<(int)R; r++) {
				const double p = (pmin+(double)r)*dp;
				const double b = (p-poff)/cost;
				Vector<double> sum(src.Dimz());
				for (int n=0; n<(int)M; n++) {
					unsigned int m = (unsigned int)round((double)n*a+b);
					if (m>=0 && m<(unsigned int)M) sum+=src(m,n);
				}
				tgt(t,r)=sum/fabs(cost);
			}		
		}
	}
//	tgt*=(1.0/T);
	return tgt;
}
#endif
