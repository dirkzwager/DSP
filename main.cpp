/*
$ g++ -o dsp main.cpp
$ rm -rf *.o
$ ./dsp

Code bestaat uit 4 secties:

1. AUXILARY PLOTTING AND ANALYSIS FUNCTIONS
4. PROJECTION ALGORITMS

The section PROJECTION ALGORITMS contains the radon/hough transform algorithms.
*/

#include <iostream>
#include <fstream>
#include <string>
#include "GImg.h"

using namespace std;

string iname = "homer127.pgm";
string lname = "lhamming"+iname;
string nname = "nhamming"+iname;
string rawname = "rawhamming"+iname;
string filtname = "filthamming"+iname;
string qname = "Qhamming"+iname;

const double NORM = 255;

/*
  AUXILARY PLOTTING AND ANALYSIS FUNCTIONS
*/

Matrix<double> logTransform(Matrix<Complex> X, const double c=1.0) {
  Matrix<double> P(X.Dimx(),X.Dimy(),1);
  for (int u=0; u<X.Dimx(); u++) 
    for (int v=0; v<X.Dimy(); v++) 
      P(u,v)(0)=c*log10(1.0+X(u,v)(0).Mod());
  return P;
}

//multiply image by -1^(x+y) shifting the tranform of the image to Dimx/2,Dimy/2 (0,0)
void centreSpectrum(Matrix<double>& f) {
  for (int x=0; x<f.Dimx(); x++) 
    for (int y=0; y<f.Dimy(); y++) 
      f(x,y)(0)*=pow(-1.0,(double)(x+y));
}

//write the real component of a signal to a text file to plot with eg gnuplot
void plotReal(Matrix<Complex>& z, const std::string& fname) {
  fstream out;
  out.open(fname.c_str(),std::ios::out|std::ios::trunc);
  if (out) {
    for(int x=0; x<z.Dimx(); x++) {
      for(int y=0; y<z.Dimy(); y++) 
	out << x << " " << y <<  " " << z(x,y)(0).Real() << endl;
      out << std::endl;
    }
    out.close();
  }
  else std::cout << "Error, unable to write file: " << fname << std:: endl;
}

//write the imaginary component of a signal to a text file to plot with eg gnuplot
void plotImag(Matrix<Complex>& z, const std::string& fname) {
  fstream out;
  out.open(fname.c_str(),std::ios::out|std::ios::trunc);
  if (out) {
    for(int x=0; x<z.Dimx(); x++) {
      for(int y=0; y<z.Dimy(); y++) 
	out << x << " " << y <<  " " << z(x,y)(0).Imag() << endl;
      out << std::endl;
    }
    out.close();
  }
  else std::cout << "Error, unable to write file: " << fname << std:: endl;
}

//plot the magnitude of a complex signal to a text file to plot with eg gnuplot
void plotMod(Matrix<Complex>& z, const std::string& fname) {
  fstream out;
  out.open(fname.c_str(),std::ios::out|std::ios::trunc);
  if (out) {
    for(int x=0; x<z.Dimx(); x++) {
      for(int y=0; y<z.Dimy(); y++) 
	out << x << " " << y <<  " " << z(x,y)(0).Mod() << endl;
      out << std::endl;
    }
    out.close();
  }
  else std::cout << "Error, unable to write file: " << fname << std:: endl;
}

//write the phase component of a signal to a text file to plot with eg gnuplot
void plotArg(Matrix<Complex>& z, const std::string& fname) {
  fstream out;
  out.open(fname.c_str(),std::ios::out|std::ios::trunc);
  if (out) {
    for(int x=0; x<z.Dimx(); x++) {
      for(int y=0; y<z.Dimy(); y++) 
	out << x << " " << y <<  " " << z(x,y)(0).Arg() << endl;
      out << std::endl;
    }
    out.close();
  }
  else std::cout << "Error, unable to write file: " << fname << std:: endl;
}

//pad a real vector with complex part of 0
Vector<Complex> vector2Complex(const Vector<double>& v) {
  Vector<Complex> z(v.Dim());
  for (int i=0; i<v.Dim(); i++) z(i)=Complex(v(i));
  return z;
}

//take the real component of a complex vector
Vector<double> vector2Real(const Vector<Complex>& v) {
  Vector<double> z(v.Dim());
  for (int i=0; i<v.Dim(); i++) z(i)=v(i).Real();
  return z;
}

//make a complex matrix out of a real matrix
Matrix<Complex> matrix2Complex(const Matrix<double>& m) {
  Matrix<Complex> z(m.Dimx(),m.Dimy(),1);
  for (int i=0; i<m.Dimx(); i++)
    for (int j=0; j<m.Dimy(); j++) z(i,j)(0)=Complex(m(i,j)(0));
  return z;
}

//take the real component of a complex matrix
Matrix<double> matrix2Real(const Matrix<Complex>& m) {
  Matrix<double> z(m.Dimx(),m.Dimy(),1);
  for (int i=0; i<m.Dimx(); i++)
    for (int j=0; j<m.Dimy(); j++) z(i,j)(0)=m(i,j)(0).Real();
  return z;
}

//take the imaginary component of a complex matrix
Matrix<double> matrix2Imag(const Matrix<Complex>& m) {
  Matrix<double> z(m.Dimx(),m.Dimy(),1);
  for (int i=0; i<m.Dimx(); i++)
    for (int j=0; j<m.Dimy(); j++) z(i,j)(0)=m(i,j)(0).Imag();
  return z;
}

//perform a 1D fourier transform on rows or cols of a matrix.
Matrix<Complex> matrixcdft(const Matrix<Complex>& f, const string& rc, bool inverse=false) {
  Matrix<Complex> F(f.Dimx(),f.Dimy(),1);
  
  if (rc=="cols") {
    for (int x=0; x<f.Dimx(); x++) 
      F.Col(cdft(f.Col(x,0),inverse),x,0);
  }
  else if (rc=="rows") {
    for (int y=0; y<f.Dimy(); y++) 
      F.Row(cdft(f.Row(y,0),inverse),y,0);
  }
  return F;
}

//write a sinogram to a textfile
void writeSinogram(const Matrix<double>& gram, const string& fname) {
  fstream out;
  out.open(fname.c_str(),ios::out);
  
  for (int y=0; y<gram.Dimy(); y++) {
    out << y << " ";
    for (int x=0; x<gram.Dimx(); x++) {
      out << gram(x,y) << " ";
    }
    out << endl;
  }
  out.close();
}

//write a complex sinogram to a file so that it can be plotted with gnuplot for example
void writeComplexSinogram(const Matrix<Complex>& gram, const string& fname) {
  fstream out;
  out.open(fname.c_str(),ios::out);
  
  for (int y=0; y<gram.Dimy(); y++) {
    out << y << " ";
    for (int x=0; x<gram.Dimx(); x++) {
      out.setf(std::ios::fixed);
      out.precision(2);
      out << gram(x,y) << " ";
    }
    out << endl;
  }
  out.close();
}

//Generate the negative frequencies of a signal that are nog generated with a dft
Vector<Complex> generateNegativeFrequencies(Vector<Complex>& rspctrm) {
  const int Nr = rspctrm.Dim();
  const int Nc = (Nr-1)*2;
  Vector<Complex> cspctrm(Nc);
  
  for (int k=Nc/2+1; k<Nc; k++) 
    cspctrm(k)=rspctrm(Nc-k).getConjugate();
  return cspctrm;
}

int rndnear(double x) {
  const double xabs = fabs(x);
  double tail = xabs - floor(xabs);
  if (tail<.5) return (int)((x/xabs)*floor(xabs));
  else return (int)((x/xabs)*ceil(xabs));
}

/*
  FILTERS & SIGNALS
*/

//Apply a hamming window as filter
void windowFilter(Vector<Complex>& h) {
  static Vector<double> w(hammingWindow(h.Dim()));
  for (int m=0; m<h.Dim(); m++) h(m)*=w(m);	
}

//describe the outline of a unit cube. So: 1,1,0.75,0.5,0.25,0,0,0,0,0.25,0.5,0.75,1,1
Vector<double> unitcube(const int K) {
  const double pi = M_PI;
  const double pi075 = .75*pi;
  const double pi025 = .25*pi;
  const double w = (2*pi)/(double)K;
  
  Vector<double> a(K);
  Vector<double> f(K);
  
  a(0)=.5;
  
  for (int ik=1; ik<K; ik++) {
    const double k =(double)ik;
    a(ik)+=(-16.0/((double)(K*K)*w*w*k*k))*(cos(k*pi075)-(cos(k*pi025))) ;
  }
  
  for (int t=0; t<K; t++) {
    for (int k=0; k<K; k++) {
      f(t)+=a(k)*cos((double)k*w*(double)t);
    }
  }

  return f;
}

//the dirac delta equation
Vector<double> dirac(int L, int inf) {
  Vector<double> d(L);
  d(inf-1)=1;
  return d;
}

//make an image with a cube in the center (2D dirac)
double cub(const double theta, const int K=360) {
  const double pi = M_PI;
  const double pi075 = .75*pi;
  const double pi025 = .25*pi;
  const double w = (2*pi)/(double)K;
  
  static Vector<double> a(K);
  static Vector<double> x(K);
  static Vector<double> y(K);
  
  a(0)=.5;
  
  for (int ik=1; ik<K; ik++) {
    const double k =(double)ik;
    a(ik)+=(-16.0/((double)(K*K)*w*w*k*k))*(cos(k*pi075)-(cos(k*pi025))) ;
  }
  
  for (int t=0; t<K; t++) {
    for (int k=0; k<K; k++) {
      x(t)+=a(k)*cos((double)k*w*(double)t);
    }
  }
  
  Vector<Complex> X = cdft(vector2Complex(x),false);
  Vector<Complex> d = vector2Complex(dirac(K,K/4));
  Vector<Complex> D = cdft(d,false);
  D.gnuPlot("D.dat");
  Vector<Complex> A = vector2Complex(a);
  Vector<Complex> F = X*D;
  Vector<double> f = vector2Real(cdft(F,true));
  f-=1.0;
  //x.gnuPlot("x.dat");
  //f.gnuPlot("f.dat");
  
  GImg<double> sqr(51,51,1);
  
  
  for (int t=0; t<K; t++) {
    const int xr = rndnear((x(t)*20)+23);
    const int yr = rndnear((f(t)*20)+26);
    sqr(xr,yr)(0)+=1.0;		 
  }
  
  sqr.Save("square.pgm");
  
  return 0;
}

//construct a ramp filter (triangle 90,45,45 with sides of length N) 
Vector<double> Ramp(const int N) {
  const double T = (double)N;
  const double off = -T+.5;
  const double dp = 1.0;
  Vector<double> r(N);
  
  for (int y=0; y<N; y++) {
    const double p = off+y;
    r(y)=fabs(p*dp);
  }	
  return r;
}


Vector<Complex> RampRec(const int W) {
  const double dp = 1.0/(SQRT2);
  const double off = -.5*((double)W)+.5;
  Vector<Complex> R(W);
  
  for (int k=0; k<W; k++) {
    const double w = (off+(double)k)*(1.0/(double)W);
    R(k).Real(fabs(w));
  }	
  return R;
}

//construct a signal which is an inverted triangle (so y,0,y)
Vector<Complex> rmpf(const int NF) {
  const double W = (double)NF;
  const double min = -.5*(W-1.0);
  const double a = 2.0/(W-1);
  
  Vector<Complex> H(NF);
  
  for (int w=0; w<NF; w++) {
    const double r = -min*fabs(((min+(double)w)*a))*100.0;
    H(w)+=r;
  }
  return H;	
}

//construcs a signal similar to a sinc function
Vector<double> rampfilter(const int N) {
  const double R = (double)N;
  const double min = -.5*(R-1.0);
  const double a = 2.0/(R-1);
  const double kmax = .5*R;
  const double cut = 2.0*PIE*kmax;
  
  Vector<double> rmp(N);
  for (int k=0; k<N; k++) {
    const double r = (-min*((min+(double)k)*a))*1.0/SQRT2;
    const double pulse = kmax*sin(cut*r)/(PIE*r);	
    const double saw   = (1.0-cos(cut*r))/(2.0*PIE*PIE*r*r);
    rmp(k)=pulse-saw;
  }

  for (int k=0; k<N; k++) {
    cout << k << " " << rmp(k) << endl;
  }
  return rmp;
}

//similar to rampfilter but complex
Vector<Complex> backramp(const int N, double off =-1.0) {
  Vector<Complex> h(N);
  if (off==-1.0) off = (double)(N/-2);
  
  const double dp = 1.0/SQRT2;
  const double dnom = 1.0/(dp*dp);
  
  for (int m=0; m<N; m++) {
    const double k = off+(double)m;
    if (k==0.0) h(m)=1.0/(8.0*dp*dp);
    else if ((int)k%2 !=0) h(m)=-1.0/(k*k*PIE*PIE*dp*dp);
    else h(m)=Complex(0.0);
  }
  return h;
}

//construct arbitrary signal
Vector<Complex> signal(int N, int srate, int t1, int t2, int t3) {
  Vector<Complex> v(N);
  double dt = 1.0/(double)srate;
  for (int i=0; i<N; i++) {
    double x = (double)i*dt;
    v(i)+=cos(x*t1*PIE)+cos(x*t2*PIE)+cos(x*t3*PIE);
  }
  return v;
}

Matrix<double> idht(const Matrix<double>& src);

//some experimental filter combining two blackman bandpass filters (not included in examples)
GImg<double> filter3(GImg<double>& a) {
  GImg<double> p = drt(a.Data()); 
  p.Resize(p.Dimx(),p.Dimy()+1);
  
  const Matrix<double> mp = p.Data();
  const int N = a.Dimx();
  const int R = mp.Dimy();
  const int T = mp.Dimx();
  
  const double lcut = 0.050;
  const double hcut = 0.005;
  const double band = .25;//fabs(hcut-lcut);
  const int    M    = (int)(4.0/band)+1;
  
  cout << endl << endl;
  cout << "N: " << N << endl
       << "R: " << R << endl
       << "T: " << T << endl
       << "M: " << M << endl;
  cout << endl << endl;
  
  const Vector<double> bp1 = bpassSinc(R,(int)M,lcut,hcut,BLACKMAN);
  const Vector<Complex> BP1 = rdft(bp1);
  
  const Vector<double>  hp1 = hpassSinc(R,(int)M,hcut,BLACKMAN);
  const Vector<Complex> HP1 = rdft(hp1);
  
  const Vector<double> rmp = rampfilter(R);
  const Vector<Complex> RMP = rdft(rmp);
  
  Matrix<double> q(T,R,1);
  
  for (int t=0; t<T; t++) {
    const Vector<Complex> Pt = rdft(mp.Col(t,0));
    const Vector<Complex> Qt = Pt*RMP;
    const Vector<double>  qt = irdft(Qt);
    q.Col(qt,t,0);
  }	
  
  p.normalize(0.0,512.0);
  p.Save("uout.pgm");
  GImg<double> out = idht(p.Data());
  out.normalize(0.0,512.0);
  out.Save("uradon.pgm");
  return GImg<double>(q);
}

/*
  DSP ALGORITMS
*/

//perform the convolution of two functions f and g
Vector<double> convolve(const Vector<double>& f, const Vector<double>& g) {
  Vector<double> y(f.Dim()+g.Dim()-1);
  
  for (int i=0; i<f.Dim(); i++) 
    for (int j=0; j<g.Dim(); j++) 
      y(i+j) += f(i)*g(j);
  return y;
}

//real discrete 2D fourier transformation
Matrix<Complex> rdft2d(const Matrix<double>& f) {
  const int N = f.Dimx();
  const int M = f.Dimy();
  const double wn = (2.0*PIE)/(double)N;
  const double wm = (2.0*PIE)/(double)M;
  
  Matrix<Complex> F(N,M,1);
  Matrix<Complex> F2(N,M,1);
  
  //cols
  for (int x=0; x<N; x++) 
    for (int k=0; k<M; k++) {
      const double kwm = (double)k*wm;
      for (int y=0; y<M; y++) {
	const double kwmy = kwm*(double)y;
	F2(x,k)(0)+=Complex::Polar(f(x,y)(0),-kwmy);
      }
    }
  //rows
  for (int y=0; y<M; y++) 
    for (int k=0; k<M; k++) {
      const double kwn = (double)k*wn;
      for (int x=0; x<N; x++) {
	const double kwnx = kwn*(double)x;
	F(k,y)(0)+=F2(x,y)(0)*Complex::Polar(1.0,-kwnx);
      }
    }
  return F;
}

//inverse real discrete 2D fourier transformation
Matrix<double> irdft2d(const Matrix<Complex>& F) {
  const int N = F.Dimx();
  const int M = F.Dimy();
  const double wn = (2.0*PIE)/(double)N;
  const double wm = (2.0*PIE)/(double)M;
  const double norm = 1.0/(double)(N*M);
  
  Matrix<double> f(N,M,1);
  Matrix<double> f2(N,M,1);
  
  //cols
  for (int u=0; u<N; u++) 
    for (int k=0; k<M; k++) {
      const double kwm = (double)k*wm;
      for (int v=0; v<M; v++) {
	const double kwmv = kwm*(double)v;
	f2(u,k)(0)+=F(u,v)(0).Real()*cos(kwmv) + F(u,v)(0).Imag()*sin(kwmv);
      }
    }
  //rows
  for (int v=0; v<M; v++) 
    for (int k=0; k<M; k++) {
      const double kwn = (double)k*wn;
      for (int u=0; u<N; u++) {
	const double kwnu = kwn*(double)u;
	f(k,v)(0)+=(f2(u,v)(0)*(cos(kwnu) + sin(kwnu)));
      }
      f(k,v)(0)*=norm;
    }
  return f;
}

/*
PROJECTION ALGORITMS DSP ALGORITMS FILTERS & SIGNALS
*/

//first attempt to filtered backprojection (not included in working examples)
GImg<double> terugprojectie(const GImg<double> f) {
  const int N = f.Dimx();
  const int R = 2*N;
  const int T = (int)ceil(PIE*((double)(N-1)));
  
  const Vector<Complex> H = RampRec(R/2+1);
  Matrix<double> c = f.Data();
  Matrix<double> p = drt(c); 
  cout << "R: " << R <<  " pR " << p.Dimy() << endl;
  p.Resize(T,R);
  
  Matrix<double> q(T,R,1);	
  
  for (int t=0; t<T; t++) {
    Vector<double> pt = p.Col(t,0);
    Vector<Complex> Pt = rdft(pt);
    Vector<Complex> Qt = H*Pt;
    Vector<double>  qt = irdft(Qt);
    q.Col(qt,t,0);
  }
  return GImg<double>(q);
}

//stack based radon transform. dwz een bin bijhouden onderwijl loopen over m en n.
Matrix<double> stackRadon(const Matrix<double> f) {
  const int N = f.Dimx();
  const int R = 2*N-1;                                     //number of beams
  const int T = (int)ceil((double)(N-1)*PIE);              //number of projections
  
  const double dt = PIE/(double)T;                         //angle between projections
  const double dr = 1.0/SQRT2;                             //distance between beams
  
  const double mmin = (double)N/-2.0+.5;                   //beam offset
  const double rmin = (((double)R-1.0)/-2.0)*dr;           //projection offset
  
  cout << endl << "stack: " << endl;
  cout << "N     : " << N << endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dp    : " << dr << endl
       << "mmin : " << mmin << endl
       << "pmin  : " << rmin << endl;
  cout << endl << endl;
  
  Matrix<double> p(T,R,1);
  
  for (int t=0; t<T; t++) {                                //for all projections
    const double theta = (double)t*dt;                     //the particular projection angle
    const double cost = cos(theta);
    const double sint = sin(theta);
    for (int m=0; m<N; m++) {                              //for all beams
      const double x = (double)m+mmin;                     //x of a cell
      for (int n=0; n<N; n++) {
	const double y = (double)n+mmin;                   //y of a cell
	const double rho = x*cost+y*sint;
	const double bin = round((rho-rmin)/dr);           //determine the bin onto which the beam projects
	p(t,(int)bin)(0)+=f(m,n)(0);                       //add the pixel value of a cell to the bin
      }
    }
  }
  return p;
}

//inverse stack radon transformatie, zie stackRadon()
Matrix<double> istackRadon(const Matrix<double> p) {
  const int R = p.Dimy();
  const int T = p.Dimx();
  const int N = (R+1)/2;
  
  const double dt = PIE/(double)T;
  const double dr = 1.0/SQRT2;
  
  const double mmin = (double)N/-2.0+.5;
  const double rmin = (((double)R-1.0)/-2.0)*dr;
  
  cout << endl << "istack: " << endl;
  cout << "N     : " << N << endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dp    : " << dr << endl
       << "mmin : " << mmin << endl
       << "rmin : " << rmin << endl;
  cout << endl << endl;
  
  Matrix<double> p_1(N,N,1);
  
  for (int t=0; t<T; t++) {
    const double theta = (double)t*dt;
    const double cost = cos(theta);
    const double sint = sin(theta);
    for (int m=0; m<N; m++) {
      const double x = (double)m+mmin;
      for (int n=0; n<N; n++) {
	const double y = (double)n+mmin;
	const double rho = x*cost+y*sint;
	const double bin = round((rho-rmin)/dr);
	p_1(m,n)(0)+=p(t,(int)bin)(0);
      }
    }
  }
  return p_1;
}

//discrete normal radon transform
Matrix<double> dnrt(const Matrix<double>& src) {
  const double M = (double)src.Dimx();
  const double T = ceil(PIE*(M-1));                           //number of projection angles
  const double R = 2*M-1;                                     //number of beams
  const double dt = PIE/T;                                    //angle between different projections
  const double dp = 1.0/sqrt(2.0);                            //distance between two beams
  const double xymin = -(M-1)/2.0;                            //min cartesian offset (assuming dimx==dimy)
  const double pmin = -(R-1)/2;                               //min polar offset

  cout << endl << "inter: " << endl;
  cout << "M     : " << M<< endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dp    : " << dp << endl
       << "mmin : " << xymin << endl
       << "rmin  : " << pmin*dp << endl;
  cout << endl << endl;

  //aux vars
  const int iT = (int)ceil(T);
  const int iR = (int)R;
  const int iM = (int)M;
  
  Matrix<double> tgt(iT,iR,src.Dimz());                        //sinogram

  for (int t=0; t<iT; t++) {                                   //for all projection-angles
    const double theta = t*dt;                                 //radians
    const double cost = cos(theta);     
    const double sint = sin(theta);
    const double poff = xymin*(cost+sint);                     //factor
    if (sint > 1.0/SQRT2) {                                    //project on x-axis
      const double a = -cost/sint;                             //x/y (slope)
      for (int r=0; r<iR; r++) {                               //for all beams
	const double p = (pmin+(double)r)*dp;                  //rho of line
	const double b = (p-poff)/sint;                        //scale rho to function with period
	Vector<double> sum(src.Dimz());                        //bin (don`t mind the vector)
	for (int m=0; m<iM; m++) {                             //all values x
	  unsigned int n = (unsigned int)round((double)m*a+b); //y=ax+b
	  if (n>=0 && n<iM)                                    //within boundaries?
	    sum+=src(m,n);
	}
	tgt(t,r)=sum;///sint;
      }
    }
    else {                                                     //project on y-axis
      const double a = -sint/cost;
      for (int r=0; r<iR; r++) {
	const double p = (pmin+(double)r)*dp;
	const double b = (p-poff)/cost;
	Vector<double> sum(src.Dimz());
	for (int n=0; n<iM; n++) {
	  unsigned int m = (unsigned int)round((double)n*a+b); //all values x
	  if (m>=0 && m<iM)
	    sum+=src(m,n);
	}
	tgt(t,r)=sum;///fabs(cost);
      }		
    }
  }
  //tgt*=(1.0/(double)T);
  return tgt;
}

//inverse discrete normal radon transform (with some nice programming dnrt==idnrt)
Matrix<double> idnrt(const Matrix<double>& src) {
  const int T = src.Dimx();                                 //number of projection angles
  const int R = src.Dimy();                                 //number of beams
  const int M = (R+1)/2;                                    //dimx and dimy of the inverse sinogram
  
  const double dt = PIE/(double)T;                          //angle sampling density
  const double dp = 1.0/sqrt(2.0);                          //distance between beams
  const double xymin = -((double)M-1.0)/2.0;                //min cartesian offset
  const double pmin = -((double)R-1.0)/2.0;                 //min polar offset
  
  cout << endl << "config: " << endl;
  cout << "M     : " << M<< endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dp    : " << dp << endl
       << "mmin : " << xymin << endl
       << "rmin  : " << pmin*dp << endl;
  cout << endl << endl;
  
  Matrix<double> tgt(M,M,1);                                  //inverse sinogram
  
  for (int t=0; t<(int)ceil(T); t++) {                        //all projection angles
    const double theta = t*dt;                                //radians
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double poff = xymin*(cost+sint);                     //factor
    if (sint > 1.0/SQRT2) {
      const double a = -cost/sint;                             //x/y (slope)
      for (int r=0; r<(int)R; r++) {                           //all beams
	const double p = (pmin+(double)r)*dp;                  //beam rho
	const double b = (p-poff)/sint;                        //scale rho to function with period
	for (int m=0; m<(int)M; m++) {
	  unsigned int n = (unsigned int)round((double)m*a+b); //y=ax+b
	  if (n>=0 && n<(int)M) 
	    tgt(m,n)+=src(t,r);
	}
      }
    }
    else {
      const double a = -sint/cost;
      for (int r=0; r<(int)R; r++) {
	const double p = (pmin+(double)r)*dp;
	const double b = (p-poff)/cost;
	for (int n=0; n<(int)M; n++) {
	  unsigned int m = (unsigned int)round((double)n*a+b);
	  if (m>=0 && m<(unsigned int)M) 
	    tgt(m,n)+=src(t,r);
	}
      }		
    }
  }
  return tgt;
}

//dit is een poging tot een gefilterde backpojectie
void backProjectGoed(const Matrix<double>& f) {
  const int M = f.Dimx();
  const int T = (int)ceil((double)(M-1)*PIE); //number of projection angles
  const int R = M*2-1;                        //number of beams
  
  const double A = (double)M;
  const double dr = 1.0/SQRT2;                //distance between two beams
  const double dt = PIE/(double)T;            //angle between two projections
  const double pmax = .5*A*SQRT2;             //max polar
  const double mnmin = .5*(A-1);              //min cartesian
  
  Matrix<Complex> p(T,R,1);                   //sinogram
  Matrix<Complex> P(T,R+1,1);                 //frequency domain sinogram
  Matrix<Complex> Q(T,R+1,1);                 //frequency domain filtered sinogram
  Matrix<Complex> q(T,R,1);                   //time domain filtered sinogram
  Matrix<Complex> g(M,M,1);                   //filtered image
  Matrix<Complex> g2(M,M,1);                  //filtered image (real)
  
  Vector<Complex> h = backramp(R,0.0);       //time domain ramp filter
  h.Resize(R+1);
  h.gnuPlot("h.dat");
  Vector<Complex> H = cdft(h,false);         //frequency domain filter
  H.gnuPlot("H.dat");
  
  //1. for each angle compute one 'column' of the sinogram and:
  //2. filter in the frequency domain with filter H
  //3. compute the time domain filtered 'column' of the sinogram
  //4. project back the filtered time domain signal
  for (int t=0; t<T; t++) {
    const double theta = (double)t*dt;
    const double cost = cos(theta);
    const double sint = sin(theta);
    Vector<Complex> pt(R+1);
    for (int m=0; m<M; m++) {
      const double x = (double)m-mnmin;
      for (int n=0; n<M; n++) {
	const double y = (double)n-mnmin;
	const double r = (x*cost+y*sint)/dr+A-1;
	const int ri = (int)round(r);
	pt(ri)+=f(m,n)(0);
	p(t,ri)(0)+=f(n,m)(0);
      }
    }
    
    //filtering
    Vector<Complex> Pt = cdft(pt,false);
    Vector<Complex> Qt = Pt*H;
    Vector<Complex> qt = cdft(Qt,true);
    qt.Resize(R);
    P.Col(Pt,t,0);
    Q.Col(Qt,t,0);
    q.Col(qt,t,0);
    
    //filter in frequency domain and project back
    for (int m=0; m<M; m++) {
      const double x = (double)m-mnmin;
      for (int n=0; n<M; n++) {
	const double y = (double)n-mnmin;
	const double r = (x*cost+y*sint)/dr+A-1;
	const int ri = (int)round(r);
	g(m,n)(0)+=p(t,ri)(0);
	if (theta > dr) 
	  g2(m,n)(0)+=q(t,ri)(0).Real()*sint;
	else
	  g2(m,n)(0)+=q(t,ri)(0).Real()*cost;
      }
    }
    
  }
  
  GImg<double> Rdn(matrix2Real(p));
  Rdn.Save("R.pgm");
  cout << "written R.pgm" << endl;
  
  GImg<double> Pimg(logTransform(P));
  Pimg.Save("P.pgm");
  cout << "written P.pgm" << endl;
  
  GImg<double> Qimg(logTransform(Q));
  Qimg.Save("Q.pgm");	
  cout << "written Q.pgm" << endl;
  
  GImg<double> qimg(matrix2Real(q));
  qimg.normalize(0.0,255.0);
  qimg.Save("q.pgm");
  cout << "written q.pgm" << endl;	
  
  GImg<double> gimg(matrix2Real(g));
  gimg.normalize(0.0,255.0);
  gimg.Save("g.pgm");
  cout << "written g.pgm" << endl;	
  
  GImg<double> g2img(matrix2Real(g2));
  g2img.normalize(0.0,255.0);
  g2img.Save("g2.pgm");	
  cout << "written g2.pgm" << endl;
}

//filtered backprojection using the dnrt algorithm and the backramp filter
Matrix<double> rampdrt(const Matrix<double>& src) {
  fstream out;
  out.open("check.dat");
  const double M = (double)src.Dimx();
  const double T = PIE*(M-1);
  const double R = 2*M-1;
  const double dt = PIE/T;
  const double dp = 1.0/sqrt(2.0);//.5*sqrt(1^2+1^2)
  const double xymin = -(M-1)/2;
  const double pmin = -(R-1)/2;
  
  Matrix<double> tgt((int)ceil(T),(int)R,src.Dimz());
  Matrix<double> g(src.Dimx(),src.Dimy(),1);
  
  Vector<Complex> h = backramp((int)R+1,0.0);
  Vector<Complex> H = cdft(h,false);
  
  int i=0;
  
  for (int t=0; t<(int)ceil(T); t++) {
    const double theta = t*dt;
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double poff = xymin*(cost+sint);
    if (sint > 1.0/SQRT2) {
      const double a = -cost/sint;//x/y
      for (int r=0; r<(int)R; r++) {
	const double p = (pmin+(double)r)*dp;
	i++;
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
    
    Vector<Complex> p = vector2Complex(tgt.Col(t,0));
    p.Resize((int)R+1);
    Vector<Complex> P = cdft(p,false);
    Vector<Complex> q = cdft(P*H,true);
    q.Resize((int)R);
    
    if (sint > 1.0/SQRT2) {
      const double a = -cost/sint;//x/y
      for (int r=0; r<(int)R; r++) {
	const double p = (pmin+(double)r)*dp;
	const double b = (p-poff)/sint;//scale rho to function with period
	Vector<double> sum(src.Dimz());
	for (int m=0; m<(int)M; m++) {
	  unsigned int n = (unsigned int)round((double)m*a+b);
	  if (n>=0 && n<(int)M) g(m,n)(0)+=q(r).Real();
	}
	//g(t,r)=sum/sint;
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
	  if (m>=0 && m<(unsigned int)M) g(m,n)+=q(r).Real();
	}
	//tgt(t,r)=sum/fabs(cost);
      }		
    }
  }
  GImg<double> gout(g);
  gout.normalize(0.0,255.0);
  gout.Save("gout.pgm");
  //	tgt*=(1.0/T);
  return tgt;
}

//filtered backprojection using discrete 'normal' radon transform (liniair interpolation) and the backramp filter
Matrix<double> lindrt(const Matrix<double>& f, int T=-1, int R=-1) {
  const int M = f.Dimx();
  const double A = (double)M;
  R = 2*M-1;
  T = 2*(int)ceil((A-1)*PIE);
  
  const double rmax = .5*(A*SQRT2);
  const double dr = (2*rmax)/(double)(R-1);
  const double dt = PIE/(double)T;
  const double mmax = .5*(A-1);
  
  cout << endl << "config: " << endl;
  cout << "M     : " << M<< endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dr    : " << dr << endl
       << "rmax : " << rmax << endl
       << "mmin : " << mmax << endl;
  cout << endl << endl;
  
  Matrix<double> p(T,R,1);
  Matrix<Complex> q(M,M,1);
  
  Vector<Complex> h = backramp(R);
  Vector<Complex> H = cdft(h,false);
  
  int i=0;
  
  for (int t=0; t<T; t++) {
    const double theta = dt*(double)t;
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double roff = -mmax*(cost+sint);
    Vector<Complex> pt(R);
    if (sint > 1.0/SQRT2) {
      const double alpha = -cost/sint;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rmax;
	const double beta = (rho-roff)/sint;
	for (int m=0; m<M; m++) {
	  const int n = (int)round(alpha*(double)m+beta);
	  if (n>=0 && n<M) {
	    pt(r)+=f(m,n)(0);
	    p(t,r)+=f(m,n);
	  }
	}
	i++;
	pt*=Complex((1.0/sint));
	p(t,r)*=(1.0/sint);
      }
      Vector<Complex> Pt = cdft(pt,false);
      Vector<Complex> qt = cdft(Pt*H,true);
      
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rmax;
	const double beta = (rho-roff)/sint;
	for (int m=0; m<M; m++) {
	  const int n = (int)round(alpha*(double)m+beta);
	  if (n>=0 && n<M) {
	    q(m,n)+=qt(r);
	  }
	}
      }
      
    }
    else {
      const double alpha = -sint/cost;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rmax;
	const double beta = (rho-roff)/cost;
	for (int n=0; n<M; n++) {
	  const int m = (int)round(alpha*(double)n+beta);
	  if (m>=0 && m<M) {
	    pt(r)+=f(m,n)(0);
	    p(t,r)+=f(m,n);
	  }
	}
	i++;
	pt*=Complex((1.0/sint));
	p(t,r)*=(1.0/fabs(cost));
      }
      
      Vector<Complex> Pt = cdft(pt,false);
      Vector<Complex> qt = cdft(Pt*H,true);
      
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rmax;
	const double beta = (rho-roff)/cost;
	for (int n=0; n<M; n++) {
	  const int m = (int)round(alpha*(double)n+beta);
	  if (m>=0 && m<M) {
	    q(m,n)+=qt(r);
	  }
	}
      }
    }
  }
  
  GImg<double> qout(matrix2Real(q));
  qout.normalize(0.0,255.0);
  qout.Save("qout.pgm");
  return p;
}

//inverse discrete normal radon transform using liniair interpolation
Matrix<double> linidrt(const Matrix<double>& p) {
  const int T = p.Dimx();
  const int R = p.Dimy();
  const int M = (R+1)/2;
  
  const double A = (double)M;
  const double rmax = .5*(A*SQRT2);
  const double dr = (2*rmax)/(double)(R-1);
  const double dt = PIE/(double)T;
  const double mmax = .5*(A-1);
  
  cout << endl << "config: " << endl;
  cout << "M     : " << M<< endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dr    : " << dr << endl
       << "rmax : " << rmax << endl
       << "mmin : " << mmax << endl;
  cout << endl << endl;
  
  Matrix<double> g(M,M,1);
  
  for (int t=0; t<T; t++) {
    const double theta = dt*(double)t;
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double roff = -mmax*(cost+sint);
    if (sint > 1.0/SQRT2) {
      const double alpha = -cost/sint;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rmax;
	const double beta = (rho-roff)/sint;
	for (int m=0; m<M; m++) {
	  const int n = (int)round(alpha*(double)m+beta);
	  if (n>=0 && n<M) {
	    g(m,n)+=p(t,r);
	  }
	}
      }
    }
    else {
      const double alpha = -sint/cost;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rmax;
	const double beta = (rho-roff)/cost;
	for (int n=0; n<M; n++) {
	  const int m = (int)round(alpha*(double)n+beta);
	  if (m>=0 && m<M) {
	    g(m,n)+=p(t,r);
	  }
	}
      }
    }
  }
  return g;
}

//discrete hough transform
Matrix<double> dht(const Matrix<double>& src) {
  const double M = (double)src.Dimx();
  const double T = PIE*(M-1);
  const double R = 2*M-1;
  const double dt = PIE/T;
  const double dp = 1.0/SQRT2;
  const double xymin = -(M-1)/2;
  const double pmin = SQRT2*xymin;
  
  Matrix<double> tgt((int)ceil(T),(int)R,src.Dimz());
  
  for (unsigned int m=0; m<(int)M; m++) {
    double x = xymin+(double)m;
    for (unsigned int n=0; n<(int)M; n++) {
      Vector<double> pixval = src(m,(int)(M-1)-n);
      if (pixval.vectorLength() > 0.0) {
	double y = xymin+((M-1)-(double)n);
	for (int t=0; t<(int)ceil(T-1); t++) {
	  double theta = (double)t*dt;
	  double p = x*cos(theta)+y*sin(theta);
	  int pidx = (int)round((p-pmin)/dp);
	  if (pidx >= 0 && pidx < (int)R)	tgt(t,pidx)+=pixval;
	}
      }
    }
  }	
  return tgt;
}

//inverse discrete hough transform
Matrix<double> idht(const Matrix<double>& src) {
  const int M = (src.Dimy()+1)/2;
  const int T = src.Dimx();
  const int R = src.Dimy();
  const double dt = PIE/(double)T;
  const double dp = 1.0/SQRT2;
  const double xymin = -.5*(double)(M-1);
  const double pmin = SQRT2*xymin;
  
  Matrix<double> tgt(M,M,src.Dimz());
  for (unsigned int m=0; m<(int)M; m++) {
    double x = xymin+(double)m;
    for (unsigned int n=0; n<(int)M; n++) {
      double y = xymin+((M-1)-(double)n);
      for (int t=0; t<(int)ceil(T-1); t++) {
	double theta = (double)t*dt;
	double p = x*cos(theta)+y*sin(theta);
	int pidx = (int)round((p-pmin)/dp);
	if (pidx >= 0 && pidx < (int)R)	tgt(m,n)+=src(t,pidx);
      }
    }
  }	
  return tgt;
}

void doTheStackBasedRadon(GImg<double> src) {
  GImg<double> out(stackRadon(src.Data()));
  out.Save("stackRadon.pgm");
  cout << "The sinogram was saved to stackRadon.pgm" << endl;
  GImg<double> tgt(istackRadon(out.Data()));
  tgt.Save("istackRadon.pgm");
  cout << "The backprojected sinogram was saved to istackRadon.pgm" << endl;
}

void doTheDiscreteNormalRadon(GImg<double> src) {
  GImg<double> out(dnrt(src.Data()));
  out.Save("dnrt.pgm");
  cout << "The sinogram was saved to dnrt.pgm" << endl;
  GImg<double> tgt(idnrt(out.Data()));
  tgt.normalize(0,255);
  tgt.Save("idnrt.pgm");
  cout << "The backprojected sinogram was saved to idnrt.pgm" << endl; 
}

void doTheLiniairDiscreteRadonTransform(GImg<double> src) {
  GImg<double> out(lindrt(src.Data()));
  out.Save("lindrt.pgm");
  cout << "The sinogram was saved to lindrt.pgm" << endl;
  GImg<double> tgt(linidrt(out.Data()));
  tgt.normalize(0,255);
  tgt.Save("linidrt.pgm");
  cout << "The backprojected sinogram was saved to linidrt.pgm" << endl;
}

void doTheBackProjectGoed(GImg<double> src) {
  backProjectGoed(src.Data());
}

void doTheDiscreteHough(GImg<double> src) {
  GImg<double> out(dht(src.Data()));
  out.Save("dht.pgm");
  cout << "The sinogram was saved to dht.pgm" << endl;
  GImg<double> tgt(idht(out.Data()));
  tgt.normalize(0,255);
  tgt.Save("idht.pgm");
  cout << "The backprojected sinogram was saved to idth.pgm" << endl;
}

int main() {
  cout << endl << endl;
  string filename = "";
  int choice = 0;
  
  cout << "Welcome\n\nPlease enter filename: ";
  cin >> filename;

  GImg<double> file(filename);

  if (!(file.Dimx()*file.Dimy()>0)) {
    cout << "Fatal error: Invalid file specified." << endl;
    return 0;
  }

  cout << endl
       << "(1) Stack based radon transform." << endl
       << "(2) Discrete Normal Radon transform." << endl
       << "(3) Discrete Radon tranform using liniair interpolation and filtered backprojection with a ramp filter." << endl
       << "(4) Filtered backprojection using a stack based radon-transform and a ramp filter." << endl
       << "(5) Discrete Hough transform." << endl
       << "(6) Exit" << endl
       << endl
       << "Please enter choice: ";

  cin >> choice;
  
  switch (choice) {
  case 1: doTheStackBasedRadon(file.Data());
    break;
  case 2: doTheDiscreteNormalRadon(file.Data());
    break;
  case 3: doTheLiniairDiscreteRadonTransform(file.Data());
    break;
  case 4: doTheBackProjectGoed(file.Data());
    break;
  case 5: doTheDiscreteHough(file.Data());
    break;
  default:
    cout << "Invalid choice, exiting...";
    break;
  }
  cout << endl << endl;
  return 0;
  
}

/*
  GImg<double> rampprojectie(const GImg<double> f) {
  const int N = f.Dimx();
  const int R = 2*N;
  const int T = (int)ceil(PIE*((double)(N-1)));
  
  const Vector<double> h = rampfilter();
  Matrix<double> c = f.Data();
  Matrix<double> p = drt(c); 
  cout << "R: " << R <<  " pR " << p.Dimy() << endl;
  p.Resize(T,R);
  
  Matrix<double> q(T,R,1);	
  
  for (int t=0; t<T; t++) {
  Vector<double> pt = p.Col(t,0);
  Vector<Complex> Pt = rdft(pt);
  Vector<Complex> Qt = H*Pt;
  Vector<double>  qt = irdft(Qt);
  q.Col(qt,t,0);
  }
  return GImg<double>(q);
  }
*/

/*

void writeComplexSinogramCode(const Matrix<Complex>& gram, const string& fname) {
  fstream out;
  out.open(fname.c_str(),ios::out);
  out << " { ";
  for (int y=0; y<gram.Dimy(); y++) {
    out << " { ";
    for (int x=0; x<gram.Dimx(); x++) {
      out << " { " << gram(x,y) << " }, ";
    }
    out << " }, " << endl;
  }
  out << " } ";
  out.close();
}

void writeSinogramCode(const Matrix<double>& gram, const string& fname) {
  fstream out;
  out.open(fname.c_str(),ios::out);
  out << " { ";
  for (int x=0; x<gram.Dimx(); x++) {
    out << " { ";
    for (int y=0; y<gram.Dimy(); y++) {
      out << gram(x,y) << ", ";
    }
    out << " }, " << endl;
  }
  out << " } ";
  out.close();
}


int main() {
  GImg<double> out(40,40,1);
  Vector<double> x(3);
  Vector<double> y(3);
  
  const double A = 10.0;
  
  const int T = 39;
  const double dt = PIE/(double)T;
  const double R = 10.0;
  for (int i=0; i<T; i++) {
  double t = (double)i*dt;
  double cost = cos(t);
  double sint = sin(t);
  double sint2 = sint*sint;
  double r = R*(cost+sint);
  double x = r*cost;
  double y = r*sint;
  out(rndnear(x+(int)R),rndnear(y+(int)R))(0)=1.0;
  }
  out.Save("uit.pgm");
  }*/


/*int main() {
  cout << endl << endl;
  
  const int N = 101;
  const int T = (int)ceil((double)(N-1)*M_PI);
  const int R = 2*N-1;
  
  const double dt = M_PI/(double)T;
  const double dim = (double)(N-1)*SQRT2;
  const double off = -.5*dim;
  
  GImg<double> scan(N,N,1);
  
  const double dr = 1.0/(double)(R-1);
  for (int t=0; t<T; t++) {
  const double theta = (double)t*dt;
  for (int r=0; r<R; r++) {
  const double rdr = (double)r*dr;
  const double rho = off + ((double)r*dr)*dim;
  cout << "rdr: " << rdr << " rho: " << rho << " 2xrho: " << 2*rho << endl;
  scan.polarLine(rho,theta);
  }
  cout << endl;
  }
  
  cout << "R: " << R << " dim " << dim << " off: " << off << " last " << scan(N-1,N-1)(0) << endl;
  
  //	for (int x=49; x<51; x++) 
  //		for (int y=49; y<51; y++) 
  //			scan(x,y)=0;
  
  scan.Save("scan.pgm");
  
  cout << endl << endl;
  return 0;
  }*/



/*int main() {
  
  void binsplot(int M) {
  const double A = (double)M;
  const double R = 2*A-1;
  const double dp = 1.0/SQRT2;
  const double pmax = A*SQRT2*.5;
  const double mmin = (double)((M-1)/2);
  
  
  const double theta = .75*PIE;
  const double cost = cos(theta);
  const double sint = sin(theta);
  
  for (int j=0; j<M; j++) {
  for (int i=0; i<M; i++) {
  const double x = (double)i-mmin;
  const double y = (double)j-mmin;
  const double p = x*cost+y*sint;
  const int ri = (int)round(p/dp+A-1);
  cout << ri << " ";
  }
  cout << endl;
  
  }
  }
  
  int main() {
  cout << endl << endl;
  
  const int N = 256;
  
  Vector<double> f = signal(N,N,8,16,32);
  Vector<Complex> z = pad2Complex(f);
  Vector<Complex> F = cdft(z,true);
  F.gnuPlot("F.dat");
  
  
  double* fr = new double[N];
  double* fi = new double[N];
  double* Fr = new double[N];
  double* Fi = new double[N];
  
  for (int i=0; i<f.Dim(); i++) {
  fr[i]=f(i);
  fi[i]=0.0;
  Fr[i]=0.0;
  Fi[i]=0.0;
  }
  
  DFT(N,false,fr,fi,Fr,Fi);
  
  fstream out;
  out.open("Dude.dat", ios::out);
  if (out) {
  for (int i=0; i<f.Dim(); i++) 
  out << i << " " << Fr[i] << " " << Fi[i] << endl;
  out.close();
  }
  
  Vector<Complex> Zz(N);
  
  for (int i=0; i<f.Dim(); i++) {
  Zz(i)=Complex::Rectangular(Fr[i],Fi[i]);
  }
  
  Zz.gnuPlot("Zz.dat");
  
  cout << endl << endl;Vector<Complex> backFilter(const int N, double off =-1.0) {
  Vector<Complex> h(N);
  if (off==-1.0) off = (double)(N/-2);
  
  const double nom = -4.0/(PIE*PIE);
  for (int m=0; m<N; m++) {
  const double k = off+(double)m;
  if (k==0.0) h(m)=Complex(1.0);
  else if ((int)k%2 !=0) h(m)=Complex(nom/(k*k));
  else h(m)=Complex(0.0);
  cout << "m: " << m << " k " << k << " h " << h(m) << endl;
  //h(m) = k==0.0 ? 1.0 : (int)k%2 != 0 ? nom/(k*k) : 0.0;
  }	
  return h;
  }
  return 0;
  }
  
  
  cout << endl << endl;
  
  GImg<double> in("brain.pgm");
  GImg<double> uit = rugprojectie(in);
  uit.normalize(0.0,512.0);
  uit.Resize(uit.Dimx(),uit.Dimy()-1);
  GImg<double> g = idht(uit.Data());
  g.normalize(0.0,512);
  g.Save("G.pgm");
  uit.Save("uit.pgm");
  
  cout << endl << endl;
  return 0;
  }*/

/*


Matrix<double> interpol(const Matrix<double>& f) {
  //	fstream out;
  //	out.open("check2.dat",ios::out);
  const int M = f.Dimx();
  const int T = (int)ceil((double)(M-1.0)*PIE);
  const int R = M*2-1;
  
  const double A = (double)M;
  const double rhomax = .05*SQRT2*A;
  const double dr = (2.0*rhomax)/(double)(R);
  const double dt = PIE/(double)T;
  const double mnmin = -.5*(A-1);
  
  
  cout << endl << "config: " << endl;
  cout << "M     : " << M<< endl
       << "T     : " << T << endl
       << "R     : " << R << endl
       << "dt    : " << dt << endl
       << "dr    : " << dr << endl
       << "mmin : " << mnmin << endl;
  cout << endl << endl;
  
  Vector<Complex> h = backramp(R+1,0.0);
  Vector<double> window = hammingWindow(R);
  //	for (int i=0; i<R; i++) h(i)*=window(i);
  //	h.Resize(R+1);
  Vector<Complex> H = cdft(h,false);
  
  h.gnuPlot("lh.dat");
  H.gnuPlot("lH.dat");
  
  Matrix<double> p(T,R+1,1);
  Matrix<double> pfiltered(T,R+1,1);
  Matrix<Complex> Q(T,R+1,1);
  Matrix<double> q(M,M,1);
  
  int i=0; 
  
  for (int t=0; t<T; t++) {
    const double theta = (double)t*dt;
    const double cost = cos(theta);
    const double sint = sin(theta);
    const double rmin = mnmin*(cost+sint);	
    Vector<Complex> pt(R+1);	
    if (sint > dr) {
      const double a = -cost/sint;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rhomax;
	const double b = (rho-rmin)/sint;
	
	for (int m=0; m<M; m++) {
	  int n = (int)round((double)m*a+b);
	  if (n>=0 && n<M) {//check if not interpolating outside the image
	    pt(r) += f(m,n)(0);
	    p(t,r)(0) += f(m,n)(0);
	  }
	}
	pt(r)*=Complex((1.0/sint));
      }
    }
    else {
      const double a = -sint/cost;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rhomax;
	const double b = (rho-rmin)/cost;
	
	for (int n=0; n<M; n++) {
	  int m = (int)round((double)n*a+b);
	  if (m>=0 && m<M) {
	    pt(r) += f(m,n)(0);
	    p(t,r)(0) += f(m,n)(0);
	  }
	}
	pt(r)*=Complex(1.0/fabs(cost));
      }
      
    }
    Vector<Complex> Pt = cdft(pt,false);
    Vector<Complex> Qt=Pt*Complex(1.0);
    Q.Col(Qt,t,0);
    Vector<Complex> qt = cdft(Qt,true);
    pfiltered.Col(vector2Real(qt),t,0);
    
    if (sint > dr) {
      const double a = -cost/sint;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rhomax;
	const double b = (rho-rmin)/sint;
	
	for (int m=0; m<M; m++) {
	  int n = (int)round((double)m*a+b);
	  if (n>=0 && n<M) q(m,n)(0) += qt(r).Real();
	}
      }
    }
    else {
      const double a = -sint/cost;
      for (int r=0; r<R; r++) {
	const double rho = (double)r*dr-rhomax;
	const double b = (rho-rmin)/cost;
	
	for (int n=0; n<M; n++) {
	  int m = (int)round((double)n*a+b);
	  if (m>=0 && m<M) q(m,n) += qt(r).Real();
	}
      }
      
    }
    
  }
  //out.close();
  
  GImg<double> b(logTransform(Q));
  GImg<double> pi(p);
  pi.normalize(0.0,255.0);
  pi.Save(rawname);
  GImg<double> pk(pfiltered);
  pk.normalize(0.0,255.0);
  pk.Save(filtname);
  
  b.Save(qname);
  return q;
}


*/
