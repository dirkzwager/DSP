#include <iostream>
#include <fstream>
#include <string>
#include "VectorMatrix.h"
#include "DSP.h"

#define PI M_PI
#define SQRT2 M_SQRT2

using namespace std;

std::string isPPM(const std::string& text) { return (text=="ppm" ? "P3" : (text=="P3") ? "ppm" : "udef" ); }
std::string isPGM(const std::string& text) { return (text=="pgm" ? "P2" : (text=="P2") ? "pgm" : isPPM(text) ); }

template <class T> class GImg;
template <class T> std::ostream& operator<<(std::ostream&, const GImg<T>&);

template <class T>
class GImg {
protected:
  void readFileHeader(std::ifstream& is, std::string& ftype, int& dimx, int& dimy, int& dimz);
  void writeFileHeader(std::ofstream& os, std::string ftype);
  void writeFileData(std::ofstream& os);
  Matrix<T> data;	
public:
  GImg(const unsigned int dimx, const unsigned int dimy, const unsigned int dimz);
  GImg(const Matrix<T>& m);
  GImg(const GImg& img);
  GImg(const std::string& fname);
  ~GImg();
  
  GImg& Resize(const int dmx, const int dmy);
  
  GImg& operator=(const GImg& img);
  Vector<T>& operator () (const unsigned int x, const unsigned int y);
  const Vector<T>& operator () (const unsigned int x, const unsigned int y) const; 
  friend std::ostream& operator<< <T>(std::ostream&, const GImg<T>&);
  GImg  operator  + (const T x);
  GImg  operator  + (const GImg& m);
  GImg& operator += (const T x);
  GImg& operator += (const GImg& m);
  GImg  operator  - (const T x);
  GImg  operator  - (const GImg& m);
  GImg& operator -= (const T x);
  GImg& operator -= (const GImg& m);
  GImg  operator *  (const T x);
  GImg  operator *  (const GImg& m);  
  GImg& operator *= (const T x);
  GImg& operator *= (const GImg& m); 
  
  GImg& normalize(const T min, const T max);
  //	GImg convolve(const Mask<T>& h, const unsigned int levels=256);
  GImg scale(const int dmx, const int dmy);
  Vector<T> average();
  Vector<T> standardDeviation();	
  
  const int Dimx() const { return data.Dimx(); }
  const int Dimy() const { return data.Dimy(); }
  const int Dimz() const { return data.Dimz(); }
  const Vector<T> maxVector() const;
  const Vector<T> minVector() const;
  const T maxValue() const;
  const T minValue() const;
  Matrix<T> Data()  { Matrix<T> _data(data); return _data; }
  const Matrix<T> Data() const { Matrix<T> _data(data); return _data; }
  
  void Save(const std::string& fname);
  void gnuPlot(const std::string& fname); 
  void polarLine(const double rho, const double theta);
};

template <class T>
inline GImg<T>::GImg(const unsigned int dimx, const unsigned int dimy, const unsigned int dimz) : data(dimx,dimy,dimz) { }

template <class T>
GImg<T>::GImg(const Matrix<T>& m) : data(m) { }

template <class T>
GImg<T>::GImg(const std::string& fname) {
  ifstream is;
  is.open(fname.c_str(),std::ios::in);
  
  int dmx, dmy, dmz;
  std::string extstr;
  if (is.is_open()) readFileHeader(is,extstr,dmx,dmy,dmz);
  else std::cout << "Error, unable to locate file" << std::endl;
  if (extstr=="udef") std::cout << "Error, undefined file header" << std::endl;
  else {
    data=Matrix<T>(dmx,dmy,dmz); 
    
    for (unsigned int y=0; y<dmy; y++)
      for (unsigned int x=0; x<dmx; x++) 
	for (unsigned int z=0; z<dmz; z++) 
	  is >> data(x,y)(z);
  }
  is.close();	
}

template <class T>
GImg<T>::GImg(const GImg& img) { data=img.data; }

template <class T>
GImg<T>::~GImg() { }

template <class T>
GImg<T>& GImg<T>::Resize(const int dmx, const int dmy) {
  data.Resize(dmx,dmy);
  return *this;
}

template <class T>
void GImg<T>::readFileHeader(std::ifstream& is, std::string& extstr, int& dimx, int& dimy, int& dimz) {
  std::string tmp;
  is >> tmp;
  extstr=isPGM(tmp);
  if (extstr=="pgm") {
    getline(is,tmp); 		
    getline(is,tmp); 
    is >> dimx >> dimy;
    getline(is,tmp); is >> tmp;		
    dimz=1;
  }
  else if (extstr=="ppm") {
    getline(is,tmp); 		
    getline(is,tmp); 
    is >> dimx >> dimy;
    getline(is,tmp); is >> tmp;			
    dimz=3;
  }	
}

template <class T>
void GImg<T>::writeFileHeader(std::ofstream& os, std::string extstr) {
  os << isPGM(extstr) << std::endl
     << "#" << std::endl 
     << Dimx() << " " 
     << Dimy() << std::endl 
     << maxValue() << std::endl;
}

template <class T>
void GImg<T>::writeFileData(std::ofstream& os) {
  for (int y=0; y<Dimy(); y++) {
    for (int x=0; x<Dimx(); x++)  {
      //		if (x>70) os << endl;
      os << data(x,y) << endl;
    }
    //		cout << endl;
  }
}

template <class T>
void GImg<T>::Save(const std::string& fname) {
  const unsigned int extpos = fname.rfind(".");
  const std::string  extstr = fname.substr(extpos+1,extpos+3);
  
  ofstream os;
  os.open(fname.c_str(),std::ios::out);
  
  if (isPGM(extstr) == "udef") {
    std::cout << "Error, incorrect file description" << std::endl;
  }
  else {
    os.setf(std::ios::fixed);
    os.precision(0);
    writeFileHeader(os,extstr);
    writeFileData(os);
    std::cout << "Saved: " << fname << std::endl;
  }		
  os.close();
}

template <class T>
GImg<T>& GImg<T>::normalize(const T min, const T max) {
  const T imgmin = minValue();
  const T imgmax = maxValue();
  const T scale  = (max-min)/(imgmax-imgmin);
  data-=imgmin;
  data*=scale;
  return *this; 
}

/*template <class T>
  GImg<T> GImg<T>::convolve(const Mask<T>& h, const unsigned int levels) {
  GImg<T> g(*this);
  for (int x=1; x<Dimx()-1; x++)
  for (int y=1; y<Dimy()-1; y++) 
  for (int m=0; m<3; m++)
  for (int n=0; n<3; n++) 
  g.data(x,y)+=h(m,n)*data(x+m-1,y+n-1);
  g.normalize(T(0),(double)levels);
  return g;	
  }*/

template <class T>
GImg<T> GImg<T>::scale(const int dmx, const int dmy) {
  const double sx = dmx/Dimx();
  const double sy = dmx/Dimy();
  const double dx = 1.0/sx;
  const double dy = 1.0/sy;
  
  
}

template <class T>
Vector<T>& GImg<T>::operator () (const unsigned int x, const unsigned int y) {
  if (x>=Dimx() || y>=Dimy()) std::cout << "Error, unbound array subscript" << std::endl;
  else return data(x,y);
}

template <class T>
const Vector<T>& GImg<T>::operator () (const unsigned int x, const unsigned int y) const {
  if (x>=Dimx() || y>=Dimy()) std::cout << "Error, unbound array subscript" << std::endl;
  else return data(x,y);
}

template <class T>
std::ostream& operator<<(std::ostream& os, const GImg<T>& img) {
  for (int y=0; y<img.Dimy(); y++) {
    for (int x=0; x<img.Dimx(); x++)
      os << "| " << img.data(x,y) << "| ";
    os << std::endl;  
  }
  return os;
}

template <class T>
GImg<T>& GImg<T>::operator=(const GImg &img) {
  if(this==&img) return *this;
  data=img.data;
  return *this;
}

template <class T>
GImg<T> GImg<T>::operator  + (const T x) {
  return GImg(data+x);
}

template <class T>
GImg<T> GImg<T>::operator  + (const GImg& m) {
  if (!(Dimx()==m.Dimx() && Dimy()==m.Dimy() && Dimz()==m.Dimz())) std::cout << "Error, incompatible dimensions" << std::endl;
  else return GImg(data+m.data);
}

template <class T>
GImg<T>& GImg<T>::operator += (const T x) {
  data+=x;
  return *this;
}

template <class T>
GImg<T>& GImg<T>::operator += (const GImg& m){
  if (!(Dimx()==m.Dimx() && Dimy()==m.Dimy() && Dimz()==m.Dimz())) std::cout << "Error, incompatible dimensions" << std::endl;
  else {
    data+=m.data;
    return *this;
  }
}

template <class T>
GImg<T> GImg<T>::operator  - (const T x) {
  return GImg(data-x);
}

template <class T>
GImg<T> GImg<T>::operator - (const GImg& m) {
  if (!(Dimx()==m.Dimx() && Dimy()==m.Dimy() && Dimz()==m.Dimz())) std::cout << "Error, incompatible dimensions" << std::endl;
  else return GImg(data-m.data);
}

template <class T>
GImg<T>& GImg<T>::operator -= (const T x) {
  data-=x;
  return *this;
}

template <class T>
GImg<T>& GImg<T>::operator -= (const GImg& m){
  if (!(Dimx()==m.Dimx() && Dimy()==m.Dimy() && Dimz()==m.Dimz())) std::cout << "Error, incompatible dimensions" << std::endl;
  else {
    data-=m.data;
    return *this;
  }
}

template <class T>
GImg<T> GImg<T>::operator * (const T x) {
  return GImg(data*x);
}

template <class T>
GImg<T> GImg<T>::operator * (const GImg& m) {
  if (!(Dimx()==m.Dimx() && Dimy()==m.Dimy() && Dimz()==m.Dimz())) std::cout << "Error, incompatible dimensions" << std::endl;
  else return GImg(data*m.data);
}

template <class T>
GImg<T>& GImg<T>::operator *= (const T x) {
  data*=x;
  return *this;
}

template <class T>
GImg<T>& GImg<T>::operator *= (const GImg& m){
  if (!(Dimx()==m.Dimx() && Dimy()==m.Dimy() && Dimz()==m.Dimz())) std::cout << "Error, incompatible dimensions" << std::endl;
  else {
    data*=m.data;
    return *this;
  }
}

template <class T>
inline const Vector<T> GImg<T>::maxVector() const {
  Vector<T> max=data(0,0);
  for (int x=0; x<Dimx(); x++) 
    for (int y=0; y<Dimy(); y++) 
      max = data(x,y) > max ? data(x,y) : max;
  return max;
}

template <class T>
inline const Vector<T> GImg<T>::minVector() const {
  Vector<T> min=data(0,0);
  for (int x=0; x<Dimx(); x++) 
    for (int y=0; y<Dimy(); y++) 
      min = data(x,y) < min ? data(x,y) : min;
  return min;
}

template <class T>
inline const T GImg<T>::maxValue() const {
  T max=data(0,0).Max();
  for (int x=0; x<Dimx(); x++) 
    for (int y=0; y<Dimy(); y++) 
      max = data(x,y).Max() > max ? data(x,y).Max() : max;
  return max;
}

template <class T>
inline const T GImg<T>::minValue() const {
  T min=data(0,0).Min();
  for (int x=0; x<Dimx(); x++) 
    for (int y=0; y<Dimy(); y++) 
      min = data(x,y).Min() < min ? data(x,y).Min() : min;
  return min;
}

template <class T>
void GImg<T>::gnuPlot(const std::string& fname) {
  data.gnuPlot(fname);
}

template <class T>
void GImg<T>::polarLine(const double rho, const double theta) {
  const double cost = cos(theta);
  const double sint = sin(theta);
  
  if (sint > 1.0/M_SQRT2) {//xas
    const double slope = -cost/sint;
    const double interc = rho/sint;
    const double max = .5*(double)(data.Dimx()-1);
    //		cout << "r: " << rho << " theta: " << theta << " slope " << slope << " interc " << interc << endl; 
    
    for (int m=0; m<Dimx(); m++) {
      const double x = (double)m-max;
      const double y = slope*x+interc;
      const int yidx = (int)round(y+max);
      if (yidx>=0 && yidx<Dimx()) data(m,yidx)(0)+=(1.0/sint);				
    }
  }
  else {
    const double slope = -sint/cost;
    const double interc = rho/cost;
    const double max = .5*(double)(data.Dimy()-1);
    
    //		cout << "r: " << rho << " theta: " << theta << " slope " << slope << " interc " << interc << endl;
    
    for (int n=0; n<Dimy(); n++) {
      const double y = (double)n-max;
      const double x = slope*y+interc;
      const int xidx = (int)round(x+max);
      if (xidx>=0 && xidx<Dimy()) data(xidx,n)(0)+=(1.0/fabs(cost));				
    }
    
  }
}



template <class T>
class Mask {
private:
	Matrix<T> M;
	Mask(const Matrix& _M) : M(_M) { }
public:
	~Mask();
	static Mask Sobel();
	static Mask hSobel();
	static Mask vSobel();
	static Mask Laplacian();
	static Mask Kirsch();
	static Mask vKirsch();
	static Mask hKirsch();
	static Mask d1Kirsch();
	static Mask d2Kirsch();
	static Mask FreiChen();
	static Mask Color();

	
	const T operator () (const unsigned int x, const unsigned int y, const unsigned int z) const { return M(x,y)(z); };
	const Vector<T>& operator () (const unsigned int x, const unsigned int y) const { return M(x,y); }
	Vector<T>& operator () (const unsigned int x, const unsigned int y) { return M(x,y); }
};

Mask<T>::~Mask() {}

//const double Mask::operator () (const unsigned int x, const unsigned int y, const unsigned int z) const {
//	return M(x,y)(z);
//}

//Vector& Mask::operator () (const unsigned int x, const unsigned int y) {
//	return M(x,y);
//}

template <class T>
inline Mask<T> Mask<T>::Sobel() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=-2.0; _M(1,0)(0)=-2.0; _M(2,0)(0)=0.0;
	_M(0,1)(0)=-2.0; _M(1,1)(0)=0.0; _M(2,1)(0)=2.0;
	_M(0,2)(0)= 0.0; _M(1,2)(0)=2.0; _M(2,2)(0)=2.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::hSobel() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=-1.0; _M(1,0)(0)=-2.0; _M(2,0)(0)=-1.0;
	_M(0,1)(0)= 0.0; _M(1,1)(0)= 0.0; _M(2,1)(0)= 0.0;
	_M(0,2)(0)= 1.0; _M(1,2)(0)= 2.0; _M(2,2)(0)= 1.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::vSobel() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=-1.0; _M(1,0)(0)=0.0; _M(2,0)(0)=1.0;
	_M(0,1)(0)=-2.0; _M(1,1)(0)=0.0; _M(2,1)(0)=2.0;
	_M(0,2)(0)=-1.0; _M(1,2)(0)=0.0; _M(2,2)(0)=1.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::Laplacian() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=0.0; _M(1,0)(0)=1.0; _M(2,0)(0)=0.0;
	_M(0,1)(0)=1.0; _M(1,1)(0)=-4.0; _M(2,1)(0)=1.0;
	_M(0,2)(0)=0.0; _M(1,2)(0)=1.0; _M(2,2)(0)=0.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::Kirsch() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)= 1.0; _M(1,0)(0)= 3.0; _M(2,0)(0)= 3.0;
	_M(0,1)(0)=-1.0; _M(1,1)(0)= 0.0; _M(2,1)(0)= 1.0;
	_M(0,2)(0)=-3.0; _M(1,2)(0)=-3.0; _M(2,2)(0)=-1.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::vKirsch() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=-1.0; _M(1,0)(0)=0.0; _M(2,0)(0)=1.0;
	_M(0,1)(0)=-1.0; _M(1,1)(0)=0.0; _M(2,1)(0)=1.0;
	_M(0,2)(0)=-1.0; _M(1,2)(0)=0.0; _M(2,2)(0)=1.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::hKirsch() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)= 1.0; _M(1,0)(0)= 1.0; _M(2,0)(0)= 1.0;
	_M(0,1)(0)= 0.0; _M(1,1)(0)= 0.0; _M(2,1)(0)= 0.0;
	_M(0,2)(0)=-1.0; _M(1,2)(0)=-1.0; _M(2,2)(0)=-1.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::d1Kirsch() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=-1.0; _M(1,0)(0)=-1.0; _M(2,0)(0)=0.0;
	_M(0,1)(0)=-1.0; _M(1,1)(0)= 0.0; _M(2,1)(0)=1.0;
	_M(0,2)(0)= 0.0; _M(1,2)(0)= 1.0; _M(2,2)(0)=1.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::d2Kirsch() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=0.0; _M(1,0)(0)=-1.0; _M(2,0)(0)=-1.0;
	_M(0,1)(0)=1.0; _M(1,1)(0)= 0.0; _M(2,1)(0)=-1.0;
	_M(0,2)(0)=1.0; _M(1,2)(0)= 1.0; _M(2,2)(0)= 0.0;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::FreiChen() {
	Matrix<T> _M(3,3,1);
	_M(0,0)(0)=1+SQRT2; _M(1,0)(0)=-1.0+SQRT2; _M(2,0)(0)=1.0+SQRT2;
	_M(0,1)(0)=-1.0+SQRT2; _M(1,1)(0)= 9.0; _M(2,1)(0)=-1.0-SQRT2;
	_M(0,2)(0)=1.0-SQRT2; _M(1,2)(0)= 3.0-SQRT2; _M(2,2)(0)= -3-SQRT2;
	return Mask(_M);
}

template <class T>
inline Mask<T> Mask<T>::Color() {
	Matrix<T> _M(3,3,3);
	_M(0,0)(0)=1; _M(1,0)(0)=0; _M(2,0)(0)=-1;
	_M(0,1)(0)=1; _M(1,1)(0)=0; _M(2,1)(0)=-1;
	_M(0,2)(0)=1; _M(1,2)(0)=0; _M(2,2)(0)=-1;
	
	_M(0,0)(1)=1; _M(1,0)(1)=1; _M(2,0)(1)= 1;
	_M(0,1)(1)=0; _M(1,1)(1)=0; _M(2,1)(1)= 0;
	_M(0,2)(1)=-1; _M(1,2)(1)=-1; _M(2,2)(1)=-1;
	
	_M(0,0)(2)=0; _M(1,0)(2)=-1; _M(2,0)(2)=-1;
	_M(0,1)(2)=1; _M(1,1)(2)=0; _M(2,1)(2)=-1;
	_M(0,2)(2)=1; _M(1,2)(2)=1; _M(2,2)(2)=0;

	return Mask(_M);
}

/*
struct ColorMap {
	unsigned int levels;
	unsigned int max;	
	Matrix map;
	inline ColorMap(const unsigned int _levels, const unsigned int _max, 
		     const double pr, const double tr, 
		     const double pg, const double tg, 
		     const double pb, const double tb) 
		     : levels(_levels), max(_max), map(_levels,1,3) {
		const double spr = (PI*pr)/(double)levels;
		const double spg = (PI*pg)/(double)levels;
		const double spb = (PI*pb)/(double)levels;
		const double str = (tr*PI)/(double)levels;
		const double stg = (tr*PI)/(double)levels;
		const double stb = (tr*PI)/(double)levels;
	
		for (int i=0; i<levels; i++) {
			map(i,0)(0) = fabs(sin((double)i*spr+str))*(double)_max;
			map(i,0)(1) = fabs(sin((double)i*spg+stg))*(double)_max;
			map(i,0)(2) = fabs(sin((double)i*spb+stb))*(double)_max;
		}
	}
	
	inline ColorMap(const unsigned int _levels, const unsigned int _max)
		     : levels(_levels), max(_max), map(_levels,1,3) {
		const double sr=.222, sg=.707, sb=.0071;
	
		for (int i=0; i<levels; i++) {
			map(i,0)(0) = (sr*(double)i)*(double)_max;
			map(i,0)(1) = (sg*(double)i)*(double)_max;
			map(i,0)(2) = (sb*(double)i)*(double)_max;
		}
	}
	
	const unsigned int operator () (const unsigned int entry, const unsigned int rgb) const {
		if (entry >= levels || rgb >= 3) std::cout << "Error, incorrect entry value 1" << std::endl;
		else {
			const double val = map(entry,0)(rgb);
			const double flrval = floor(val);

			return (unsigned int)(val-flrval < .5 ? flrval : ceil(val));
		}
	}
	
	void operator () (const unsigned int entry, double& R, double& G, double& B) const {
		if (entry >= levels) std::cout << "Error, incorrect entry value 2" << std::endl;
		R = map(entry,0)(0);
		G = map(entry,0)(1);
		B = map(entry,0)(2);
	}
	
	friend std::ostream& operator<< (std::ostream& os, const ColorMap& cm) {
		for(int i=0; i<cm.levels; i++) {
			os  << i << " " << cm(i,0)
			   << " " << cm(i,1)
			   << " " << cm(i,2) << endl; 
		}
		return os;
	}
};*/
