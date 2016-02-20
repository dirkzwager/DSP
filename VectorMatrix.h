#ifndef _VECTORMATRIX_H_
#define _VECTORMATRIX_H_

#include <iostream>
#include <cassert>
#include <math.h>

using namespace std;

template <class T> class Vector;
template <class T>  std::ostream& operator<<(std::ostream&, const Vector<T>&);

template <class T> 
class Vector {
private:
	int dim;
	T* data;
public:
	Vector(int _dim=0);		
	Vector(const T* _data, int _dim);
	Vector(const Vector&);			
	~Vector();				
	void init(const T* _data, int _dim);
	Vector& Resize(const int _dim);		
	inline int Dim() const { return dim; }
	T& operator() (int idx);
	const T& operator()(int idx) const;
	Vector& operator=(const Vector<T>&);
	const int operator == (const Vector& v) const;
	const int operator != (const Vector& v) const;
	const int operator <  (const Vector& v) const;
	const int operator <= (const Vector& v) const;
	const int operator >  (const Vector& v) const;
	const int operator >= (const Vector& v) const;
	Vector  operator  + (const T x) 		const;
  	Vector  operator  + (const Vector& v) 	const;
  	Vector& operator += (const T x);
  	Vector& operator += (const Vector& v);
  	Vector  operator  - (const T x) 		const;
  	Vector  operator  - (const Vector& v) 	const;
  	Vector& operator -= (const T x);
  	Vector& operator -= (const Vector& v);
  	Vector  operator *  (const T x) 		const;
  	Vector  operator *  (const Vector& z) 	const;  
  	Vector& operator *= (const T x);
  	Vector& operator *= (const Vector& v);
  	Vector  operator /  (const T x) 		const;
  	Vector  operator /  (const Vector& z)	const;  
  	Vector& operator /= (const T x);
  	Vector& operator /= (const Vector& v);
  	friend ostream& operator << <T>(std::ostream&, const Vector<T>&);
  	
  	T vectorLength() const;
  	T Max() const;
  	T Min() const;
  	Vector& unitVector();
  	void gnuPlot(const std::string& fname) const;
};

template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
  for (int i=0; i<v.dim; i++) {
  	os << v.data[i] << " ";
//  	os << i << " " << v.data[i] << endl;
  }
  return os;
}

template <class T>
void Vector<T>::init(const T *_data, int _dim) {
	dim=_dim;
	data=new T[dim];
	assert(data!=0);
	for(int i=0; i<dim; i++) data[i]= (_data!=0) ? _data[i] : T(0);
}

template <class T>
Vector<T>::Vector(int _dim) { init(0,_dim); }

template <class T>
Vector<T>::Vector(const Vector& v) { init(v.data,v.dim); }

template <class T>
Vector<T>::Vector(const T* _data, int _dim) { init(data,_dim); }

template <class T>
Vector<T>::~Vector() { delete [] data; }

template <class T>
Vector<T>& Vector<T>::Resize(const int _dim) {
	const int mindim = _dim < dim ? _dim : dim;
	T* cdata = new T[dim];
	for (int i=0; i<dim; i++) cdata[i]=data[i];
	delete[] data; data = NULL; 
	data=new T[_dim];
	
	for (int k=0; k<mindim; k++) data[k]=cdata[k];	
	for (int k=mindim; k<_dim; k++) data[k]=T(0);
	dim=_dim;
	return *this;	
}

template <class T>
T& Vector<T>::operator () (int index) {
	assert(index>=0 && index<dim);
	return data[index];
}

template <class T>
const T& Vector<T>::operator () (int index) const {
	assert(index>=0 && index<dim);
	return data[index];
}

template <class T>
Vector<T>& Vector<T>::operator=(const Vector &v) {
	if(this==&v) return *this;//*this==v??
	delete [] data;
	init(v.data,v.dim);
	return *this;
}

template <class T>
const int Vector<T>::operator == (const Vector& v) const {
	unsigned int i=0;
	while (i<dim) {
		if (data[i]!=v(i)) return 0;
		else i++;
	}
	return 1;
}

template <class T>
const int Vector<T>::operator != (const Vector& v) const {
	return !(*this==v);
}

template <class T>
const int Vector<T>::operator < (const Vector& v) const {
	return vectorLength() < v.vectorLength();
}

template <class T>
const int Vector<T>::operator <= (const Vector& v) const {
	return vectorLength() <= v.vectorLength();
}

template <class T>
const int Vector<T>::operator > (const Vector& v) const {
	return vectorLength() > v.vectorLength();
}

template <class T>
const int Vector<T>::operator >= (const Vector& v) const {
	return vectorLength() >= v.vectorLength();
}

template <class T>
Vector<T> Vector<T>::operator + (const T x) const {
  Vector res(*this);
  return res+=x;
}

template <class T>
Vector<T>  Vector<T>::operator + (const Vector& v) const {
	if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  	else {
	  	Vector res(*this); 
	  	return res+=v;
  	}
}

template <class T>
Vector<T>  Vector<T>::operator - (const T x) const {
  Vector res(*this);
  return res-=x;
}

template <class T>
Vector<T>  Vector<T>::operator - (const Vector& v) const {
	if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  	else {
	  	Vector res(*this); 
	  	return res-=v;
  	}
}

template <class T>
Vector<T>& Vector<T>::operator += (const T x){
  for (int i=0; i<dim; i++) data[i]+=x;
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator += (const Vector& v){
  if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  else {
	  for (int i=0; i<dim; i++) data[i]+=v(i);
	  return *this;
  }
}

template <class T>
Vector<T>& Vector<T>::operator -= (const T x){
  for (int i=0; i<dim; i++) data[i]-=x;
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator -= (const Vector& v){
  if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  else {
	  for (int i=0; i<dim; i++) data[i]-=v(i);
	  return *this;
  }
}

template <class T>
Vector<T> Vector<T>::operator * (const T x) const {
  Vector res(dim);
  for (int i=0; i<dim; i++) res(i)=data[i]*x;
  return res;
}

template <class T>
Vector<T> Vector<T>::operator * (const Vector& v) const {
  if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
//  else {
	  Vector res(dim);
	  for (int i=0; i<dim; i++) res(i)=data[i]*v(i);
	  return res;
//  }
}

template <class T>
Vector<T>& Vector<T>::operator *= (const T x){
  for (int i=0; i<dim; i++) data[i]*=x;
  return *this;
}

template <class T>
Vector<T>& Vector<T>::operator *= (const Vector& v){
  if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  else {
	  for (int i=0; i<dim; i++) data[i]*=v(i);
	  return *this;
  }
}

template <class T>
Vector<T> Vector<T>::operator / (const T x) const{
	if (x!=0.0) {
 	 	Vector res(dim);
  		for (int i=0; i<dim; i++) res(i)=data[i]/x;
  		return res;
  	}
  	else std::cout << "Error, division by 0!\n" << std::endl;
}

template <class T>
Vector<T> Vector<T>::operator / (const Vector& v) const{
  if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  else {
	  Vector res(dim);
	  for (int i=0; i<dim; i++) {
	  	if (v(i) != 0.0) res(i)=data[i]/v(i);
	  	else std::cout << "Error, division by 0!\n" << std::endl;
	  }
	  return res;
  }
}

template <class T>
Vector<T>& Vector<T>::operator /= (const T x) {
	if (x!=0.0) {
  		for (int i=0; i<dim; i++) data[i]/=x;
  		return *this;
  	}
  	else std::cout << "Error, division by 0!\n" << std::endl;
}

template <class T>
Vector<T>& Vector<T>::operator /= (const Vector& v){
  if (dim != v.dim) std::cout << "Error, incompatible array dimensions!\n" << std::endl;
  else {
	  for (int i=0; i<dim; i++) {
	  	if (v(i) != T(0)) data[i]/=v(i);
	  	else std::cout << "Error, division by 0!\n" << std::endl;
	  }
	  return *this;
  }
}

template <class T>
T Vector<T>::vectorLength() const {
	T square=T(0);
	for (int i=0; i<dim; i++) square+=(data[i]*data[i]);
	return sqrt(square);
}


template <class T>
Vector<T>& Vector<T>::unitVector() {
	T length=vectorLength();
	if (length!=T(0)) (*this)/=length;
	return *this;
}

template <class T>
T Vector<T>::Max() const {
	T max=data[0];
	for (int i=1; i<dim; i++) max = data[i] > max ? data[i] : max;
	return max;
}

template <class T>
T Vector<T>::Min() const {
	T min=data[0];
	for (int i=1; i<dim; i++) min = data[i] < min ? data[i] : min;
	return min;
}

template <class T> 
void Vector<T>::gnuPlot(const std::string& fname) const {
	std::fstream out;
	out.open(fname.c_str(),std::ios::out|std::ios::trunc);
	if (out) {
		for (int i=0; i<dim; i++) out << i << " " << data[i] << endl;
		out.close();
	}
	else std::cout << "Error, unable to write file: " << fname << std:: endl;
}

template <class T> class Matrix;
template <class T> std::ostream& operator<<(std::ostream&, const Matrix<T>&);

template <class T> 
class Matrix {
protected:
	Vector<T>* data;
	int dimx,dimy,dimz;
	void init(const Vector<T>* _data, int _dimx, int _dimy, int _dimz);
public:
	Matrix(int _dimx=0, int _dimy=0, int _dimz=1);		
	Matrix(const Vector<T>* _data, int _dimx, int _dimy, int _dimz);
	Matrix(const Matrix&);			
	~Matrix();			
	Matrix& Resize(const unsigned int _dimx, const unsigned int _dimy);						
	Matrix& Rescale(const unsigned int _dimx, const unsigned int _dimy);
	inline int Dimx() {return dimx;}
	inline int Dimy() {return dimy;}
	inline int Dimz() {return dimz;}
	
	Vector<T>& operator()(int, int);
	const Vector<T>& operator()(int, int) 		const;
	Matrix& operator=(const Matrix&);
  	Matrix  operator  + (const T x) 			const;
  	Matrix  operator  + (const Matrix& m) 		const;
  	Matrix& operator += (const T x);
  	Matrix& operator += (const Matrix& m);
  	Matrix  operator  - (const T x) 			const;
  	Matrix  operator  - (const Matrix& m)		const;
  	Matrix& operator -= (const T x);
  	Matrix& operator -= (const Matrix& m);
  	Matrix  operator *  (const T x)				const;
  	Matrix  operator *  (const Matrix& m)		const;  
  	Matrix& operator *= (const T x);
  	Matrix& operator *= (const Matrix& m);
  	Matrix  operator *  (const Vector<T>& v)	const;  
  	Matrix& operator *= (const Vector<T>& v);
   	Matrix  operator /  (const T x)				const;
  	Matrix  operator /  (const Matrix& m)		const;  
  	Matrix& operator /= (const T x);
  	Matrix& operator /= (const Matrix& m);
  	friend std::ostream& operator<< <T>(std::ostream&, const Matrix&);
  	
  	const int Dimx() const { return dimx; }
  	const int Dimy() const { return dimy; }
  	const int Dimz() const { return dimz; }

  	Vector<T> Row(int y, int z) const;
  	Vector<T> Col(int x, int z) const;
  	Matrix<T>& Row(const Vector<T>& r, int y, int z);
  	Matrix<T>& Col(const Vector<T>& c, int x, int z);
  	void gnuPlot(const std::string& fname);
  	
/*  	Vector<Vector<T> > Row(int idx) const;
  	Vector<Vector<T> > Col(int idx) const;
  	void Row(const Vector<Vector<T> >& r, int y);
  	void Col(const Vector<Vector<T> >& c, int x);*/
  	
  	Matrix& unitVectorMatrix();
};

template <class T> 
std::ostream& operator<<(std::ostream& os,Matrix<T>& v) {
	for(int y=0; y<v.Dimy(); y++) {
	  	for(int x=0; x<v.Dimx(); x++) 
			os << v(x,y) << " ";
		os << std::endl;
	}
	return os;
}

template <class T> 
void Matrix<T>::gnuPlot(const std::string& fname) {
	fstream out;
	out.open(fname.c_str(),std::ios::out|std::ios::trunc);
	if (out) {
		for(int y=0; y<dimy; y++) {
		  	for(int x=0; x<dimx; x++) 
				out << data[y*dimx+x] << " " << std::endl;
			out << std::endl;
		}
		out.close();
	}
	else std::cout << "Error, unable to write file: " << fname << std:: endl;
}

template <class T> 
void Matrix<T>::init(const Vector<T>* _data, int _dimx, int _dimy, int _dimz) {
	dimx=_dimx; dimy=_dimy; dimz=_dimz;
	data=new Vector<T>[dimx*dimy];
	assert(data!=0);
	for(int i=0; i<dimx*dimy; i++) data[i]= _data!=0 ? _data[i] : Vector<T>(dimz);
}

template <class T> 
Matrix<T>::Matrix(int _dimx, int _dimy, int _dimz) { init(0,_dimx,_dimy,_dimz); }

template <class T> 
Matrix<T>::Matrix(const Matrix& v) { init(v.data,v.dimx,v.dimy,v.dimz); }

template <class T> 
Matrix<T>::Matrix(const Vector<T>* _data, int _dimx, int _dimy, int _dimz) { init(data,_dimx,_dimy,_dimz); }

template <class T> 
Matrix<T>::~Matrix() { delete [] data; }

template <class T> 
Matrix<T>& Matrix<T>::Resize(const unsigned int _dimx, const unsigned int _dimy) {
	const unsigned int minx = _dimx < dimx ? _dimx : dimx;
	const unsigned int miny = _dimy < dimy ? _dimy : dimy;
	Vector<T>* cdata = new Vector<T>[dimx*dimy];
	for (int i=0; i<dimx*dimy; i++) cdata[i]=data[i];
	
	delete[] data; data=NULL;
	data = new Vector<T>[_dimx*_dimy];
	
	for (int i=0; i<_dimx*_dimy; i++) data[i]=Vector<T>(dimz);
		
	for (int x=0; x<minx; x++)
		for (int y=0; y<miny; y++)
			data[y*_dimx+x] = cdata[y*dimx+x];
			
	dimx=_dimx;
	dimy=_dimy;
	return *this;
}

template <class T> 
Vector<T>& Matrix<T>::operator () (int x, int y) {
	assert(x>=0 && x<dimx);
	return data[y*dimx+x];
}

template <class T> 
const Vector<T>& Matrix<T>::operator () (int x, int y) const {
	assert(x>=0 && x<dimx);
	return data[y*dimx+x];
}

template <class T> 
Matrix<T>& Matrix<T>::operator=(const Matrix &v) {
	if(this==&v) return *this;//*this==v??
	delete [] data;
	init(v.data,v.dimx,v.dimy,v.dimz);
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator  + (const T x) const {
	Matrix M(*this);
	return M+=x;
}

template <class T> 
Matrix<T>& Matrix<T>::operator += (const T x) {
	for (unsigned int i=0; i<dimx*dimy; i++) data[i]+=x;
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator + (const Matrix& m) const {
	if (!(dimx==m.dimx && dimy==m.dimy && dimz==m.dimz)) std::cout << "Error, incompatible dimensions" << std::endl;
	else {
		Matrix M=*this;
		return M+=m;	
	}
}

template <class T> 
Matrix<T>& Matrix<T>::operator += (const Matrix& m) {
	if (!(dimx==m.dimx && dimy==m.dimy && dimz==m.dimz)) std::cout << "Error, incompatible dimensions" << std::endl;
	else for (int i=0; i<m.dimx*m.dimy; i++) data[i]+=m.data[i];
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator - (const T x) const {
	Matrix M(*this);
	return M-=x;
}

template <class T> 
Matrix<T>& Matrix<T>::operator -= (const T x) {
	for (unsigned int i=0; i<dimx*dimy; i++) data[i]-=x;
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator - (const Matrix& m) const {
	if (!(dimx==m.dimx && dimy==m.dimy && dimz==m.dimz)) std::cout << "Error, incompatible dimensions" << std::endl;
	else {
		Matrix M(*this);
		return M-=m;	
	}
}

template <class T> 
Matrix<T>& Matrix<T>::operator -= (const Matrix& m) {
	if (!(dimx==m.dimx && dimy==m.dimy && dimz==m.dimz)) std::cout << "Error, incompatible dimensions" << std::endl;
	else for (int i=0; i<dimx*dimy; i++) data[i]-=m.data[i];
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator *  (const T x) const {
	Matrix M(*this);
	return M*=x;
}

template <class T> 
Matrix<T>& Matrix<T>::operator *= (const T x) {
	for (unsigned int i=0; i<dimx*dimy; i++) data[i]*=x;
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator * (const Matrix& m) const {
	Matrix M(*this);
	return M*=m;
}

template <class T> 
Matrix<T>& Matrix<T>::operator *= (const Matrix& M) {
	const int thisx = dimx;
	const int thisy = dimy;
	const int Mx = M.dimx;
	
	Matrix C(*this);
	delete[] data;
	init(0,Mx,thisy,dimz);
	for (int i=0; i<thisy; ++i) {
		for (int j=0; j<Mx; ++j) {
			Vector<T> sum(dimz);
			for (int k=0; k<thisx; ++k) 
			 	sum+=C(k,i)*M(j,k);
			data[i*Mx+j]=sum;				 
		}
	}
	return *this;
}

template <class T> 
Matrix<T>  Matrix<T>::operator * (const Vector<T>& v) const {
	if (v.Dim() == dimy) {
		Matrix M(*this);
		for (int x=0; x<dimx; x++) {
			for (int y=0; y<dimy; y++) 
				M.data[y*dimx+x]*=v(y);
		} 
		return M;
	}
	else std::cout << "\nError, vector dimensions different from matrix" << std::endl;
}

template <class T> 
Matrix<T>& Matrix<T>::operator *= (const Vector<T>& v) {
	if (v.Dim() == dimy) {
		for (int x=0; x<dimx; x++) {
			for (int y=0; y<dimy; y++) 
				data[y*dimx+x]*=v(y);
		} 
		return *this;
	}
	else std::cout << "\nError, vector dimensions different from matrix" << std::endl;
}

template <class T> 
Matrix<T>  Matrix<T>::operator / (const T x) const {
	if (x != T(0)) {
		Matrix M(*this);
		return M*(T(1.0)/x);
	}
	else {
	  	std::cout << "Error, division by 0!\n" << std::endl;
	}
}

template <class T> 
Matrix<T>& Matrix<T>::operator /= (const T x) {
	if (x != T(0)) {
		for (int i=0; i<Dimx()*Dimy(); i++) data[i]/=x;
		return *this;
	}
	else {
	  	std::cout << "Error, division by 0!\n" << std::endl;
	}
}

template <class T> 
Matrix<T>  Matrix<T>::operator / (const Matrix& m) const {
	Matrix M(*this);
	return M/=m;
}

template <class T> 
Matrix<T>& Matrix<T>::operator /= (const Matrix& M) {
	const int thisx = dimx;
	const int thisy = dimy;
	const int Mx = M.dimx;
	
	Matrix C(*this);
	delete[] data;
	init(0,Mx,thisy,dimz);
	for (int i=0; i<thisy; ++i) {
		for (int j=0; j<Mx; ++j) {
			Vector<T> sum(dimz);
			for (int k=0; k<thisx; ++k) {
				sum+=C(k,i)/M(j,k);
			}
			data[i*Mx+j]=sum;				 
		}
	}
	return *this;
}

template <class T> 
Matrix<T>& Matrix<T>::unitVectorMatrix() {
	for (unsigned int i=0; i<dimx*dimy; i++) data[i].unitVector();
	return *this;
}


template <class T> 
Vector<T> Matrix<T>::Col(int x, int z)  const {
	Vector<T> c(Dimy());
	for (int y=0; y<Dimy(); y++) 
		c(y)=data[y*Dimx()+x](z);
	return c;
}

template <class T> 
Vector<T> Matrix<T>::Row(int y, int z) const {
	Vector<T> r(Dimx());
	for (int x=0; x<Dimx(); x++) 
		r(x)=data[y*Dimx()+x](z);
	return r;
}

template <class T> 
Matrix<T>& Matrix<T>::Col(const Vector<T>& c, int x, int z) {
	if (!(c.Dim() == Dimy())) std::cout << "Error, incompatible Col dimensions." << std::endl; 
	else {
		for (int y=0; y<Dimy(); y++) 
			data[y*Dimx()+x](z)=c(y);
	}
	return *this;
}

template <class T> 
Matrix<T>& Matrix<T>::Row(const Vector<T>& r, int y, int z) {
	if (!(r.Dim() == Dimx())) std::cout << "Error, incompatible Row dimensions." << std::endl; 
	else {
		for (int x=0; x<Dimx(); x++) 
			data[y*Dimx()+x](z)=r(x);
	}
	return *this;
}



/*

template <class T> 
void Matrix<T>::Col(const Vector<Vector<T> >& c, int x) {
	if (c.Dim() == Dimx()) {
		for (int y=0; y<Dimy(); y++) 
			data[y*Dimx()+x]=c(y);
	}
}

template <class T> 
void Matrix<T>::Row(const Vector<Vector<T> >& r, int y) {
	if (r.Dim() == Dimy()) {
		for (int x=0; x<Dimx(); x++) 
			data[y*Dimx()+x]=r(x);
	}
}


template <class T> 
Vector<Vector<T> > Matrix<T>::Col(int x)  const {
	Vector<Vector<T> > c(Dimy());
	for (int y=0; y<Dimy(); y++) 
		c(y)=data[y*Dimx()+x];
	return c;
}

template <class T> 
Vector<Vector<T> > Matrix<T>::Row(int y) const {
	Vector<Vector<T> > r(Dimx());
	for (int x=0; x<Dimx(); x++) 
		r(x)=data[y*Dimx()+x];
	return r;
}

template <class T> 
void Matrix<T>::Col(const Vector<Vector<T> >& c, int x) {
	if (c.Dim() == Dimx()) {
		for (int y=0; y<Dimy(); y++) 
			data[y*Dimx()+x]=c(y);
	}
}

template <class T> 
void Matrix<T>::Row(const Vector<Vector<T> >& r, int y) {
	if (r.Dim() == Dimy()) {
		for (int x=0; x<Dimx(); x++) 
			data[y*Dimx()+x]=r(x);
	}
}*/

#endif
