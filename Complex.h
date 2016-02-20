#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include <math.h>
#include <iostream>

class Complex;
std::ostream& operator<<(std::ostream&, const Complex&);  

class Complex {
protected:
  double mod;
  double arg;
  double imag;
  double real;

  Complex(const double _mod, const double _arg, const double _real, const double _imag);
  void updateRectangular(const double _real, const double _imag);
  void updatePolar(const double _mod, const double _arg);
public:
  Complex();
  Complex(const double _real);
  Complex(const Complex& z);
   ~Complex();
  static Complex Polar(const double _mod, const double _arg);
  static Complex Rectangular(const double _real, const double _imag);

  Complex& operator = (const Complex& z);
  const int operator == (const Complex& z) const;
  const int operator <  (const Complex& z) const;
  const int operator <= (const Complex& z) const;
  const int operator >  (const Complex& z) const;
  const int operator >= (const Complex& z) const;
  const int operator != (const Complex& z) const;
  Complex  operator  +  (const double x) const;
  Complex  operator  +  (const Complex& z) const;
  Complex& operator +=  (const double x);
  Complex& operator +=  (const Complex& z);
  Complex  operator  -  (const double x) const;
  Complex  operator  -  (const Complex& z) const;
  Complex& operator -=  (const double x);
  Complex& operator -=  (const Complex& z);
  Complex  operator *   (const double x) const;
  Complex  operator *   (const Complex& z) const;  
  Complex& operator *=  (const double x);
  Complex& operator *=  (const Complex& z);

  Complex nthSqrt(const int nth) const;
  Complex nthPow(const int nth) const;
  Complex& Conjugate();
  Complex getConjugate();
  
  static double computeModulus(const double _real, const double _imag);
  static double computeArgument(const double _real, const double _imag);
  static double computeReal(const double _mod, const double _arg);
  static double computeImaginary(const double _mod, const double _arg);
  
//  Complex& Real(const double _real);
//  Complex& Imag(const double _imag);
//  Complex& Mod(const double _mod);
//  Complex& Arg(const double _arg);
  void Real(const double _real);
  void Imag(const double _imag);
  void Mod(const double _mod);
  void Arg(const double _arg);
  const double Real() const;
  const double Imag() const;
  const double Mod() const;
  const double Arg() const;  
  
  friend std::ostream& operator<<(std::ostream&, const Complex&);  
};

inline Complex::Complex(const double _mod, const double _arg, const double _real, const double _imag) :
  mod(_mod), arg(_arg), real(_real), imag(_imag) {}

inline Complex Complex::Polar(const double _mod, const double _arg) {
  return Complex(_mod,_arg,computeReal(_mod,_arg),computeImaginary(_mod,_arg));
}

inline Complex Complex::Rectangular(const double _real, const double _imag) {
  return Complex(computeModulus(_real,_imag), computeArgument(_real,_imag), _real, _imag);
}

inline double Complex::computeModulus(const double real, const double imag) {
	return sqrt(real*real + imag*imag);
}

inline double Complex::computeArgument(const double real, const double imag) {
  if (real > 0.0) {
    if (imag >= 0.0) return atan2(imag,real);
    else if (imag < 0.0) return atan2(imag,real) + 2*M_PI;
  }
  else if (real < 0.0) {
    if (imag >= 0.0) return atan2(imag,real) + M_PI;
    else if (imag < 0.0) return atan2(imag,real) + M_PI;
  }
  else {
    if (imag < 0.0) return 1.5*M_PI;
    else if (imag > 0.0) return .5*M_PI;
    else return 0;//undef
  }
  return 0;
}

inline double Complex::computeReal(const double _mod, const double _arg) {
  return _mod*cos(_arg);
}

inline double Complex::computeImaginary(const double _mod, const double _arg) {
  return _mod*sin(_arg);
}

Complex::Complex() : mod(0.0), arg(0.0), real(0.0), imag(0.0) {}

Complex::Complex(const double _real) : mod(_real), arg(0.0), real(_real), imag(0.0) {}

Complex::Complex(const Complex& z) : mod(z.mod), arg(z.arg), real(z.real), imag(z.imag) {}

Complex::~Complex() {}

Complex& Complex::operator=(const Complex& z) {
	mod=z.mod; arg=z.arg; real=z.real; imag=z.imag;
}

const int Complex::operator == (const Complex& z) const {
	return mod==z.real && arg==z.imag;
}

const int Complex::operator != (const Complex& z) const {
	return !(*this==z);
}

const int Complex::operator < (const Complex& z) const {
	return mod<z.mod;
}

const int Complex::operator <= (const Complex& z) const {
	return mod<=z.mod;
}

const int Complex::operator > (const Complex& z) const {
	return mod>z.mod;
}

const int Complex::operator >= (const Complex& z) const {
	return mod>=z.mod;
}

Complex  Complex::operator  + (const double x) const {
  return Rectangular(real+x,imag);
}

Complex  Complex::operator  + (const Complex& z) const {
  return Rectangular(real+z.real, imag+z.imag);
}

Complex& Complex::operator += (const double x) {
  updateRectangular(real+x,imag);
  return *this;
}

Complex& Complex::operator += (const Complex& z) {
  updateRectangular(real+z.real,imag+z.imag);
  return *this;
}

Complex  Complex::operator  - (const double x) const {
  return Rectangular(real-x, imag);
}

Complex  Complex::operator  - (const Complex& z) const {
  return Rectangular(real-z.real, imag-z.imag);
}

Complex& Complex::operator -= (const double x) {
  updateRectangular(real-x,imag);
  return *this;
}

Complex& Complex::operator -= (const Complex& z) {
  updateRectangular(real-z.real,imag-z.imag);
  return *this;
}

Complex  Complex::operator * (const double x) const {
  double _real = real*x - imag*0.0;
  double _imag = real*0.0 + imag*x;
  return Rectangular(_real,_imag);
}

Complex  Complex::operator *  (const Complex& z) const {
  double _real = real*z.real - imag*z.imag;
  double _imag = real*z.imag + imag*z.real;
  return Rectangular(_real,_imag);
}

Complex& Complex::operator *= (const double x) {
  double _real = real*x - imag*0.0;
  double _imag = real*0.0 + imag*x;
  updateRectangular(_real,_imag);
  return *this;
}

Complex& Complex::operator *= (const Complex& z) {
  double _real = real*z.real - imag*z.imag;
  double _imag = real*z.imag + imag*z.real;
  updateRectangular(_real,_imag);
  return *this;
}

std::ostream& operator<<(std::ostream& os, const Complex& z) {
  return os << z.Real() << " " << z.Imag() << " " << z.Mod() << " " << z.Arg();
}

Complex Complex::nthSqrt(const int nth) const {
  double n = (double)nth;
  return Polar(pow(mod,(1.0/n)),(arg/n));
}

Complex Complex::nthPow(const int nth) const {
  double n = (double)nth;
  return Polar(pow(mod,n),(arg*n));
}

Complex& Complex::Conjugate() {
	updateRectangular(real,-1.0*imag);
	return *this;
}

Complex Complex::getConjugate() {
	return Complex::Rectangular(real,-1.0*imag);
}

void Complex::Real(const double _real) { updateRectangular(_real,imag); }

void Complex::Imag(const double _imag) { updateRectangular(real,_imag); }

void Complex::Mod(const double _mod) { updatePolar(_mod,arg); }

void Complex::Arg(const double _arg) { updatePolar(mod,_arg); }

const double Complex::Real() const { return real; }

const double Complex::Imag() const { return imag; }

const double Complex::Mod() const { return mod; }

const double Complex::Arg() const { return arg; } 

void Complex::updatePolar(const double _mod, const double _arg) {
  mod = _mod;
  arg = _arg;
  real=computeReal(mod,arg);
  imag=computeImaginary(mod,arg);
}

void Complex::updateRectangular(const double _real, const double _imag) {
  real = _real;
  imag = _imag;
  mod=computeModulus(real,imag);
  arg=computeArgument(real,imag);
}

#endif

