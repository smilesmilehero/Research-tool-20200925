// complex define & complex function propotype

#include 	<math.h>
#include 	<stdio.h>
#include	<stdlib.h>

#ifndef COMPLEX
#define COMPLEX
typedef struct{double re,im;} complex;
#define		M_PI	3.14159265358979
#define		M_PI_2	M_PI/2.0
#endif

void errmess(char*);
complex wcomplex(double,double);
complex coper(complex,char,complex);
complex conj(complex);
double norm(complex);
double arg(complex);
complex polar(double,double);
void fft(complex*,short);
void ifft(complex*,short);
void co_to_mp(complex*,short);
void mp_to_co(complex*,short);
void energy(complex*,short);
double mod_twopi(double);
double mod_pi(double);
void symmetry(complex*,short,complex);
