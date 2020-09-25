#include "wcomplex.h"

complex wcomplex(double x,double y)
{
	complex z;
	z.re=x; z.im=y;
	return(z);
}

complex coper(complex x,char op,complex y)
{
	complex z;
	double ft;

	switch(op) {
		case '+' :
			z.re=x.re+y.re;
			z.im=x.im+y.im;
			return(z);
		case '-' :
			z.re=x.re-y.re;
			z.im=x.im-y.im;
			return(z);
		case '*' :
			z.re=x.re*y.re-x.im*y.im;
			z.im=x.re*y.im+x.im*y.re;
			return(z);
		case '/' :
			ft=y.re*y.re+y.im*y.im;
			z.re=(x.re*y.re+x.im*y.im)/ft;
			z.im=(x.im*y.re-x.re*y.im)/ft;
			return(z);
		default :
			printf("\nUnknown complex operator");
			exit(0);z.re=z.im=0.0;
			return(z);
	}
}

complex conj(complex x)
{
	complex y;
	y.re=x.re;
	y.im=-x.im;
	return(y);
}
double norm(complex x)
{
	double y;
	y=x.re*x.re+x.im*x.im;
	return(y);
}
double arg(complex x)
{
	double y;
	if(x.re==0.0){
		if(x.im>0.0)return(M_PI_2);
		else if(x.im<0.0) return(-M_PI_2);
		else return(0.0);
//		else {printf("\nmagnitude and phase are 0 !"); exit(0);}
	}
	y=atan2(x.im,x.re);
	return(y);
}
complex polar(double mag,double ang)
{
	complex x;
	x.re=mag*cos(ang);
	x.im=mag*sin(ang);
	return(x);
}

void symmetry(complex x[],short N,complex cx)
{
	short i;

	x[N>>1]=cx;
	for(i=1;i<(N>>1);i++) x[N-i]=conj(x[i]);
	return;
}

void co_to_mp(complex x[],short N)
{
	short i;

	for(i=0;i<(N>>1)+1;i++)
		x[i] = wcomplex(sqrt(norm(x[i])),arg(x[i]));
	symmetry(x,N,x[N>>1]);
	return;
}

void energy(complex x[], short N)
{
	short i;

	for(i=0;i<(N>>1)+1;i++)	//x[i].re = sqrt(norm(x[i]));
		x[i].re = norm(x[i]);
//	symmetry(x,N,x[N>>1]);
	return;
}

void mp_to_co(complex x[],short N)
{
	short i;

	for(i=0;i<(N>>1)+1;i++) x[i] = polar(x[i].re,x[i].im);
	symmetry(x,N,x[N>>1]);
	return;
}

double mod_twopi(double x)
{
	x = fmod(x,(2.0*M_PI));
	while(x>M_PI) { x -= 2.0*M_PI; }
	while(x<-M_PI) { x+= 2.0*M_PI; }
	return(x);
}
double mod_pi(double x)
{
	x = fmod(x,M_PI);
	while(x>M_PI_2) { x -= M_PI; }
	while(x<-M_PI_2) { x += M_PI; }
	return(x);
}

void time_shift(complex x[],short N,short shift)
{
	short i;
	double ang;

	ang=2.0*M_PI*(double)shift/(double)N;
	for(i=1;i<=N>>1;i++) {
		x[i].im=mod_twopi(x[i].im+ang*i);
		x[N-i].im=-x[i].im;
	}
	return;
}

void errmess(char x[])
{
//	printf("\nInsufficient memory for '%s' arrary.",x);
	printf("\nModify '%s' function.",x);
	exit (1);
}

void bit_reverse(complex x[],short Dim)
{
	short i,j,k;
	complex temp;

	for(i=1;i<1<<Dim;i++) {
		k = 0;
		for(j=0;j<Dim;j++) k += ((i>>j)%2)<<(Dim-1-j);
		if(k>i) { temp=x[i]; x[i]=x[k]; x[k]=temp; }
	}
	return;
}

void DIT(complex x[],short Dim, short pn)
{
	short N=1<<Dim,i,j,k;//	double NORM=sqrt((double)N);
	complex y[1024],w[512],temp;

	if (N>1024) errmess("FFT");
	for(i=0;i<Dim;i++) { //begin fft
		for(j=0;j<1<<i;j++) {
			if(j==0 || (N>>i+1)==0) temp=wcomplex(1.0,0.0);
			else temp =	wcomplex(cos(-2.0*M_PI*pn*j*(N>>i+1)/N),
						sin(-2.0*M_PI*pn*j*(N>>i+1)/N));
			for(k=0;k<(N>>i+1);k++) w[k+j*(N>>i+1)] = temp;
		}
		for(j=0;j<(N>>1);j++) {
			if(i!=0) x[(j<<1)+1] = coper(x[(j<<1)+1],'*',w[j]);
			y[j] = coper(x[j<<1],'+',x[(j<<1)+1]);
			y[j+(N>>1)] = coper(x[j<<1],'-',x[(j<<1)+1]);
		}
		for(j=0;j<N;j++) x[j]=y[j];  // in-place
	}//end fft for loop
	if(pn==-1) for(i=0;i<N;i++) {x[i].re/=N;x[i].im/=N;}
	return;
}

void fft(complex x[],short N)
{
	short i,Dim=0;

	for(i=1;i<14;i++) if(N==1<<i) {Dim=i; i=15;}
	if(Dim==0) { printf("\nfft dimension = %d!!",Dim); exit(1);}
	bit_reverse(x,Dim);
	DIT(x,Dim,1);
	return;
}

void ifft(complex x[],short N)
{
	short i,Dim=0;

	for(i=1;i<14;i++) if(N==1<<i) {Dim=i; i=15;}
	if(Dim==0) {printf("\nfft dimension = %d!!",Dim); exit(1);}
	bit_reverse(x,Dim);
	DIT(x,Dim,-1);
	return;
}
