#ifndef _ipicmath_h_
#define _ipicmath_h_
#include "assert.h"
#include "math.h"
#include "stdlib.h" // for rand
#include "debug.h"

// valid if roundup power is representable.
inline int
pow2roundup (int x)
{
    assert(x>=0);
    //if (x < 0)
    //    return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

// does not work if highest non-sign bit is set
inline int
pow2rounddown (int x)
{
    assert(x>=0);
    //if (x < 0)
    //    return 0;

    // set all bits below highest bit
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    // set the bit higher than the highest bit
    x++;
    // shift it down and return it
    return (x >> 1);
}

// return integer ceiling of n over m
inline int ceiling_of_ratio(int n, int m)
{
  // probably this way is cheapest
  return ((n-1)/m+1);
  //return (n+m-1)/m;
  //const int out = ceil(n/double(m));
  //return out;
}

// round n up to next multiple of m
inline int roundup_to_multiple(int n, int m)
{
  //return ((n-1)/m+1)*m;
  return (n+m-1)/m*m;
}

// sample from clopen unit interval (0,1]
inline double sample_clopen_u_double()
{
  // old way (retained for sake of bit-wise code agreement)
  const double harvest = rand() / (double) RAND_MAX;
  return 1.0 - .999999 * harvest;
  // better way
  const double max_inv = 1./(double(RAND_MAX)+1);
  return (double(rand())+1)*max_inv;
}

// sample from open unit interval (0,1)
static inline double sample_open_u_double()
{
  const double max_inv = 1./(double(RAND_MAX)+2);
  return (double(rand())+1)*max_inv;
}

// sample from unit interval [0,1]
inline double sample_u_double()
{
  // old way
  return rand()/double(RAND_MAX);
  // faster way
  const double max_inv = 1./(double(RAND_MAX));
  return double(rand())*max_inv;
}

inline void sample_standard_maxwellian(double& u)
{
  // we sample a single component by pretending that it
  // is part of a two-dimensional joint distribution.
  const double prob = sqrt(-2.0 * log(sample_clopen_u_double()));
  const double theta = 2.0 * M_PI * sample_u_double();
  u = prob * cos(theta);
}

inline void sample_standard_maxwellian(double& u, double& v)
{
  // the distribution of the magnitude of (u,v)
  // can be integrated analytically
  const double prob = sqrt(-2.0 * log(sample_clopen_u_double()));
  const double theta = 2.0 * M_PI * sample_u_double();
  u = prob * cos(theta);
  v = prob * sin(theta);
}

inline void sample_standard_maxwellian(double& u, double& v, double& w)
{
  sample_standard_maxwellian(u,v);
  sample_standard_maxwellian(w);
}

inline void sample_maxwellian(double& u, double ut)
{
  sample_standard_maxwellian(u);
  u *= ut;
}

inline void sample_maxwellian(double& u, double& v, double& w,
  double ut, double vt, double wt)
{
  sample_standard_maxwellian(u,v,w);
  u *= ut;
  v *= vt;
  w *= wt;
}

inline void sample_maxwellian(
  double& u, double& v, double& w,
  double ut, double vt, double wt,
  double u0, double v0, double w0)
{
  sample_standard_maxwellian(u,v,w);
  u = u0 + ut*u;
  v = v0 + vt*v;
  w = w0 + wt*w;
}


// GLD mar 15
inline double fexp(double va, double ve, double lam, double fmax)
{
    double ff;
    ff=fmax*exp(-sqrt(pow(va,2)+pow(ve,2))/lam);
    return ff;
}

// GLD mar 15 - assumes that parallel is z
inline double sample_ramscb(
  double& u, double& v, double& w,
  double dvpa, double dvpe, double df,
  double vpa[798], double vpe[798],double f[798][798],
  double x[100001], double fun[100001],double lam, double fmax,double C)
{
  int N=798;
  int i, j, ip, ind, ind1, ind2;
  double b, R, xi, vec1, vec2, vec3, vec4;
  double w1, w2, fint, ffint;  
  int m=1;
  double vec0[1];

  //  cout << endl;
  //cout <<"calling sample_ramscb"<< endl;
  //cout <<"---------------------"<< endl;
  // Rejection-method

  ip=0;
  do
  {
    ip=ip+1;
   
    b=sample_u_double();
    ind=floor(b/df);
    R=lam*(x[ind]+(x[ind+1]-x[ind])/df*(b-fun[ind]));       

    xi=0.5*M_PI*sample_u_double();
    vec1=R*sin(xi); // parallel
    vec2=R*cos(xi); // perpendicular
    vec3=sample_u_double();

    ind1=floor((vec1-vpa[0])/dvpa);
    ind2=floor((vec2-vpe[0])/dvpe);
       
    if (((ind1>=0)&&(ind2>=0))&&((ind1<N-1)&&(ind2<N-1)))
    {
      w1=(vec1-vpa[ind1])/dvpa;
      w2=(vec2-vpe[ind2])/dvpe;
           
      fint=f[ind1+1][ind2+1]*w1*w2+f[ind1][ind2]*(1.0-w1)*(1.0-w2)+f[ind1+1][ind2]*w1*(1.0-w2)+f[ind1][ind2+1]*(1.0-w1)*w2;
      ffint=fexp(vec1,vec2,lam,fmax);
      
      if (vec3<fint/(C*ffint)) //accept
      {
	vec4=sample_u_double();
	if (vec4>0.5)
	  {
	    w=vec1;
	  }
	else
	  {
	    w=-vec1;
	  }
	vec3=2.0*M_PI*sample_u_double();
	u=vec2*cos(vec3);
	v=vec2*sin(vec3);
      }
      else
	ip=ip-1;
    }
    else
            ip=ip-1;
  } while (ip<1);
}

// add or subtract multiples of L until x lies
// between 0 and L.
//
/** RIFAI QUESTA PARTE questo e' la tomba delle performance*/
//inline void MODULO(double *x, double L)
//{
//  *x = *x - floor(*x / L) * L;
//}
// version of previous method that assumes Linv = 1/L
// (faster if 1/L is precomputed)
inline double modulo(double x, double L, double Linv)
{
  return x - floor(x * Linv) * L;
}

#endif
