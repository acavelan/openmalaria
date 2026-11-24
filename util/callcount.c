/* LD_PRELOAD call-count shim: counts calls to libm/GSL/BLAS entry points.
 * Build: gcc -shared -fPIC -O2 -o callcount.so callcount.c -ldl -lm
 * Run:   LD_PRELOAD=./callcount.so openMalaria ... ; counts print at exit.
 * Used in PERF round 6 to prove library call counts are identical across trees.
 * NOTE: LD_PRELOAD interposition must run with LD_DEBUG=bindings if counts read 0. */
/* Call-count shim: LD_PRELOAD this, counts calls to libm/GSL/BLAS entry points. */
#define _GNU_SOURCE
#include <stdio.h>
#include <dlfcn.h>
#include <math.h>

static long n_exp, n_pow, n_log, n_log1p, n_expm1, n_ln, n_gsl_ln, n_gsl_P, n_gsl_Pinv, n_gsl_flat, n_gsl_pois, n_dgemm, n_dgemv, n_erf;
double (*r_exp)(double); double (*r_pow)(double,double); double (*r_log)(double);
double (*r_log1p)(double); double (*r_expm1)(double);
double (*r_gsl_ln)(const void*,double,double); double (*r_gsl_P)(double,double,double);
double (*r_gsl_Pinv)(double,double,double); double (*r_gsl_flat)(const void*,double,double);
long (*r_gsl_pois)(const void*,double);
void (*r_dgemm_)(const char*,const char*,const int*,const int*,const int*,const double*,const double*,const int*,const double*,const int*,const double*,double*,const int*);
void (*r_dgemv_)(const char*,const int*,const int*,const double*,const double*,const int*,const double*,const int*,const double*,double*,const int*);
double (*r_erf)(double);

#define RES(v,s) if(!v){ v=dlsym(RTLD_NEXT,s); if(!v) fprintf(stderr,"SHIM: missing %s\n",s); }
double exp(double x){ RES(r_exp,"exp"); __sync_fetch_and_add(&n_exp,1); return r_exp(x); }
double pow(double a,double b){ RES(r_pow,"pow"); __sync_fetch_and_add(&n_pow,1); return r_pow(a,b); }
double log(double x){ RES(r_log,"log"); __sync_fetch_and_add(&n_log,1); return r_log(x); }
double log1p(double x){ RES(r_log1p,"log1p"); __sync_fetch_and_add(&n_log1p,1); return r_log1p(x); }
double expm1(double x){ RES(r_expm1,"expm1"); __sync_fetch_and_add(&n_expm1,1); return r_expm1(x); }
double erf(double x){ RES(r_erf,"erf"); __sync_fetch_and_add(&n_erf,1); return r_erf(x); }
double gsl_ran_lognormal(const void*s,double m,double z){ RES(r_gsl_ln,"gsl_ran_lognormal"); __sync_fetch_and_add(&n_gsl_ln,1); return r_gsl_ln(s,m,z); }
double gsl_cdf_ugaussian_P(double x,double m,double z){ RES(r_gsl_P,"gsl_cdf_ugaussian_P"); __sync_fetch_and_add(&n_gsl_P,1); return r_gsl_P(x,m,z); }
double gsl_cdf_ugaussian_Pinv(double p,double m,double z){ RES(r_gsl_Pinv,"gsl_cdf_ugaussian_Pinv"); __sync_fetch_and_add(&n_gsl_Pinv,1); return r_gsl_Pinv(p,m,z); }
double gsl_ran_flat(const void*s,double a,double b){ RES(r_gsl_flat,"gsl_ran_flat"); __sync_fetch_and_add(&n_gsl_flat,1); return r_gsl_flat(s,a,b); }
long gsl_ran_poisson(const void*s,double m){ RES(r_gsl_pois,"gsl_ran_poisson"); __sync_fetch_and_add(&n_gsl_pois,1); return r_gsl_pois(s,m); }
void dgemm_(const char*a,const char*b,const int*m,const int*n,const int*k,const double*al,const double*x,const int*ldx,const double*y,const int*ldy,const double*be,double*c,const int*ldc){ RES(r_dgemm_,"dgemm_"); __sync_fetch_and_add(&n_dgemm,1); r_dgemm_(a,b,m,n,k,al,x,ldx,y,ldy,be,c,ldc); }
void dgemv_(const char*m,const int*r,const int*c,const double*a,const double*x,const int*ix,const double*y,const int*iy,const double*b,double*z,const int*iz){ RES(r_dgemv_,"dgemv_"); __sync_fetch_and_add(&n_dgemv,1); r_dgemv_(m,r,c,a,x,ix,y,iy,b,z,iz); }

__attribute__((destructor)) static void dump(void){
    fprintf(stderr,"CALLCOUNTS exp=%ld pow=%ld log=%ld log1p=%ld expm1=%ld erf=%ld dgemm=%ld dgemv=%ld gsl_lognormal=%ld gsl_P=%ld gsl_Pinv=%ld gsl_flat=%ld gsl_poisson=%ld\n",
        n_exp,n_pow,n_log,n_log1p,n_expm1,n_erf,n_dgemm,n_dgemv,n_gsl_ln,n_gsl_P,n_gsl_Pinv,n_gsl_flat,n_gsl_pois);
}
