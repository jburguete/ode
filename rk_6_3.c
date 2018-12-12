/*
ODE: a program to get optime Runge-Kutta and multi-steps methods.

Copyright 2011-2018, Javier Burguete Tolosa.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright notice,
		this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright notice,
		this list of conditions and the following disclaimer in the
		documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Javier Burguete Tolosa ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL Javier Burguete Tolosa OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * \file rk_6_4.c
 * \brief Source file to optimize Runge-Kutta 6 steps 4th order methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include "utils.h"
#include "rk.h"
#include "rk_6_2.h"
#include "rk_6_3.h"

long double value_rk_6_3[32] = {
  0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L
};

long double interval_rk_6_3[32] = {
  0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L
};

int imax_rk_6_3[32] =
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1,
  0, 1, 0, 1, 0, 1, 0
};

void
rk_make_random_6_3 (long double *rk, long double *vo)
{
  vo[0] = t1 (rk);
/*
	vo[1]=t3(rk);
	vo[2]=t4(rk);
	vo[3]=t2(rk);
	vo[4]=a54(rk);
	vo[5]=t5(rk);
	vo[6]=a65(rk);
	vo[7]=a43(rk);
	vo[8]=a32(rk);
	vo[9]=a21(rk);
	vo[10]=c63(rk);
	vo[11]=c62(rk);
	vo[12]=c61(rk);
	vo[13]=c21(rk);
	vo[14]=c32(rk);
	vo[15]=c43(rk);
	vo[16]=c54(rk);
	vo[17]=c65(rk);
	vo[18]=c41(rk);
	vo[19]=a41(rk);
	vo[20]=c51(rk);
	vo[21]=a51(rk);
	vo[22]=c64(rk);
	vo[23]=a64(rk);
	vo[24]=c53(rk);
	vo[25]=a53(rk);
	vo[26]=c42(rk);
	vo[27]=a42(rk);
	vo[28]=c52(rk);
	vo[29]=a52(rk);
	vo[30]=c31(rk);
	vo[31]=a31(rk);
*/
}

void
rk_print_maxima_6_3 (long double *rk, FILE * file)
{
  rk_print_maxima_6_2 (rk, file);
  fprintf (file, "b62*b21*t1+b63*(b31*t1+b32*t2)+b64*(b41*t1+b42*t2+b43*t3)"
           "+b65*(b51*t1+b52*t2+b53*t3+b54*t4)-1/6;\n");
  fprintf (file, "b61*t1^2+b62*t2^2+b63*t3^2+b64*t4^2+b65*t5^2-1/3;\n");
}

void
rk_6_3 (long double *data)
{
  int thread;
  long long int i;
  long double A[3], B[3], C[3], D[3], k, k2, rk[2 * size];
  gsl_rng *r;
  thread = (size_t) data;
  r = rng[thread];
  k2 = kmin;
  t6 (rk) = 1.L;
  for (i = cell[thread]; i < cell[thread + 1]; ++i)
    {
      t1 (rk) = vm[0] + vi[0] * gsl_rng_uniform (r);
/*
		t3(rk)=vm[1]+vi[1]*gsl_rng_uniform(r);
		t4(rk)=vm[2]+vi[2]*gsl_rng_uniform(r);
		t2(rk)=vm[3]+vi[3]*gsl_rng_uniform(r);
		a54(rk)=vm[4]+vi[4]*gsl_rng_uniform(r);
		t5(rk)=vm[5]+vi[5]*gsl_rng_uniform(r);
		a65(rk)=vm[6]+vi[6]*gsl_rng_uniform(r);
		a43(rk)=vm[7]+vi[7]*gsl_rng_uniform(r);
		a32(rk)=vm[8]+vi[8]*gsl_rng_uniform(r);
		a21(rk)=vm[9]+vi[9]*gsl_rng_uniform(r);
		c63(rk)=vm[10]+vi[10]*gsl_rng_uniform(r);
		c62(rk)=vm[11]+vi[11]*gsl_rng_uniform(r);
		c61(rk)=vm[12]+vi[12]*gsl_rng_uniform(r);
		c21(rk)=vm[13]+vi[13]*gsl_rng_uniform(r);
		c32(rk)=vm[14]+vi[14]*gsl_rng_uniform(r);
		c43(rk)=vm[15]+vi[15]*gsl_rng_uniform(r);
		c54(rk)=vm[16]+vi[16]*gsl_rng_uniform(r);
		c65(rk)=vm[17]+vi[17]*gsl_rng_uniform(r);
		c41(rk)=vm[18]+vi[18]*gsl_rng_uniform(r);
		a41(rk)=vm[19]+vi[19]*random_zero(r);
		c51(rk)=vm[20]+vi[20]*gsl_rng_uniform(r);
		a51(rk)=vm[21]+vi[21]*random_zero(r);
		c64(rk)=vm[22]+vi[22]*gsl_rng_uniform(r);
		a64(rk)=vm[23]+vi[23]*random_zero(r);
		c53(rk)=vm[24]+vi[24]*gsl_rng_uniform(r);
		a53(rk)=vm[25]+vi[25]*random_zero(r);
		c42(rk)=vm[26]+vi[26]*gsl_rng_uniform(r);
		a42(rk)=vm[27]+vi[27]*random_zero(r);
		c52(rk)=vm[28]+vi[28]*gsl_rng_uniform(r);
		a52(rk)=vm[29]+vi[29]*random_zero(r);
		c31(rk)=vm[30]+vi[30]*gsl_rng_uniform(r);
		a31(rk)=vm[31]+vi[31]*random_zero(r);
		A[0]=t1(rk);
		B[0]=t2(rk);
		C[0]=t3(rk);
		D[0]=0.5L-b64(rk)*t4(rk)-b65(rk)*t5(rk);
		A[1]=sqr(t1(rk));
		B[1]=sqr(t2(rk));
		C[1]=sqr(t3(rk));
		D[1]=1.L/3.L-b64(rk)*sqr(t4(rk))-b65(rk)*sqr(t5(rk));
		A[2]=0.L;
		B[2]=b21(rk)*t1(rk);
		C[2]=b31(rk)*t1(rk)+b32(rk)*t2(rk);
		D[2]=1.L/6.L-b64(rk)*(b41(rk)*t1(rk)+b42(rk)*t2(rk)+b43(rk)*t3(rk))
			-b65(rk)*(b51(rk)*t1(rk)+b52(rk)*t2(rk)+b53(rk)*t3(rk)
			+b54(rk)*t4(rk));
		solve_3(A,B,C,D);
		a63(rk)=(D[2]-a64(rk)*b43(rk)-a65(rk)*b53(rk))/c63(rk);
		if (a63(rk)<0.L || isnan(a63(rk))) continue;
		a62(rk)=(D[1]-a63(rk)*b32(rk)-a64(rk)*b42(rk)-a65(rk)*b52(rk))/c62(rk);
		if (a62(rk)<0.L || isnan(a62(rk))) continue;
		a61(rk)=(D[0]-a62(rk)*b21(rk)-a63(rk)*b31(rk)-a64(rk)*b41(rk)
			-a65(rk)*b51(rk))/c61(rk);
		if (a61(rk)<0.L || isnan(a61(rk))) continue;
*/
      c65 (rk) = c54 (rk) = c43 (rk) = c32 (rk) = c21 (rk) = t1 (rk);
      t5 (rk) = 1.L - c65 (rk);
      t2 (rk) = t1 (rk) + c21 (rk);
      t3 (rk) = t2 (rk) + c32 (rk);
      a32 (rk) = a21 (rk) = a65 (rk) = 1.L;
      c63 (rk) = a63 (rk) = c62 (rk) = a62 (rk) = c61 (rk) = a61 (rk) =
        c41 (rk) = a41 (rk) = c51 (rk) = a41 (rk) = c64 (rk) = a64 (rk) =
        c53 (rk) = a53 (rk) = c42 (rk) = a42 (rk) = c52 (rk) = a52 (rk) =
        c31 (rk) = a31 (rk) = 0.L;
      A[0] = c21 (rk) * t1 (rk) + c32 (rk) * t2 (rk) + c43 (rk) * t3 (rk);
      B[0] = 1.L;
      C[0] = 0.L;
      D[0] = 0.5L - b65 (rk) * t5 (rk);
      A[1] =
        c21 (rk) * sqr (t1 (rk)) + c32 (rk) * sqr (t2 (rk)) +
        c43 (rk) * sqr (t3 (rk));
      B[1] = 0.L;
      C[1] = 1.L;
      D[1] = 1.L / 3.L - b65 (rk) * sqr (t5 (rk));
      A[2] = c32 (rk) * c21 (rk) * t1 (rk) + (c65 (rk) + c54 (rk) + c43 (rk))
        * (c32 (rk) * t2 (rk) + c21 (rk) * t1 (rk)) + (c65 (rk) +
                                                       c54 (rk)) * c43 (rk) *
        t3 (rk);
      B[2] = c65 (rk);
      C[2] = 0.L;
      D[2] = 1.L / 6.L;
      solve_3 (A, B, C, D);
      t4 (rk) = D[2] / D[1];
      if (isnan (t4 (rk)))
        continue;
      a54 (rk) = D[1] / (c54 (rk) * t4 (rk));
      if (a54 (rk) <= 0.L || isnan (a54 (rk)))
        continue;
      a43 (rk) = D[0] / a54 (rk);
      if (a43 (rk) < 0.L || isnan (a43 (rk)))
        continue;
/*
		c65(rk)=c54(rk)=c43(rk)=c32(rk)=c21(rk)=t1(rk);
		t5(rk)=1.L-c65(rk);
		t2(rk)=t1(rk)+c21(rk);
		a21(rk)=a65(rk)=1.L;
		c63(rk)=a63(rk)=c62(rk)=a62(rk)=c61(rk)=a61(rk)=c41(rk)=a41(rk)
			=c51(rk)=a41(rk)=c64(rk)=a64(rk)=c53(rk)=a53(rk)=c42(rk)=a42(rk)
			=c52(rk)=a52(rk)=c31(rk)=a31(rk)=0.L;
		A[0]=c21(rk)*t1(rk)+c32(rk)*t2(rk);
		B[0]=c43(rk)*t3(rk);
		C[0]=c54(rk)*t4(rk);
		D[0]=0.5L-b65(rk)*t5(rk);
		A[1]=c21(rk)*sqr(t1(rk))+c32(rk)*sqr(t2(rk));
		B[1]=c43(rk)*sqr(t3(rk));
		C[1]=c54(rk)*sqr(t4(rk));
		D[1]=1.L/3.L-b65(rk)*sqr(t5(rk));
		A[2]=c32(rk)*c21(rk)*t1(rk)+(c65(rk)+c54(rk)+c43(rk))
			*(c32(rk)*t2(rk)+c21(rk)*t1(rk));
		B[2]=(c65(rk)+c54(rk))*c43(rk)*t3(rk);
		C[2]=c65(rk)*c54(rk)*t4(rk);
		D[2]=1.L/6.L;
		solve_3(A,B,C,D);
		a54(rk)=D[2];
		if (a54(rk)<=0.L || isnan(a54(rk))) continue;
		a43(rk)=D[1]/D[2];
		if (a43(rk)<=0.L || isnan(a43(rk))) continue;
		a32(rk)=D[0]/D[1];
		if (a32(rk)<0.L || isnan(a32(rk))) continue;
*/
/*
		c65(rk)=c54(rk)=c43(rk)=c32(rk)=c21(rk)=t1(rk);
		t5(rk)=1.L-c65(rk);
		a65(rk)=1.L;
		c63(rk)=a63(rk)=c62(rk)=a62(rk)=c61(rk)=a61(rk)=c41(rk)=a41(rk)
			=c51(rk)=a41(rk)=c64(rk)=a64(rk)=c53(rk)=a53(rk)=c42(rk)=a42(rk)
			=c52(rk)=a52(rk)=c31(rk)=a31(rk)=0.L;
		A[0]=t1(rk);
		B[0]=t2(rk);
		C[0]=t3(rk);
		D[0]=0.5L-b64(rk)*t4(rk)-b65(rk)*t5(rk);
		A[1]=sqr(t1(rk));
		B[1]=sqr(t2(rk));
		C[1]=sqr(t3(rk));
		D[1]=1.L/3.L-b64(rk)*sqr(t4(rk))-b65(rk)*sqr(t5(rk));
		A[2]=(c65(rk)+c54(rk)+c43(rk)+c32(rk))*t1(rk);
		B[2]=(c65(rk)+c54(rk)+c43(rk))*t2(rk);
		C[2]=(c65(rk)+c54(rk))*t3(rk);
		D[2]=1.L/6.L-b65(rk)*b54(rk)*t4(rk);
		solve_3(A,B,C,D);
		a43(rk)=D[2]/(a65(rk)*a54(rk)*c43(rk));
		if (a43(rk)<0.L || isnan(a43(rk))) continue;
		a32(rk)=D[1]/(a65(rk)*a54(rk)*a43(rk)*c32(rk));
		if (a32(rk)<0.L || isnan(a32(rk))) continue;
		a21(rk)=D[0]/(a65(rk)*a54(rk)*a43(rk)*a32(rk)*c21(rk));
		if (a21(rk)<0.L || isnan(a21(rk))) continue;
*/
      if (!rk_tvd_6 (rk))
        continue;
      k = rk_max_6 (rk);
      if (k < k2)
        {
          k2 = k;
          memcpy (rk + size, rk, size * sizeof (long double));
          printf ("rank=%d thread=%d k2=%.19Le t=%Lg\n", rank, thread, k2,
                  (100.L * (i - cell[thread])) / (cell[thread + 1] -
                                                  cell[thread]));
        }
    }
  if (k2 < kmin)
    {
      g_mutex_lock (mutex);
      kmin = k2;
      memcpy (rko, rk + size, size * sizeof (long double));
      g_mutex_unlock (mutex);
      rk_print_6 (rko, stdout);
      printf ("rank=%d thread=%d kmin=%.19Le\n", rank, thread, kmin);
    }
}
