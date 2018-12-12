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
 * \file rk_6_2.c
 * \brief Source file to optimize Runge-Kutta 6 steps 2nd order methods.
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

long double value_rk_6_2[34] = {
  0.5L, 1.0L, 1.5L, 2.0L, 2.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L, 0.5L, 0.5L
};

long double interval_rk_6_2[34] = {
  0.5L, 1.0L, 1.5L, 2.0L, 2.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L, 0.5L,
  0.5L, 0.5L,
  0.5L, 0.5L, 0.5L, 0.5L
};

int imax_rk_6_2[34] =
  { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
  0, 1, 0, 1, 0, 1, 0, 1, 0
};

void
rk_make_random_6_2 (long double *rk, long double *vo)
{
  vo[0] = t1 (rk);
  vo[1] = t2 (rk);
  vo[2] = t3 (rk);
  vo[3] = t4 (rk);
  vo[4] = t5 (rk);
  vo[5] = c61 (rk);
  vo[6] = a65 (rk);
  vo[7] = a54 (rk);
/*
	vo[8]=a43(rk);
	vo[9]=a32(rk);
	vo[10]=a21(rk);
	vo[11]=c21(rk);
	vo[12]=c32(rk);
	vo[13]=c43(rk);
	vo[14]=c54(rk);
	vo[15]=c65(rk);
	vo[16]=c64(rk);
	vo[17]=a64(rk);
	vo[18]=c53(rk);
	vo[19]=a53(rk);
	vo[20]=c42(rk);
	vo[21]=a42(rk);
	vo[22]=c41(rk);
	vo[23]=a41(rk);
	vo[24]=c51(rk);
	vo[25]=a51(rk);
	vo[26]=c63(rk);
	vo[27]=a63(rk);
	vo[28]=c62(rk);
	vo[29]=a62(rk);
	vo[30]=c52(rk);
	vo[31]=a52(rk);
	vo[32]=c31(rk);
	vo[33]=a31(rk);
*/
}

void
rk_print_maxima_6_2 (long double *rk, FILE * file)
{
  fprintf (file, "a60+a61+a62+a63+a64+a65-1;\n");
  fprintf (file, "b61*t1+b62*t2+b63*t3+b64*t4+b65*t5-1/2;\n");
}

static inline double
prueba (double x)
{
  return fmaxl (0., x + x - 1.);
}

void
rk_6_2 (long double *data)
{
  int thread;
  long long int i;
  long double k, k2, rk[2 * size];
  gsl_rng *r;
  thread = (size_t) data;
  r = rng[thread];
  k2 = kmin;
  t6 (rk) = 1.L;
  for (i = cell[thread]; i < cell[thread + 1]; ++i)
    {
      t1 (rk) = vm[0] + vi[0] * gsl_rng_uniform (r);
      t2 (rk) = vm[1] + vi[1] * gsl_rng_uniform (r);
      t3 (rk) = vm[2] + vi[2] * gsl_rng_uniform (r);
      t4 (rk) = vm[3] + vi[3] * gsl_rng_uniform (r);
      t5 (rk) = vm[4] + vi[4] * gsl_rng_uniform (r);
      c61 (rk) = vm[5] + vi[5] * gsl_rng_uniform (r);
      a65 (rk) = vm[6] + vi[6] * gsl_rng_uniform (r);
      a54 (rk) = vm[7] + vi[7] * gsl_rng_uniform (r);
/*
		a43(rk)=vm[8]+vi[8]*gsl_rng_uniform(r);
		a32(rk)=vm[9]+vi[9]*gsl_rng_uniform(r);
		a21(rk)=vm[10]+vi[10]*gsl_rng_uniform(r);
		c21(rk)=vm[5]+vi[5]*gsl_rng_uniform(r);
		c32(rk)=vm[6]+vi[6]*gsl_rng_uniform(r);
		c43(rk)=vm[7]+vi[7]*gsl_rng_uniform(r);
		c54(rk)=vm[8]+vi[8]*gsl_rng_uniform(r);
		c65(rk)=vm[10]+vi[10]*gsl_rng_uniform(r);
		c64(rk)=vm[16]+vi[16]*gsl_rng_uniform(r);
		a64(rk)=vm[17]+vi[17]*prueba(gsl_rng_uniform(r));
		c53(rk)=vm[18]+vi[18]*gsl_rng_uniform(r);
		a53(rk)=vm[19]+vi[19]*prueba(gsl_rng_uniform(r));
		c42(rk)=vm[20]+vi[20]*gsl_rng_uniform(r);
		a42(rk)=vm[21]+vi[21]*prueba(gsl_rng_uniform(r));
		c41(rk)=vm[22]+vi[22]*gsl_rng_uniform(r);
		a41(rk)=vm[23]+vi[23]*prueba(gsl_rng_uniform(r));
		c51(rk)=vm[24]+vi[24]*gsl_rng_uniform(r);
		a51(rk)=vm[25]+vi[25]*prueba(gsl_rng_uniform(r));
		c63(rk)=vm[26]+vi[26]*gsl_rng_uniform(r);
		a63(rk)=vm[27]+vi[27]*prueba(gsl_rng_uniform(r));
		c62(rk)=vm[28]+vi[28]*gsl_rng_uniform(r);
		a62(rk)=vm[29]+vi[29]*prueba(gsl_rng_uniform(r));
		c52(rk)=vm[30]+vi[30]*gsl_rng_uniform(r);
		a52(rk)=vm[31]+vi[31]*prueba(gsl_rng_uniform(r));
		c31(rk)=vm[32]+vi[32]*gsl_rng_uniform(r);
		a31(rk)=vm[33]+vi[33]*prueba(gsl_rng_uniform(r));
		a61(rk)=((0.5L-b62(rk)*t2(rk)-b63(rk)*t3(rk)-b64(rk)*t4(rk)
			-b65(rk)*t5(rk))/t1(rk)-a62(rk)*b21(rk)-a63(rk)*b31(rk)
			-a64(rk)*b41(rk)-a65(rk)*b51(rk))/c61(rk);
		if (a61(rk)<0.L || isnan(a61(rk))) continue;
*/
      c21 (rk) = c32 (rk) = c43 (rk) = c54 (rk) = c65 (rk) = t1 (rk);
      t2 (rk) = t1 (rk) + c21 (rk);
      t3 (rk) = t2 (rk) + c32 (rk);
      a32 (rk) = a21 (rk) = 1.L;
      c61 (rk) = a61 (rk) = c64 (rk) = a64 (rk) = c53 (rk) = a53 (rk) =
        c51 (rk) = a51 (rk) = c63 (rk) = a63 (rk) = c62 (rk) = a62 (rk) =
        c52 (rk) = a52 (rk) = c31 (rk) = a31 (rk) = 0.L;
      a43 (rk) =
        (((0.5L / a65 (rk) - c65 (rk) * t5 (rk)) / a54 (rk) -
          c54 (rk) * t4 (rk))) / (c43 (rk) * t3 (rk) + c32 (rk) * t2 (rk) +
                                  c21 (rk) * t1 (rk));
      if (a43 (rk) < 0.L || isnan (a43 (rk)))
        continue;
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
