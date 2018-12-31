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
#include <libxml/parser.h>
#include <glib.h>
#include <libintl.h>
#include <gsl/gsl_rng.h>
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "rk.h"
#include "rk_6_2.h"

#define DEBUG_RK_6_2 0          ///< macro to debug.

/**
 * Function to print a maxima format file to check the accuracy order of a 6
 * steps 2nd order Runge-Kutta method.
 */
void
tb_print_maxima_6_2 (FILE * file,       ///< file.
                     unsigned int nsteps __attribute__ ((unused)),
                     ///< steps number.
                     unsigned int order __attribute__ ((unused)))
  ///< accuracy order.
{
  fprintf (file, "b60+b61+b62+b63+b64+b65-1;\n");
  fprintf (file, "b61*t1+b62*t2+b63*t3+b64*t4+b65*t5-1/2;\n");
  fprintf (file, "b61*t1^2+b62*t2^2+b63*t3^2+b64*t4^2+b65*t5^2-1/3;\n");
}

/**
 * Function to print a maxima format file to check the accuracy order of a 6
 * steps 2nd order Runge-Kutta method.
 */
void
rk_print_maxima_6_2 (FILE * file,       ///< file.
                     unsigned int nsteps,       ///< steps number.
                     unsigned int order)        ///< accuracy order.
{
  tb_print_maxima_6_2 (file, nsteps, order);
  rk_print_maxima_6 (file);
}

/**
 * Function to obtain the coefficients of a 6 steps 2nd order Runge-Kutta 
 * method.
 */
int
rk_tb_6_2 (Optimize * optimize) ///< Optimize struct.
{
  long double *tb, *r;
#if DEBUG_RK_6_2
  fprintf (stderr, "rk_tb_6_2: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t6 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = r[1];
  b21 (tb) = r[2];
  t3 (tb) = r[3];
  b31 (tb) = r[4];
  b32 (tb) = r[5];
  t4 (tb) = r[6];
  b41 (tb) = r[7];
  b42 (tb) = r[8];
  b43 (tb) = r[9];
  t5 (tb) = r[10];
  b51 (tb) = r[11];
  b52 (tb) = r[12];
  b53 (tb) = r[13];
  b54 (tb) = r[14];
  b62 (tb) = r[15];
  b63 (tb) = r[16];
  b64 (tb) = r[17];
  b65 (tb) = r[18];
  b61 (tb) = (0.5L - b62 (tb) * t2 (tb) - b63 (tb) * t3 (tb)
              - b64 (tb) * t4 (tb) - b65 (tb) * t5 (tb)) / t1 (tb);
  rk_b_6 (tb);
#if DEBUG_RK_6_2
  rk_print_tb_6 (tb, "rk_tb_6_2", stderr);
  fprintf (stderr, "rk_tb_6_2: end\n");
#endif
  if (isnan (b61 (tb)))
    return 0;
  return 1;
}

/**
 * Function to obtain the coefficients of a 6 steps 2nd order, 3rd order in
 * equations depending only in time, Runge-Kutta method.
 */
int
rk_tb_6_2t (Optimize * optimize)        ///< Optimize struct.
{
  long double *tb, *r;
#if DEBUG_RK_6_2
  fprintf (stderr, "rk_tb_6_2t: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t6 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = r[1];
  b21 (tb) = r[2];
  t3 (tb) = r[3];
  b31 (tb) = r[4];
  b32 (tb) = r[5];
  t4 (tb) = r[6];
  b41 (tb) = r[7];
  b42 (tb) = r[8];
  b43 (tb) = r[9];
  t5 (tb) = r[10];
  b51 (tb) = r[11];
  b52 (tb) = r[12];
  b53 (tb) = r[13];
  b54 (tb) = r[14];
  b61 (tb) = r[15];
  b62 (tb) = r[16];
  b63 (tb) = r[17];
  b64 (tb) = (1.L / 3.L - 0.5L * t5 (tb)
              - b61 (tb) * t1 (tb) * (t1 (tb) - t5 (tb))
              - b62 (tb) * t2 (tb) * (t2 (tb) - t5 (tb))
              - b63 (tb) * t3 (tb) * (t3 (tb) - t5 (tb)))
    / (t4 (tb) * (t4 (tb) - t5 (tb)));
  b65 (tb) = (0.5L - b61 (tb) * t1 (tb) - b62 (tb) * t2 (tb)
              - b63 (tb) * t3 (tb) - b64 (tb) * t4 (tb)) / t5 (tb);
  rk_b_6 (tb);
#if DEBUG_RK_6_2
  rk_print_tb_6 (tb, "rk_tb_6_2t", stderr);
  fprintf (stderr, "rk_tb_6_2t: end\n");
#endif
  if (isnan (b65 (tb)) || isnan (b64 (tb)))
    return 0;
  return 1;
}

/**
 * Function to calculate the objective function of a 6 steps 2nd order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_2 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_6_2
  fprintf (stderr, "rk_objective_tb_6_2: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b50 (tb) < 0.L)
    o += b50 (tb);
  if (b60 (tb) < 0.L)
    o += b60 (tb);
  if (b61 (tb) < 0.L)
    o += b61 (tb);
  if (o < 0.L)
    {
      o = 40.L - o;
      goto end;
    }
  o = 30.L
    + fmaxl (1.L,
             fmaxl (t1 (tb),
                    fmaxl (t2 (tb),
                           fmaxl (t3 (tb), fmaxl (t4 (tb), t5 (tb))))));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_6_2
  fprintf (stderr, "rk_objective_tb_6_2: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_2: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 6 steps 2nd order, third
 * order in equations depending only on time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_2t (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_6_2
  fprintf (stderr, "rk_objective_tb_6_2t: start\n");
#endif
  tb = rk->tb->coefficient;
#if DEBUG_RK_6_2
  rk_print_tb_6 (tb, "rk_objective_tb_6_2t", stderr);
#endif
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b50 (tb) < 0.L)
    o += b50 (tb);
  if (b60 (tb) < 0.L)
    o += b60 (tb);
  if (b64 (tb) < 0.L)
    o += b64 (tb);
  if (b65 (tb) < 0.L)
    o += b65 (tb);
  if (o < 0.L)
    {
      o = 40.L - o;
      goto end;
    }
  o = 30.L
    + fmaxl (1.L,
             fmaxl (t1 (tb),
                    fmaxl (t2 (tb),
                           fmaxl (t3 (tb), fmaxl (t4 (tb), t5 (tb))))));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_6_2
  fprintf (stderr, "rk_objective_tb_6_2t: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_2t: end\n");
#endif
  return o;
}
