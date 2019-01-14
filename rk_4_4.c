/*
Multisteps: a program to get optime Runge-Kutta and multi-steps methods.

Copyright 2011-2019, Javier Burguete Tolosa.

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
 * \file rk_4_4.c
 * \brief Source file to optimize Runge-Kutta 5 steps 4th order methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2019.
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
#include "rk_4_4.h"

#define DEBUG_RK_4_4 0          ///< macro to debug.

/**
 * Function to obtain the coefficients of a 4 steps 4th order Runge-Kutta 
 * method.
 */
int
rk_tb_4_4 (Optimize * optimize) ///< Optimize struct.
{
  long double *tb, *r;
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_tb_4_4: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t4 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = r[1];
  t3 (tb) = 1.L;
  b43 (tb) = (0.25L - 1.L / 3.L * t1 (tb)
              - (1.L / 3.L - 0.5L * t1 (tb)) * t2 (tb))
    / (t3 (tb) * (t3 (tb) - t2 (tb)) * (t3 (tb) - t1 (tb)));
  if (isnan (b43 (tb)))
    return 0;
  b42 (tb) = (1.L / 3.L - 0.5L * t1 (tb)
              - b43 (tb) * t3 (tb) * (t3 (tb) - t1 (tb)))
    / (t2 (tb) * (t2 (tb) - t1 (tb)));
  if (isnan (b42 (tb)))
    return 0;
  b41 (tb) = (0.5L - b42 (tb) * t2 (tb) - b43 (tb) * t3 (tb)) / t1 (tb);
  if (isnan (b41 (tb)))
    return 0;
  b32 (tb) = (1.L / 12.L - 1.L / 6.L * t1 (tb))
    / (b43 (tb) * t2 (tb) * (t2 (tb) - t1 (tb)));
  if (isnan (b32 (tb)))
    return 0;
  b31 (tb) = ((0.125L - 1.L / 6.L * t2 (tb)) / (b43 (tb) * (t3 (tb) - t2 (tb)))
              - b32 (tb) * t2 (tb)) / t1 (tb);
  if (isnan (b31 (tb)))
    return 0;
  b21 (tb) = 1.L / 24.L / (t1 (tb) * b43 (tb) * b32 (tb));
  if (isnan (b21 (tb)))
    return 0;
  rk_b_4 (tb);
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_tb_4_4: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 4 steps 4th order, 5th order in
 * equations depending only in time, Runge-Kutta method.
 */
int
rk_tb_4_4t (Optimize * optimize)        ///< Optimize struct.
{
  long double *tb, *r;
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_tb_4_4t: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t4 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = 0.5L * (t1 (tb) - 0.6L) / (t1 (tb) - 0.5L);
  t3 (tb) = 1.L;
  b43 (tb) = (0.25L - 1.L / 3.L * t1 (tb)
              - (1.L / 3.L - 0.5L * t1 (tb)) * t2 (tb))
    / (t3 (tb) * (t3 (tb) - t2 (tb)) * (t3 (tb) - t1 (tb)));
  if (isnan (b43 (tb)))
    return 0;
  b42 (tb) = (1.L / 3.L - 0.5L * t1 (tb)
              - b43 (tb) * t3 (tb) * (t3 (tb) - t1 (tb)))
    / (t2 (tb) * (t2 (tb) - t1 (tb)));
  if (isnan (b42 (tb)))
    return 0;
  b41 (tb) = (0.5L - b42 (tb) * t2 (tb) - b43 (tb) * t3 (tb)) / t1 (tb);
  if (isnan (b41 (tb)))
    return 0;
  b32 (tb) = (1.L / 12.L - 1.L / 6.L * t1 (tb))
    / (b43 (tb) * t2 (tb) * (t2 (tb) - t1 (tb)));
  if (isnan (b32 (tb)))
    return 0;
  b31 (tb) = ((0.125L - 1.L / 6.L * t2 (tb)) / (b43 (tb) * (t3 (tb) - t2 (tb)))
              - b32 (tb) * t2 (tb)) / t1 (tb);
  if (isnan (b31 (tb)))
    return 0;
  b21 (tb) = 1.L / 24.L / (t1 (tb) * b43 (tb) * b32 (tb));
  if (isnan (b21 (tb)))
    return 0;
  rk_b_4 (tb);
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_tb_4_4t: end\n");
#endif
  return 1;
}

/**
 * Function to calculate the objective function of a 4 steps 4th order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_4_4 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_objective_tb_4_4: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b21 (tb) < 0.L)
    o += b21 (tb);
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b31 (tb) < 0.L)
    o += b31 (tb);
  if (b32 (tb) < 0.L)
    o += b32 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b41 (tb) < 0.L)
    o += b41 (tb);
  if (b42 (tb) < 0.L)
    o += b42 (tb);
  if (b43 (tb) < 0.L)
    o += b43 (tb);
  if (o < 0.L)
    {
      o = 40.L - o;
      goto end;
    }
  o = 30.L + fmaxl (1.L, fmaxl (t1 (tb), fmaxl (t2 (tb), t3 (tb))));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_objective_tb_4_4: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_4_4: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 4 steps 4th order, 5th
 * order in equations depending only in time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_4_4t (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_objective_tb_4_4: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b21 (tb) < 0.L)
    o += b21 (tb);
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b31 (tb) < 0.L)
    o += b31 (tb);
  if (b32 (tb) < 0.L)
    o += b32 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b41 (tb) < 0.L)
    o += b41 (tb);
  if (b42 (tb) < 0.L)
    o += b42 (tb);
  if (b43 (tb) < 0.L)
    o += b43 (tb);
  if (o < 0.L)
    {
      o = 40.L - o;
      goto end;
    }
  o = 30.L + fmaxl (1.L, fmaxl (t1 (tb), fmaxl (t2 (tb), t3 (tb))));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_4_4
  fprintf (stderr, "rk_objective_tb_4_4: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_4_4: end\n");
#endif
  return o;
}
