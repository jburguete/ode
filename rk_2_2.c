/*
ODE: a program to get optime Runge-Kutta and multi-steps methods.

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
 * \file rk_2_2.c
 * \brief Source file to optimize Runge-Kutta 3 steps 2nd order methods.
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
#include "rk_2_2.h"

#define DEBUG_RK_2_2 0          ///< macro to debug.

/**
 * Function to obtain the coefficients of a 2 steps 2nd order Runge-Kutta 
 * method.
 */
int
rk_tb_2_2 (Optimize * optimize) ///< Optimize struct.
{
  long double *tb, *r;
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_tb_2_2: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t2 (tb) = 1.L;
  t1 (tb) = r[0];
  b21 (tb) = 0.5L / t1 (tb);
  rk_b_2 (tb);
#if DEBUG_RK_2_2
  rk_print_tb (optimize, "rk_tb_2_2", stderr);
  fprintf (stderr, "rk_tb_2_2: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 2 steps 2nd order, 3rd order in
 * equations depending only in time, Runge-Kutta method.
 */
int
rk_tb_2_2t (Optimize * optimize)        ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_tb_2_2t: start\n");
#endif
  tb = optimize->coefficient;
  t2 (tb) = 1.L;
  t1 (tb) = 2.L / 3.L;
  b21 (tb) = 0.5L / t1 (tb);
  rk_b_2 (tb);
#if DEBUG_RK_2_2
  rk_print_tb (optimize, "rk_tb_2_2t", stderr);
  fprintf (stderr, "rk_tb_2_2t: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 2 steps 1st-2nd order Runge-Kutta 
 * pair.
 */
int
rk_tb_2_2p (Optimize * optimize)        ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_tb_2_2p: start\n");
#endif
  if (!rk_tb_2_2 (optimize))
    return 0;
  tb = optimize->coefficient;
  e20 (tb) = 1.L;
#if DEBUG_RK_2_2
  rk_print_tb (optimize, "rk_tb_2_2p", stderr);
  fprintf (stderr, "rk_tb_2_2p: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 2 steps 1st-2nd order, 1st-3rd order
 * in equations depending only on time, Runge-Kutta pair.
 */
int
rk_tb_2_2tp (Optimize * optimize)       ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_tb_2_2tp: start\n");
#endif
  if (!rk_tb_2_2t (optimize))
    return 0;
  tb = optimize->coefficient;
  e20 (tb) = 1.L;
#if DEBUG_RK_2_2
  rk_print_tb (optimize, "rk_tb_2_2tp", stderr);
  fprintf (stderr, "rk_tb_2_2tp: end\n");
#endif
  return 1;
}

/**
 * Function to calculate the objective function of a 2 steps 2nd order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_2_2 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_objective_tb_2_2: start\n");
#endif
  tb = rk->tb->coefficient;
  if (b20 (tb) < 0.L)
    {
      o = 40.L - b20 (tb);
      goto end;
    }
  o = 30.L + fmaxl (1.L, t1 (tb));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_objective_tb_2_2: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_2_2: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 2 steps 2nd order, 3rd 
 * order in equations depending only in time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_2_2t (RK * rk)  ///< RK struct.
{
  long double o;
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_objective_tb_2_2t: start\n");
#endif
  o = 31.L;
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
#if DEBUG_RK_2_2
  fprintf (stderr, "rk_objective_tb_2_2t: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_2_2t: end\n");
#endif
  return o;
}
