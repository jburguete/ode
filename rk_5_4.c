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
 * \file rk_5_4.c
 * \brief Source file to optimize Runge-Kutta 5 steps 4th order methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "rk.h"
#include "rk_5_3.h"
#include "rk_5_4.h"

#define DEBUG_RK_5_4 0          ///< macro to debug.

///> array of minimum freedom degree values for the t-b coefficients of the 5
///> steps 4th order Runge-Kutta method.
const long double minimum_tb_5_4[7] = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of minimum freedom degree intervals for the t-b coefficients of the 5
///> steps 4th order Runge-Kutta method.
const long double interval_tb_5_4[7] = { 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L };

///> array of freedom degree random function types for the t-b coefficients of
///> the 5 steps 4th order Runge-Kutta method.
const unsigned int random_tb_5_4[7] = { 0, 0, 0, 0, 0, 0, 0 };

/**
 * Function to print a maxima format file to check the accuracy order of a 5
 * steps 4th order Runge-Kutta method.
 */
void
tb_print_maxima_5_4 (FILE * file,       ///< file.
                     unsigned int nsteps,       ///< steps number.
                     unsigned int order)        ///< accuracy order.
{
  tb_print_maxima_5_3 (file, nsteps, order);
  fprintf (file, "b53*b32*b21*t1+b54*(b42*b21*t1+b43*(b31*t1+b32*t2))-1/24;\n");
  fprintf (file, "b52*b21*t1^2+b53*(b31*t1^2+b32*t2^2)"
           "+b54*(b41*t1^2+b42*t2^2+b43*t3^2)-1/12;\n");
  fprintf (file, "b52*t2*b21*t1+b53*t3*(b31*t1+b32*t2)"
           "+b54*t4*(b41*t1+b42*t2+b43*t3)-1/8;\n");
  fprintf (file, "b51*t1^3+b52*t2^3+b53*t3^3+b54*t4^3-1/4;\n");
  fprintf (file, "b51*t1^4+b52*t2^4+b53*t3^4+b54*t4^4-1/5;\n");
}

/**
 * Function to print a maxima format file to check the accuracy order of a 5
 * steps 4th order Runge-Kutta method.
 */
void
rk_print_maxima_5_4 (FILE * file,       ///< file.
                     unsigned int nsteps,       ///< steps number.
                     unsigned int order)        ///< accuracy order.
{
  tb_print_maxima_5_4 (file, nsteps, order);
  rk_print_maxima_5 (file);
}

/**
 * Function to obtain the coefficients of a 5 steps 4th order Runge-Kutta 
 * method.
 */
void
rk_tb_5_4 (Optimize * optimize) ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *tb, *r;
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_tb_5_4: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t5 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = r[1];
  b21 (tb) = r[2];
  t3 (tb) = r[3];
  b31 (tb) = r[4];
  b32 (tb) = r[5];
  t4 (tb) = r[6];
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = t4 (tb);
  E[0] = 0.5L;
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = D[0] * t4 (tb);
  E[1] = 1.L / 3.L;
  A[2] = A[1] * t1 (tb);
  B[2] = B[1] * t2 (tb);
  C[2] = C[1] * t3 (tb);
  D[2] = D[1] * t4 (tb);
  E[2] = 0.25L;
  A[3] = D[3] = 0.L;
  B[3] = b21 (tb) * t1 (tb) * (t2 (tb) - t4 (tb));
  C[3] = (b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb)) * (t3 (tb) - t4 (tb));
  E[3] = 0.125L - 1.L / 6.L * t4 (tb);
  solve_4 (A, B, C, D, E);
  b54 (tb) = E[3];
  b53 (tb) = E[2];
  b52 (tb) = E[1];
  b51 (tb) = E[0];
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = 1.L / 6.L - b52 (tb) * b21 (tb) * t1 (tb)
    - b53 (tb) * (b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb));
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = 1.L / 12.L - b52 (tb) * b21 (tb) * sqr (t1 (tb))
    - b53 (tb) * (b31 (tb) * sqr (t1 (tb)) + b32 (tb) * sqr (t2 (tb)));
  A[2] = 0.L;
  B[2] = b21 (tb) * t1 (tb);
  C[2] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[2] = 1.L / 24.L - b53 (tb) * b32 (tb) * b21 (tb) * t1 (tb);
  solve_3 (A, B, C, D);
  b43 (tb) = D[2] / b54 (tb);
  b42 (tb) = D[1] / b54 (tb);
  b41 (tb) = D[0] / b54 (tb);
  rk_b_5 (tb);
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_tb_5_4: end\n");
#endif
}

/**
 * Function to obtain the coefficients of a 5 steps 4th order, 5th order in
 * equations depending only on time, Runge-Kutta method.
 */
void
rk_tb_5_4t (Optimize * optimize) ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *tb, *r;
#if RK_PAIR
	long double e52, e53;
#endif
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_tb_5_4t: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t5 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = r[1];
  b21 (tb) = r[2];
  t3 (tb) = r[3];
  t4 (tb) = r[4];
#if !RK_PAIR
	b31 (tb) = r[5];
#endif
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = t4 (tb);
  E[0] = 0.5L;
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = D[0] * t4 (tb);
  E[1] = 1.L / 3.L;
  A[2] = A[1] * t1 (tb);
  B[2] = B[1] * t2 (tb);
  C[2] = C[1] * t3 (tb);
  D[2] = D[1] * t4 (tb);
  E[2] = 0.25L;
  A[3] = A[2] * t1 (tb);
  B[3] = B[2] * t2 (tb);
  C[3] = C[2] * t3 (tb);
  D[3] = D[2] * t4 (tb);
  E[3] = 0.2L;
  solve_4 (A, B, C, D, E);
  b54 (tb) = E[3];
  b53 (tb) = E[2];
  b52 (tb) = E[1];
  b51 (tb) = E[0];
	b32 (tb) = (1.L / 6.L * t4 (tb) - 0.125L 
			        - t1 (tb) * (b52 (tb) * b21 (tb) * (t4 (tb) - t2 (tb))
			                     + b53 (tb) * b31 (tb) * (t4 (tb) - t3 (tb))))
		/ (b53 (tb) * t2 (tb) * (t4 (tb) - t3 (tb)));
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = 1.L / 6.L - b52 (tb) * b21 (tb) * t1 (tb)
    - b53 (tb) * (b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb));
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = 1.L / 12.L - b52 (tb) * b21 (tb) * sqr (t1 (tb))
    - b53 (tb) * (b31 (tb) * sqr (t1 (tb)) + b32 (tb) * sqr (t2 (tb)));
  A[2] = 0.L;
  B[2] = b21 (tb) * t1 (tb);
  C[2] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[2] = 1.L / 24.L - b53 (tb) * b32 (tb) * b21 (tb) * t1 (tb);
  solve_3 (A, B, C, D);
  b43 (tb) = D[2] / b54 (tb);
  b42 (tb) = D[1] / b54 (tb);
  b41 (tb) = D[0] / b54 (tb);
#if RK_PAIR
  e53 = (0.25L - 1.L / 3.L * t1 (tb) - (1.L / 3.L - 0.5L * t1 (tb) )* t2 (tb))
		/ (t3 (tb) * (t3 (tb) - t2 (tb)) * (t2 (tb) - t1 (tb)));
  e52 = (1.L / 3.L - 0.5L * t1 (tb) - t3 (tb) * (t3 (tb) - t1 (tb)) * e53)
		/ (t2 (tb) * (t2 (tb) - t1 (tb)));
	b31 (tb) = ((1.L / 6.L - e52 * b21 (tb) * t1 (tb)) / e53 - b32 (tb) * t2 (tb))
		/ t1 (tb);
#endif
  rk_b_5 (tb);
#if DEBUG_RK_5_4
	rk_print_tb_5 (tb, "rk_tb_5_4t", stderr);
  fprintf (stderr, "rk_tb_5_4t: end\n");
#endif
}

/**
 * Function to calculate the objective function of a 5 steps 4th order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_5_4 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_objective_tb_5_4: start\n");
#endif
  tb = rk->tb->coefficient;
  if (isnan (b41 (tb)) || isnan (b42 (tb)) || isnan (b43 (tb))
      || isnan (b51 (tb)) || isnan (b52 (tb)) || isnan (b53 (tb))
      || isnan (b54 (tb)))
    {
      o = INFINITY;
      goto end;
    }
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b41 (tb) < 0.L)
    o += b41 (tb);
  if (b42 (tb) < 0.L)
    o += b42 (tb);
  if (b43 (tb) < 0.L)
    o += b43 (tb);
  if (b50 (tb) < 0.L)
    o += b50 (tb);
  if (b51 (tb) < 0.L)
    o += b51 (tb);
  if (b52 (tb) < 0.L)
    o += b52 (tb);
  if (b53 (tb) < 0.L)
    o += b53 (tb);
  if (b54 (tb) < 0.L)
    o += b54 (tb);
  if (o < 0.L)
    {
      o = 40.L - o;
      goto end;
    }
  o = 30.L
    + fmaxl (1.L, fmaxl (t1 (tb), fmaxl (t2 (tb), fmaxl (t3 (tb), t4 (tb)))));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_objective_tb_5_4: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_5_4: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 5 steps 4th order, 5th
 * order in equations depending only on time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_5_4t (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_objective_tb_5_4t: start\n");
#endif
  tb = rk->tb->coefficient;
  if (isnan (b31 (tb)) || isnan (b32 (tb)) || isnan (b41 (tb)) || 
			isnan (b42 (tb)) || isnan (b43 (tb)) || isnan (b51 (tb)) || 
			isnan (b52 (tb)) || isnan (b53 (tb)) || isnan (b54 (tb)))
    {
      o = INFINITY;
      goto end;
    }
  o = fminl (0.L, b20 (tb));
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
  if (b50 (tb) < 0.L)
    o += b50 (tb);
  if (b51 (tb) < 0.L)
    o += b51 (tb);
  if (b52 (tb) < 0.L)
    o += b52 (tb);
  if (b53 (tb) < 0.L)
    o += b53 (tb);
  if (b54 (tb) < 0.L)
    o += b54 (tb);
  if (o < 0.L)
    {
      o = 40.L - o;
      goto end;
    }
  o = 30.L
    + fmaxl (1.L, fmaxl (t1 (tb), fmaxl (t2 (tb), fmaxl (t3 (tb), t4 (tb)))));
  if (rk->strong)
    {
      rk_bucle_ac (rk);
      o = fminl (o, *rk->ac0->optimal);
    }
end:
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_objective_tb_5_4t: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_5_4t: end\n");
#endif
  return o;
}
