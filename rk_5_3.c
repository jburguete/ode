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
 * \file rk_5_3.c
 * \brief Source file to optimize Runge-Kutta 5 steps 3rd order methods.
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
#include "rk_5_3.h"

#define DEBUG_RK_5_3 0          ///< macro to debug.

/**
 * Function to obtain the coefficients of a 5 steps 3rd order Runge-Kutta 
 * method.
 */
int
rk_tb_5_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *tb, *r;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3: start\n");
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
  b41 (tb) = r[7];
  b42 (tb) = r[8];
  b43 (tb) = r[9];
  b54 (tb) = r[10];
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = 0.5L - b54 (tb) * t4 (tb);
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = 1.L / 3.L - b54 (tb) * sqr (t4 (tb));
  A[2] = 0.L;
  B[2] = b21 (tb) * t1 (tb);
  C[2] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[2] = 1.L / 6.L - b54 (tb) * (b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb)
                                 + b43 (tb) * t3 (tb));
  solve_3 (A, B, C, D);
  if (isnan (D[0]) || isnan (D[1]) || isnan (D[2]))
    return 0;
  b53 (tb) = D[2];
  b52 (tb) = D[1];
  b51 (tb) = D[0];
  rk_b_5 (tb);
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 5 steps 3rd order, 4th order in
 * equations depending only in time, Runge-Kutta method.
 */
int
rk_tb_5_3t (Optimize * optimize)        ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *tb, *r;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3t: start\n");
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
  b41 (tb) = r[7];
  b42 (tb) = r[8];
  b43 (tb) = r[9];
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
  A[3] = 0.L;
  B[3] = b21 (tb) * t1 (tb);
  C[3] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[3] = b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb) + b43 (tb) * t3 (tb);
  E[3] = 1.L / 6.L;
  solve_4 (A, B, C, D, E);
  if (isnan (E[0]) || isnan (E[1]) || isnan (E[2]) || isnan (E[3]))
    return 0;
  b54 (tb) = E[3];
  b53 (tb) = E[2];
  b52 (tb) = E[1];
  b51 (tb) = E[0];
  rk_b_5 (tb);
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3t: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 5 steps 2nd-3rd order Runge-Kutta 
 * pair.
 */
int
rk_tb_5_3p (Optimize * optimize)        ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3p: start\n");
#endif
  if (!rk_tb_5_3 (optimize))
    return 0;
  tb = optimize->coefficient;
  e51 (tb) = 0.5L / t1 (tb);
  e52 (tb) = e53 (tb) = 0.L;
  rk_e_5 (tb);
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3p: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 5 steps 2nd-3rd order, 3rd-4th order
 * in equations depending only in time, Runge-Kutta pair.
 */
int
rk_tb_5_3tp (Optimize * optimize)       ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3tp: start\n");
#endif
  if (!rk_tb_5_3t (optimize))
    return 0;
  tb = optimize->coefficient;
  e53 (tb) = 0.L;
  e52 (tb) = (1.L / 3.L - 0.5L * t1 (tb)) / (t2 (tb) * (t2 (tb) - t1 (tb)));
  if (isnan (e52 (tb)))
    return 0;
  e51 (tb) = (0.5L - e52 (tb) * t2 (tb)) / t1 (tb);
  if (isnan (e51 (tb)))
    return 0;
  rk_e_5 (tb);
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_tb_5_3tp: end\n");
#endif
  return 1;
}

/**
 * Function to calculate the objective function of a 5 steps 3rd order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_5_3 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b50 (tb) < 0.L)
    o += b50 (tb);
  if (b51 (tb) < 0.L)
    o += b51 (tb);
  if (b52 (tb) < 0.L)
    o += b52 (tb);
  if (b53 (tb) < 0.L)
    o += b53 (tb);
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
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_5_3: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 5 steps 3rd order, 4th 
 * order in equations depending only in time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_5_3t (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3t: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
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
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3t: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_5_3t: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 5 steps 2nd-3rd order 
 * Runge-Kutta pair.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_5_3p (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3p: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
  if (b50 (tb) < 0.L)
    o += b50 (tb);
  if (b51 (tb) < 0.L)
    o += b51 (tb);
  if (b52 (tb) < 0.L)
    o += b52 (tb);
  if (b53 (tb) < 0.L)
    o += b53 (tb);
  if (e50 (tb) < 0.L)
    o += e50 (tb);
  if (e51 (tb) < 0.L)
    o += e51 (tb);
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
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3p: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_5_3p: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 5 steps 2nd-3rd order, 
 * 3rd-4th order in equations depending only in time, Runge-Kutta pair.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_5_3tp (RK * rk) ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3tp: start\n");
#endif
  tb = rk->tb->coefficient;
  o = fminl (0.L, b20 (tb));
  if (b30 (tb) < 0.L)
    o += b30 (tb);
  if (b40 (tb) < 0.L)
    o += b40 (tb);
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
  if (e50 (tb) < 0.L)
    o += e50 (tb);
  if (e51 (tb) < 0.L)
    o += e51 (tb);
  if (e52 (tb) < 0.L)
    o += e52 (tb);
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
#if DEBUG_RK_5_3
  fprintf (stderr, "rk_objective_tb_5_3tp: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_5_3tp: end\n");
#endif
  return o;
}
