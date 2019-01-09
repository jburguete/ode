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
 * \file rk_6_4.c
 * \brief Source file to optimize Runge-Kutta 6 steps 4th order methods.
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
#include "rk_6_3.h"
#include "rk_6_4.h"

#define DEBUG_RK_6_4 0          ///< macro to debug.

/**
 * Function to obtain the coefficients of a 6 steps 4th order Runge-Kutta 
 * method.
 */
int
rk_tb_6_4 (Optimize * optimize) ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *tb, *r;
#if DEBUG_RK_6_4
  fprintf (stderr, "rk_tb_6_4: start\n");
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
  b54 (tb) = r[11];
  b65 (tb) = r[12];
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = t4 (tb);
  E[0] = 0.5L - b65 (tb) * t5 (tb);
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = D[0] * t4 (tb);
  E[1] = 1.L / 3.L - b65 (tb) * sqr (t5 (tb));
  A[2] = A[1] * t1 (tb);
  B[2] = B[1] * t2 (tb);
  C[2] = C[1] * t3 (tb);
  D[2] = D[1] * t4 (tb);
  E[2] = 0.25L - b65 (tb) * sqr (t5 (tb)) * t5 (tb);
  A[3] = 0.L;
  B[3] = b21 (tb) * t1 (tb) * (t2 (tb) - t5 (tb));
  C[3] = (b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb)) * (t3 (tb) - t5 (tb));
  D[3] = (b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb) + b43 (tb) * t3 (tb))
    * (t4 (tb) - t5 (tb));
  E[3] = 0.125L - 1.L / 6.L * t5 (tb);
  solve_4 (A, B, C, D, E);
  b64 (tb) = E[3];
  b63 (tb) = E[2];
  b62 (tb) = E[1];
  b61 (tb) = E[0];
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = (1.L / 6.L - b62 (tb) * b21 (tb) * t1 (tb)
          - b63 (tb) * (b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb))
          - b64 (tb) * (b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb)
                        + b43 (tb) * t3 (tb))) / b65 (tb) - b54 (tb) * t4 (tb);
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = (1.L / 12.L - b62 (tb) * b21 (tb) * sqr (t1 (tb))
          - b63 (tb) * (b31 (tb) * sqr (t1 (tb)) + b32 (tb) * sqr (t2 (tb)))
          - b64 (tb) * (b41 (tb) * sqr (t1 (tb)) + b42 (tb) * sqr (t2 (tb))
                        + b43 (tb) * sqr (t3 (tb)))) / b65 (tb)
    - b54 (tb) * sqr (t4 (tb));
  A[2] = 0.L;
  B[2] = b21 (tb) * t1 (tb);
  C[2] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[2] = (1.L / 24.L - b63 (tb) * b32 (tb) * b21 (tb) * t1 (tb)
          - b64 (tb) * (b42 (tb) * b21 (tb) * t1 (tb)
                        + b43 (tb) * (b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb))))
    / b65 (tb) - b54 (tb) * (b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb)
                             + b43 (tb) * t3 (tb));
  solve_3 (A, B, C, D);
  b53 (tb) = D[2];
  b52 (tb) = D[1];
  b51 (tb) = D[0];
  rk_b_6 (tb);
#if DEBUG_RK_6_4
  fprintf (stderr, "rk_tb_6_4: end\n");
#endif
  if (isnan (b51 (tb)) || isnan (b52 (tb)) || isnan (b53 (tb))
      || isnan (b61 (tb)) || isnan (b62 (tb)) || isnan (b63 (tb))
      || isnan (b64 (tb)))
    return 0;
  return 1;
}

/**
 * Function to obtain the coefficients of a 6 steps 4th order, 5th order in
 * equations depending only on time, Runge-Kutta method.
 */
int
rk_tb_6_4t (Optimize * optimize)        ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *tb, *r;
#if RK_PAIR
  long double e52, e53;
#endif
#if DEBUG_RK_6_4
  fprintf (stderr, "rk_tb_6_4t: start\n");
#endif
  tb = optimize->coefficient;
  r = optimize->random_data;
  t6 (tb) = 1.L;
  t1 (tb) = r[0];
  t2 (tb) = r[1];
  t3 (tb) = r[2];
  t4 (tb) = r[3];
  t5 (tb) = r[4];
  b31 (tb) = r[4];
#if !RK_PAIR
  b21 (tb) = r[5];
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
#if RK_PAIR
  e53 = (0.25L - 1.L / 3.L * t1 (tb) - (1.L / 3.L - 0.5L * t1 (tb)) * t2 (tb))
    / (t3 (tb) * (t3 (tb) - t2 (tb)) * (t3 (tb) - t1 (tb)));
  e52 = (1.L / 3.L - 0.5L * t1 (tb) - t3 (tb) * (t3 (tb) - t1 (tb)) * e53)
    / (t2 (tb) * (t2 (tb) - t1 (tb)));
  b21 (tb) = (1.L / 6.L * b53 (tb) * (t4 (tb) - t3 (tb))
              + e53 * (0.125L - 1.L / 6.L * t4 (tb)))
    / (t1 (tb) * (e52 * b53 (tb) * (t4 (tb) - t3 (tb))
                  - e53 * b52 (tb) * (t4 (tb) - t2 (tb))));
#endif
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
  rk_b_5 (tb);
#if DEBUG_RK_6_4
  rk_print_tb (optimize, "rk_tb_6_4t", stderr);
  fprintf (stderr, "rk_tb_6_4t: end\n");
#endif
#if RK_PAIR
  if (isnan (e51) || isnan (e52) || isnan (e53))
    return 0;
#endif
  if (isnan (b41 (tb)) || isnan (b42 (tb)) || isnan (b43 (tb)) ||
      isnan (b32 (tb)) || isnan (b21 (tb)) || isnan (b51 (tb)) ||
      isnan (b52 (tb)) || isnan (b53 (tb)) || isnan (b54 (tb)))
    return 0;
  return 1;
}

/**
 * Function to calculate the objective function of a 6 steps 4th order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_4 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_objective_tb_5_4: start\n");
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
  if (b60 (tb) < 0.L)
    o += b60 (tb);
  if (b61 (tb) < 0.L)
    o += b61 (tb);
  if (b62 (tb) < 0.L)
    o += b62 (tb);
  if (b63 (tb) < 0.L)
    o += b63 (tb);
  if (b64 (tb) < 0.L)
    o += b64 (tb);
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
#if DEBUG_RK_6_4
  fprintf (stderr, "rk_objective_tb_6_4: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_4: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 6 steps 4th order, 5th
 * order in equations depending only on time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_4t (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if RK_PAIR
  long double e50, e51, e52, e53, e54;
#endif
#if DEBUG_RK_5_4
  fprintf (stderr, "rk_objective_tb_5_4t: start\n");
#endif
  tb = rk->tb->coefficient;
#if DEBUG_RK_5_4
  rk_print_tb (optimize, "rk_objective_tb_5_4t", stderr);
#endif
  o = fminl (0.L, b20 (tb));
  if (b21 (tb) < 0.L)
    o += b21 (tb);
  if (b30 (tb) < 0.L)
    o += b30 (tb);
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
#if RK_PAIR
  e53 = (0.25L - 1.L / 3.L * t1 (tb) - (1.L / 3.L - 0.5L * t1 (tb)) * t2 (tb))
    / (t3 (tb) * (t3 (tb) - t2 (tb)) * (t3 (tb) - t1 (tb)));
  e52 = (1.L / 3.L - 0.5L * t1 (tb) - t3 (tb) * (t3 (tb) - t1 (tb)) * e53)
    / (t2 (tb) * (t2 (tb) - t1 (tb)));
  e51 = (0.5L - t2 (tb) * e52 - t3 (tb) * e53) / t1 (tb);
  e50 = 1.L - e51 - e52 - e53;
  if (e50 < 0.L)
    o += e50;
  if (e51 < 0.L)
    o += e51;
  if (e52 < 0.L)
    o += e52;
  if (e53 < 0.L)
    o += e53;
  if (e54 < 0.L)
    o += e54;
#endif
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
