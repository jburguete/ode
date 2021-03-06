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

#define DEBUG_RK_6_3 0          ///< macro to debug.

/**
 * Function to obtain the coefficients of a 6 steps 3rd order Runge-Kutta 
 * method.
 */
int
rk_tb_6_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *tb, *r;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3: start\n");
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
  b65 (tb) = r[15];
  b64 (tb) = r[16];
  A[0] = t1 (tb);
  B[0] = t2 (tb);
  C[0] = t3 (tb);
  D[0] = 0.5L - b64 (tb) * t4 (tb) - b65 (tb) * t5 (tb);
  A[1] = A[0] * t1 (tb);
  B[1] = B[0] * t2 (tb);
  C[1] = C[0] * t3 (tb);
  D[1] = 1.L / 3.L - b64 (tb) * sqr (t4 (tb)) - b65 (tb) * sqr (t5 (tb));
  A[2] = 0.L;
  B[2] = b21 (tb) * t1 (tb);
  C[2] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[2] = 1.L / 6.L - b64 (tb) * (b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb)
                                 + b43 (tb) * t3 (tb))
    - b65 (tb) * (b51 (tb) * t1 (tb) + b52 (tb) * t2 (tb)
                  + b53 (tb) * t3 (tb) + b54 (tb) * t4 (tb));
  solve_3 (A, B, C, D);
  if (isnan (D[0]) || isnan (D[1]) || isnan (D[2]))
    return 0;
  b63 (tb) = D[2];
  b62 (tb) = D[1];
  b61 (tb) = D[0];
  rk_b_6 (tb);
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 6 steps 3rd order, 4th order in
 * equations depending only in time, Runge-Kutta method.
 */
int
rk_tb_6_3t (Optimize * optimize)        ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *tb, *r;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3t: start\n");
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
  b65 (tb) = r[15];
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
  B[3] = b21 (tb) * t1 (tb);
  C[3] = b31 (tb) * t1 (tb) + b32 (tb) * t2 (tb);
  D[3] = b41 (tb) * t1 (tb) + b42 (tb) * t2 (tb) + b43 (tb) * t3 (tb);
  E[3] = 1.L / 6.L - b65 (tb) * (b51 (tb) * t1 (tb) + b52 (tb) * t2 (tb)
                                 + b53 (tb) * t3 (tb) + b54 (tb) * t4 (tb));
  solve_4 (A, B, C, D, E);
  if (isnan (E[0]) || isnan (E[1]) || isnan (E[2]) || isnan (E[3]))
    return 0;
  b64 (tb) = E[3];
  b63 (tb) = E[2];
  b62 (tb) = E[1];
  b61 (tb) = E[0];
  rk_b_6 (tb);
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3t: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 6 steps 2nd-3rd order Runge-Kutta 
 * pair.
 */
int
rk_tb_6_3p (Optimize * optimize)        ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3p: start\n");
#endif
  if (!rk_tb_6_3 (optimize))
    return 0;
  tb = optimize->coefficient;
  e51 (tb) = 0.5L / t1 (tb);
  e52 (tb) = e53 (tb) = 0.L;
  rk_e_6 (tb);
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3p: end\n");
#endif
  return 1;
}

/**
 * Function to obtain the coefficients of a 6 steps 2nd-3rd order, 3rd-4th order
 * in equations depending only in time, Runge-Kutta pair.
 */
int
rk_tb_6_3tp (Optimize * optimize)       ///< Optimize struct.
{
  long double *tb;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3tp: start\n");
#endif
  if (!rk_tb_6_3t (optimize))
    return 0;
  tb = optimize->coefficient;
  e63 (tb) = e64 (tb) = 0.L;
  e62 (tb) = (1.L / 3.L - 0.5L * t1 (tb)) / (t2 (tb) * (t2 (tb) - t1 (tb)));
  if (isnan (e62 (tb)))
    return 0;
  e61 (tb) = (0.5L - e52 (tb) * t2 (tb)) / t1 (tb);
  if (isnan (e61 (tb)))
    return 0;
  rk_e_6 (tb);
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_tb_6_3tp: end\n");
#endif
  return 1;
}

/**
 * Function to calculate the objective function of a 6 steps 3rd order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_3 (RK * rk)   ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3: start\n");
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
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_3: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 6 steps 3rd order, 4th 
 * order in equations depending only in time, Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_3t (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3t: start\n");
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
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3t: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_3t: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 6 steps 2nd-3rd order 
 * Runge-Kutta pair.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_3p (RK * rk)  ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3p: start\n");
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
  if (b62 (tb) < 0.L)
    o += b62 (tb);
  if (b63 (tb) < 0.L)
    o += b63 (tb);
  if (b64 (tb) < 0.L)
    o += b64 (tb);
  if (e60 (tb) < 0.L)
    o += e60 (tb);
  if (e61 (tb) < 0.L)
    o += e61 (tb);
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
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3p: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_3p: end\n");
#endif
  return o;
}

/**
 * Function to calculate the objective function of a 6 steps 2nd-3rd order, 
 * 3rd-4th order in equations depending only in time, Runge-Kutta pair.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_6_3tp (RK * rk) ///< RK struct.
{
  long double *tb;
  long double o;
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3tp: start\n");
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
  if (b62 (tb) < 0.L)
    o += b62 (tb);
  if (b63 (tb) < 0.L)
    o += b63 (tb);
  if (b64 (tb) < 0.L)
    o += b64 (tb);
  if (e60 (tb) < 0.L)
    o += e60 (tb);
  if (e61 (tb) < 0.L)
    o += e61 (tb);
  if (e62 (tb) < 0.L)
    o += e62 (tb);
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
#if DEBUG_RK_6_3
  fprintf (stderr, "rk_objective_tb_6_3tp: optimal=%Lg\n", o);
  fprintf (stderr, "rk_objective_tb_6_3tp: end\n");
#endif
  return o;
}
