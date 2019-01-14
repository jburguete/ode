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
 * \file rk.c
 * \brief Source file with common variables and functions to optimize
 *   Runge-Kutta methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2019.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <libxml/parser.h>
#include <glib.h>
#include <libintl.h>
#include <gsl/gsl_rng.h>
#if HAVE_MPI
#include <mpi.h>
#endif
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "rk.h"
#include "rk_2_2.h"
#include "rk_3_2.h"
#include "rk_3_3.h"
#include "rk_4_2.h"
#include "rk_4_3.h"
#include "rk_4_4.h"
#include "rk_5_2.h"
#include "rk_5_3.h"
#include "rk_5_4.h"
#include "rk_6_2.h"
#include "rk_6_3.h"
#include "rk_6_4.h"

#define DEBUG_RK 0              ///< macro to debug.

/**
 * Function to print the t-b Runge-Kutta coefficients.
 */
void
rk_print_tb (Optimize * tb,     ///< Optimize struct.
             char *label,       ///< label.
             FILE * file)       ///< file.
{
  long double *x;
  unsigned int i, j, k;
  x = tb->coefficient;
  fprintf (file, "%s: t1=%.19Le\n", label, x[0]);
  for (i = 2, k = 0; i <= tb->nsteps; ++i)
    {
      fprintf (file, "%s: t%u=%.19Le\n", label, i, x[++k]);
      for (j = 0; j < i; ++j)
        fprintf (file, "%s: b%u%u=%.19Le\n", label, i, j, x[++k]);
    }
}

/**
 * Function to print the e Runge-Kutta coefficients.
 */
void
rk_print_e (Optimize * tb,      ///< Optimize struct.
            char *label,        ///< label.
            FILE * file)        ///< file.
{
  long double *x;
  unsigned int i, k, nsteps;
  x = tb->coefficient;
  nsteps = tb->nsteps;
  k = (nsteps + 2) * (nsteps + 1) / 2 - 2;
  for (i = 0; i < nsteps - 1; ++i)
    fprintf (file, "%s: e%u%u=%.19Le\n", label, nsteps, i, x[k++]);
}

/**
 * Function to print in a maxima file the Runge-Kutta coefficients.
 */
static void
rk_print (RK * rk,              ///< RK struct.
          FILE * file)          ///< file.
{
  Optimize *tb, *ac;
  long double *x, *y;
  unsigned int i, j, k, l, nsteps;
  tb = rk->tb;
  ac = rk->ac;
  x = tb->coefficient;
  y = ac->coefficient;
  fprintf (file, "t1:%.19Le;\n", x[0]);
  nsteps = tb->nsteps;
  for (i = 2, k = l = 0; i <= nsteps; ++i)
    {
      fprintf (file, "t%u:%.19Le;\n", i, x[++k]);
      for (j = 0; j < i; ++j)
        fprintf (file, "b%u%u:%.19Le;\n", i, j, x[++k]);
      if (!rk->strong)
        continue;
      for (j = 0; j < i; ++j)
        fprintf (file, "a%u%u:%.19Le;\n", i, j, y[l++]);
      for (j = 0; j < i; ++j)
        fprintf (file, "c%u%u:%.19Le;\n", i, j, y[l++]);
    }
  if (rk->pair)
    for (i = 0; i < nsteps - 1; ++i)
      fprintf (file, "e%u%u:%.19Le;\n", nsteps, i, x[++k]);
}

/**
 * Function to print a maxima format file to check the accuracy order of the
 * Runge-Kutta simple stable methods.
 */
static void
rk_print_maxima (FILE * file,   ///< file.
                 unsigned int nsteps,   ///< steps number.
                 unsigned int ncoefficients,    ///< coefficients number.
                 unsigned int order,    ///< accuracy order.
                 char label)    ///< coefficient label.
{
  unsigned int i, j, k, l;
  // b_{ij}=1 (1st order)
  for (i = 0; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u+", label, nsteps, i);
  fprintf (file, "-1;\n");
  // b_{ij}t_j=1/2 (2nd order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u+", label, nsteps, i, i);
  fprintf (file, "-1/2;\n");
  if (order < 2)
    return;
  // b_{ij}t_j^2=1/3 (3rd order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^2+", label, nsteps, i, i);
  fprintf (file, "-1/3;\n");
  if (order < 3)
    return;
  // b_{ij}b_{jk}t_k=1/6 (3rd order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/6;\n");
  // b_{ij}t_j^3=1/4 (4th order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^3+", label, nsteps, i, i);
  fprintf (file, "-1/4;\n");
  if (order < 4)
    return;
  // b_{ij}b_{jk}b_{kl}t_l=1/24 (4th order)
  for (i = 3; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 2; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 1; k < j; ++k)
            fprintf (file, "b%u%u*t%u+", j, k, k);
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/24;\n");
  // b_{ij}b_{jk}t_k^2=1/12 (4th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u^2+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/12;\n");
  // b_{ij}t_jb_{jk}t_k=1/8 (4th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u*(", label, nsteps, i, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/8;\n");
  // b_{ij}t_j^4=1/5 (5th order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^4+", label, nsteps, i, i);
  fprintf (file, "-1/5;\n");
  if (order < 5)
    return;
  // b_{ij}t_j^2b_{jk}t_k^2=1/10 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u^2*(", label, nsteps, i, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/10;\n");
  // b_{ij}b_{jk}t_k^3=1/20 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u^3+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/20;\n");
  // b_{ij}(b_{jk}t_k)^2=1/20 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)^2+");
    }
  fprintf (file, "-1/20;\n");
  // b_{ij}t_jb_{jk}t_k^2=1/15 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u*(", label, nsteps, i, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u^2+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/8;\n");
  // b_{ij}t_jb_{jk}b_{kl}t_l=1/24 (5th order)
  for (i = 3; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u*(", label, nsteps, i, i);
      for (j = 2; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 1; k < j; ++k)
            fprintf (file, "b%u%u*t%u+", j, k, k);
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-7/120;\n");
  // b_{ij}b_{jk}b_{kl}t_l^2=1/60 (5th order)
  for (i = 3; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 2; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 1; k < j; ++k)
            fprintf (file, "b%u%u*t%u^2+", j, k, k);
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/60;\n");
  // b_{ij}b_{jk}b_{kl}b_{lm}t_m=1/120 (5th order)
  for (i = 4; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 3; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 2; k < j; ++k)
            {
              fprintf (file, "b%u%u*(", j, k);
              for (l = 1; l < k; ++l)
                fprintf (file, "b%u%u*t%u+", k, l, l);
              fprintf (file, "0)+");
            }
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/120;\n");
  // b_{ij}t_j^5=1/6 (6th order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^5+", label, nsteps, i, i);
  fprintf (file, "-1/6;\n");
}

/**
 * Function to print in a maxima file the a-c Runge-Kutta equations.
 */
static void
ac_print_maxima (FILE * file,   ///< file.
                 unsigned int nsteps)   ///< steps number.
{
  unsigned int i, j, k;
  for (i = 2; i <= nsteps; ++i)
    {
      // a_{i,j}=1
      for (j = 0; j < i; ++j)
        fprintf (file, "a%u%u+", i, j);
      fprintf (file, "-1;\n");
      // a_{i0}c_{i0}+a_{ij}b_{j0}=b_{i0}
      fprintf (file, "a%u0*c%u0+a%u1*t1+", i, i, i);
      for (j = 2; j < i; ++j)
        fprintf (file, "a%u%u*b%u0+", i, j, j);
      fprintf (file, "-b%u0;\n", i);
      // a_{ij}c_{ij}+a_{ik}b_{kj}=b_{ij}
      for (j = 1; j < i; ++j)
        {
          fprintf (file, "a%u%u*c%u%u+", i, j, i, j);
          for (k = j; ++k < i;)
            fprintf (file, "a%u%u*b%u%u+", i, k, k, j);
          fprintf (file, "-b%u%u;\n", i, j);
        }
    }
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 2nd step.
 */
static int
rk_ac_2 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_2: start\n");
#endif
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  r = rk->ac->random_data;
  c21 (ac) = r[0];
  a21 (ac) = b21 (tb) / c21 (ac);
  a20 (ac) = 1.L - a21 (ac);
  ac0 = b20 (tb) - a21 (ac) * t1 (tb);
  c20 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / a20 (ac);
#if DEBUG_RK
  fprintf (stderr, "rk_ac_2: a20=%Lg c20=%Lg\n", a20 (ac), c20 (ac));
  fprintf (stderr, "rk_ac_2: a21=%Lg c21=%Lg\n", a21 (ac), c21 (ac));
  fprintf (stderr, "rk_ac_2: end\n");
#endif
  if (isnan (c20 (ac)) || isnan (a21 (ac)))
    return 0;
  return 1;
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 3rd step.
 */
static int
rk_ac_3 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_3: start\n");
#endif
  if (!rk_ac_2 (rk))
    {
#if DEBUG_RK
      fprintf (stderr, "rk_ac_3: end\n");
#endif
      return 0;
    }
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  r = rk->ac->random_data;
  c31 (ac) = r[1];
  c32 (ac) = r[2];
  a32 (ac) = b32 (tb) / c32 (ac);
  ac0 = b31 (tb) - a32 (ac) * b21 (tb);
  a31 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c31 (ac);
  a30 (ac) = 1.L - a31 (ac) - a32 (ac);
  ac0 = b30 (tb) - a31 (ac) * t1 (tb) - a32 (ac) * b20 (tb);
  c30 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / a30 (ac);
#if DEBUG_RK
  fprintf (stderr, "rk_ac_3: a30=%Lg c30=%Lg\n", a30 (ac), c30 (ac));
  fprintf (stderr, "rk_ac_3: a31=%Lg c31=%Lg\n", a31 (ac), c31 (ac));
  fprintf (stderr, "rk_ac_3: a32=%Lg c32=%Lg\n", a32 (ac), c32 (ac));
  fprintf (stderr, "rk_ac_3: end\n");
#endif
  if (isnan (c30 (ac)) || isnan (a31 (ac)) || isnan (a32 (ac)))
    return 0;
  return 1;
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 4th step.
 */
static int
rk_ac_4 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_4: start\n");
#endif
  if (!rk_ac_3 (rk))
    {
#if DEBUG_RK
      fprintf (stderr, "rk_ac_4: end\n");
#endif
      return 0;
    }
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  r = rk->ac->random_data;
  c41 (ac) = r[3];
  c42 (ac) = r[4];
  c43 (ac) = r[5];
  a43 (ac) = b43 (tb) / c43 (ac);
  ac0 = b42 (tb) - a43 (ac) * b32 (tb);
  a42 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c42 (ac);
  ac0 = b41 (tb) - a42 (ac) * b21 (tb) - a43 (ac) * b31 (tb);
  a41 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c41 (ac);
  a40 (ac) = 1.L - a41 (ac) - a42 (ac) - a43 (ac);
  ac0 = b40 (tb) - a41 (ac) * t1 (tb) - a42 (ac) * b20 (tb)
    - a43 (ac) * b30 (tb);
  c40 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / a40 (ac);
#if DEBUG_RK
  fprintf (stderr, "rk_ac_4: a40=%Lg c40=%Lg\n", a40 (ac), c40 (ac));
  fprintf (stderr, "rk_ac_4: a41=%Lg c41=%Lg\n", a41 (ac), c41 (ac));
  fprintf (stderr, "rk_ac_4: a42=%Lg c42=%Lg\n", a42 (ac), c42 (ac));
  fprintf (stderr, "rk_ac_4: a43=%Lg c43=%Lg\n", a43 (ac), c43 (ac));
  fprintf (stderr, "rk_ac_4: end\n");
#endif
  if (isnan (c40 (ac)) || isnan (a41 (ac)) || isnan (a42 (ac))
      || isnan (a43 (ac)))
    return 0;
  return 1;
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 5th step.
 */
static int
rk_ac_5 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_5: start\n");
#endif
  if (!rk_ac_4 (rk))
    {
#if DEBUG_RK
      fprintf (stderr, "rk_ac_5: end\n");
#endif
      return 0;
    }
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  r = rk->ac->random_data;
  c51 (ac) = r[6];
  c52 (ac) = r[7];
  c53 (ac) = r[8];
  c54 (ac) = r[9];
  a54 (ac) = b54 (tb) / c54 (ac);
  ac0 = b53 (tb) - a54 (ac) * b43 (tb);
  a53 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c53 (ac);
  ac0 = b52 (tb) - a53 (ac) * b32 (tb) - a54 (ac) * b42 (tb);
  a52 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c52 (ac);
  ac0 = b51 (tb) - a52 (ac) * b21 (tb) - a53 (ac) * b31 (tb)
    - a54 (ac) * b41 (tb);
  a51 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c51 (ac);
  a50 (ac) = 1.L - a51 (ac) - a52 (ac) - a53 (ac) - a54 (ac);
  ac0 = b50 (tb) - a51 (ac) * t1 (tb) - a52 (ac) * b20 (tb)
    - a53 (ac) * b30 (tb) - a54 (ac) * b40 (tb);
  c50 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / a50 (ac);
#if DEBUG_RK
  fprintf (stderr, "rk_ac_5: a50=%Lg c50=%Lg\n", a50 (ac), c50 (ac));
  fprintf (stderr, "rk_ac_5: a51=%Lg c51=%Lg\n", a51 (ac), c51 (ac));
  fprintf (stderr, "rk_ac_5: a52=%Lg c52=%Lg\n", a52 (ac), c52 (ac));
  fprintf (stderr, "rk_ac_5: a53=%Lg c53=%Lg\n", a53 (ac), c53 (ac));
  fprintf (stderr, "rk_ac_5: a54=%Lg c54=%Lg\n", a54 (ac), c54 (ac));
  fprintf (stderr, "rk_ac_5: end\n");
#endif
  if (isnan (c50 (ac)) || isnan (a51 (ac)) || isnan (a52 (ac))
      || isnan (a53 (ac)) || isnan (a54 (ac)))
    return 0;
  return 1;
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 6th step.
 */
static int
rk_ac_6 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_6: start\n");
#endif
  if (!rk_ac_5 (rk))
    {
#if DEBUG_RK
      fprintf (stderr, "rk_ac_6: end\n");
#endif
      return 0;
    }
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  r = rk->ac->random_data;
  c61 (ac) = r[10];
  c62 (ac) = r[11];
  c63 (ac) = r[12];
  c64 (ac) = r[13];
  c65 (ac) = r[14];
  a65 (ac) = b65 (tb) / c65 (ac);
  ac0 = b64 (tb) - a65 (ac) * b54 (tb);
  a64 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c64 (ac);
  ac0 = b63 (tb) - a64 (ac) * b43 (tb) - a65 (ac) * b53 (tb);
  a63 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c63 (ac);
  ac0 = b62 (tb) - a63 (ac) * b32 (tb) - a64 (ac) * b42 (tb)
    - a65 (ac) * b52 (tb);
  a62 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c62 (ac);
  ac0 = b61 (tb) - a62 (ac) * b21 (tb) - a63 (ac) * b31 (tb)
    - a64 (ac) * b41 (tb) - a65 (ac) * b51 (tb);
  a61 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / c61 (ac);
  a60 (ac) = 1.L - a61 (ac) - a62 (ac) - a63 (ac) - a64 (ac) - a65 (ac);
  ac0 = b60 (tb) - a61 (ac) * t1 (tb) - a62 (ac) * b20 (tb)
    - a63 (ac) * b30 (tb) - a64 (ac) * b40 (tb) - a65 (ac) * b50 (tb);
  c60 (ac) = fabsl (ac0) < LDBL_EPSILON ? 0.L : ac0 / a60 (ac);
#if DEBUG_RK
  fprintf (stderr, "rk_ac_6: a60=%Lg c60=%Lg\n", a60 (ac), c60 (ac));
  fprintf (stderr, "rk_ac_6: a61=%Lg c61=%Lg\n", a61 (ac), c61 (ac));
  fprintf (stderr, "rk_ac_6: a62=%Lg c62=%Lg\n", a62 (ac), c62 (ac));
  fprintf (stderr, "rk_ac_6: a63=%Lg c63=%Lg\n", a63 (ac), c63 (ac));
  fprintf (stderr, "rk_ac_6: a64=%Lg c64=%Lg\n", a64 (ac), c64 (ac));
  fprintf (stderr, "rk_ac_6: a65=%Lg c65=%Lg\n", a65 (ac), c65 (ac));
  fprintf (stderr, "rk_ac_6: end\n");
#endif
  if (isnan (c60 (ac)) || isnan (a61 (ac)) || isnan (a62 (ac))
      || isnan (a63 (ac)) || isnan (a64 (ac)) || isnan (a65 (ac)))
    return 0;
  return 1;
}

/**
 * Function to get the objective function of 2 steps Runge-Kutta methods.
 *
 * \return objective function value.
 */
static long double
rk_objective_ac_2 (RK * rk)     ///< RK struct.
{
  long double *tb, *ac;
  long double k;
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_2: start\n");
#endif
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  k = fminl (0.L, a20 (ac));
  if (a21 (ac) < 0.L)
    k += a21 (ac);
  if (k < 0.L)
    {
      k = 20.L - k;
      goto end;
    }
  k = fminl (0.L, c20 (ac));
  if (c21 (ac) < 0.L)
    k += c21 (ac);
  if (k < 0.L)
    {
      k = 10.L - k;
      goto end;
    }
  k = 1.L / rk_cfl_2 (tb, ac);
end:
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_2: objective=%Lg\n", k);
  fprintf (stderr, "rk_objective ac_2: end\n");
#endif
  return k;
}

/**
 * Function to get the objective function of 3 steps Runge-Kutta methods.
 *
 * \return objective function value.
 */
static long double
rk_objective_ac_3 (RK * rk)     ///< RK struct.
{
  long double *tb, *ac;
  long double k;
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_3: start\n");
#endif
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  k = fminl (0.L, a20 (ac));
  if (a21 (ac) < 0.L)
    k += a21 (ac);
  if (a30 (ac) < 0.L)
    k += a30 (ac);
  if (a31 (ac) < 0.L)
    k += a31 (ac);
  if (a32 (ac) < 0.L)
    k += a32 (ac);
  if (k < 0.L)
    {
      k = 20.L - k;
      goto end;
    }
  k = fminl (0.L, c20 (ac));
  if (c21 (ac) < 0.L)
    k += c21 (ac);
  if (c30 (ac) < 0.L)
    k += c30 (ac);
  if (c31 (ac) < 0.L)
    k += c31 (ac);
  if (c32 (ac) < 0.L)
    k += c32 (ac);
  if (k < 0.L)
    {
      k = 10.L - k;
      goto end;
    }
  k = 1.L / rk_cfl_3 (tb, ac);
end:
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_3: objective=%Lg\n", k);
  fprintf (stderr, "rk_objective ac_3: end\n");
#endif
  return k;
}

/**
 * Function to get the objective function of 4 steps Runge-Kutta methods.
 *
 * \return objective function value.
 */
static long double
rk_objective_ac_4 (RK * rk)     ///< RK struct.
{
  long double *tb, *ac;
  long double k;
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_4: start\n");
#endif
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  k = fminl (0.L, a20 (ac));
  if (a21 (ac) < 0.L)
    k += a21 (ac);
  if (a30 (ac) < 0.L)
    k += a30 (ac);
  if (a31 (ac) < 0.L)
    k += a31 (ac);
  if (a32 (ac) < 0.L)
    k += a32 (ac);
  if (a40 (ac) < 0.L)
    k += a40 (ac);
  if (a41 (ac) < 0.L)
    k += a41 (ac);
  if (a42 (ac) < 0.L)
    k += a42 (ac);
  if (a43 (ac) < 0.L)
    k += a43 (ac);
  if (k < 0.L)
    {
      k = 20.L - k;
      goto end;
    }
  k = fminl (0.L, c20 (ac));
  if (c21 (ac) < 0.L)
    k += c21 (ac);
  if (c30 (ac) < 0.L)
    k += c30 (ac);
  if (c31 (ac) < 0.L)
    k += c31 (ac);
  if (c32 (ac) < 0.L)
    k += c32 (ac);
  if (c40 (ac) < 0.L)
    k += c40 (ac);
  if (c41 (ac) < 0.L)
    k += c41 (ac);
  if (c42 (ac) < 0.L)
    k += c42 (ac);
  if (c43 (ac) < 0.L)
    k += c43 (ac);
  if (k < 0.L)
    {
      k = 10.L - k;
      goto end;
    }
  k = 1.L / rk_cfl_4 (tb, ac);
end:
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_4: objective=%Lg\n", k);
  fprintf (stderr, "rk_objective ac_4: end\n");
#endif
  return k;
}

/**
 * Function to get the objective function of 5 steps Runge-Kutta methods.
 *
 * \return objective function value.
 */
static long double
rk_objective_ac_5 (RK * rk)     ///< RK struct.
{
  long double *tb, *ac;
  long double k;
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_5: start\n");
#endif
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  k = fminl (0.L, a20 (ac));
  if (a21 (ac) < 0.L)
    k += a21 (ac);
  if (a30 (ac) < 0.L)
    k += a30 (ac);
  if (a31 (ac) < 0.L)
    k += a31 (ac);
  if (a32 (ac) < 0.L)
    k += a32 (ac);
  if (a40 (ac) < 0.L)
    k += a40 (ac);
  if (a41 (ac) < 0.L)
    k += a41 (ac);
  if (a42 (ac) < 0.L)
    k += a42 (ac);
  if (a43 (ac) < 0.L)
    k += a43 (ac);
  if (a50 (ac) < 0.L)
    k += a50 (ac);
  if (a51 (ac) < 0.L)
    k += a51 (ac);
  if (a52 (ac) < 0.L)
    k += a52 (ac);
  if (a53 (ac) < 0.L)
    k += a53 (ac);
  if (a54 (ac) < 0.L)
    k += a54 (ac);
  if (k < 0.L)
    {
      k = 20.L - k;
      goto end;
    }
  k = fminl (0.L, c20 (ac));
  if (c21 (ac) < 0.L)
    k += c21 (ac);
  if (c30 (ac) < 0.L)
    k += c30 (ac);
  if (c31 (ac) < 0.L)
    k += c31 (ac);
  if (c32 (ac) < 0.L)
    k += c32 (ac);
  if (c40 (ac) < 0.L)
    k += c40 (ac);
  if (c41 (ac) < 0.L)
    k += c41 (ac);
  if (c42 (ac) < 0.L)
    k += c42 (ac);
  if (c43 (ac) < 0.L)
    k += c43 (ac);
  if (c50 (ac) < 0.L)
    k += c50 (ac);
  if (c51 (ac) < 0.L)
    k += c51 (ac);
  if (c52 (ac) < 0.L)
    k += c52 (ac);
  if (c53 (ac) < 0.L)
    k += c53 (ac);
  if (c54 (ac) < 0.L)
    k += c54 (ac);
  if (k < 0.L)
    {
      k = 10.L - k;
      goto end;
    }
  k = 1.L / rk_cfl_5 (tb, ac);
end:
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_5: objective=%Lg\n", k);
  fprintf (stderr, "rk_objective ac_5: end\n");
#endif
  return k;
}

/**
 * Function to get the objective function of 6 steps Runge-Kutta methods.
 *
 * \return objective function value.
 */
static long double
rk_objective_ac_6 (RK * rk)     ///< RK struct.
{
  long double *tb, *ac;
  long double k;
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_6: start\n");
#endif
  tb = rk->tb->coefficient;
  ac = rk->ac->coefficient;
  k = fminl (0.L, a20 (ac));
  if (a21 (ac) < 0.L)
    k += a21 (ac);
  if (a30 (ac) < 0.L)
    k += a30 (ac);
  if (a31 (ac) < 0.L)
    k += a31 (ac);
  if (a32 (ac) < 0.L)
    k += a32 (ac);
  if (a40 (ac) < 0.L)
    k += a40 (ac);
  if (a41 (ac) < 0.L)
    k += a41 (ac);
  if (a42 (ac) < 0.L)
    k += a42 (ac);
  if (a43 (ac) < 0.L)
    k += a43 (ac);
  if (a50 (ac) < 0.L)
    k += a50 (ac);
  if (a51 (ac) < 0.L)
    k += a51 (ac);
  if (a52 (ac) < 0.L)
    k += a52 (ac);
  if (a53 (ac) < 0.L)
    k += a53 (ac);
  if (a54 (ac) < 0.L)
    k += a54 (ac);
  if (a60 (ac) < 0.L)
    k += a60 (ac);
  if (a61 (ac) < 0.L)
    k += a61 (ac);
  if (a62 (ac) < 0.L)
    k += a62 (ac);
  if (a63 (ac) < 0.L)
    k += a63 (ac);
  if (a64 (ac) < 0.L)
    k += a64 (ac);
  if (a65 (ac) < 0.L)
    k += a65 (ac);
  if (k < 0.L)
    {
      k = 20.L - k;
      goto end;
    }
  k = fminl (0.L, c20 (ac));
  if (c21 (ac) < 0.L)
    k += c21 (ac);
  if (c30 (ac) < 0.L)
    k += c30 (ac);
  if (c31 (ac) < 0.L)
    k += c31 (ac);
  if (c32 (ac) < 0.L)
    k += c32 (ac);
  if (c40 (ac) < 0.L)
    k += c40 (ac);
  if (c41 (ac) < 0.L)
    k += c41 (ac);
  if (c42 (ac) < 0.L)
    k += c42 (ac);
  if (c43 (ac) < 0.L)
    k += c43 (ac);
  if (c50 (ac) < 0.L)
    k += c50 (ac);
  if (c51 (ac) < 0.L)
    k += c51 (ac);
  if (c52 (ac) < 0.L)
    k += c52 (ac);
  if (c53 (ac) < 0.L)
    k += c53 (ac);
  if (c54 (ac) < 0.L)
    k += c54 (ac);
  if (c60 (ac) < 0.L)
    k += c60 (ac);
  if (c61 (ac) < 0.L)
    k += c61 (ac);
  if (c62 (ac) < 0.L)
    k += c62 (ac);
  if (c63 (ac) < 0.L)
    k += c63 (ac);
  if (c64 (ac) < 0.L)
    k += c64 (ac);
  if (c65 (ac) < 0.L)
    k += c65 (ac);
  if (k < 0.L)
    {
      k = 10.L - k;
      goto end;
    }
  k = 1.L / rk_cfl_6 (tb, ac);
end:
#if DEBUG_RK
  fprintf (stderr, "rk_objective ac_6: objective=%Lg\n", k);
  fprintf (stderr, "rk_objective ac_6: end\n");
#endif
  return k;
}

/**
 * Function to init required variables on a RK struct data.
 */
static inline void
rk_init (RK * rk,               ///< RK struct.
         gsl_rng * rng,         ///< GSL pseudo-random number generator struct.
         unsigned int thread)   ///< thread number.
{
  optimize_init (rk->tb, rng, thread);
  if (rk->strong)
    optimize_init (rk->ac0, rng, 0);
}

/**
 * Function to free the memory allocated by a RK struct.
 */
static inline void
rk_delete (RK * rk)             ///< RK struct.
{
  if (rk->strong)
    optimize_delete (rk->ac0);
  optimize_delete (rk->tb);
}

/**
 * Function to perform every optimization step for the a-c Runge-Kutta 
 * coefficients.
 */
static inline void
rk_step_ac (RK * rk)            ///< RK struct.
{
  Optimize *tb, *ac;
  long double *is, *vo, *vo2;
  long double o, o2, v, f;
  unsigned long long int ii, nsimulations;
  unsigned int i, j, k, n, nfree;

#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: start\n");
#endif

  // save optimal values
#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: save optimal values\n");
#endif
  tb = rk->tb;
  ac = rk->ac;
  nfree = ac->nfree;
  o2 = INFINITY;
  vo = (long double *) alloca (nfree * sizeof (long double));
  vo2 = (long double *) alloca (nfree * sizeof (long double));
  memcpy (vo, ac->value_optimal, nfree * sizeof (long double));

  // optimzation algorithm sampling
#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: optimization algorithm sampling\n");
  fprintf (stderr, "rk_step_ac: nsimulations=%Lu nclimbings=%u nfree=%u\n",
           ac->nsimulations, ac->nclimbings, ac->nfree);
#endif
  nsimulations = ac->nsimulations;
  for (ii = 0L; ii < nsimulations; ++ii)
    {

      // random freedom degrees
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: random freedom degrees\n");
#endif
      optimize_generate_freedom (ac, ii);

      // method coefficients
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: method coefficients\n");
#endif
      if (!ac->method ((Optimize *) rk))
        o = INFINITY;
      else
        o = ac->objective ((Optimize *) rk);
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: objective=%Lg o2=%Lg\n", o, o2);
#endif
      if (o < o2)
        {
          o2 = o;
          memcpy (vo, ac->random_data, nfree * sizeof (long double));
        }
      if (file_variables)
        {
          g_mutex_lock (mutex);
          print_variables (tb->random_data, tb->nfree, file_variables);
          print_variables (ac->random_data, nfree, file_variables);
          fprintf (file_variables, "%.19Le\n", o);
          g_mutex_unlock (mutex);
        }
    }

  // array of intervals to climb around the optimal
#if DEBUG_RK
  fprintf (stderr,
           "rk_step_ac: array of intervals to climb around the optimal\n");
  fprintf (stderr, "rk_step_ac: nclimbings=%u climbing_factor=%Lg\n",
           ac->nclimbings, ac->climbing_factor);
#endif
  is = (long double *) alloca (nfree * sizeof (long double));
  for (j = 0; j < nfree; ++j)
    is[j] = ac->interval0[j] * ac->climbing_factor;
#if DEBUG_RK
  for (j = 0; j < nfree; ++j)
    fprintf (stderr, "rk_step_ac: i=%u is=%Lg\n", j, is[j]);
#endif

  // hill climbing algorithm bucle
#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: hill climbing algorithm bucle\n");
#endif
  memcpy (vo2, vo, nfree * sizeof (long double));
  memcpy (ac->random_data, vo, nfree * sizeof (long double));
  n = ac->nclimbings;
  for (i = 0; i < n; ++i)
    {
#if DEBUG_RK
      for (j = 0; j < nfree; ++j)
        fprintf (stderr, "rk_step_ac: j=%u is=%Lg\n", j, is[j]);
#endif
      for (j = k = 0; j < nfree; ++j)
        {
          v = vo[j];
          ac->random_data[j] = v + is[j];
#if DEBUG_RK
          fprintf (stderr, "rk_step_ac: j=%u random=%Lg\n", j,
                   ac->random_data[j]);
#endif
          if (!ac->method ((Optimize *) rk))
            o = INFINITY;
          else
            o = ac->objective ((Optimize *) rk);
#if DEBUG_RK
          fprintf (stderr, "rk_step_ac: k=%u objective=%Lg o2=%Lg\n", k, o, o2);
#endif
          if (o < o2)
            {
              k = 1;
              o2 = o;
              memcpy (vo2, ac->random_data, nfree * sizeof (long double));
            }
          if (file_variables)
            {
              g_mutex_lock (mutex);
              print_variables (tb->random_data, tb->nfree, file_variables);
              print_variables (ac->random_data, nfree, file_variables);
              fprintf (file_variables, "%.19Le\n", o);
              g_mutex_unlock (mutex);
            }
          ac->random_data[j] = fmaxl (0.L, v - is[j]);
          if (!ac->method ((Optimize *) rk))
            o = INFINITY;
          else
            o = ac->objective ((Optimize *) rk);
#if DEBUG_RK
          fprintf (stderr, "rk_step_ac: k=%u objective=%Lg o2=%Lg\n", k, o, o2);
#endif
          if (o < o2)
            {
              k = 1;
              o2 = o;
              memcpy (vo2, ac->random_data, nfree * sizeof (long double));
            }
          if (file_variables)
            {
              g_mutex_lock (mutex);
              print_variables (tb->random_data, tb->nfree, file_variables);
              print_variables (ac->random_data, nfree, file_variables);
              fprintf (file_variables, "%.19Le\n", o);
              g_mutex_unlock (mutex);
            }
          ac->random_data[j] = v;
        }


      // update optimal values and increase or reduce intervals if converging or
      // not
      if (!k)
        f = 0.5L;
      else
        {
          f = 1.2L;
          memcpy (vo, vo2, nfree * sizeof (long double));
        }
      for (j = 0; j < nfree; ++j)
        is[j] *= f;
    }

  // update optimal values
#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: update optimal values\n");
  fprintf (stderr, "rk_step_ac: optimal=%Lg o2=%Lg\n", *ac->optimal, o2);
#endif
  if (o2 < *ac->optimal)
    {
      *ac->optimal = o2;
      memcpy (ac->value_optimal, vo2, nfree * sizeof (long double));
    }

#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: end\n");
#endif
}

/**
 * Function to do the optimization bucle for the a-c Runge-Kutta coefficients.
 */
void
rk_bucle_ac (RK * rk)           ///< RK struct.
{
  Optimize *tb, *ac, *ac0;
  long double *vo;
  long double optimal;
  unsigned int i, nfree;

#if DEBUG_RK
  fprintf (stderr, "rk_bucle_ac: start\n");
#endif

  tb = rk->tb;
  ac = rk->ac;
  ac0 = rk->ac0;
  nfree = ac0->nfree;
  vo = (long double *) alloca (nfree * sizeof (long double));

  // Init some parameters
#if DEBUG_RK
  fprintf (stderr, "rk_bucle_ac: nfree=%u optimal=%Lg\n", nfree, *tb->optimal);
#endif
  *ac0->optimal = optimal = *tb->optimal;
  for (i = 0; i < nfree; ++i)
    vo[i] = ac0->minimum[i] + 0.5L * ac0->interval[i];
  memcpy (ac, ac0, sizeof (Optimize));
  ac->optimal = &optimal;
  ac->value_optimal = vo;
  optimize_init (ac, ac0->rng, 0);
#if DEBUG_RK
  for (i = 0; i < tb->nfree; ++i)
    fprintf (stderr, "rk_bucle_ac: i=%u random=%Lg\n", i, tb->random_data[i]);
  for (i = 0; i < nfree; ++i)
    fprintf (stderr, "rk_bucle_ac: i=%u minimum=%Lg interval=%Lg type=%u\n",
             i, ac->minimum[i], ac->interval[i], ac->random_type[i]);
#endif

  // Iterate
#if DEBUG_RK
  fprintf (stderr, "rk_bucle_ac: iterate\n");
#endif
  for (i = 0; i < ac->niterations; ++i)
    {

      // Optimization step
      rk_step_ac (rk);

      // Updating coefficient intervals to converge
      optimize_converge (ac);

      // Iterate
#if DEBUG_RK
      fprintf (stderr, "Iteration ac %u\n", i);
#endif
    }

  // Check and save optimal
  if (optimal < *ac0->optimal)
    {
#if DEBUG_RK
      fprintf (stderr, "rk_bucle_ac: optimal=%Lg\n", *ac0->optimal);
#endif
      g_mutex_lock (mutex);
      *ac0->optimal = optimal;
      memcpy (ac0->value_optimal, vo, nfree * sizeof (long double));
      g_mutex_unlock (mutex);
#if DEBUG_RK
      fprintf (stderr, "rk_bucle_ac: optimal=%Lg\n", *ac0->optimal);
      for (i = 0; i < ac0->nfree; ++i)
        fprintf (stderr, "rk_bucle_ac: vo%u=%Lg\n", i, ac0->value_optimal[i]);
#endif
    }

  // Free memory
  optimize_delete (ac);

#if DEBUG_RK
  fprintf (stderr, "rk_bucle_ac: end\n");
#endif
}

/**
 * Function to perform every optimization step for the t-b Runge-Kutta 
 * coefficients.
 */
static inline void
rk_step_tb (RK * rk)            ///< RK struct.
{
  Optimize *tb;
  long double *is, *vo;
  long double o, v, f;
  unsigned long long int ii, nrandom;
  unsigned int b, i, j, k, n, nfree;

#if DEBUG_RK
  fprintf (stderr, "rk_step_tb: start\n");
#endif

  // save optimal values
#if DEBUG_RK
  fprintf (stderr, "rk_step_tb: save optimal values\n");
#endif
  tb = rk->tb;
  nfree = tb->nfree;
  vo = (long double *) alloca (nfree * sizeof (long double));
  b = (file_variables && !rk->strong) ? 1 : 0;

  // optimization algorithm sampling
#if DEBUG_RK
  fprintf (stderr, "rk_step_tb: optimization algorithm sampling\n");
  fprintf (stderr, "rk_step_tb: nsimulations=%Lu nclimbings=%u\n",
           tb->nsimulations, tb->nclimbings);
#endif
  ii = tb->nsimulations * (rank * nthreads + tb->thread) / (nnodes * nthreads);
  nrandom = tb->nsimulations * (rank * nthreads + tb->thread + 1)
    / (nnodes * nthreads);
  for (; ii < nrandom; ++ii)
    {

      // random freedom degrees
#if DEBUG_RK
      fprintf (stderr, "rk_step_tb: random freedom degrees\n");
#endif
      optimize_generate_freedom (tb, ii);

      // method coefficients
#if DEBUG_RK
      fprintf (stderr, "rk_step_tb: method coefficients\n");
#endif
      if (!tb->method (tb))
        o = INFINITY;
      else
        o = tb->objective (tb);
      if (o < *tb->optimal)
        {
          g_mutex_lock (mutex);
          *tb->optimal = o;
          memcpy (tb->value_optimal, tb->random_data,
                  nfree * sizeof (long double));
          g_mutex_unlock (mutex);
        }
      if (b)
        {
          g_mutex_lock (mutex);
          print_variables (tb->random_data, nfree, file_variables);
          fprintf (file_variables, "%.19Le\n", o);
          g_mutex_unlock (mutex);
        }
    }

  // array of intervals to climb around the optimal
#if DEBUG_RK
  fprintf (stderr,
           "rk_step_tb: array of intervals to climb around the optimal\n");
#endif
  is = (long double *) alloca (nfree * sizeof (long double));
  for (j = 0; j < nfree; ++j)
    is[j] = tb->interval0[j] * tb->climbing_factor;

  // hill climbing algorithm bucle
#if DEBUG_RK
  fprintf (stderr, "rk_step_tb: hill climbing algorithm bucle\n");
#endif
  memcpy (tb->random_data, tb->value_optimal, nfree * sizeof (long double));
  memcpy (vo, tb->value_optimal, nfree * sizeof (long double));
  n = tb->nclimbings;
  for (i = 0; i < n; ++i)
    {
      for (j = k = 0; j < nfree; ++j)
        {
          v = vo[j];
          tb->random_data[j] = v + is[j];
          if (!tb->method (tb))
            o = INFINITY;
          else
            o = tb->objective (tb);
          if (o < *tb->optimal)
            {
              k = 1;
              g_mutex_lock (mutex);
              *tb->optimal = o;
              memcpy (tb->value_optimal, tb->random_data,
                      nfree * sizeof (long double));
              g_mutex_unlock (mutex);
            }
          if (b)
            {
              g_mutex_lock (mutex);
              print_variables (tb->random_data, nfree, file_variables);
              fprintf (file_variables, "%.19Le\n", o);
              g_mutex_unlock (mutex);
            }
          tb->random_data[j] = fmaxl (0.L, v - is[j]);
          if (!tb->method (tb))
            o = INFINITY;
          else
            o = tb->objective (tb);
          if (o < *tb->optimal)
            {
              k = 1;
              g_mutex_lock (mutex);
              *tb->optimal = o;
              memcpy (tb->value_optimal, tb->random_data,
                      nfree * sizeof (long double));
              g_mutex_unlock (mutex);
            }
          if (b)
            {
              g_mutex_lock (mutex);
              print_variables (tb->random_data, nfree, file_variables);
              fprintf (file_variables, "%.19Le\n", o);
              g_mutex_unlock (mutex);
            }
          tb->random_data[j] = v;
        }

      // increase or reduce intervals if converging or not
      if (!k)
        f = 0.5L;
      else
        {
          f = 1.2L;
          memcpy (vo, tb->value_optimal, nfree * sizeof (long double));
        }
      for (j = 0; j < nfree; ++j)
        is[j] *= f;
    }
#if DEBUG_RK
  fprintf (stderr, "rk_step_tb: end\n");
#endif
}

/**
 * Function to do the optimization bucle.
 */
static inline void
rk_bucle_tb (RK * rk)           ///< RK struct.
{
  GThread *thread[nthreads];
  Optimize *tb, *ac;
#if HAVE_MPI
  long double *vo;
  MPI_Status status;
#endif
  unsigned int i, j, nfree, nfree2, strong;

#if DEBUG_RK
  fprintf (stderr, "rk_bucle_tb: start\n");
#endif

  // Allocate local array of optimal values
#if DEBUG_RK
  fprintf (stderr, "rk_bucle_tb: allocate local array of optimal values\n");
#endif
  tb = rk->tb;
  nfree = tb->nfree;
  strong = rk->strong;
  if (strong)
    {
      ac = rk->ac0;
      nfree2 = ac->nfree;
    }
  else
    nfree2 = 0;
#if HAVE_MPI
  vo = (long double *) alloca ((1 + nfree + nfree2) * sizeof (long double));
#endif

  // Init some parameters
#if DEBUG_RK
  fprintf (stderr, "rk_bucle_tb: init some parameters\n");
  fprintf (stderr, "rk_bucle_tb: nfree=%u\n", nfree);
#endif
  *tb->optimal = INFINITY;
  for (i = 0; i < nfree; ++i)
    tb->value_optimal[i] = tb->minimum[i] + 0.5L * tb->interval[i];
  if (strong)
    for (i = 0; i < nfree2; ++i)
      ac->value_optimal[i] = ac->minimum[i] + 0.5L * ac->interval[i];

  // Iterate
#if DEBUG_RK
  fprintf (stderr, "rk_bucle_tb: iterate\n");
#endif
  for (i = 0; i < tb->niterations; ++i)
    {

      // Optimization step parallelized for every node by GThreads
      if (nthreads > 1)
        {
          for (j = 0; j < nthreads; ++j)
            thread[j]
              = g_thread_new (NULL,
                              (GThreadFunc) (void (*)(void)) rk_step_tb,
                              (void *) (rk + j));
          for (j = 0; j < nthreads; ++j)
            g_thread_join (thread[j]);
        }
      else
        rk_step_tb (rk);

#if HAVE_MPI
      if (rank > 0)
        {

          // Secondary nodes send the optimal coefficients to the master node
          vo[0] = *tb->optimal;
          memcpy (vo + 1, tb->value_optimal, nfree * sizeof (long double));
          if (strong)
            memcpy (vo + 1 + nfree, ac->value_optimal,
                    nfree2 * sizeof (long double));
          MPI_Send (vo, 1 + nfree + nfree2, MPI_LONG_DOUBLE, 0, 1,
                    MPI_COMM_WORLD);

          // Secondary nodes receive the optimal coefficients
          MPI_Recv (vo, 1 + nfree + nfree2, MPI_LONG_DOUBLE, 0, 1,
                    MPI_COMM_WORLD, &status);
          *tb->optimal = *ac->optimal = vo[0];
          memcpy (tb->value_optimal, vo + 1, nfree * sizeof (long double));
          if (strong)
            memcpy (ac->value_optimal, vo + 1 + nfree,
                    nfree2 * sizeof (long double));
        }
      else
        {
          printf ("rank=%d optimal=%.19Le\n", rank, *tb->optimal);

          for (j = 1; j < nnodes; ++j)
            {

              // Master node receives the optimal coefficients obtained by
              // secondary nodes
              MPI_Recv (vo, 1 + nfree + nfree2, MPI_LONG_DOUBLE, j, 1,
                        MPI_COMM_WORLD, &status);

              // Master node selects the optimal coefficients
              if (vo[0] < *tb->optimal)
                {
                  *tb->optimal = *ac->optimal = vo[0];
                  memcpy (tb->value_optimal, vo + 1,
                          nfree * sizeof (long double));
                  if (strong)
                    memcpy (ac->value_optimal, vo + 1 + nfree,
                            nfree2 * sizeof (long double));
                }
            }

          // Master node sends the optimal coefficients to secondary nodes
          vo[0] = *tb->optimal;
          memcpy (vo + 1, tb->value_optimal, nfree * sizeof (long double));
          if (strong)
            memcpy (vo + 1 + nfree, ac->value_optimal,
                    nfree2 * sizeof (long double));
          for (j = 1; j < nnodes; ++i)
            MPI_Send (vo, 1 + nfree + nfree2, MPI_LONG_DOUBLE, j, 1,
                      MPI_COMM_WORLD);
        }

#endif

      // Print the optimal coefficients
#if DEBUG_RK
      optimize_print_random (tb, stderr);
      if (strong)
        optimize_print_random (ac, stderr);
      fprintf (stderr, "optimal=%.19Le\n", *tb->optimal);
#endif

      // Updating coefficient intervals to converge
      optimize_converge (tb);

      // Iterate
      printf ("Iteration %u Optimal %.19Le\n", i, *tb->optimal);
    }

#if DEBUG_RK
  fprintf (stderr, "rk_bucle_tb: end\n");
#endif
}

/**
 * Function to select the Runge-Kutta method.
 *
 * \return 1 on success, 0 on error.
 */
static inline int
rk_select (RK * rk,             ///< RK struct.
           unsigned int nsteps, ///< steps number.
           unsigned int order)  ///< accuracy order.
{
  static int (*tb_method[7][6]) (Optimize *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_2_2, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_3_2, &rk_tb_3_3, NULL, NULL},
    {
    NULL, NULL, &rk_tb_4_2, &rk_tb_4_3, &rk_tb_4_4, NULL},
    {
    NULL, NULL, &rk_tb_5_2, &rk_tb_5_3, &rk_tb_5_4, NULL},
    {
    NULL, NULL, &rk_tb_6_2, &rk_tb_6_3, &rk_tb_6_4, NULL}
  };
  static int (*tb_method_t[7][6]) (Optimize *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_2_2t, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_3_2t, &rk_tb_3_3t, NULL, NULL},
    {
    NULL, NULL, &rk_tb_4_2t, &rk_tb_4_3t, &rk_tb_4_4t, NULL},
    {
    NULL, NULL, &rk_tb_5_2t, &rk_tb_5_3t, &rk_tb_5_4t, NULL},
    {
    NULL, NULL, &rk_tb_6_2t, &rk_tb_6_3t, &rk_tb_6_4t, NULL}
  };
  static int (*tb_method_p[7][6]) (Optimize *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_2_2p, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_3_2p, &rk_tb_3_3p, NULL, NULL},
    {
    NULL, NULL, &rk_tb_4_2p, &rk_tb_4_3p, NULL, NULL},
    {
    NULL, NULL, &rk_tb_5_2p, &rk_tb_5_3p, &rk_tb_5_4p, NULL},
    {
    NULL, NULL, &rk_tb_6_2p, &rk_tb_6_3p, &rk_tb_6_4p, NULL}
  };
  static int (*tb_method_tp[7][6]) (Optimize *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_2_2tp, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_tb_3_2tp, &rk_tb_3_3tp, NULL, NULL},
    {
    NULL, NULL, &rk_tb_4_2tp, &rk_tb_4_3tp, NULL, NULL},
    {
    NULL, NULL, &rk_tb_5_2tp, &rk_tb_5_3tp, &rk_tb_5_4tp, NULL},
    {
    NULL, NULL, &rk_tb_6_2tp, &rk_tb_6_3tp, &rk_tb_6_4tp, NULL}
  };
  static long double (*tb_objective[7][6]) (RK *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_2_2, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_2_2, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_3_2, &rk_objective_tb_3_3, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_4_2, &rk_objective_tb_4_3,
        &rk_objective_tb_4_4, NULL},
    {
    NULL, NULL, &rk_objective_tb_5_2, &rk_objective_tb_5_3,
        &rk_objective_tb_5_4, NULL},
    {
    NULL, NULL, &rk_objective_tb_6_2, &rk_objective_tb_6_3,
        &rk_objective_tb_6_4, NULL}
  };
  static long double (*tb_objective_t[7][6]) (RK *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_2_2t, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_2_2t, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_3_2t, &rk_objective_tb_3_3t, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_4_2t, &rk_objective_tb_4_3t,
        &rk_objective_tb_4_4t, NULL},
    {
    NULL, NULL, &rk_objective_tb_5_2t, &rk_objective_tb_5_3t,
        &rk_objective_tb_5_4t, NULL},
    {
    NULL, NULL, &rk_objective_tb_6_2t, &rk_objective_tb_6_3t,
        &rk_objective_tb_6_4t, NULL}
  };
  static long double (*tb_objective_p[7][6]) (RK *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_2_2, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_3_2, &rk_objective_tb_3_3p, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_4_2, &rk_objective_tb_4_3p, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_5_2, &rk_objective_tb_5_3p,
        &rk_objective_tb_5_4p, NULL},
    {
    NULL, NULL, &rk_objective_tb_6_2, &rk_objective_tb_6_3p,
        &rk_objective_tb_6_4p, NULL}
  };
  static long double (*tb_objective_tp[7][6]) (RK *) =
  {
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, NULL, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_2_2t, NULL, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_3_2t, &rk_objective_tb_3_3tp, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_4_2t, &rk_objective_tb_4_3tp, NULL, NULL},
    {
    NULL, NULL, &rk_objective_tb_5_2t, &rk_objective_tb_5_3tp,
        &rk_objective_tb_5_4tp, NULL},
    {
    NULL, NULL, &rk_objective_tb_6_2t, &rk_objective_tb_6_3tp,
        &rk_objective_tb_6_4tp, NULL}
  };
  static int (*ac_method[7]) (RK *) =
  {
  NULL, NULL, &rk_ac_2, &rk_ac_3, &rk_ac_4, &rk_ac_5, &rk_ac_6};
  static long double (*ac_objective[7]) (RK * rk) =
  {
  NULL, NULL, &rk_objective_ac_2, &rk_objective_ac_3, &rk_objective_ac_4,
      &rk_objective_ac_5, &rk_objective_ac_6};
  const unsigned int nequations[6] = { 0, 1, 2, 4, 8, 16 };
  Optimize *tb, *ac;
#if DEBUG_RK
  fprintf (stderr, "rk_select: start\n");
#endif
  tb = rk->tb;
  ac = rk->ac;
  tb->nsteps = nsteps;
  tb->order = order;
  tb->size = nsteps * (nsteps + 3) / 2 - 1;
  tb->nfree = tb->size - nsteps + 1 - nequations[order];
  if (rk->pair)
    {
      tb->size += nsteps - 1;
      if (rk->time_accuracy)
        {
          tb->method = tb_method_tp[nsteps][order];
          tb->objective = (OptimizeObjective) tb_objective_tp[nsteps][order];
        }
      else
        {
          tb->method = tb_method_p[nsteps][order];
          tb->objective = (OptimizeObjective) tb_objective_p[nsteps][order];
        }
    }
  else
    {
      if (rk->time_accuracy)
        {
          tb->method = tb_method_t[nsteps][order];
          tb->objective = (OptimizeObjective) tb_objective_t[nsteps][order];
        }
      else
        {
          tb->method = tb_method[nsteps][order];
          tb->objective = (OptimizeObjective) tb_objective[nsteps][order];
        }
    }
  if (!tb->method)
    goto exit_on_error;
  switch (nsteps)
    {
    case 5:
      switch (order)
        {
        case 4:
          if (rk->time_accuracy && rk->pair)
            --tb->nfree;
        }
      break;
    }
  if (rk->time_accuracy)
    --tb->nfree;
  tb->minimum0
    = (long double *) g_slice_alloc (tb->nfree * sizeof (long double));
  tb->interval0
    = (long double *) g_slice_alloc (tb->nfree * sizeof (long double));
  tb->random_type
    = (unsigned int *) g_slice_alloc (tb->nfree * sizeof (unsigned int));
  if (rk->strong)
    {
      ac = rk->ac0;
      ac->size = nsteps * (nsteps + 1) - 2;
      ac->nfree = nsteps * (nsteps - 1) / 2;
      ac->minimum0
        = (long double *) g_slice_alloc (ac->nfree * sizeof (long double));
      ac->interval0
        = (long double *) g_slice_alloc (ac->nfree * sizeof (long double));
      ac->random_type
        = (unsigned int *) g_slice_alloc (ac->nfree * sizeof (unsigned int));
      ac->method = (OptimizeMethod) ac_method[nsteps];
      if (!ac->method)
        goto exit_on_error;
      ac->objective = (OptimizeObjective) ac_objective[nsteps];
    }

#if DEBUG_RK
  fprintf (stderr, "rk_select: end\n");
#endif
  return 1;

exit_on_error:
  error_message = g_strdup (_("Unknown method"));
#if DEBUG_RK
  fprintf (stderr, "rk_select: end\n");
#endif
  return 0;
}

/**
 * Function to read the Runge-Kutta method data on a XML node.
 *
 * \return 1 on success, 0 on error.
 */
int
rk_run (xmlNode * node,         ///< XML node.
        gsl_rng ** rng)         ///< array of gsl_rng structs.
{
  RK rk[nthreads];
  char filename[64];
  Optimize *tb, *ac;
  gchar *buffer;
  xmlChar *prop;
  FILE *file;
  long double *value_optimal, *value_optimal2;
  long double optimal, optimal2;
  int code;
  unsigned int i, j, nsteps, order, nfree, nfree2;

#if DEBUG_RK
  fprintf (stderr, "rk_run: start\n");
#endif

  tb = rk->tb;
  nsteps = xml_node_get_uint (node, XML_STEPS, &code);
  if (code)
    {
      error_message = g_strdup (_("Bad steps number"));
      goto exit_on_error;
    }
  order = xml_node_get_uint (node, XML_ORDER, &code);
  if (code)
    {
      error_message = g_strdup (_("Bad order"));
      goto exit_on_error;
    }
  prop = xmlGetProp (node, XML_STRONG);
  if (!prop || !xmlStrcmp (prop, XML_NO))
    rk->strong = 0;
  else if (!xmlStrcmp (prop, XML_YES))
    rk->strong = 1;
  else
    {
      error_message = g_strdup (_("Bad strong stability"));
      goto exit_on_error;
    }
  xmlFree (prop);
  prop = xmlGetProp (node, XML_PAIR);
  if (!prop || !xmlStrcmp (prop, XML_NO))
    rk->pair = 0;
  else if (!xmlStrcmp (prop, XML_YES))
    rk->pair = 1;
  else
    {
      error_message = g_strdup (_("Bad pair"));
      goto exit_on_error;
    }
  xmlFree (prop);
  prop = xmlGetProp (node, XML_TIME_ACCURACY);
  if (!prop || !xmlStrcmp (prop, XML_NO))
    rk->time_accuracy = 0;
  else if (!xmlStrcmp (prop, XML_YES))
    rk->time_accuracy = 1;
  else
    {
      error_message = g_strdup (_("Bad time accuracy"));
      goto exit_on_error;
    }
  xmlFree (prop);
  if (!rk_select (rk, nsteps, order))
    goto exit_on_error;
  if (!optimize_read (tb, node))
    goto exit_on_error;
  nfree = tb->nfree;
  value_optimal = (long double *) g_slice_alloc (nfree * sizeof (long double));
  optimize_create (tb, &optimal, value_optimal);
  node = node->children;
  for (i = 0; i < nfree; ++i, node = node->next)
    if (!read_variable (node, tb->minimum0, tb->interval0, tb->random_type, i))
      goto exit_on_error;
  if (rk->strong)
    {
      ac = rk->ac0;
      if (!node)
        {
          error_message = g_strdup (_("No a-c coefficients data"));
          goto exit_on_error;
        }
      if (xmlStrcmp (node->name, XML_AC))
        {
          error_message = g_strdup (_("Bad a-c coefficients XML node"));
          goto exit_on_error;
        }
      if (!optimize_read (ac, node))
        {
          buffer = error_message;
          error_message
            = g_strconcat (_("a-c coefficients"), ":\n", error_message, NULL);
          g_free (buffer);
          goto exit_on_error;
        }
      nfree2 = ac->nfree;
      value_optimal2
        = (long double *) g_slice_alloc (nfree2 * sizeof (long double));
      optimize_create (ac, &optimal2, value_optimal2);
      for (i = 0; i < nfree2; ++i)
        {
          node = node->next;
          if (!read_variable (node, ac->minimum0, ac->interval0,
                              ac->random_type, i))
            goto exit_on_error;
        }
    }
  for (i = 1; i < nthreads; ++i)
    memcpy (rk + i, rk, sizeof (RK));
  j = rank * nthreads;
  for (i = 0; i < nthreads; ++i)
    rk_init (rk + i, rng[j + i], i);

  // Method bucle
  printf ("Optimize bucle\n");
  rk_bucle_tb (rk);

  // Print the optimal coefficients
  printf ("Print the optimal coefficients\n");
  memcpy (tb->random_data, tb->value_optimal, nfree * sizeof (long double));
  code = tb->method (tb);
  if (rk->strong)
    {
      memcpy (ac->random_data, ac->value_optimal,
              nfree2 * sizeof (long double));
      memcpy (rk->ac, ac, sizeof (Optimize));
      code = ac->method ((Optimize *) rk);
    }
  snprintf (filename, 64, "rk-%u-%u-%u-%u-%u.mc",
            nsteps, order, rk->time_accuracy, rk->pair, rk->strong);
  file = fopen (filename, "w");
  print_maxima_precision (file);
  rk_print (rk, file);
  rk_print_maxima (file, nsteps, nsteps, order, 'b');
  if (rk->pair)
    rk_print_maxima (file, nsteps, nsteps - 1, order - 1, 'e');
  if (rk->strong)
    ac_print_maxima (file, nsteps);
  fclose (file);
  snprintf (filename, 64, "sed -i 's/e+/b+/g' rk-%u-%u-%u-%u-%u.mc",
            nsteps, order, rk->time_accuracy, rk->pair, rk->strong);
  code = system (filename);
  snprintf (filename, 64, "sed -i 's/e-/b-/g' rk-%u-%u-%u-%u-%u.mc",
            nsteps, order, rk->time_accuracy, rk->pair, rk->strong);
  code = system (filename);

  // Free memory
  if (rk->strong)
    {
      g_slice_free1 (nfree2 * sizeof (unsigned int), ac->random_type);
      g_slice_free1 (nfree2 * sizeof (long double), ac->interval0);
      g_slice_free1 (nfree2 * sizeof (long double), ac->minimum0);
      g_slice_free1 (nfree2 * sizeof (long double), value_optimal2);
    }
  for (i = 0; i < nthreads; ++i)
    rk_delete (rk + i);
  g_slice_free1 (nfree * sizeof (unsigned int), tb->random_type);
  g_slice_free1 (nfree * sizeof (long double), tb->interval0);
  g_slice_free1 (nfree * sizeof (long double), tb->minimum0);
  g_slice_free1 (nfree * sizeof (long double), value_optimal);


#if DEBUG_RK
  fprintf (stderr, "rk_run: end\n");
#endif
  return 1;

exit_on_error:
  buffer = error_message;
  error_message = g_strconcat ("Runge-Kutta:\n", buffer, NULL);
  g_free (buffer);
#if DEBUG_RK
  fprintf (stderr, "rk_run: end\n");
#endif
  return 0;
}
