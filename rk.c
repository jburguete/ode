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
 * \file rk.c
 * \brief Source file with common variables and functions to optimize
 *   Runge-Kutta methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
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

///> array of minimum freedom degree values for the a-c coefficients of the 2
///> steps Runge-Kutta methods.
const long double minimum_ac_rk_2[1] = { 0.L };

///> array of minimum freedom degree intervals for the a-c coefficients of the 2
///> steps Runge-Kutta methods.
const long double interval_ac_rk_2[1] = { 2.L };

///> array of freedom degree random function types for the a-c coefficients of
///> the 2 steps Runge-Kutta methods.
const unsigned int random_ac_rk_2[1] = { 0 };

///> array of minimum freedom degree values for the a-c coefficients of the 3
///> steps Runge-Kutta methods.
const long double minimum_ac_rk_3[3] = { 0.L, 0.L, 0.L };

///> array of minimum freedom degree intervals for the a-c coefficients of the 3
///> steps Runge-Kutta methods.
const long double interval_ac_rk_3[3] = { 2.L, 2.L, 2.L };

///> array of freedom degree random function types for the a-c coefficients of
///> the 3 steps Runge-Kutta methods.
const unsigned int random_ac_rk_3[3] = { 0, 0, 0 };

///> array of minimum freedom degree values for the a-c coefficients of the 4
///> steps Runge-Kutta methods.
const long double minimum_ac_rk_4[6] = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of minimum freedom degree intervals for the a-c coefficients of the 4
///> steps Runge-Kutta methods.
const long double interval_ac_rk_4[6] = { 2.L, 2.L, 2.L, 2.L, 2.L, 2.L };

///> array of freedom degree random function types for the a-c coefficients of
///> the 4 steps Runge-Kutta methods.
const unsigned int random_ac_rk_4[6] = { 0, 0, 0, 0, 0, 0 };

///> array of minimum freedom degree values for the a-c coefficients of the 5
///> steps Runge-Kutta methods.
const long double minimum_ac_rk_5[10]
  = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of minimum freedom degree intervals for the a-c coefficients of the 5
///> steps Runge-Kutta methods.
const long double interval_ac_rk_5[10]
  = { 2.L, 2.L, 2.L, 2.L, 2.L, 2.L, 2.L, 2.L, 2.L, 2.L };

///> array of freedom degree random function types for the a-c coefficients of
///> the 5 steps Runge-Kutta methods.
const unsigned int random_ac_rk_5[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

///> array of minimum freedom degree values for the a-c coefficients of the 6
///> steps Runge-Kutta methods.
const long double minimum_ac_rk_6[15]
  =
  { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of minimum freedom degree intervals for the a-c coefficients of the 6
///> steps Runge-Kutta methods.
const long double interval_ac_rk_6[15]
  =
  { 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L, 1.L };

///> array of freedom degree random function types for the a-c coefficients of
///> the 6 steps Runge-Kutta methods.
const unsigned int random_ac_rk_6[15] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/**
 * Function to print the 2nd step t-b Runge-Kutta coefficients.
 */
void
rk_print_tb_2 (long double *tb, ///< array of t-b Runge-Kutta coefficients.
               char *label,     ///< label.
               FILE * file)     ///< file.
{
  fprintf (file, "%s: t1=%.19Le\n", label, t1 (tb));
  fprintf (file, "%s: t2=%.19Le\n", label, t2 (tb));
  fprintf (file, "%s: b20=%.19Le\n", label, b20 (tb));
  fprintf (file, "%s: b21=%.19Le\n", label, b21 (tb));
}

/**
 * Function to print the 3nd step t-b Runge-Kutta coefficients.
 */
void
rk_print_tb_3 (long double *tb, ///< array of t-b Runge-Kutta coefficients.
               char *label,     ///< label.
               FILE * file)     ///< file.
{
  rk_print_tb_2 (tb, label, file);
  fprintf (file, "%s: t3=%.19Le\n", label, t3 (tb));
  fprintf (file, "%s: b30=%.19Le\n", label, b30 (tb));
  fprintf (file, "%s: b31=%.19Le\n", label, b31 (tb));
  fprintf (file, "%s: b32=%.19Le\n", label, b32 (tb));
}

/**
 * Function to print the 4th step t-b Runge-Kutta coefficients.
 */
void
rk_print_tb_4 (long double *tb, ///< array of t-b Runge-Kutta coefficients.
               char *label,     ///< label.
               FILE * file)     ///< file.
{
  rk_print_tb_3 (tb, label, file);
  fprintf (file, "%s: t4=%.19Le\n", label, t4 (tb));
  fprintf (file, "%s: b40=%.19Le\n", label, b40 (tb));
  fprintf (file, "%s: b41=%.19Le\n", label, b41 (tb));
  fprintf (file, "%s: b42=%.19Le\n", label, b42 (tb));
  fprintf (file, "%s: b43=%.19Le\n", label, b43 (tb));
}

/**
 * Function to print the 5th step t-b Runge-Kutta coefficients.
 */
void
rk_print_tb_5 (long double *tb, ///< array of t-b Runge-Kutta coefficients.
               char *label,     ///< label.
               FILE * file)     ///< file.
{
  rk_print_tb_4 (tb, label, file);
  fprintf (file, "%s: t5=%.19Le\n", label, t5 (tb));
  fprintf (file, "%s: b50=%.19Le\n", label, b50 (tb));
  fprintf (file, "%s: b51=%.19Le\n", label, b51 (tb));
  fprintf (file, "%s: b52=%.19Le\n", label, b52 (tb));
  fprintf (file, "%s: b53=%.19Le\n", label, b53 (tb));
  fprintf (file, "%s: b54=%.19Le\n", label, b54 (tb));
}

/**
 * Function to print the 6th step t-b Runge-Kutta coefficients.
 */
void
rk_print_tb_6 (long double *tb, ///< array of t-b Runge-Kutta coefficients.
               char *label,     ///< label.
               FILE * file)     ///< file.
{
  rk_print_tb_5 (tb, label, file);
  fprintf (file, "%s: t6=%.19Le\n", label, t6 (tb));
  fprintf (file, "%s: b60=%.19Le\n", label, b60 (tb));
  fprintf (file, "%s: b61=%.19Le\n", label, b61 (tb));
  fprintf (file, "%s: b62=%.19Le\n", label, b62 (tb));
  fprintf (file, "%s: b63=%.19Le\n", label, b63 (tb));
  fprintf (file, "%s: b64=%.19Le\n", label, b64 (tb));
  fprintf (file, "%s: b65=%.19Le\n", label, b65 (tb));
}

/**
 * Function to print in a maxima file the 2nd step Runge-Kutta coefficients.
 */
void
rk_print_2 (RK * rk,            ///< RK struct.
            FILE * file)        ///< file.
{
  long double *tb, *ac;
  tb = rk->tb->coefficient;
  fprintf (file, "t1:%.19Le;\n", t1 (tb));
  fprintf (file, "t2:%.19Le;\n", t2 (tb));
  fprintf (file, "b20:%.19Le;\n", b20 (tb));
  fprintf (file, "b21:%.19Le;\n", b21 (tb));
  if (!rk->strong)
    return;
  ac = rk->ac->coefficient;
  fprintf (file, "a20:%.19Le;\n", a20 (ac));
  fprintf (file, "a21:%.19Le;\n", a21 (ac));
  fprintf (file, "c20:%.19Le;\n", c20 (ac));
  fprintf (file, "c21:%.19Le;\n", c21 (ac));
}

/**
 * Function to print in a maxima file the 3rd step Runge-Kutta coefficients.
 */
void
rk_print_3 (RK * rk,            ///< RK struct.
            FILE * file)        ///< file.
{
  long double *tb, *ac;
  rk_print_2 (rk, file);
  tb = rk->tb->coefficient;
  fprintf (file, "t3:%.19Le;\n", t3 (tb));
  fprintf (file, "b30:%.19Le;\n", b30 (tb));
  fprintf (file, "b31:%.19Le;\n", b31 (tb));
  fprintf (file, "b32:%.19Le;\n", b32 (tb));
  if (!rk->strong)
    return;
  ac = rk->ac->coefficient;
  fprintf (file, "a30:%.19Le;\n", a30 (ac));
  fprintf (file, "a31:%.19Le;\n", a31 (ac));
  fprintf (file, "a32:%.19Le;\n", a32 (ac));
  fprintf (file, "c30:%.19Le;\n", c30 (ac));
  fprintf (file, "c31:%.19Le;\n", c31 (ac));
  fprintf (file, "c32:%.19Le;\n", c32 (ac));
}

/**
 * Function to print in a maxima file the 4th step Runge-Kutta coefficients.
 */
void
rk_print_4 (RK * rk,            ///< RK struct.
            FILE * file)        ///< file.
{
  long double *tb, *ac;
  rk_print_3 (rk, file);
  tb = rk->tb->coefficient;
  fprintf (file, "t4:%.19Le;\n", t4 (tb));
  fprintf (file, "b40:%.19Le;\n", b40 (tb));
  fprintf (file, "b41:%.19Le;\n", b41 (tb));
  fprintf (file, "b42:%.19Le;\n", b42 (tb));
  fprintf (file, "b43:%.19Le;\n", b43 (tb));
  if (!rk->strong)
    return;
  ac = rk->ac->coefficient;
  fprintf (file, "a40:%.19Le;\n", a40 (ac));
  fprintf (file, "a41:%.19Le;\n", a41 (ac));
  fprintf (file, "a42:%.19Le;\n", a42 (ac));
  fprintf (file, "a43:%.19Le;\n", a43 (ac));
  fprintf (file, "c40:%.19Le;\n", c40 (ac));
  fprintf (file, "c41:%.19Le;\n", c41 (ac));
  fprintf (file, "c42:%.19Le;\n", c42 (ac));
  fprintf (file, "c43:%.19Le;\n", c43 (ac));
}

/**
 * Function to print in a maxima file the 5th step Runge-Kutta coefficients.
 */
void
rk_print_5 (RK * rk,            ///< RK struct.
            FILE * file)        ///< file.
{
  long double *tb, *ac;
  rk_print_4 (rk, file);
  tb = rk->tb->coefficient;
  fprintf (file, "t5:%.19Le;\n", t5 (tb));
  fprintf (file, "b50:%.19Le;\n", b50 (tb));
  fprintf (file, "b51:%.19Le;\n", b51 (tb));
  fprintf (file, "b52:%.19Le;\n", b52 (tb));
  fprintf (file, "b53:%.19Le;\n", b53 (tb));
  fprintf (file, "b54:%.19Le;\n", b54 (tb));
  if (!rk->strong)
    return;
  ac = rk->ac->coefficient;
  fprintf (file, "a50:%.19Le;\n", a50 (ac));
  fprintf (file, "a51:%.19Le;\n", a51 (ac));
  fprintf (file, "a52:%.19Le;\n", a52 (ac));
  fprintf (file, "a53:%.19Le;\n", a53 (ac));
  fprintf (file, "a54:%.19Le;\n", a54 (ac));
  fprintf (file, "c50:%.19Le;\n", c50 (ac));
  fprintf (file, "c51:%.19Le;\n", c51 (ac));
  fprintf (file, "c52:%.19Le;\n", c52 (ac));
  fprintf (file, "c53:%.19Le;\n", c53 (ac));
  fprintf (file, "c54:%.19Le;\n", c54 (ac));
}

/**
 * Function to print in a maxima file the 6th step Runge-Kutta coefficients.
 */
void
rk_print_6 (RK * rk,            ///< RK struct.
            FILE * file)        ///< file.
{
  long double *tb, *ac;
  rk_print_5 (rk, file);
  tb = rk->tb->coefficient;
  fprintf (file, "t6:%.19Le;\n", t6 (tb));
  fprintf (file, "b60:%.19Le;\n", b60 (tb));
  fprintf (file, "b61:%.19Le;\n", b61 (tb));
  fprintf (file, "b62:%.19Le;\n", b62 (tb));
  fprintf (file, "b63:%.19Le;\n", b63 (tb));
  fprintf (file, "b64:%.19Le;\n", b64 (tb));
  fprintf (file, "b65:%.19Le;\n", b65 (tb));
  if (!rk->strong)
    return;
  ac = rk->ac->coefficient;
  fprintf (file, "a60:%.19Le;\n", a60 (ac));
  fprintf (file, "a61:%.19Le;\n", a61 (ac));
  fprintf (file, "a62:%.19Le;\n", a62 (ac));
  fprintf (file, "a63:%.19Le;\n", a63 (ac));
  fprintf (file, "a64:%.19Le;\n", a64 (ac));
  fprintf (file, "a65:%.19Le;\n", a65 (ac));
  fprintf (file, "c60:%.19Le;\n", c60 (ac));
  fprintf (file, "c61:%.19Le;\n", c61 (ac));
  fprintf (file, "c62:%.19Le;\n", c62 (ac));
  fprintf (file, "c63:%.19Le;\n", c63 (ac));
  fprintf (file, "c64:%.19Le;\n", c64 (ac));
  fprintf (file, "c65:%.19Le;\n", c65 (ac));
}

/**
 * Function to print in a maxima file the 2nd step Runge-Kutta equations.
 */
void
rk_print_maxima_2 (FILE * file) ///< file.
{
  fprintf (file, "a20+a21-1;\n");
  fprintf (file, "a20*c20+a21*t1-b20;\n");
  fprintf (file, "a21*c21-b21;\n");
}

/**
 * Function to print in a maxima file the 3rd step Runge-Kutta equations.
 */
void
rk_print_maxima_3 (FILE * file) ///< file.
{
  rk_print_maxima_2 (file);
  fprintf (file, "a30+a31+a32-1;\n");
  fprintf (file, "a30*c30+a31*t1+a32*b20-b30;\n");
  fprintf (file, "a31*c31+a32*b21-b31;\n");
  fprintf (file, "a32*c32-b32;\n");
}

/**
 * Function to print in a maxima file the 4th step Runge-Kutta equations.
 */
void
rk_print_maxima_4 (FILE * file) ///< file.
{
  rk_print_maxima_3 (file);
  fprintf (file, "a40+a41+a42+a43-1;\n");
  fprintf (file, "a40*c40+a41*t1+a42*b20+a43*b30-b40;\n");
  fprintf (file, "a41*c41+a42*b21+a43*b31-b41;\n");
  fprintf (file, "a42*c42+a43*b32-b42;\n");
  fprintf (file, "a43*c43-b43;\n");
}

/**
 * Function to print in a maxima file the 5th step Runge-Kutta equations.
 */
void
rk_print_maxima_5 (FILE * file) ///< file.
{
  rk_print_maxima_4 (file);
  fprintf (file, "a50+a51+a52+a53+a54-1;\n");
  fprintf (file, "a50*c50+a51*t1+a52*b20+a53*b30+a54*b40-b50;\n");
  fprintf (file, "a51*c51+a52*b21+a53*b31+a54*b41-b51;\n");
  fprintf (file, "a52*c52+a53*b32+a54*b42-b52;\n");
  fprintf (file, "a53*c53+a54*b43-b53;\n");
  fprintf (file, "a54*c54-b54;\n");
}

/**
 * Function to print in a maxima file the 6th step Runge-Kutta equations.
 */
void
rk_print_maxima_6 (FILE * file) ///< file.
{
  rk_print_maxima_5 (file);
  fprintf (file, "a60+a61+a62+a63+a64+a65-1;\n");
  fprintf (file, "a60*c60+a61*t1+a62*b20+a63*b30+a64*b40+a65*b50-b60;\n");
  fprintf (file, "a61*c61+a62*b21+a63*b31+a64*b41+a65*b51-b61;\n");
  fprintf (file, "a62*c62+a63*b32+a64*b42+a65*b52-b62;\n");
  fprintf (file, "a63*c63+a64*b43+a65*b53-b63;\n");
  fprintf (file, "a64*c64+a65*b54-b64;\n");
  fprintf (file, "a65*c65-b65;\n");
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 2nd step.
 */
void
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
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 3rd step.
 */
void
rk_ac_3 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_3: start\n");
#endif
  rk_ac_2 (rk);
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
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 4th step.
 */
void
rk_ac_4 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_4: start\n");
#endif
  rk_ac_3 (rk);
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
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 5th step.
 */
void
rk_ac_5 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_5: start\n");
#endif
  rk_ac_4 (rk);
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
}

/**
 * Function to get \f$a_{ij}\f$ and \f$c_{ij}\f$ coefficients of the Runge-Kutta
 * 6th step.
 */
void
rk_ac_6 (RK * rk)               ///< RK struct.
{
  long double *tb, *ac, *r;
  register long double ac0;
#if DEBUG_RK
  fprintf (stderr, "rk_ac_6: start\n");
#endif
  rk_ac_5 (rk);
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
}

/**
 * Function to get the objective function of 2 steps Runge-Kutta methods.
 *
 * \return objective function value.
 */
long double
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
long double
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
long double
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
long double
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
long double
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
void
rk_init (RK * rk,               ///< RK struct.
         gsl_rng * rng)         ///< GSL pseudo-random number generator struct.
{
  optimize_init (rk->tb, rng);
  if (rk->strong)
    optimize_init (rk->ac0, rng);
}

/**
 * Function to free the memory allocated by a RK struct.
 */
void
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
void
rk_step_ac (RK * rk)            ///< RK struct.
{
  Optimize *ac;
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
  ac = rk->ac;
  nfree = ac->nfree;
  o2 = *ac->optimal;
  vo = (long double *) alloca (nfree * sizeof (long double));
  vo2 = (long double *) alloca (nfree * sizeof (long double));
  memcpy (vo, ac->value_optimal, nfree * sizeof (long double));

  // Monte-Carlo algorithm
#if DEBUG_RK
  fprintf (stderr, "rk_step_ac: Monte-Carlo algorithm\n");
  fprintf (stderr, "rk_step_ac: nsimulations=%Lu nsearch=%u\n",
           ac->nsimulations, ac->nsearch);
#endif
  nsimulations = ac->nsimulations;
  for (ii = 0L; ii < nsimulations; ++ii)
    {

      // random freedom degrees
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: random freedom degrees\n");
#endif
      optimize_generate_random (ac);
#if PRINT_RANDOM
      print_random (ac->random_data, nfree, file_random2);
#endif

      // method coefficients
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: method coefficients\n");
#endif
      ac->method ((Optimize *) rk);

      // objective function
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: objective function\n");
#endif
      o = ac->objective ((Optimize *) rk);
#if DEBUG_RK
      fprintf (stderr, "rk_step_ac: objective=%Lg o2=%Lg\n", o, o2);
#endif
      if (o < o2)
        {
          o2 = o;
          memcpy (vo, ac->random_data, nfree * sizeof (long double));
        }
    }

  // array of intervals to search around the optimal
#if DEBUG_RK
  fprintf (stderr,
           "rk_step_ac: array of intervals to search around the optimal\n");
  fprintf (stderr, "rk_step_ac: nsearch=%u search_factor=%Lg\n",
           ac->nsearch, ac->search_factor);
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
  n = ac->nclimbings;
  for (i = 0; i < n; ++i)
    {
#if DEBUG_RK
      for (j = 0; j < nfree; ++j)
        fprintf (stderr, "rk_step_ac: j=%u is=%Lg\n", j, is[j]);
#endif
      memcpy (ac->random_data, vo, nfree * sizeof (long double));
      for (j = k = 0; j < nfree; ++j)
        {
          v = vo[j];
          ac->random_data[j] = v + is[j];
#if DEBUG_RK
          fprintf (stderr, "rk_step_ac: j=%u random=%Lg\n", j,
                   ac->random_data[j]);
#endif
          ac->method ((Optimize *) rk);
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
          ac->random_data[j] = fmaxl(0.L, v - is[j]);
          ac->method ((Optimize *) rk);
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
          ac->random_data[j] = v;
        }

      // increase or reduce intervals if converging or not
      if (!k)
				f = 0.5L;
			else
				f = 1.2L;
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
  optimize_init (ac, ac0->rng);
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
  for (i = 0; i < tb->niterations; ++i)
    {

      // Optimization step
      rk_step_ac (rk);

      // Updating coefficient intervals to converge
      optimize_converge (tb);

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
void
rk_step_tb (RK * rk)            ///< RK struct.
{
  Optimize *tb;
  long double *is, *vo;
  long double o, v, f;
  unsigned long long int ii, nrandom;
  unsigned int i, j, k, n, nfree;

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

  // Monte-Carlo algorithm
#if DEBUG_RK
  fprintf (stderr, "rk_step_tb: Monte-Carlo algorithm\n");
  fprintf (stderr, "rk_step_tb: nrandom=%Lu nsearch=%u\n",
           tb->nrandom, tb->nsearch);
#endif
  nrandom = tb->nrandom;
  for (ii = 0L; ii < nrandom; ++ii)
    {

      // random freedom degrees
#if DEBUG_RK
      fprintf (stderr, "rk_step_tb: random freedom degrees\n");
#endif
      optimize_generate_random (tb);
#if PRINT_RANDOM
      print_random (tb->random_data, nfree, file_random);
#endif
#if DEBUG_RK
      print_random (tb->random_data, nfree, stderr);
#endif

      // method coefficients
#if DEBUG_RK
      fprintf (stderr, "rk_step_tb: method coefficients\n");
#endif
      tb->method (tb);

      // objective function
#if DEBUG_RK
      fprintf (stderr, "rk_step_tb: objective function\n");
#endif
      o = tb->objective (tb);
      if (o < *tb->optimal)
        {
          g_mutex_lock (mutex);
          *tb->optimal = o;
          memcpy (tb->value_optimal, tb->random_data,
                  nfree * sizeof (long double));
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
  n = tb->nclimbings;
  for (i = 0; i < n; ++i)
    {
      memcpy (vo, tb->value_optimal, nfree * sizeof (long double));
      for (j = k = 0; j < nfree; ++j)
        {
          v = vo[j];
          tb->random_data[j] = v + is[j];
          tb->method (tb);
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
          tb->random_data[j] = fmaxl (0.L, v - is[j]);
          tb->method (tb);
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
          tb->random_data[j] = v;
        }

      // increase or reduce intervals if converging or not
      if (!k)
				f = 0.5L;
			else
				f = 1.2L;
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
void
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
int
rk_select (RK * rk,             ///< RK struct.
           unsigned int nsteps, ///< steps number.
           unsigned int order)  ///< accuracy order.
{
	const char *message[] = { _("Bad order"), _("Bad steps number") };
  Optimize *tb, *ac;
  unsigned int strong, code;
#if DEBUG_RK
  fprintf (stderr, "rk_select: start\n");
#endif
  strong = rk->strong;
  tb = rk->tb;
  tb->size = nsteps * (nsteps + 3) / 2 - 1;
  ac = rk->ac0;
  ac->size = nsteps * (nsteps + 1) - 2;
  ac->nfree = nsteps * (nsteps - 1) / 2;
  switch (nsteps)
    {
    case 2:
      tb->print = (OptimizePrint) rk_print_2;
      ac->minimum0 = minimum_ac_rk_2;
      ac->interval0 = interval_ac_rk_2;
      ac->random_type = random_ac_rk_2;
      ac->method = (OptimizeMethod) rk_ac_2;
      ac->objective = (OptimizeObjective) rk_objective_ac_2;
      switch (order)
        {
        case 2:
          tb->nfree = 1;
          if (strong)
            tb->print_maxima = rk_print_maxima_2_2;
          else
            tb->print_maxima = tb_print_maxima_2_2;
          tb->minimum0 = minimum_tb_2_2;
          tb->interval0 = interval_tb_2_2;
          tb->random_type = random_tb_2_2;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_2_2t;
              tb->objective = (OptimizeObjective) rk_objective_tb_2_2t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_2_2;
              tb->objective = (OptimizeObjective) rk_objective_tb_2_2;
            }
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 3:
      tb->print = (OptimizePrint) rk_print_3;
      ac->minimum0 = minimum_ac_rk_3;
      ac->interval0 = interval_ac_rk_3;
      ac->random_type = random_ac_rk_3;
      ac->method = (OptimizeMethod) rk_ac_3;
      ac->objective = (OptimizeObjective) rk_objective_ac_3;
      switch (order)
        {
        case 2:
          tb->nfree = 4;
          if (strong)
            tb->print_maxima = rk_print_maxima_3_2;
          else
            tb->print_maxima = tb_print_maxima_3_2;
          tb->minimum0 = minimum_tb_3_2;
          tb->interval0 = interval_tb_3_2;
          tb->random_type = random_tb_3_2;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_3_2t;
              tb->objective = (OptimizeObjective) rk_objective_tb_3_2t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_3_2;
              tb->objective = (OptimizeObjective) rk_objective_tb_3_2;
            }
          break;
        case 3:
          tb->nfree = 2;
          if (strong)
            tb->print_maxima = rk_print_maxima_3_3;
          else
            tb->print_maxima = tb_print_maxima_3_3;
          tb->minimum0 = minimum_tb_3_3;
          tb->interval0 = interval_tb_3_3;
          tb->random_type = random_tb_3_3;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_3_3t;
              tb->objective = (OptimizeObjective) rk_objective_tb_3_3t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_3_3;
              tb->objective = (OptimizeObjective) rk_objective_tb_3_3;
            }
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 4:
      tb->print = (OptimizePrint) rk_print_4;
      ac->minimum0 = minimum_ac_rk_4;
      ac->interval0 = interval_ac_rk_4;
      ac->random_type = random_ac_rk_4;
      ac->method = (OptimizeMethod) rk_ac_4;
      ac->objective = (OptimizeObjective) rk_objective_ac_4;
      switch (order)
        {
        case 2:
          tb->nfree = 8;
          if (strong)
            tb->print_maxima = rk_print_maxima_4_2;
          else
            tb->print_maxima = tb_print_maxima_4_2;
          tb->minimum0 = minimum_tb_4_2;
          tb->interval0 = interval_tb_4_2;
          tb->random_type = random_tb_4_2;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_4_2t;
              tb->objective = (OptimizeObjective) rk_objective_tb_4_2t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_4_2;
              tb->objective = (OptimizeObjective) rk_objective_tb_4_2;
            }
          break;
        case 3:
          tb->nfree = 6;
          if (strong)
            tb->print_maxima = rk_print_maxima_4_3;
          else
            tb->print_maxima = tb_print_maxima_4_3;
          tb->minimum0 = minimum_tb_4_3;
          tb->interval0 = interval_tb_4_3;
          tb->random_type = random_tb_4_3;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_4_3t;
              tb->objective = (OptimizeObjective) rk_objective_tb_4_3t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_4_3;
              tb->objective = (OptimizeObjective) rk_objective_tb_4_3;
            }
          break;
        case 4:
          tb->nfree = 2;
          if (strong)
            tb->print_maxima = rk_print_maxima_4_4;
          else
            tb->print_maxima = tb_print_maxima_4_4;
          tb->minimum0 = minimum_tb_4_4;
          tb->interval0 = interval_tb_4_4;
          tb->random_type = random_tb_4_4;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_4_4t;
              tb->objective = (OptimizeObjective) rk_objective_tb_4_4t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_4_4;
              tb->objective = (OptimizeObjective) rk_objective_tb_4_4;
            }
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 5:
      tb->print = (OptimizePrint) rk_print_5;
      ac->minimum0 = minimum_ac_rk_5;
      ac->interval0 = interval_ac_rk_5;
      ac->random_type = random_ac_rk_5;
      ac->method = (OptimizeMethod) rk_ac_5;
      ac->objective = (OptimizeObjective) rk_objective_ac_5;
      switch (order)
        {
        case 2:
          tb->nfree = 13;
          if (strong)
            tb->print_maxima = rk_print_maxima_5_2;
          else
            tb->print_maxima = tb_print_maxima_5_2;
          tb->minimum0 = minimum_tb_5_2;
          tb->interval0 = interval_tb_5_2;
          tb->random_type = random_tb_5_2;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_5_2t;
              tb->objective = (OptimizeObjective) rk_objective_tb_5_2t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_5_2;
              tb->objective = (OptimizeObjective) rk_objective_tb_5_2;
            }
          break;
        case 3:
          tb->nfree = 11;
          if (strong)
            tb->print_maxima = rk_print_maxima_5_3;
          else
            tb->print_maxima = tb_print_maxima_5_3;
          tb->minimum0 = minimum_tb_5_3;
          tb->interval0 = interval_tb_5_3;
          tb->random_type = random_tb_5_3;
          if (rk->time_accuracy)
            {
              tb->method = (OptimizeMethod) rk_tb_5_3t;
              tb->objective = (OptimizeObjective) rk_objective_tb_5_3t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_5_3;
              tb->objective = (OptimizeObjective) rk_objective_tb_5_3;
            }
          break;
        case 4:
          tb->nfree = 7;
          if (strong)
            tb->print_maxima = rk_print_maxima_5_4;
          else
            tb->print_maxima = tb_print_maxima_5_4;
          tb->minimum0 = minimum_tb_5_4;
          tb->interval0 = interval_tb_5_4;
          tb->random_type = random_tb_5_4;
          if (rk->time_accuracy)
            {
							if (rk->pair)
                --tb->nfree;
              tb->method = (OptimizeMethod) rk_tb_5_4t;
              tb->objective = (OptimizeObjective) rk_objective_tb_5_4t;
            }
          else
            {
              tb->method = (OptimizeMethod) rk_tb_5_4;
              tb->objective = (OptimizeObjective) rk_objective_tb_5_4;
            }
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 6:
      tb->print = (OptimizePrint) rk_print_6;
      ac->minimum0 = minimum_ac_rk_6;
      ac->interval0 = interval_ac_rk_6;
      ac->random_type = random_ac_rk_6;
      ac->method = (OptimizeMethod) rk_ac_6;
      ac->objective = (OptimizeObjective) rk_objective_ac_6;
      switch (order)
        {
        case 2:
          tb->nfree = 19;
          if (strong)
            tb->print_maxima = rk_print_maxima_6_2;
          else
            tb->print_maxima = tb_print_maxima_6_2;
          tb->minimum0 = minimum_tb_6_2;
          tb->interval0 = interval_tb_6_2;
          tb->random_type = random_tb_6_2;
	     	  if (rk->time_accuracy)
			      {
              tb->method = (OptimizeMethod) rk_tb_6_2t;
              tb->objective = (OptimizeObjective) rk_objective_tb_6_2t;
			      }
		      else
			      {
              tb->method = (OptimizeMethod) rk_tb_6_2;
              tb->objective = (OptimizeObjective) rk_objective_tb_6_2;
			      }
          break;
        case 3:
          tb->nfree = 17;
          if (strong)
            tb->print_maxima = rk_print_maxima_6_3;
          else
            tb->print_maxima = tb_print_maxima_6_3;
          tb->minimum0 = minimum_tb_6_3;
          tb->interval0 = interval_tb_6_3;
          tb->random_type = random_tb_6_3;
	     	  if (rk->time_accuracy)
			      {
              tb->method = (OptimizeMethod) rk_tb_6_3t;
              tb->objective = (OptimizeObjective) rk_objective_tb_6_3t;
			      }
		      else
			      {
              tb->method = (OptimizeMethod) rk_tb_6_3;
              tb->objective = (OptimizeObjective) rk_objective_tb_6_3;
			      }
          break;
        case 4:
          tb->nfree = 13;
          if (strong)
            tb->print_maxima = rk_print_maxima_6_4;
          else
            tb->print_maxima = tb_print_maxima_6_4;
          tb->minimum0 = minimum_tb_6_4;
          tb->interval0 = interval_tb_6_4;
          tb->random_type = random_tb_6_4;
	     	  if (rk->time_accuracy)
			      {
              tb->method = (OptimizeMethod) rk_tb_6_4t;
              tb->objective = (OptimizeObjective) rk_objective_tb_6_4t;
			      }
		      else
			      {
              tb->method = (OptimizeMethod) rk_tb_6_4;
              tb->objective = (OptimizeObjective) rk_objective_tb_6_4;
			      }
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    default:
      code = 1;
      goto exit_on_error;
    }
  if (rk->time_accuracy)
    --tb->nfree;

#if DEBUG_RK
  fprintf (stderr, "rk_select: end\n");
#endif
  return 1;

exit_on_error:
	error_message = g_strdup (message[code]);
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
rk_run (xmlNode * node, ///< XML node.
		    gsl_rng ** rng) ///< array of gsl_rng structs.
{
  RK rk[nthreads];
	char filename[32];
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
  value_optimal
    = (long double *) g_slice_alloc (nfree * sizeof (long double));
  optimize_create (tb, &optimal, value_optimal);
	if (rk->strong)
	  {
      ac = rk->ac0;
			node = node->children;
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
		}
  for (i = 1; i < nthreads; ++i)
    memcpy (rk + i, rk, sizeof (RK));
  j = rank * nthreads;
  for (i = 0; i < nthreads; ++i)
    rk_init (rk + i, rng[j + i]);

  // Method bucle
  printf ("Optimize bucle\n");
  rk_bucle_tb (rk);

  // Print the optimal coefficients
  printf ("Print the optimal coefficients\n");
  memcpy (tb->random_data, tb->value_optimal, nfree * sizeof (long double));
  tb->method (tb);
	if (rk->strong)
	  {
      memcpy (ac->random_data, ac->value_optimal,
              nfree2 * sizeof (long double));
      memcpy (rk->ac, ac, sizeof (Optimize));
      ac->method ((Optimize *) rk);
		}
  snprintf (filename, 32, "rk-%u-%u-%u-%u-%u.mc",
			      nsteps, order, rk->time_accuracy, rk->pair, rk->strong);
  file = fopen (filename, "w");
  tb->print ((Optimize *) rk, file);
  tb->print_maxima (file, nsteps, order);
  fclose (file);

#if PRINT_RANDOM
  fclose (file_random2);
#endif

  // Free memory
	if (rk->strong)
    g_slice_free1 (nfree2 * sizeof (long double), value_optimal2);
  for (i = 0; i < nthreads; ++i)
    rk_delete (rk + i);
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
