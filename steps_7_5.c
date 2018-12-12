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
 * \file steps_7_5.c
 * \brief Source file with common variables and functions to optimize the 7
 *   steps 5th order multi-steps method.
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
#include "steps.h"
#include "steps_7_4.h"
#include "steps_7_5.h"

///> array of minimum freedom degree values for the 7 steps 4th order 
///> multi-steps method.
const long double steps_minimum_7_5[8] = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of freedom degree intervals for the 7 steps 4th order multi-steps
///> method.
const long double steps_interval_7_5[8] = { 8.L, 8.L, 8.L, 8.L, 8.L, 8.L, 8.L, 1.L };

///> array of freedom degree random function types for the 7 steps 4th order
///> multi-steps method.
const unsigned int steps_random_7_5[8] = { 0, 1, 1, 1, 1, 1, 1, 0 };

/**
 * Function to get the coefficients on a 7 steps 4th order multi-steps method.
 */
void
steps_7_5 (Optimize * optimize) ///< Optimize struct.
{
  long double A[5], B[5], C[5], D[5], E[5], F[5];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  c4 (x) = r[4];
  c5 (x) = r[5];
  c6 (x) = r[6];
  a6 (x) = r[7];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = 1.L + a6 (x) * (6.L - c6 (x)) - c0 (x) * (1.L - a6 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 1.L - a6 (x) * (36.L - 12.L * c6 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = 1.L + a6 (x) * (216.L - 108.L * c6 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1.L - a6 (x) * (1296.L - 864.L * c6 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = 1.L + a6 (x) * (7776.L - 6480.L * c6 (x));
  solve_5 (A, B, C, D, E, F);
  a5 (x) = F[4];
  a4 (x) = F[3];
  a3 (x) = F[2];
  a2 (x) = F[1];
  a1 (x) = F[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x);
}
