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
 * \file steps_9_6.c
 * \brief Source file with common variables and functions to optimize the 9
 *   steps 6th order multi-steps method.
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
#include "steps_9_5.h"
#include "steps_9_6.h"

///> array of minimum freedom degree values for the 9 steps 6th order 
///> multi-steps method.
const long double steps_minimum_9_6[11]
  = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of freedom degree intervals for the 9 steps 6th order multi-steps
///> method.
const long double steps_interval_9_6[11]
  = { 10.L, 10.L, 10.L, 10.L, 10.L, 10.L, 10.L, 10.L, 1.L };

///> array of freedom degree random function types for the 9 steps 6th order
///> multi-steps method.
const unsigned int steps_random_9_6[11] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 };

/**
 * Function to get the coefficients on a 8 steps 5th order multi-steps method.
 */
void
steps_9_6 (Optimize * optimize) ///< Optimize struct.
{
  long double A[6], B[6], C[6], D[6], E[6], F[6], G[6];
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
  c7 (x) = r[7];
  c8 (x) = r[8];
  a8 (x) = r[9];
  a7 (x) = r[10];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = -6.L + c6 (x) - c0 (x);
  G[0] = 1.L + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
		     - c0 (x) * (1.L - a7 (x) - a8 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 36.L - 12.L * c6 (x);
  G[1] = 1.L - a7 (x) * (49.L - 14.L * c7 (x)) 
		     - a8 (x) * (64.L - 16.L * c8 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = -216.L + 108.L * c6 (x);
  G[2] = 1.L + a7 (x) * (343.L - 147.L * c7 (x))
		     + a8 (x) * (512.L - 192.L * c8 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1296.L - 864.L * c6 (x);
  G[3] = 1.L - a7 (x) * (2401.L - 1372.L * c7 (x))
		     - a8 (x) * (4096.L - 2048.L * c8 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = -7776.L + 6480.L * c6 (x);
  G[4] = 1.L + a7 (x) * (16807.L - 12005.L * c7 (x))
		     + a8 (x) * (32768.L - 20480.L * c8 (x));
  A[5] = 1.L - 6.L * c1 (x);
  B[5] = 64.L - 192.L * c2 (x);
  C[5] = 729.L - 1458.L * c3 (x);
  D[5] = 4096.L - 6144.L * c4 (x);
  E[5] = 15625.L - 18750.L * c5 (x);
  F[5] = 46656.L - 46656.L * c6 (x);
  G[5] = 1.L - a7 (x) * (117649.L - 100842.L * c7 (x))
		     - a8 (x) * (262144.L - 196608.L * c8 (x));
  solve_6 (A, B, C, D, E, F, G);
  a6 (x) = G[5];
  a5 (x) = G[4];
  a4 (x) = G[3];
  a3 (x) = G[2];
  a2 (x) = G[1];
  a1 (x) = G[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x);
}

