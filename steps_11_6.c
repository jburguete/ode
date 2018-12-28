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
 * \file steps_11_6.c
 * \brief Source file with common variables and functions to optimize the 11
 *   steps 6th order multi-steps method.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
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
#include "steps.h"
#include "steps_11_5.h"
#include "steps_11_6.h"

///> array of minimum freedom degree values for the 11 steps 6th order 
///> multi-steps method.
#if OPTIMIZE_STEPS_11_6
const long double steps_minimum_11_6[10]
  = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };
#else
const long double steps_minimum_11_6[15]
  = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 
		  0.L };
#endif

///> array of freedom degree intervals for the 11 steps 6th order multi-steps
///> method.
#if OPTIMIZE_STEPS_11_6
const long double steps_interval_11_6[10]
  = { 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 1.L, 1.L };
#else
const long double steps_interval_11_6[15]
  = { 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 12.L, 1.L,
		  1.L, 1.L };
#endif

///> array of freedom degree random function types for the 11 steps 6th order
///> multi-steps method.
#if OPTIMIZE_STEPS_11_6
const unsigned int steps_random_11_6[10] 
  = { 0, 1, 1, 1, 1, 1, 1, 1, 0, 1 };
#else
const unsigned int steps_random_11_6[15] 
  = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1 };
#endif

/**
 * Function to get the coefficients on a 11 steps 6th order multi-steps method.
 */
void
steps_11_6 (Optimize * optimize) ///< Optimize struct.
{
  long double A[6], B[6], C[6], D[6], E[6], F[6], G[6];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
#if OPTIMIZE_STEPS_11_6
	a2 (x) = a7 (x) = c2 (x) = c7 (x) = c10 (x) = 0.L;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c3 (x) = r[2];
  c4 (x) = r[3];
  c5 (x) = r[4];
  c6 (x) = r[5];
  c8 (x) = r[6];
  c9 (x) = r[7];
  a10 (x) = r[8];
  a9 (x) = r[9];
  a8 (x) = r[10];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -3.L + c3 (x) - c0 (x);
  C[0] = -4.L + c4 (x) - c0 (x);
  D[0] = -5.L + c5 (x) - c0 (x);
  E[0] = -6.L + c6 (x) - c0 (x);
  F[0] = -8.L + c8 (x) - c0 (x);
  G[0] = 1.L 
		     + a9 (x) * (9.L - c9 (x)) + a10 (x) * 10.L
         - c0 (x) * (1.L - a9 (x) - a10 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 9.L - 6.L * c3 (x);
  C[1] = 16.L - 8.L * c4 (x);
  D[1] = 25.L - 10.L * c5 (x);
  E[1] = 36.L - 12.L * c6 (x);
  F[1] = 64.L - 16.L * c8 (x);
  G[1] = 1.L
	 	     - a9 (x) * (81.L - 18.L * c9 (x)) - a10 (x) * 100.L;
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -27.L + 27.L * c3 (x);
  C[2] = -64.L + 48.L * c4 (x);
  D[2] = -125.L + 75.L * c5 (x);
  E[2] = -216.L + 108.L * c6 (x);
  F[2] = -512.L + 192.L * c8 (x);
  G[2] = 1.L
				 + a9 (x) * (729.L - 243.L * c9 (x))
				 + a10 (x) * 1000.L;
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 81.L - 108.L * c3 (x);
  C[3] = 256.L - 256.L * c4 (x);
  D[3] = 625.L - 500.L * c5 (x);
  E[3] = 1296.L - 864.L * c6 (x);
  F[3] = 4096.L - 2048.L * c8 (x);
  G[3] = 1.L 
				 - a9 (x) * (6561.L - 2916.L * c9 (x))
				 - a10 (x) * 10000.L;
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -243.L + 405.L * c3 (x);
  C[4] = -1024.L + 1280.L * c4 (x);
  D[4] = -3125.L + 3125.L * c5 (x);
  E[4] = -7776.L + 6480.L * c6 (x);
  F[4] = -32768.L + 20480.L * c8 (x);
  G[4] = 1.L
				 + a9 (x) * (59049.L - 32805.L * c9 (x))
				 + a10 (x) * 100000.L;
  A[5] = 1.L - 6.L * c1 (x);
  B[5] = 729.L - 1458.L * c3 (x);
  C[5] = 4096.L - 6144.L * c4 (x);
  D[5] = 15625.L - 18750.L * c5 (x);
  E[5] = 46656.L - 46656.L * c6 (x);
  F[5] = 262144.L - 196608.L * c8 (x);
  G[5] = 1.L
				 - a9 (x) * (531441.L - 354294.L * c9 (x))
				 - a10 (x) * 1000000.L;
  solve_6 (A, B, C, D, E, F, G);
  a8 (x) = G[5];
  a6 (x) = G[4];
  a5 (x) = G[3];
  a4 (x) = G[2];
  a3 (x) = G[1];
  a1 (x) = G[0];
  a0 (x) = 1.L - a1 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x)
		       - a8 (x) - a9 (x) - a10 (x);
#else
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  c4 (x) = r[4];
  c5 (x) = r[5];
  c6 (x) = r[6];
  c7 (x) = r[7];
  c8 (x) = r[8];
  c9 (x) = r[9];
  c10 (x) = r[10];
  a10 (x) = r[11];
  a9 (x) = r[12];
  a8 (x) = r[13];
  a7 (x) = r[14];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = -6.L + c6 (x) - c0 (x);
  G[0] = 1.L + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x)) 
		     + a9 (x) * (9.L - c9 (x)) + a10 (x) * (10.L - c10 (x))
         - c0 (x) * (1.L - a7 (x) - a8 (x) - a9 (x) - a10 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 36.L - 12.L * c6 (x);
  G[1] = 1.L - a7 (x) * (49.L - 14.L * c7 (x)) - a8 (x) * (64.L - 16.L * c8 (x))
	 	     - a9 (x) * (81.L - 18.L * c9 (x)) - a10 (x) * (100.L - 20.L * c10 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = -216.L + 108.L * c6 (x);
  G[2] = 1.L + a7 (x) * (343.L - 147.L * c7 (x)) 
		     + a8 (x) * (512.L - 192.L * c8 (x))
				 + a9 (x) * (729.L - 243.L * c9 (x))
				 + a10 (x) * (1000.L - 300.L * c10 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1296.L - 864.L * c6 (x);
  G[3] = 1.L - a7 (x) * (2401.L - 1372.L * c7 (x))
				 - a8 (x) * (4096.L - 2048.L * c8 (x))
				 - a9 (x) * (6561.L - 2916.L * c9 (x))
				 - a10 (x) * (10000.L - 4000.L * c10 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = -7776.L + 6480.L * c6 (x);
  G[4] = 1.L + a7 (x) * (16807.L - 12005.L * c7 (x))
				 + a8 (x) * (32768.L - 20480.L * c8 (x))
				 + a9 (x) * (59049.L - 32805.L * c9 (x))
				 + a10 (x) * (100000.L - 50000.L * c10 (x));
  A[5] = 1.L - 6.L * c1 (x);
  B[5] = 64.L - 192.L * c2 (x);
  C[5] = 729.L - 1458.L * c3 (x);
  D[5] = 4096.L - 6144.L * c4 (x);
  E[5] = 15625.L - 18750.L * c5 (x);
  F[5] = 46656.L - 46656.L * c6 (x);
  G[5] = 1.L - a7 (x) * (117649.L - 100842.L * c7 (x))
				 - a8 (x) * (262144.L - 196608.L * c8 (x))
				 - a9 (x) * (531441.L - 354294.L * c9 (x))
				 - a10 (x) * (1000000.L - 600000.L * c10 (x));
  solve_6 (A, B, C, D, E, F, G);
  a6 (x) = G[5];
  a5 (x) = G[4];
  a4 (x) = G[3];
  a3 (x) = G[2];
  a2 (x) = G[1];
  a1 (x) = G[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
		       - a8 (x) - a9 (x) - a10 (x);
#endif
}
