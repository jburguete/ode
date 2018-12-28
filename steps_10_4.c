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
 * \file steps_10_4.c
 * \brief Source file with common variables and functions to optimize the 10
 *   steps 4th order multi-steps method.
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
#include "steps_10_3.h"
#include "steps_10_4.h"

///> array of minimum freedom degree values for the 10 steps 4th order 
///> multi-steps method.
const long double steps_minimum_10_4[15]
  = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L, 0.L,
      0.L	};

///> array of freedom degree intervals for the 10 steps 4th order multi-steps
///> method.
const long double steps_interval_10_4[15]
  = { 11.L, 11.L, 11.L, 11.L, 11.L, 11.L, 11.L, 11.L, 11.L, 11.L, 1.L, 1.L, 1.L,
	 	  1.L };

///> array of freedom degree random function types for the 10 steps 4th order
///> multi-steps method.
const unsigned int steps_random_10_4[15] 
  = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };

/**
 * Function to get the coefficients on a 10 steps 4th order multi-steps method.
 */
void
steps_10_4 (Optimize * optimize) ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
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
  c9 (x) = r[9];
  a9 (x) = r[10];
  a8 (x) = r[11];
  a7 (x) = r[12];
  a6 (x) = r[13];
  a5 (x) = r[14];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L + a5 (x) * (5.L - c5 (x))  + a6 (x) * (6.L - c6 (x)) 
		     + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x)) 
				 + a9 (x) * (9.L - c9 (x))
         - c0 (x) * (1.L - a5 (x) - a6 (x) - a7 (x) - a8 (x) - a9 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 1.L - a5 (x) * (25.L - 10.L * c5 (x)) - a6 (x) * (36.L - 12.L * c6 (x))
	 	     - a7 (x) * (49.L - 14.L * c7 (x)) - a8 (x) * (64.L - 16.L * c8 (x)) 
				 - a9 (x) * (81.L - 18.L * c9 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = 1.L + a5 (x) * (125.L - 75.L * c5 (x)) 
		     + a6 (x) * (216.L - 108.L * c6 (x)) + a7 (x) * (343.L - 147.L * c7 (x))
				 + a8 (x) * (512.L - 192.L * c8 (x))
				 + a9 (x) * (729.L - 243.L * c9 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x)) 
		     - a6 (x) * (1296.L - 864.L * c6 (x)) 
				 - a7 (x) * (2401.L - 1372.L * c7 (x))
				 - a8 (x) * (4096.L - 2048.L * c8 (x))
				 - a9 (x) * (6561.L - 2916.L * c9 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
		       - a8 (x) - a9 (x);
}
