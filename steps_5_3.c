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
 * \file steps_5_3.c
 * \brief Source file with common variables and functions to optimize the 5
 *   steps 3rd order multi-steps method.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <float.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "steps.h"
#include "steps_5_2.h"
#include "steps_5_3.h"

///> array of minimum freedom degree values for the 5 steps 3rd order 
///> multi-steps method.
const long double steps_minimum_5_3[6] = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L };

///> array of freedom degree intervals for the 5 steps 3rd order multi-steps
///> method.
const long double steps_interval_5_3[6] = { 6.L, 6.L, 6.L, 6.L, 6.L, 1.L };

///> array of freedom degree random function types for the 5 steps 3rd order
///> multi-steps method.
const unsigned int steps_random_5_3[6] = { 0, 1, 1, 1, 1, 0 };

/**
 * Function to get the coefficients on a 5 steps 3rd order multi-steps method.
 */
void
steps_5_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  c4 (x) = r[4];
  a4 (x) = r[5];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) - c0 (x) * (1.L - a4 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L + a4 (x) * (64.L - 48.L * c4 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
}
