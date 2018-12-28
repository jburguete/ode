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
 * \file steps_3_2.c
 * \brief Source file with common variables and functions to optimize the 3
 *   steps 2nd order multi-steps method.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <float.h>
#include <libxml/parser.h>
#include <glib.h>
#include <libintl.h>
#include <gsl/gsl_rng.h>
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "steps.h"
#include "steps_3_2.h"

///> array of minimum freedom degree values for the 3 steps 2nd order 
///> multi-steps method.
const long double steps_minimum_3_2[3] = { 0.L, 0.L, 0.L };

///> array of freedom degree intervals for the 3 steps 2nd order multi-steps
///> method.
const long double steps_interval_3_2[3] = { 4.L, 4.L, 4.L };

///> array of freedom degree random function types for the 3 steps 2nd order
///> multi-steps method.
const unsigned int steps_random_3_2[3] = { 0, 1, 1 };

/**
 * Function to get the coefficients on a 3 steps 2nd order multi-steps method.
 */
void
steps_3_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2 + c2 (x) - c0 (x);
  C[0] = 1.L - c0 (x);
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L;
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x);
}
