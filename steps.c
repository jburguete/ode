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
 * \file steps.c
 * \brief Source file with common variables and functions to optimize
 *   multi-steps methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
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

#define c0(x) x[0]
///< c0 multi-steps coefficient.
#define a0(x) x[1]
///< a0 multi-steps coefficient.
#define c1(x) x[2]
///< c1 multi-steps coefficient.
#define a1(x) x[3]
///< a1 multi-steps coefficient.
#define c2(x) x[4]
///< c2 multi-steps coefficient.
#define a2(x) x[5]
///< a2 multi-steps coefficient.
#define c3(x) x[6]
///< c3 multi-steps coefficient.
#define a3(x) x[7]
///< a3 multi-steps coefficient.
#define c4(x) x[8]
///< c4 multi-steps coefficient.
#define a4(x) x[9]
///< a4 multi-steps coefficient.
#define c5(x) x[10]
///< c5 multi-steps coefficient.
#define a5(x) x[11]
///< a5 multi-steps coefficient.
#define c6(x) x[12]
///< c6 multi-steps coefficient.
#define a6(x) x[13]
///< a6 multi-steps coefficient.
#define c7(x) x[14]
///< c7 multi-steps coefficient.
#define a7(x) x[15]
///< a7 multi-steps coefficient.
#define c8(x) x[16]
///< c8 multi-steps coefficient.
#define a8(x) x[17]
///< a8 multi-steps coefficient.
#define c9(x) x[18]
///< c9 multi-steps coefficient.
#define a9(x) x[19]
///< a9 multi-steps coefficient.
#define c10(x) x[20]
///< c10 multi-steps coefficient.
#define a10(x) x[21]
///< a10 multi-steps coefficient.
#define b0(x) (a0(x)*c0(x))
///< b0 multi-steps coefficient.
#define b1(x) (a1(x)*c1(x))
///< b1 multi-steps coefficient.
#define b2(x) (a2(x)*c2(x))
///< b2 multi-steps coefficient.
#define b3(x) (a3(x)*c3(x))
///< b3 multi-steps coefficient.
#define b4(x) (a4(x)*c4(x))
///< b4 multi-steps coefficient.
#define b5(x) (a5(x)*c5(x))
///< b5 multi-steps coefficient.
#define b6(x) (a6(x)*c6(x))
///< b6 multi-steps coefficient.
#define b7(x) (a7(x)*c7(x))
///< b7 multi-steps coefficient.
#define b8(x) (a8(x)*c8(x))
///< b8 multi-steps coefficient.
#define b9(x) (a9(x)*c9(x))
///< b9 multi-steps coefficient.
#define b10(x) (a10(x)*c10(x))
///< b10 multi-steps coefficient.

/**
 * Function to get the coefficients on a 3 steps 2nd order multi-steps method.
 */
static int
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
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 3 steps 3rd order multi-steps method.
 */
static int
steps_3_3 (Optimize * optimize) ///< Optimize struct.
{
  long double *x, *r;
  long double A0, A1, B0, B1;
  x = optimize->coefficient;
  r = optimize->random_data;
  c1 (x) = r[0];
  c2 (x) = r[1];
  A0 = 1.L - 2.L * c1 (x);
  B0 = 4.L - 4.L * c2 (x);
  A1 = -1.L + 3.L * c1 (x);
  B1 = -8.L + 12.L * c2 (x);
  a2 (x) = (A0 - A1) / (A0 * B1 - A1 * B0);
  a1 (x) = (1.L - B0 * a2 (x)) / A0;
  a0 (x) = 1.L - a1 (x) - a2 (x);
  c0 (x) = (1.L + a1 (x) * (1.L - c1 (x)) + a2 (x) * (2.L - c2 (x))) / a0 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (c0 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 4 steps 2nd order multi-steps method.
 */
static int
steps_4_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  a3 (x) = r[4];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2 + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) - c0 (x) * (1.L - a3 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 4 steps 3rd order multi-steps method.
 */
static int
steps_4_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.l + c3 (x) - c0 (x);
  D[0] = 1.L - c0 (x);
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L;
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L;
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 5 steps 2nd order multi-steps method.
 */
static int
steps_5_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  c4 (x) = r[4];
  a4 (x) = r[5];
  a3 (x) = r[6];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    - c0 (x) * (1.L - a3 (x) - a4 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x)) - a4 (x) * (16.L - 8.L * c4 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 5 steps 3rd order multi-steps method.
 */
static int
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
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 5 steps 4th order multi-steps method.
 */
static int
steps_5_4 (Optimize * optimize) ///< Optimize struct.
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
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L - c0 (x);
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 1.L;
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = 1.L;
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L;
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 6 steps 2nd order multi-steps method.
 */
static int
steps_6_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  c0 (x) = r[0];
  c1 (x) = r[1];
  c2 (x) = r[2];
  c3 (x) = r[3];
  c4 (x) = r[4];
  c5 (x) = r[5];
  a5 (x) = r[6];
  a4 (x) = r[7];
  a3 (x) = r[8];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    + a5 (x) * (5.L - c5 (x)) - c0 (x) * (1.L - a3 (x) - a4 (x) - a5 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x))
    - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 6 steps 3rd order multi-steps method.
 */
static int
steps_6_3 (Optimize * optimize) ///< Optimize struct.
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
  c5 (x) = r[5];
  a5 (x) = r[6];
  a4 (x) = r[7];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) + a5 (x) * (5.L - c5 (x))
    - c0 (x) * (1.L - a4 (x) - a5 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] =
    1.L + a4 (x) * (64.L - 48.L * c4 (x)) + a5 (x) * (125.L - 75.L * c5 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 6 steps 4th order multi-steps method.
 */
static int
steps_6_4 (Optimize * optimize) ///< Optimize struct.
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
  a5 (x) = r[6];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L + a5 (x) * (5.L - c5 (x)) - c0 (x) * (1.L - a5 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 1.L - a5 (x) * (25.L - 10.L * c5 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = 1.L + a5 (x) * (125.L - 75.L * c5 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 6 steps 5th order multi-steps method.
 */
static int
steps_6_5 (Optimize * optimize) ///< Optimize struct.
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
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = 1.L - c0 (x);
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 1.L;
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = 1.L;
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1.L;
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = 1.L;
  solve_5 (A, B, C, D, E, F);
  a5 (x) = F[4];
  a4 (x) = F[3];
  a3 (x) = F[2];
  a2 (x) = F[1];
  a1 (x) = F[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 7 steps 2nd order multi-steps method.
 */
static int
steps_7_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
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
  a5 (x) = r[8];
  a4 (x) = r[9];
  a3 (x) = r[10];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    - c0 (x) * (1.L - a3 (x) - a4 (x) - a5 (x) - a6 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x))
    - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36 - 12.L * c6 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 7 steps 3rd order multi-steps method.
 */
static int
steps_7_3 (Optimize * optimize) ///< Optimize struct.
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
  c5 (x) = r[5];
  c6 (x) = r[6];
  a6 (x) = r[7];
  a5 (x) = r[8];
  a4 (x) = r[9];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) + a5 (x) * (5.L - c5 (x))
    + a6 (x) * (6.L - c6 (x)) - c0 (x) * (1.L - a4 (x) - a5 (x) - a6 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x))
    - a5 (x) * (25.L - 10.L * c5 (x)) - a6 (x) * (36.L - 12.L * c6 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L + a4 (x) * (64.L - 48.L * c4 (x))
    + a5 (x) * (125.L - 75.L * c5 (x)) + a6 (x) * (216.L - 108.L * c6 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 7 steps 4th order multi-steps method.
 */
static int
steps_7_4 (Optimize * optimize) ///< Optimize struct.
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
  a6 (x) = r[7];
  a5 (x) = r[8];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    - c0 (x) * (1.L - a5 (x) - a6 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] =
    1.L - a5 (x) * (25.L - 10.L * c5 (x)) - a6 (x) * (36.L - 12.L * c6 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] =
    1.L + a5 (x) * (125.L - 75.L * c5 (x)) + a6 (x) * (216.L - 108.L * c6 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x))
    - a6 (x) * (1296.L - 864.L * c6 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 7 steps 5th order multi-steps method.
 */
static int
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
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 8 steps 2nd order multi-steps method.
 */
static int
steps_8_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
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
  a7 (x) = r[8];
  a6 (x) = r[9];
  a5 (x) = r[10];
  a4 (x) = r[11];
  a3 (x) = r[12];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x))
    - c0 (x) * (1.L - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x))
    - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x)) -
    a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 8 steps 3rd order multi-steps method.
 */
static int
steps_8_3 (Optimize * optimize) ///< Optimize struct.
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
  c5 (x) = r[5];
  c6 (x) = r[6];
  c7 (x) = r[7];
  a7 (x) = r[8];
  a6 (x) = r[9];
  a5 (x) = r[10];
  a4 (x) = r[11];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) + a5 (x) * (5.L - c5 (x))
    + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    - c0 (x) * (1.L - a4 (x) - a5 (x) - a6 (x) - a7 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x))
    - a5 (x) * (25.L - 10.L * c5 (x)) - a6 (x) * (36.L - 12.L * c6 (x))
    - a7 (x) * (49.L - 14.L * c7 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L + a4 (x) * (64.L - 48.L * c4 (x))
    + a5 (x) * (125.L - 75.L * c5 (x)) + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 8 steps 4th order multi-steps method.
 */
static int
steps_8_4 (Optimize * optimize) ///< Optimize struct.
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
  a7 (x) = r[8];
  a6 (x) = r[9];
  a5 (x) = r[10];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x)) - c0 (x) * (1.L - a5 (x) - a6 (x) - a7 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 1.L - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = 1.L + a5 (x) * (125.L - 75.L * c5 (x))
    + a6 (x) * (216.L - 108.L * c6 (x)) + a7 (x) * (343.L - 147.L * c7 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x))
    - a6 (x) * (1296.L - 864.L * c6 (x)) - a7 (x) * (2401.L - 1372.L * c7 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 8 steps 5th order multi-steps method.
 */
static int
steps_8_5 (Optimize * optimize) ///< Optimize struct.
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
  c7 (x) = r[7];
  a7 (x) = r[8];
  a6 (x) = r[9];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = 1.L + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    - c0 (x) * (1.L - a6 (x) - a7 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] =
    1.L - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] =
    1.L + a6 (x) * (216.L - 108.L * c6 (x)) + a7 (x) * (343.L - 147.L * c7 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1.L - a6 (x) * (1296.L - 864.L * c6 (x))
    - a7 (x) * (2401.L - 1372.L * c7 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = 1.L + a6 (x) * (7776.L - 6480.L * c6 (x))
    + a7 (x) * (16807.L - 12005.L * c7 (x));
  solve_5 (A, B, C, D, E, F);
  a5 (x) = F[4];
  a4 (x) = F[3];
  a3 (x) = F[2];
  a2 (x) = F[1];
  a1 (x) = F[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 8 steps 6th order multi-steps method.
 */
static int
steps_8_6 (Optimize * optimize) ///< Optimize struct.
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
  a7 (x) = r[8];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = -6.L + c6 (x) - c0 (x);
  G[0] = 1.L + a7 (x) * (7.L - c7 (x)) - c0 (x) * (1.L - a7 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 36.L - 12.L * c6 (x);
  G[1] = 1.L - a7 (x) * (49.L - 14.L * c7 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = -216.L + 108.L * c6 (x);
  G[2] = 1.L + a7 (x) * (343.L - 147.L * c7 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1296.L - 864.L * c6 (x);
  G[3] = 1.L - a7 (x) * (2401.L - 1372.L * c7 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = -7776.L + 6480.L * c6 (x);
  G[4] = 1.L + a7 (x) * (16807.L - 12005.L * c7 (x));
  A[5] = 1.L - 6.L * c1 (x);
  B[5] = 64.L - 192.L * c2 (x);
  C[5] = 729.L - 1458.L * c3 (x);
  D[5] = 4096.L - 6144.L * c4 (x);
  E[5] = 15625.L - 18750.L * c5 (x);
  F[5] = 46656.L - 46656.L * c6 (x);
  G[5] = 1.L - a7 (x) * (117649.L - 100842.L * c7 (x));
  solve_6 (A, B, C, D, E, F, G);
  a6 (x) = G[5];
  a5 (x) = G[4];
  a4 (x) = G[3];
  a3 (x) = G[2];
  a2 (x) = G[1];
  a1 (x) = G[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)) || isnan (a6 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 9 steps 2nd order multi-steps method.
 */
static int
steps_9_2 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
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
  a6 (x) = r[11];
  a5 (x) = r[12];
  a4 (x) = r[13];
  a3 (x) = r[14];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
    - c0 (x) * (1.L - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x) - a8 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x))
    - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 9 steps 3rd order multi-steps method.
 */
static int
steps_9_3 (Optimize * optimize) ///< Optimize struct.
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
  c5 (x) = r[5];
  c6 (x) = r[6];
  c7 (x) = r[7];
  c8 (x) = r[8];
  a8 (x) = r[9];
  a7 (x) = r[10];
  a6 (x) = r[11];
  a5 (x) = r[12];
  a4 (x) = r[13];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) + a5 (x) * (5.L - c5 (x))
    + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    + a8 (x) * (8.L - c8 (x))
    - c0 (x) * (1.L - a4 (x) - a5 (x) - a6 (x) - a7 (x) - a8 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x))
    - a5 (x) * (25.L - 10.L * c5 (x)) - a6 (x) * (36.L - 12.L * c6 (x))
    - a7 (x) * (49.L - 14.L * c7 (x)) - a8 (x) * (64.L - 16.L * c8 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L + a4 (x) * (64.L - 48.L * c4 (x))
    + a5 (x) * (125.L - 75.L * c5 (x)) + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x)) + a8 (x) * (512.L - 192.L * c8 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 9 steps 4th order multi-steps method.
 */
static int
steps_9_4 (Optimize * optimize) ///< Optimize struct.
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
  a8 (x) = r[9];
  a7 (x) = r[10];
  a6 (x) = r[11];
  a5 (x) = r[12];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
    - c0 (x) * (1.L - a5 (x) - a6 (x) - a7 (x) - a8 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 1.L - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = 1.L + a5 (x) * (125.L - 75.L * c5 (x))
    + a6 (x) * (216.L - 108.L * c6 (x)) + a7 (x) * (343.L - 147.L * c7 (x))
    + a8 (x) * (512.L - 192.L * c8 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x))
    - a6 (x) * (1296.L - 864.L * c6 (x)) - a7 (x) * (2401.L - 1372.L * c7 (x))
    - a8 (x) * (4096.L - 2048.L * c8 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 9 steps 5th order multi-steps method.
 */
static int
steps_9_5 (Optimize * optimize) ///< Optimize struct.
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
  c7 (x) = r[7];
  c8 (x) = r[8];
  a8 (x) = r[9];
  a7 (x) = r[10];
  a6 (x) = r[11];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = 1.L + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    + a8 (x) * (8.L - c8 (x)) - c0 (x) * (1.L - a6 (x) - a7 (x) - a8 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 1.L - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = 1.L + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x)) + a8 (x) * (512.L - 192.L * c8 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1.L - a6 (x) * (1296.L - 864.L * c6 (x))
    - a7 (x) * (2401.L - 1372.L * c7 (x)) - a8 (x) * (4096.L - 2048.L * c8 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = 1.L + a6 (x) * (7776.L - 6480.L * c6 (x))
    + a7 (x) * (16807.L - 12005.L * c7 (x))
    + a8 (x) * (32768.L - 20480.L * c8 (x));
  solve_5 (A, B, C, D, E, F);
  a5 (x) = F[4];
  a4 (x) = F[3];
  a3 (x) = F[2];
  a2 (x) = F[1];
  a1 (x) = F[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 9 steps 6th order multi-steps method.
 */
static int
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
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)) || isnan (a6 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 10 steps 2nd order multi-steps method.
 */
static int
steps_10_2 (Optimize * optimize)        ///< Optimize struct.
{
  long double A[2], B[2], C[2];
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
  a4 (x) = r[15];
  a3 (x) = r[16];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
    + a9 (x) * (9.L - c9 (x))
    - c0 (x) * (1.L - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x) - a8 (x)
                - a9 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x))
    - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x)) - a9 (x) * (81.L - 18.L * c9 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 10 steps 3rd order multi-steps method.
 */
static int
steps_10_3 (Optimize * optimize)        ///< Optimize struct.
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
  a4 (x) = r[15];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) + a5 (x) * (5.L - c5 (x))
    + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    + a8 (x) * (8.L - c8 (x)) + a9 (x) * (9.L - c9 (x))
    - c0 (x) * (1.L - a4 (x) - a5 (x) - a6 (x) - a7 (x) - a8 (x) - a9 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x)) - a9 (x) * (81.L - 18.L * c9 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L + a4 (x) * (64.L - 48.L * c4 (x))
    + a5 (x) * (125.L - 75.L * c5 (x)) + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x)) + a8 (x) * (512.L - 192.L * c8 (x))
    + a9 (x) * (729.L - 243.L * c9 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 10 steps 4th order multi-steps method.
 */
static int
steps_10_4 (Optimize * optimize)        ///< Optimize struct.
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
  E[0] = 1.L + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
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
    + a8 (x) * (512.L - 192.L * c8 (x)) + a9 (x) * (729.L - 243.L * c9 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x))
    - a6 (x) * (1296.L - 864.L * c6 (x))
    - a7 (x) * (2401.L - 1372.L * c7 (x))
    - a8 (x) * (4096.L - 2048.L * c8 (x)) - a9 (x) * (6561.L - 2916.L * c9 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 10 steps 5th order multi-steps method.
 */
static int
steps_10_5 (Optimize * optimize)        ///< Optimize struct.
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
  c7 (x) = r[7];
  c8 (x) = r[8];
  c9 (x) = r[9];
  a9 (x) = r[10];
  a8 (x) = r[11];
  a7 (x) = r[12];
  a6 (x) = r[13];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = 1.L + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    + a8 (x) * (8.L - c8 (x)) + a9 (x) * (9.L - c9 (x))
    - c0 (x) * (1.L - a6 (x) - a7 (x) - a8 (x) - a9 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 1.L - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x)) - a9 (x) * (81.L - 18.L * c9 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = 1.L + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x)) + a8 (x) * (512.L - 192.L * c8 (x))
    + a9 (x) * (729.L - 243.L * c9 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1.L - a6 (x) * (1296.L - 864.L * c6 (x))
    - a7 (x) * (2401.L - 1372.L * c7 (x))
    - a8 (x) * (4096.L - 2048.L * c8 (x)) - a9 (x) * (6561.L - 2916.L * c9 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = 1.L + a6 (x) * (7776.L - 6480.L * c6 (x))
    + a7 (x) * (16807.L - 12005.L * c7 (x))
    + a8 (x) * (32768.L - 20480.L * c8 (x))
    + a9 (x) * (59049.L - 32805.L * c9 (x));
  solve_5 (A, B, C, D, E, F);
  a5 (x) = F[4];
  a4 (x) = F[3];
  a3 (x) = F[2];
  a2 (x) = F[1];
  a1 (x) = F[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 10 steps 6th order multi-steps method.
 */
static int
steps_10_6 (Optimize * optimize)        ///< Optimize struct.
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
  c9 (x) = r[9];
  a9 (x) = r[10];
  a8 (x) = r[11];
  a7 (x) = r[12];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = -6.L + c6 (x) - c0 (x);
  G[0] = 1.L + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
    + a9 (x) * (9.L - c9 (x)) - c0 (x) * (1.L - a7 (x) - a8 (x) - a9 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 36.L - 12.L * c6 (x);
  G[1] = 1.L - a7 (x) * (49.L - 14.L * c7 (x)) - a8 (x) * (64.L - 16.L * c8 (x))
    - a9 (x) * (81.L - 18.L * c9 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = -216.L + 108.L * c6 (x);
  G[2] = 1.L + a7 (x) * (343.L - 147.L * c7 (x))
    + a8 (x) * (512.L - 192.L * c8 (x)) + a9 (x) * (729.L - 243.L * c9 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1296.L - 864.L * c6 (x);
  G[3] = 1.L - a7 (x) * (2401.L - 1372.L * c7 (x))
    - a8 (x) * (4096.L - 2048.L * c8 (x)) - a9 (x) * (6561.L - 2916.L * c9 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = -7776.L + 6480.L * c6 (x);
  G[4] = 1.L + a7 (x) * (16807.L - 12005.L * c7 (x))
    + a8 (x) * (32768.L - 20480.L * c8 (x))
    + a9 (x) * (59049.L - 32805.L * c9 (x));
  A[5] = 1.L - 6.L * c1 (x);
  B[5] = 64.L - 192.L * c2 (x);
  C[5] = 729.L - 1458.L * c3 (x);
  D[5] = 4096.L - 6144.L * c4 (x);
  E[5] = 15625.L - 18750.L * c5 (x);
  F[5] = 46656.L - 46656.L * c6 (x);
  G[5] = 1.L - a7 (x) * (117649.L - 100842.L * c7 (x))
    - a8 (x) * (262144.L - 196608.L * c8 (x))
    - a9 (x) * (531441.L - 354294.L * c9 (x));
  solve_6 (A, B, C, D, E, F, G);
  a6 (x) = G[5];
  a5 (x) = G[4];
  a4 (x) = G[3];
  a3 (x) = G[2];
  a2 (x) = G[1];
  a1 (x) = G[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)) || isnan (a6 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 11 steps 2nd order multi-steps method.
 */
static int
steps_11_2 (Optimize * optimize)        ///< Optimize struct.
{
  long double A[2], B[2], C[2];
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
  c10 (x) = r[10];
  a10 (x) = r[11];
  a9 (x) = r[12];
  a8 (x) = r[13];
  a7 (x) = r[14];
  a6 (x) = r[15];
  a5 (x) = r[16];
  a4 (x) = r[17];
  a3 (x) = r[18];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = 1.L + a3 (x) * (3.L - c3 (x)) + a4 (x) * (4.L - c4 (x))
    + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
    + a9 (x) * (9.L - c9 (x)) + a10 (x) * (10.L - c10 (x))
    - c0 (x) * (1.L - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x) - a8 (x)
                - a9 (x) - a10 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 1.L - a3 (x) * (9.L - 6.L * c3 (x))
    - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x)) - a9 (x) * (81.L - 18.L * c9 (x))
    - a10 (x) * (100.L - 20.L * c10 (x));
  solve_2 (A, B, C);
  a2 (x) = C[1];
  a1 (x) = C[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x) - a10 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 11 steps 3rd order multi-steps method.
 */
static int
steps_11_3 (Optimize * optimize)        ///< Optimize struct.
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
  a6 (x) = r[15];
  a5 (x) = r[16];
  a4 (x) = r[17];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = 1.L + a4 (x) * (4.L - c4 (x)) + a5 (x) * (5.L - c5 (x))
    + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    + a8 (x) * (8.L - c8 (x)) + a9 (x) * (9.L - c9 (x))
    + a10 (x) * (10.L - c8 (x))
    - c0 (x) * (1.L - a4 (x) - a5 (x) - a6 (x) - a7 (x) - a8 (x) - a9 (x)
                - a10 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 1.L - a4 (x) * (16.L - 8.L * c4 (x)) - a5 (x) * (25.L - 10.L * c5 (x))
    - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x)) - a9 (x) * (81.L - 18.L * c9 (x))
    - a10 (x) * (100.L - 20.L * c10 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = 1.L + a4 (x) * (64.L - 48.L * c4 (x))
    + a5 (x) * (125.L - 75.L * c5 (x)) + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x)) + a8 (x) * (512.L - 192.L * c8 (x))
    + a9 (x) * (729.L - 243.L * c9 (x)) + a10 (x) * (1000.L - 300.L * c10 (x));
  solve_3 (A, B, C, D);
  a3 (x) = D[2];
  a2 (x) = D[1];
  a1 (x) = D[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x) - a10 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 11 steps 4th order multi-steps method.
 */
static int
steps_11_4 (Optimize * optimize)        ///< Optimize struct.
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
  c10 (x) = r[10];
  a10 (x) = r[11];
  a9 (x) = r[12];
  a8 (x) = r[13];
  a7 (x) = r[14];
  a6 (x) = r[15];
  a5 (x) = r[16];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = 1.L + a5 (x) * (5.L - c5 (x)) + a6 (x) * (6.L - c6 (x))
    + a7 (x) * (7.L - c7 (x)) + a8 (x) * (8.L - c8 (x))
    + a9 (x) * (9.L - c9 (x)) + a10 (x) * (10.L - c10 (x))
    - c0 (x) * (1.L - a5 (x) - a6 (x) - a7 (x) - a8 (x) - a9 (x) - a10 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 1.L - a5 (x) * (25.L - 10.L * c5 (x)) - a6 (x) * (36.L - 12.L * c6 (x))
    - a7 (x) * (49.L - 14.L * c7 (x)) - a8 (x) * (64.L - 16.L * c8 (x))
    - a9 (x) * (81.L - 18.L * c9 (x)) - a10 (x) * (100.L - 20.L * c10 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = 1.L + a5 (x) * (125.L - 75.L * c5 (x))
    + a6 (x) * (216.L - 108.L * c6 (x)) + a7 (x) * (343.L - 147.L * c7 (x))
    + a8 (x) * (512.L - 192.L * c8 (x))
    + a9 (x) * (729.L - 243.L * c9 (x)) + a10 (x) * (1000.L - 300.L * c9 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 1.L - a5 (x) * (625.L - 500.L * c5 (x))
    - a6 (x) * (1296.L - 864.L * c6 (x))
    - a7 (x) * (2401.L - 1372.L * c7 (x))
    - a8 (x) * (4096.L - 2048.L * c8 (x))
    - a9 (x) * (6561.L - 2916.L * c9 (x))
    - a10 (x) * (10000.L - 4000.L * c10 (x));
  solve_4 (A, B, C, D, E);
  a4 (x) = E[3];
  a3 (x) = E[2];
  a2 (x) = E[1];
  a1 (x) = E[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x) - a10 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 11 steps 5th order multi-steps method.
 */
static int
steps_11_5 (Optimize * optimize)        ///< Optimize struct.
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
  c7 (x) = r[7];
  c8 (x) = r[8];
  c9 (x) = r[9];
  c10 (x) = r[10];
  a10 (x) = r[11];
  a9 (x) = r[12];
  a8 (x) = r[13];
  a7 (x) = r[14];
  a6 (x) = r[15];
  A[0] = -1.L + c1 (x) - c0 (x);
  B[0] = -2.L + c2 (x) - c0 (x);
  C[0] = -3.L + c3 (x) - c0 (x);
  D[0] = -4.L + c4 (x) - c0 (x);
  E[0] = -5.L + c5 (x) - c0 (x);
  F[0] = 1.L + a6 (x) * (6.L - c6 (x)) + a7 (x) * (7.L - c7 (x))
    + a8 (x) * (8.L - c8 (x)) + a9 (x) * (9.L - c9 (x))
    + a10 (x) * (10.L - c10 (x))
    - c0 (x) * (1.L - a6 (x) - a7 (x) - a8 (x) - a9 (x) - a10 (x));
  A[1] = 1.L - 2.L * c1 (x);
  B[1] = 4.L - 4.L * c2 (x);
  C[1] = 9.L - 6.L * c3 (x);
  D[1] = 16.L - 8.L * c4 (x);
  E[1] = 25.L - 10.L * c5 (x);
  F[1] = 1.L - a6 (x) * (36.L - 12.L * c6 (x)) - a7 (x) * (49.L - 14.L * c7 (x))
    - a8 (x) * (64.L - 16.L * c8 (x)) - a9 (x) * (81.L - 18.L * c9 (x))
    - a10 (x) * (100.L - 20.L * c10 (x));
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -8.L + 12.L * c2 (x);
  C[2] = -27.L + 27.L * c3 (x);
  D[2] = -64.L + 48.L * c4 (x);
  E[2] = -125.L + 75.L * c5 (x);
  F[2] = 1.L + a6 (x) * (216.L - 108.L * c6 (x))
    + a7 (x) * (343.L - 147.L * c7 (x)) + a8 (x) * (512.L - 192.L * c8 (x))
    + a9 (x) * (729.L - 243.L * c9 (x)) + a10 (x) * (1000.L - 300.L * c10 (x));
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 16.L - 32.L * c2 (x);
  C[3] = 81.L - 108.L * c3 (x);
  D[3] = 256.L - 256.L * c4 (x);
  E[3] = 625.L - 500.L * c5 (x);
  F[3] = 1.L - a6 (x) * (1296.L - 864.L * c6 (x))
    - a7 (x) * (2401.L - 1372.L * c7 (x))
    - a8 (x) * (4096.L - 2048.L * c8 (x))
    - a9 (x) * (6561.L - 2916.L * c9 (x))
    - a10 (x) * (10000.L - 4000.L * c10 (x));
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -32.L + 80.L * c2 (x);
  C[4] = -243.L + 405.L * c3 (x);
  D[4] = -1024.L + 1280.L * c4 (x);
  E[4] = -3125.L + 3125.L * c5 (x);
  F[4] = 1.L + a6 (x) * (7776.L - 6480.L * c6 (x))
    + a7 (x) * (16807.L - 12005.L * c7 (x))
    + a8 (x) * (32768.L - 20480.L * c8 (x))
    + a9 (x) * (59049.L - 32805.L * c9 (x))
    + a10 (x) * (100000.L - 50000.L * c10 (x));
  solve_5 (A, B, C, D, E, F);
  a5 (x) = F[4];
  a4 (x) = F[3];
  a3 (x) = F[2];
  a2 (x) = F[1];
  a1 (x) = F[0];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x) - a10 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)))
    return 0;
	return 1;
}

/**
 * Function to get the coefficients on a 11 steps 6th order multi-steps method.
 */
static int
steps_11_6 (Optimize * optimize)        ///< Optimize struct.
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
  G[1] = 1.L - a9 (x) * (81.L - 18.L * c9 (x)) - a10 (x) * 100.L;
  A[2] = -1.L + 3.L * c1 (x);
  B[2] = -27.L + 27.L * c3 (x);
  C[2] = -64.L + 48.L * c4 (x);
  D[2] = -125.L + 75.L * c5 (x);
  E[2] = -216.L + 108.L * c6 (x);
  F[2] = -512.L + 192.L * c8 (x);
  G[2] = 1.L + a9 (x) * (729.L - 243.L * c9 (x)) + a10 (x) * 1000.L;
  A[3] = 1.L - 4.L * c1 (x);
  B[3] = 81.L - 108.L * c3 (x);
  C[3] = 256.L - 256.L * c4 (x);
  D[3] = 625.L - 500.L * c5 (x);
  E[3] = 1296.L - 864.L * c6 (x);
  F[3] = 4096.L - 2048.L * c8 (x);
  G[3] = 1.L - a9 (x) * (6561.L - 2916.L * c9 (x)) - a10 (x) * 10000.L;
  A[4] = -1.L + 5.L * c1 (x);
  B[4] = -243.L + 405.L * c3 (x);
  C[4] = -1024.L + 1280.L * c4 (x);
  D[4] = -3125.L + 3125.L * c5 (x);
  E[4] = -7776.L + 6480.L * c6 (x);
  F[4] = -32768.L + 20480.L * c8 (x);
  G[4] = 1.L + a9 (x) * (59049.L - 32805.L * c9 (x)) + a10 (x) * 100000.L;
  A[5] = 1.L - 6.L * c1 (x);
  B[5] = 729.L - 1458.L * c3 (x);
  C[5] = 4096.L - 6144.L * c4 (x);
  D[5] = 15625.L - 18750.L * c5 (x);
  E[5] = 46656.L - 46656.L * c6 (x);
  F[5] = 262144.L - 196608.L * c8 (x);
  G[5] = 1.L - a9 (x) * (531441.L - 354294.L * c9 (x)) - a10 (x) * 1000000.L;
  solve_6 (A, B, C, D, E, F, G);
  a8 (x) = G[5];
  a6 (x) = G[4];
  a5 (x) = G[3];
  a4 (x) = G[2];
  a3 (x) = G[1];
  a1 (x) = G[0];
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
    + a9 (x) * (729.L - 243.L * c9 (x)) + a10 (x) * (1000.L - 300.L * c10 (x));
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
#endif
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x) - a6 (x) - a7 (x)
    - a8 (x) - a9 (x) - a10 (x);
  if (isnan (a0 (x)) || isnan (a1 (x)) || isnan (a2 (x)) || isnan (a3 (x))
			|| isnan (a4 (x)) || isnan (a5 (x)) || isnan (a6 (x)))
    return 0;
	return 1;
}

/**
 * Function to print a maxima format file to check the accuracy order of a
 * multi-steps method.
 */
static void
steps_print_maxima (FILE * file,        ///< file.
                    unsigned int nsteps,        ///< steps number.
                    unsigned int order) ///< accuracy order.
{
  int m;
  unsigned int i, j, k, l;

  // 0th order
  fprintf (file, "a0");
  for (i = 1; i < nsteps; ++i)
    fprintf (file, "+a%u", i);
  fprintf (file, "-1;\n");

  // 1st order
  fprintf (file, "b0");
  for (i = 1; i < nsteps; ++i)
    fprintf (file, "+b%u", i);
  for (i = 1; i < nsteps; ++i)
    fprintf (file, "-%u*a%u", i, i);
  fprintf (file, "-1;\n");

  // high order
  for (j = 2, m = 1; j <= order; ++j, m = -m)
    {
      for (i = 1; i < nsteps; ++i)
        {
          for (k = 1, l = i; k < j; ++k)
            l *= i;
          fprintf (file, "-%u*a%u", l, i);
        }
      for (i = 1; i < nsteps; ++i)
        {
          for (k = 2, l = i * j; k < j; ++k)
            l *= i;
          fprintf (file, "+%u*b%u", l, i);
        }
      fprintf (file, "+%d;\n", m);
    }
}

/**
 * Function to print on a file the coefficients of the multi-steps 1st step.
 */
static void
steps_print_1 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  x = optimize->coefficient;
  fprintf (file, "a0:%.19Le;\n", a0 (x));
  fprintf (file, "b0:%.19Le;\n", b0 (x));
  fprintf (file, "c0:%.19Le;\n", c0 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 2nd step.
 */
static void
steps_print_2 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_1 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a1:%.19Le;\n", a1 (x));
  fprintf (file, "b1:%.19Le;\n", b1 (x));
  fprintf (file, "c1:%.19Le;\n", c1 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 3rd step.
 */
static void
steps_print_3 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_2 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a2:%.19Le;\n", a2 (x));
  fprintf (file, "b2:%.19Le;\n", b2 (x));
  fprintf (file, "c2:%.19Le;\n", c2 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 4th step.
 */
static void
steps_print_4 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_3 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a3:%.19Le;\n", a3 (x));
  fprintf (file, "b3:%.19Le;\n", b3 (x));
  fprintf (file, "c3:%.19Le;\n", c3 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 5th step.
 */
static void
steps_print_5 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_4 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a4:%.19Le;\n", a4 (x));
  fprintf (file, "b4:%.19Le;\n", b4 (x));
  fprintf (file, "c4:%.19Le;\n", c4 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 6th step.
 */
static void
steps_print_6 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_5 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a5:%.19Le;\n", a5 (x));
  fprintf (file, "b5:%.19Le;\n", b5 (x));
  fprintf (file, "c5:%.19Le;\n", c5 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 7th step.
 */
static void
steps_print_7 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_6 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a6:%.19Le;\n", a6 (x));
  fprintf (file, "b6:%.19Le;\n", b6 (x));
  fprintf (file, "c6:%.19Le;\n", c6 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 8th step.
 */
static void
steps_print_8 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_7 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a7:%.19Le;\n", a7 (x));
  fprintf (file, "b7:%.19Le;\n", b7 (x));
  fprintf (file, "c7:%.19Le;\n", c7 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 9th step.
 */
static void
steps_print_9 (Optimize * optimize,     ///< Optimize struct.
               FILE * file)     ///< file
{
  long double *x;
  steps_print_8 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a8:%.19Le;\n", a8 (x));
  fprintf (file, "b8:%.19Le;\n", b8 (x));
  fprintf (file, "c8:%.19Le;\n", c8 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 10th step.
 */
static void
steps_print_10 (Optimize * optimize,    ///< Optimize struct.
                FILE * file)    ///< file
{
  long double *x;
  steps_print_9 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a9:%.19Le;\n", a9 (x));
  fprintf (file, "b9:%.19Le;\n", b9 (x));
  fprintf (file, "c9:%.19Le;\n", c9 (x));
}

/**
 * Function to print on a file the coefficients of the multi-steps 11th step.
 */
static void
steps_print_11 (Optimize * optimize,    ///< Optimize struct.
                FILE * file)    ///< file
{
  long double *x;
  steps_print_10 (optimize, file);
  x = optimize->coefficient;
  fprintf (file, "a10:%.19Le;\n", a10 (x));
  fprintf (file, "b10:%.19Le;\n", b10 (x));
  fprintf (file, "c10:%.19Le;\n", c10 (x));
}

/**
 * Function to get the objective function of a 3 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_3 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (k < 0.L)
    return 20.L - k;
  return fmaxl (c0 (x), fmaxl (c1 (x), c2 (x)));
}

/**
 * Function to get the objective function of a 4 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_4 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (k < 0.L)
    return 20.L - k;
  return fmaxl (c0 (x), fmaxl (c1 (x), fmaxl (c2 (x), c3 (x))));
}

/**
 * Function to get the objective function of a 5 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_5 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (k < 0.L)
    return 20.L - k;
  return fmaxl (c0 (x), fmaxl (c1 (x), fmaxl (c2 (x), fmaxl (c3 (x), c4 (x)))));
}

/**
 * Function to get the objective function of a 6 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_6 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)) || isnan (a5 (x)) || isnan (c5 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (a5 (x) < 0.L)
    k += a5 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (b5 (x) < 0.L)
    k += b5 (x);
  if (k < 0.L)
    return 20.L - k;
  return
    fmaxl (c0 (x), fmaxl (c1 (x),
                          fmaxl (c2 (x),
                                 fmaxl (c3 (x), fmaxl (c4 (x), c5 (x))))));
}

/**
 * Function to get the objective function of a 7 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_7 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)) || isnan (a5 (x)) || isnan (c5 (x))
      || isnan (a6 (x)) || isnan (c6 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (a5 (x) < 0.L)
    k += a5 (x);
  if (a6 (x) < 0.L)
    k += a6 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (b5 (x) < 0.L)
    k += b5 (x);
  if (b6 (x) < 0.L)
    k += b6 (x);
  if (k < 0.L)
    return 20.L - k;
  return
    fmaxl (c0 (x),
           fmaxl (c1 (x),
                  fmaxl (c2 (x),
                         fmaxl (c3 (x),
                                fmaxl (c4 (x), fmaxl (c5 (x), c6 (x)))))));
}

/**
 * Function to get the objective function of a 8 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_8 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)) || isnan (a5 (x)) || isnan (c5 (x))
      || isnan (a6 (x)) || isnan (c6 (x)) || isnan (a7 (x)) || isnan (c7 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (a5 (x) < 0.L)
    k += a5 (x);
  if (a6 (x) < 0.L)
    k += a6 (x);
  if (a7 (x) < 0.L)
    k += a7 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (b5 (x) < 0.L)
    k += b5 (x);
  if (b6 (x) < 0.L)
    k += b6 (x);
  if (b7 (x) < 0.L)
    k += b7 (x);
  if (k < 0.L)
    return 20.L - k;
  return
    fmaxl (c0 (x),
           fmaxl (c1 (x),
                  fmaxl (c2 (x),
                         fmaxl (c3 (x),
                                fmaxl (c4 (x),
                                       fmaxl (c5 (x),
                                              fmaxl (c6 (x), c7 (x))))))));
}

/**
 * Function to get the objective function of a 9 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_9 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)) || isnan (a5 (x)) || isnan (c5 (x))
      || isnan (a6 (x)) || isnan (c6 (x)) || isnan (a7 (x)) || isnan (c7 (x))
      || isnan (a8 (x)) || isnan (c8 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (a5 (x) < 0.L)
    k += a5 (x);
  if (a6 (x) < 0.L)
    k += a6 (x);
  if (a7 (x) < 0.L)
    k += a7 (x);
  if (a8 (x) < 0.L)
    k += a8 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (b5 (x) < 0.L)
    k += b5 (x);
  if (b6 (x) < 0.L)
    k += b6 (x);
  if (b7 (x) < 0.L)
    k += b7 (x);
  if (b8 (x) < 0.L)
    k += b8 (x);
  if (k < 0.L)
    return 20.L - k;
  return
    fmaxl (c0 (x),
           fmaxl (c1 (x),
                  fmaxl (c2 (x),
                         fmaxl (c3 (x),
                                fmaxl (c4 (x),
                                       fmaxl (c5 (x),
                                              fmaxl (c6 (x),
                                                     fmaxl (c7 (x),
                                                            c8 (x)))))))));
}

/**
 * Function to get the objective function of a 10 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_10 (Optimize * optimize)        ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)) || isnan (a5 (x)) || isnan (c5 (x))
      || isnan (a6 (x)) || isnan (c6 (x)) || isnan (a7 (x)) || isnan (c7 (x))
      || isnan (a8 (x)) || isnan (c8 (x)) || isnan (a9 (x)) || isnan (c9 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (a5 (x) < 0.L)
    k += a5 (x);
  if (a6 (x) < 0.L)
    k += a6 (x);
  if (a7 (x) < 0.L)
    k += a7 (x);
  if (a8 (x) < 0.L)
    k += a8 (x);
  if (a9 (x) < 0.L)
    k += a9 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (b5 (x) < 0.L)
    k += b5 (x);
  if (b6 (x) < 0.L)
    k += b6 (x);
  if (b7 (x) < 0.L)
    k += b7 (x);
  if (b8 (x) < 0.L)
    k += b8 (x);
  if (b9 (x) < 0.L)
    k += b9 (x);
  if (k < 0.L)
    return 20.L - k;
  k = fmaxl (c8 (x), c9 (x));
  return
    fmaxl (c0 (x),
           fmaxl (c1 (x),
                  fmaxl (c2 (x),
                         fmaxl (c3 (x),
                                fmaxl (c4 (x),
                                       fmaxl (c5 (x),
                                              fmaxl (c6 (x),
                                                     fmaxl (c7 (x), k))))))));
}

/**
 * Function to get the objective function of a 11 steps mult-steps method.
 * 
 * \return objective function value.
 */
static long double
steps_objective_11 (Optimize * optimize)        ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x))
      || isnan (a2 (x)) || isnan (c2 (x)) || isnan (a3 (x)) || isnan (c3 (x))
      || isnan (a4 (x)) || isnan (c4 (x)) || isnan (a5 (x)) || isnan (c5 (x))
      || isnan (a6 (x)) || isnan (c6 (x)) || isnan (a7 (x)) || isnan (c7 (x))
      || isnan (a8 (x)) || isnan (c8 (x)) || isnan (a9 (x)) || isnan (c9 (x))
      || isnan (a10 (x)) || isnan (c10 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (a5 (x) < 0.L)
    k += a5 (x);
  if (a6 (x) < 0.L)
    k += a6 (x);
  if (a7 (x) < 0.L)
    k += a7 (x);
  if (a8 (x) < 0.L)
    k += a8 (x);
  if (a9 (x) < 0.L)
    k += a9 (x);
  if (a10 (x) < 0.L)
    k += a10 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (b4 (x) < 0.L)
    k += b4 (x);
  if (b5 (x) < 0.L)
    k += b5 (x);
  if (b6 (x) < 0.L)
    k += b6 (x);
  if (b7 (x) < 0.L)
    k += b7 (x);
  if (b8 (x) < 0.L)
    k += b8 (x);
  if (b9 (x) < 0.L)
    k += b9 (x);
  if (b10 (x) < 0.L)
    k += b10 (x);
  if (k < 0.L)
    return 20.L - k;
  k = fmaxl (c9 (x), c10 (x));
  return
    fmaxl (c0 (x),
           fmaxl (c1 (x),
                  fmaxl (c2 (x),
                         fmaxl (c3 (x),
                                fmaxl (c4 (x),
                                       fmaxl (c5 (x),
                                              fmaxl (c6 (x),
                                                     fmaxl (c7 (x),
                                                            fmaxl (c8 (x),
                                                                   k)))))))));
}

/**
 * Function to select the multi-steps method.
 *
 * \return 1 on success, 0 on error.
 */
static inline int
steps_select (Optimize * optimize,      ///< Optimize struct.
              unsigned int nsteps,      ///< number of steps.
              unsigned int order)       ///< order of accuracy.
{
  const char *message[] = { _("Bad order"), _("Bad steps number") };
  unsigned int code;
#if DEBUG_STEPS
  fprintf (stderr, "steps_run: start\n");
#endif
  optimize->size = 2 * nsteps;
  optimize->nfree = optimize->size - order - 1;
  optimize->minimum0
    = (long double *) g_slice_alloc (optimize->nfree * sizeof (long double));
  optimize->interval0
    = (long double *) g_slice_alloc (optimize->nfree * sizeof (long double));
  optimize->random_type
    = (unsigned int *) g_slice_alloc (optimize->nfree * sizeof (unsigned int));
  optimize->data = NULL;
  optimize->print_maxima = steps_print_maxima;
  switch (nsteps)
    {
    case 3:
      optimize->print = steps_print_3;
      optimize->objective = steps_objective_3;
      switch (order)
        {
        case 2:
          optimize->method = steps_3_2;
          break;
        case 3:
          optimize->method = steps_3_3;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 4:
      optimize->print = steps_print_4;
      optimize->objective = steps_objective_4;
      switch (order)
        {
        case 2:
          optimize->method = steps_4_2;
          break;
        case 3:
          optimize->method = steps_4_3;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 5:
      optimize->print = steps_print_5;
      optimize->objective = steps_objective_5;
      switch (order)
        {
        case 2:
          optimize->method = steps_5_2;
          break;
        case 3:
          optimize->method = steps_5_3;
          break;
        case 4:
          optimize->method = steps_5_4;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 6:
      optimize->print = steps_print_6;
      optimize->objective = steps_objective_6;
      switch (order)
        {
        case 2:
          optimize->method = steps_6_2;
          break;
        case 3:
          optimize->method = steps_6_3;
          break;
        case 4:
          optimize->method = steps_6_4;
          break;
        case 5:
          optimize->method = steps_6_5;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 7:
      optimize->print = steps_print_7;
      optimize->objective = steps_objective_7;
      switch (order)
        {
        case 2:
          optimize->method = steps_7_2;
          break;
        case 3:
          optimize->method = steps_7_3;
          break;
        case 4:
          optimize->method = steps_7_4;
          break;
        case 5:
          optimize->method = steps_7_5;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 8:
      optimize->print = steps_print_8;
      optimize->objective = steps_objective_8;
      switch (order)
        {
        case 2:
          optimize->method = steps_8_2;
          break;
        case 3:
          optimize->method = steps_8_3;
          break;
        case 4:
          optimize->method = steps_8_4;
          break;
        case 5:
          optimize->method = steps_8_5;
          break;
        case 6:
          optimize->method = steps_8_6;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 9:
      optimize->print = steps_print_9;
      optimize->objective = steps_objective_9;
      switch (order)
        {
        case 2:
          optimize->method = steps_9_2;
          break;
        case 3:
          optimize->method = steps_9_3;
          break;
        case 4:
          optimize->method = steps_9_4;
          break;
        case 5:
          optimize->method = steps_9_5;
          break;
        case 6:
          optimize->method = steps_9_6;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 10:
      optimize->print = steps_print_10;
      optimize->objective = steps_objective_10;
      switch (order)
        {
        case 2:
          optimize->method = steps_10_2;
          break;
        case 3:
          optimize->method = steps_10_3;
          break;
        case 4:
          optimize->method = steps_10_4;
          break;
        case 5:
          optimize->method = steps_10_5;
          break;
        case 6:
          optimize->method = steps_10_6;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
    case 11:
      optimize->print = steps_print_11;
      optimize->objective = steps_objective_11;
      switch (order)
        {
        case 2:
          optimize->method = steps_11_2;
          break;
        case 3:
          optimize->method = steps_11_3;
          break;
        case 4:
          optimize->method = steps_11_4;
          break;
        case 5:
          optimize->method = steps_11_5;
          break;
        case 6:
#if OPTIMIZE_STEPS_11_6 == 1
          optimize->nfree = 10;
#elif OPTIMIZE_STEPS_11_6 == 2
          optimize->nfree = 4;
#endif
          optimize->method = steps_11_6;
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
#if DEBUG_STEPS
  fprintf (stderr, "steps_select: end\n");
#endif
  return 1;

exit_on_error:
  error_message = g_strdup (message[code]);
#if DEBUG_STEPS
  fprintf (stderr, "steps_select: end\n");
#endif
  return 0;
}

/**
 * Function to read the multi-steps method data on a XML node.
 *
 * \return 1 on success, 0 on error.
 */
int
steps_run (xmlNode * node,      ///< XML node.
           gsl_rng ** rng)      ///< array of gsl_rng structs.
{
  Optimize s[nthreads];
  char filename[32];
  gchar *buffer;
  FILE *file;
  long double *value_optimal;
  long double optimal;
  int code;
  unsigned int i, j, nsteps, order, nfree;

#if DEBUG_STEPS
  fprintf (stderr, "steps_run: start\n");
#endif

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
  if (!steps_select (s, nsteps, order))
    goto exit_on_error;
  if (!optimize_read (s, node))
    goto exit_on_error;
  nfree = s->nfree;
  value_optimal = (long double *) g_slice_alloc (nfree * sizeof (long double));
  optimize_create (s, &optimal, value_optimal);
  node = node->children;
  for (i = 0; i < nfree; ++i, node = node->next)
    if (!read_variable (node, s->minimum0, s->interval0, s->random_type, i))
      goto exit_on_error;
  for (i = 1; i < nthreads; ++i)
    memcpy (s + i, s, sizeof (Optimize));
  j = rank * nthreads;
  for (i = 0; i < nthreads; ++i)
    optimize_init (s + i, rng[j + i]);

  // Method bucle
  printf ("Optimize bucle\n");
  optimize_bucle (s);

  // Print the optimal coefficients
  printf ("Print the optimal coefficients\n");
  memcpy (s->random_data, s->value_optimal, nfree * sizeof (long double));
  code = s->method (s);
  snprintf (filename, 32, "steps-%u-%u.mc", nsteps, order);
  file = fopen (filename, "w");
  s->print (s, file);
  s->print_maxima (file, nsteps, order);
  fclose (file);

  // Free memory
  g_slice_free1 (nfree * sizeof (unsigned int), s->random_type);
  g_slice_free1 (nfree * sizeof (long double), s->interval0);
  g_slice_free1 (nfree * sizeof (long double), s->minimum0);
  for (i = 0; i < nthreads; ++i)
    optimize_delete (s + i);
  g_slice_free1 (nfree * sizeof (long double), value_optimal);

#if DEBUG_STEPS
  fprintf (stderr, "steps_run: end\n");
#endif
  return 1;

exit_on_error:
  buffer = error_message;
  error_message = g_strconcat ("Multi-steps:\n", buffer, NULL);
  g_free (buffer);
#if DEBUG_STEPS
  fprintf (stderr, "steps_run: end\n");
#endif
  return 0;
}
