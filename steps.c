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

#define a0(x) x[0]
///< a0 multi-steps coefficient.
#define b0(x) x[1]
///< b0 multi-steps coefficient.
#define a1(x) x[2]
///< a1 multi-steps coefficient.
#define b1(x) x[3]
///< b1 multi-steps coefficient.
#define a2(x) x[4]
///< a2 multi-steps coefficient.
#define b2(x) x[5]
///< b2 multi-steps coefficient.
#define a3(x) x[6]
///< a3 multi-steps coefficient.
#define b3(x) x[7]
///< b3 multi-steps coefficient.
#define a4(x) x[8]
///< a4 multi-steps coefficient.
#define b4(x) x[9]
///< b4 multi-steps coefficient.
#define a5(x) x[10]
///< a5 multi-steps coefficient.
#define b5(x) x[11]
///< b5 multi-steps coefficient.
#define a6(x) x[12]
///< a6 multi-steps coefficient.
#define b6(x) x[13]
///< b6 multi-steps coefficient.
#define a7(x) x[14]
///< a7 multi-steps coefficient.
#define b7(x) x[15]
///< b7 multi-steps coefficient.
#define a8(x) x[16]
///< a8 multi-steps coefficient.
#define b8(x) x[17]
///< b8 multi-steps coefficient.
#define a9(x) x[18]
///< a9 multi-steps coefficient.
#define b9(x) x[19]
///< b9 multi-steps coefficient.
#define a10(x) x[20]
///< a10 multi-steps coefficient.
#define b10(x) x[21]
///< b10 multi-steps coefficient.
#define c(a, b) (b / a)
#define c0(x) (c(a0(x), b0(x)))
///< c0 multi-steps coefficient.
#define c1(x) (c(a1(x), b1(x)))
///< c1 multi-steps coefficient.
#define c2(x) (c(a2(x), b2(x)))
///< c2 multi-steps coefficient.
#define c3(x) (c(a3(x), b3(x)))
///< c3 multi-steps coefficient.
#define c4(x) (c(a4(x), b4(x)))
///< c4 multi-steps coefficient.
#define c5(x) (c(a5(x), b5(x)))
///< c5 multi-steps coefficient.
#define c6(x) (c(a6(x), b6(x)))
///< c6 multi-steps coefficient.
#define c7(x) (c(a7(x), b7(x)))
///< c7 multi-steps coefficient.
#define c8(x) (c(a8(x), b8(x)))
///< c8 multi-steps coefficient.
#define c9(x) (c(a9(x), b9(x)))
///< c9 multi-steps coefficient.
#define c10(x) (c(a10(x), b10(x)))
///< c10 multi-steps coefficient.

/**
 * Function to get the coefficients on a 3 steps 2nd order multi-steps method.
 */
static int
steps_3_2 (Optimize * optimize) ///< Optimize struct.
{
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  b2 (x) = r[2];
  a0 (x) = 1.L - a1 (x) - a2 (x);
  b1 (x) = 0.5L * (a1 (x) + 4.L * (a2 (x) - b2 (x)) - 1.L);
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) - b1 (x) - b2 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 3 steps 3th order multi-steps method.
 */
static int
steps_3_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a0 (x) = 1.L - a1 (x) - a2 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = -1.L + a1 (x) + 4.L * a2 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 1.L + a1 (x) + 8.L * a2 (x);
  solve_2 (A, B, C);
  b2 (x) = C[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = C[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) - b1 (x) - b2 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 4 steps 2nd order multi-steps method.
 */
static int
steps_4_2 (Optimize * optimize) ///< Optimize struct.
{
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  b3 (x) = r[3];
  b2 (x) = r[4];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x);
  b1 (x) =
    0.5L * (a1 (x) + 4.L * (a2 (x) - b2 (x)) + 9.L * a3 (x) - 6.L * b3 (x) -
            1.L);
  b0 (x) =
    1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) - b1 (x) - b2 (x) - b3 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 4 steps 3th order multi-steps method.
 */
static int
steps_4_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  b3 (x) = r[3];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) - 6.L * b3 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * (a3 (x) - b3 (x));
  solve_2 (A, B, C);
  b2 (x) = C[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = C[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) - b1 (x) - b2 (x)
    - b3 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 4 steps 4th order multi-steps method.
 */
static int
steps_4_4 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = 6.L;
  D[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 27.L;
  D[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * a3 (x);
  A[2] = 4.L;
  B[2] = 32.L;
  C[2] = 108.L;
  D[2] = -1.L + a1 (x) + 16.L * a2 (x) + 81.L * a3 (x);
  solve_3 (A, B, C, D);
  b3 (x) = D[2];
  if (isnan (b3 (x)))
    return 0;
  b2 (x) = D[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = D[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) - b1 (x) - b2 (x)
    - b3 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 5 steps 2nd order multi-steps method.
 */
static int
steps_5_2 (Optimize * optimize) ///< Optimize struct.
{
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  b4 (x) = r[4];
  b3 (x) = r[5];
  b2 (x) = r[6];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
  b1 (x) = 0.5L * (a1 (x) + 4.L * (a2 (x) - b2 (x)) + 9.L * a3 (x)
                   + 16.L * a4 (x) - 6.L * b3 (x) - 8.L * b4 (x) - 1.L);
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x) - b1 (x)
    - b2 (x) - b3 (x) - b4 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 5 steps 3th order multi-steps method.
 */
static int
steps_5_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  b4 (x) = r[4];
  b3 (x) = r[5];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x)
    - 6.L * b3 (x) - 8.L * b4 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * (a3 (x) - b3 (x)) + 64.L * a4 (x)
    - 48.L * b4 (x);
  solve_2 (A, B, C);
  b2 (x) = C[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = C[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x) - b1 (x)
    - b2 (x) - b3 (x) - b4 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 5 steps 4th order multi-steps method.
 */
static int
steps_5_4 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  b4 (x) = r[4];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = 6.L;
  D[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x)
    - 8.L * b4 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 27.L;
  D[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * a3 (x) + 64.L * a4 (x)
    - 48.L * b4 (x);
  A[2] = 4.L;
  B[2] = 32.L;
  C[2] = 108.L;
  D[2] = -1.L + a1 (x) + 16.L * a2 (x) + 81.L * a3 (x)
    + 256.L * (a4 (x) - b4 (x));
  solve_3 (A, B, C, D);
  b3 (x) = D[2];
  if (isnan (b3 (x)))
    return 0;
  b2 (x) = D[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = D[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x) - b1 (x)
    - b2 (x) - b3 (x) - b4 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 5 steps 5th order multi-steps method.
 */
static int
steps_5_5 (Optimize * optimize) ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = 6.L;
  D[0] = 8.L;
  E[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 27.L;
  D[1] = 48.L;
  E[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * a3 (x) + 64.L * a4 (x);
  A[2] = 4.L;
  B[2] = 32.L;
  C[2] = 108.L;
  D[2] = 256.L;
  E[2] = -1.L + a1 (x) + 16.L * a2 (x) + 81.L * a3 (x) + 256.L * a4 (x);
  A[3] = 5.L;
  B[3] = 80.L;
  C[3] = 405.L;
  D[3] = 1280.L;
  E[3] = 1.L + a1 (x) + 32.L * a2 (x) + 243.L * a3 (x) + 1024.L * a4 (x);
  solve_4 (A, B, C, D, E);
  b4 (x) = E[3];
  if (isnan (b4 (x)))
    return 0;
  b3 (x) = E[2];
  if (isnan (b3 (x)))
    return 0;
  b2 (x) = E[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = E[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x) - b1 (x)
    - b2 (x) - b3 (x) - b4 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 6 steps 2nd order multi-steps method.
 */
static int
steps_6_2 (Optimize * optimize) ///< Optimize struct.
{
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  a5 (x) = r[4];
  b5 (x) = r[5];
  b4 (x) = r[6];
  b3 (x) = r[7];
  b2 (x) = r[8];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  b1 (x) = 0.5L * (a1 (x) + 4.L * (a2 (x) - b2 (x)) + 9.L * a3 (x)
                   + 16.L * a4 (x) + 25.L * a5 (x) - 6.L * b3 (x) - 8.L * b4 (x)
                   - 10.L * b5 (x) - 1.L);
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x)
    + 5.L * a5 (x) - b1 (x) - b2 (x) - b3 (x) - b4 (x) - b5 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 6 steps 3th order multi-steps method.
 */
static int
steps_6_3 (Optimize * optimize) ///< Optimize struct.
{
  long double A[2], B[2], C[2];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  a5 (x) = r[4];
  b5 (x) = r[5];
  b4 (x) = r[6];
  b3 (x) = r[7];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x)
    + 25.L * a5 (x) - 6.L * b3 (x) - 8.L * b4 (x) - 10.L * b5 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * (a3 (x) - b3 (x)) + 64.L * a4 (x)
    + 125.L * a5 (x) - 48.L * b4 (x) - 75.L * b5 (x);
  solve_2 (A, B, C);
  b2 (x) = C[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = C[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x)
    + 5.L * a5 (x) - b1 (x) - b2 (x) - b3 (x) - b4 (x) - b5 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 6 steps 4th order multi-steps method.
 */
static int
steps_6_4 (Optimize * optimize) ///< Optimize struct.
{
  long double A[3], B[3], C[3], D[3];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  a5 (x) = r[4];
  b5 (x) = r[5];
  b4 (x) = r[6];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = 6.L;
  D[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x)
    + 25.L * a5 (x) - 8.L * b4 (x) - 10.L * b5 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 27.L;
  D[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * a3 (x) + 64.L * a4 (x)
    + 125.L * a5 (x) - 48.L * b4 (x) - 75.L * b5 (x);
  A[2] = 4.L;
  B[2] = 32.L;
  C[2] = 108.L;
  D[2] = -1.L + a1 (x) + 16.L * a2 (x) + 81.L * a3 (x)
    + 256.L * (a4 (x) - b4 (x)) + 625.L * a5 (x) - 500.L * b5 (x);
  solve_3 (A, B, C, D);
  b3 (x) = D[2];
  if (isnan (b3 (x)))
    return 0;
  b2 (x) = D[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = D[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x)
    + 5.L * a5 (x) - b1 (x) - b2 (x) - b3 (x) - b4 (x) - b5 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 6 steps 5th order multi-steps method.
 */
static int
steps_6_5 (Optimize * optimize) ///< Optimize struct.
{
  long double A[4], B[4], C[4], D[4], E[4];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  a5 (x) = r[4];
  b5 (x) = r[5];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = 6.L;
  D[0] = 8.L;
  E[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x)
    + 25.L * a5 (x) - 10.L * b5 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 27.L;
  D[1] = 48.L;
  E[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * a3 (x) + 64.L * a4 (x)
    + 125.L * a5 (x) - 75.L * b5 (x);
  A[2] = 4.L;
  B[2] = 32.L;
  C[2] = 108.L;
  D[2] = 256.L;
  E[2] = -1.L + a1 (x) + 16.L * a2 (x) + 81.L * a3 (x) + 256.L * a4 (x)
    + 625.L * a5 (x) - 500.L * b5 (x);
  A[3] = 5.L;
  B[3] = 80.L;
  C[3] = 405.L;
  D[3] = 1280.L;
  E[3] = 1.L + a1 (x) + 32.L * a2 (x) + 243.L * a3 (x) + 1024.L * a4 (x)
    + 3125.L * (a5 (x) - b5 (x));
  solve_4 (A, B, C, D, E);
  b4 (x) = E[3];
  if (isnan (b4 (x)))
    return 0;
  b3 (x) = E[2];
  if (isnan (b3 (x)))
    return 0;
  b2 (x) = E[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = E[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x)
    + 5.L * a5 (x) - b1 (x) - b2 (x) - b3 (x) - b4 (x) - b5 (x);
  return 1;
}

/**
 * Function to get the coefficients on a 6 steps 6th order multi-steps method.
 */
static int
steps_6_6 (Optimize * optimize) ///< Optimize struct.
{
  long double A[5], B[5], C[5], D[5], E[5], F[5];
  long double *x, *r;
  x = optimize->coefficient;
  r = optimize->random_data;
  a1 (x) = r[0];
  a2 (x) = r[1];
  a3 (x) = r[2];
  a4 (x) = r[3];
  a5 (x) = r[4];
  a0 (x) = 1.L - a1 (x) - a2 (x) - a3 (x) - a4 (x) - a5 (x);
  A[0] = 2.L;
  B[0] = 4.L;
  C[0] = 6.L;
  D[0] = 8.L;
  E[0] = 10.L;
  F[0] = -1.L + a1 (x) + 4.L * a2 (x) + 9.L * a3 (x) + 16.L * a4 (x)
    + 25.L * a5 (x);
  A[1] = 3.L;
  B[1] = 12.L;
  C[1] = 27.L;
  D[1] = 48.L;
  E[1] = 75.L;
  F[1] = 1.L + a1 (x) + 8.L * a2 (x) + 27.L * a3 (x) + 64.L * a4 (x)
    + 125.L * a5 (x);
  A[2] = 4.L;
  B[2] = 32.L;
  C[2] = 108.L;
  D[2] = 256.L;
  E[2] = 500.L;
  F[2] = -1.L + a1 (x) + 16.L * a2 (x) + 81.L * a3 (x) + 256.L * a4 (x)
    + 625.L * a5 (x);
  A[3] = 5.L;
  B[3] = 80.L;
  C[3] = 405.L;
  D[3] = 1280.L;
  E[3] = 3125.L;
  F[3] = 1.L + a1 (x) + 32.L * a2 (x) + 243.L * a3 (x) + 1024.L * a4 (x)
    + 3125.L * a5 (x);
  A[4] = 6.L;
  B[4] = 192.L;
  C[4] = 1458.L;
  D[4] = 6144.L;
  E[4] = 18750.L;
  F[4] = -1.L + a1 (x) + 64.L * a2 (x) + 729.L * a3 (x) + 4096.L * a4 (x)
    + 15625.L * a5 (x);
  solve_5 (A, B, C, D, E, F);
  b5 (x) = F[4];
  if (isnan (b5 (x)))
    return 0;
  b4 (x) = F[3];
  if (isnan (b4 (x)))
    return 0;
  b3 (x) = F[2];
  if (isnan (b3 (x)))
    return 0;
  b2 (x) = F[1];
  if (isnan (b2 (x)))
    return 0;
  b1 (x) = F[0];
  if (isnan (b1 (x)))
    return 0;
  b0 (x) = 1.L + a1 (x) + 2.L * a2 (x) + 3.L * a3 (x) + 4.L * a4 (x)
    + 5.L * a5 (x) - b1 (x) - b2 (x) - b3 (x) - b4 (x) - b5 (x);
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
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (k < 0.L)
    return 20.L - k;
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)))
    return 20.L;
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
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (b2 (x) < 0.L)
    k += b2 (x);
  if (b3 (x) < 0.L)
    k += b3 (x);
  if (k < 0.L)
    return 20.L - k;
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x)))
    return 20.L;
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
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (a2 (x) < 0.L)
    k += a2 (x);
  if (a3 (x) < 0.L)
    k += a3 (x);
  if (a4 (x) < 0.L)
    k += a4 (x);
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)))
    return 20.L;
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
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)) || isnan (c5 (x)))
    return 20.L;
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
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)) || isnan (c5 (x)) || isnan (c6 (x)))
    return 20.L;
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
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)) || isnan (c5 (x)) || isnan (c6 (x)) || isnan (c7 (x)))
    return 20.L;
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
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)) || isnan (c5 (x)) || isnan (c6 (x)) || isnan (c7 (x))
      || isnan (c8 (x)))
    return 20.L;
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
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)) || isnan (c5 (x)) || isnan (c6 (x)) || isnan (c7 (x))
      || isnan (c8 (x)) || isnan (c9 (x)))
    return 20.L;
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
  if (k < 0.L)
    return 30.L - k;
  k = fminl (0.L, b0 (x));
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
  if (isnan (c0 (x)) || isnan (c1 (x)) || isnan (c2 (x)) || isnan (c3 (x))
      || isnan (c4 (x)) || isnan (c5 (x)) || isnan (c6 (x)) || isnan (c7 (x))
      || isnan (c8 (x)) || isnan (c9 (x)) || isnan (c10 (x)))
    return 20.L;
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
        case 4:
          optimize->method = steps_4_4;
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
        case 5:
          optimize->method = steps_5_5;
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
        case 6:
          optimize->method = steps_6_6;
          break;
        default:
          code = 0;
          goto exit_on_error;
        }
      break;
/*
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
*/
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
