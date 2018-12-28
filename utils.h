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
 * \file utils.h
 * \brief Header file with util functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#ifndef UTILS__H
#define UTILS__H 1

#define _(x) (gettext(x)) ///< macro to translate messages.

extern GMutex mutex[1];
extern gchar *error_message;

void show_error (const char *message);
int xml_node_get_int (xmlNode * node, const xmlChar * prop, int *error_code);
unsigned int
xml_node_get_uint (xmlNode * node, const xmlChar * prop, int *error_code);
unsigned int
xml_node_get_uint_with_default (xmlNode * node, const xmlChar * prop, 
		                            unsigned int default_value, int *error_code);
long double
xml_node_get_float (xmlNode * node, const xmlChar * prop, int *error_code);

/**
 * Function to print the random variables on a file.
 */
static inline void
print_variables (long double *r,   ///< random variables.
                 unsigned int nfree,       ///< number of freedom degrees.
                 FILE * file)      ///< file.
{
  unsigned int i;
  g_mutex_lock (mutex);
  for (i = 0; i < nfree; ++i)
    fprintf (file, "%.19Le ", r[i]);
  g_mutex_unlock (mutex);
}

/**
 * Function to calculate a random number between [0,1] being 0 and 1 the fifty
 * percent.
 *
 * \return result.
 */
static inline long double
random_zero (gsl_rng * rng)     ///< GSL random number generator struct data.
{
  register double r;
  r = gsl_rng_uniform (rng);
  if (r <= 0.25)
    return 0.L;
  if (r >= 0.75)
    return 1.L;
  return 2.L * (r - 0.25L);
}

/**
 * Function to calculate a random number between [0,1] being 1 the fifty
 * percent.
 *
 * \return result.
 */
static inline long double
random_one (gsl_rng * rng)      ///< GSL random number generator struct data.
{
  register double r;
  r = gsl_rng_uniform (rng);
  if (r >= 0.5)
    return 1.L;
  return 2.L * r;
}

/**
 * Function to calculate the square of a number.
 *
 * \return result.
 */
static inline long double
sqr (long double x)             ///< number.
{
  return x * x;
}

/**
 * Function to solve a system of two linear equations.
 */
static inline void
solve_2 (long double *A,
         ///< First matrix column modified by the algorithm.
         long double *B,
         ///< Second matrix column modified by the algorithm.
         long double *C)
 ///< Third matrix column modified by the algorithm to contain the solutions.
{
  B[1] = A[0] * B[1] - A[1] * B[0];
  C[1] = (A[0] * C[1] - A[1] * C[0]) / B[1];
#if EPSILON
  if (fabsl (C[1]) < LDBL_EPSILON)
    C[1] = 0.L;
#endif
  C[0] = (C[0] - B[0] * C[1]) / A[0];
#if EPSILON
  if (fabsl (C[0]) < LDBL_EPSILON)
    C[0] = 0.L;
#endif
}

/**
 * Function to solve a system of three linear equations.
 */
static inline void
solve_3 (long double *A,
         ///< First matrix column modified by the algorithm:
         long double *B,
         ///< Second matrix column modified by the algorithm.
         long double *C,
         ///< Third matrix column modified by the algorithm.
         long double *D)
 ///< Fourth matrix column modified by the algorithm to contain the solutions.
{
  B[1] = A[0] * B[1] - A[1] * B[0];
  C[1] = A[0] * C[1] - A[1] * C[0];
  D[1] = A[0] * D[1] - A[1] * D[0];
  B[2] = A[0] * B[2] - A[2] * B[0];
  C[2] = A[0] * C[2] - A[2] * C[0];
  D[2] = A[0] * D[2] - A[2] * D[0];
  solve_2 (B + 1, C + 1, D + 1);
  D[0] = (D[0] - B[0] * D[1] - C[0] * D[2]) / A[0];
#if EPSILON
  if (fabsl (D[0]) < LDBL_EPSILON)
    D[0] = 0.L;
#endif
}

/**
 * Function to solve a system of fourth linear equations.
 */
static inline void
solve_4 (long double *A,
         ///< First matrix column modified by the algorithm:
         long double *B,
         ///< Second matrix column modified by the algorithm.
         long double *C,
         ///< Third matrix column modified by the algorithm.
         long double *D,
         ///< Fourth matrix column modified by the algorithm.
         long double *E)
 ///< Fifth matrix column modified by the algorithm to contain the solutions.
{
  B[1] = A[0] * B[1] - A[1] * B[0];
  C[1] = A[0] * C[1] - A[1] * C[0];
  D[1] = A[0] * D[1] - A[1] * D[0];
  E[1] = A[0] * E[1] - A[1] * E[0];
  B[2] = A[0] * B[2] - A[2] * B[0];
  C[2] = A[0] * C[2] - A[2] * C[0];
  D[2] = A[0] * D[2] - A[2] * D[0];
  E[2] = A[0] * E[2] - A[2] * E[0];
  B[3] = A[0] * B[3] - A[3] * B[0];
  C[3] = A[0] * C[3] - A[3] * C[0];
  D[3] = A[0] * D[3] - A[3] * D[0];
  E[3] = A[0] * E[3] - A[3] * E[0];
  solve_3 (B + 1, C + 1, D + 1, E + 1);
  E[0] = (E[0] - B[0] * E[1] - C[0] * E[2] - D[0] * E[3]) / A[0];
#if EPSILON
  if (fabsl (E[0]) < LDBL_EPSILON)
    E[0] = 0.L;
#endif
}

/**
 * Function to solve a system of five linear equations.
 */
static inline void
solve_5 (long double *A,
         ///< First matrix column modified by the algorithm:
         long double *B,
         ///< Second matrix column modified by the algorithm.
         long double *C,
         ///< Third matrix column modified by the algorithm.
         long double *D,
         ///< Fourth matrix column modified by the algorithm.
         long double *E,
         ///< Fifh matrix column modified by the algorithm.
         long double *F)
 ///< Sixth matrix column modified by the algorithm to contain the solutions.
{
  B[1] = A[0] * B[1] - A[1] * B[0];
  C[1] = A[0] * C[1] - A[1] * C[0];
  D[1] = A[0] * D[1] - A[1] * D[0];
  E[1] = A[0] * E[1] - A[1] * E[0];
  F[1] = A[0] * F[1] - A[1] * F[0];
  B[2] = A[0] * B[2] - A[2] * B[0];
  C[2] = A[0] * C[2] - A[2] * C[0];
  D[2] = A[0] * D[2] - A[2] * D[0];
  E[2] = A[0] * E[2] - A[2] * E[0];
  F[2] = A[0] * F[2] - A[2] * F[0];
  B[3] = A[0] * B[3] - A[3] * B[0];
  C[3] = A[0] * C[3] - A[3] * C[0];
  D[3] = A[0] * D[3] - A[3] * D[0];
  E[3] = A[0] * E[3] - A[3] * E[0];
  F[3] = A[0] * F[3] - A[3] * F[0];
  B[4] = A[0] * B[4] - A[4] * B[0];
  C[4] = A[0] * C[4] - A[4] * C[0];
  D[4] = A[0] * D[4] - A[4] * D[0];
  E[4] = A[0] * E[4] - A[4] * E[0];
  F[4] = A[0] * F[4] - A[4] * F[0];
  solve_4 (B + 1, C + 1, D + 1, E + 1, F + 1);
  F[0] = (F[0] - B[0] * F[1] - C[0] * F[2] - D[0] * F[3] - E[0] * F[4]) / A[0];
#if EPSILON
  if (fabsl (F[0]) < LDBL_EPSILON)
    F[0] = 0.L;
#endif
}

/**
 * Function to solve a system of six linear equations.
 */
static inline void
solve_6 (long double *A,
         ///< First matrix column modified by the algorithm:
         long double *B,
         ///< Second matrix column modified by the algorithm.
         long double *C,
         ///< Third matrix column modified by the algorithm.
         long double *D,
         ///< Fourth matrix column modified by the algorithm.
         long double *E,
         ///< Fifh matrix column modified by the algorithm.
         long double *F,
         ///< Sixth matrix column modified by the algorithm.
         long double *G)
 ///< Seventh matrix column modified by the algorithm to contain the solutions.
{
  B[1] = A[0] * B[1] - A[1] * B[0];
  C[1] = A[0] * C[1] - A[1] * C[0];
  D[1] = A[0] * D[1] - A[1] * D[0];
  E[1] = A[0] * E[1] - A[1] * E[0];
  F[1] = A[0] * F[1] - A[1] * F[0];
  G[1] = A[0] * G[1] - A[1] * G[0];
  B[2] = A[0] * B[2] - A[2] * B[0];
  C[2] = A[0] * C[2] - A[2] * C[0];
  D[2] = A[0] * D[2] - A[2] * D[0];
  E[2] = A[0] * E[2] - A[2] * E[0];
  F[2] = A[0] * F[2] - A[2] * F[0];
  G[2] = A[0] * G[2] - A[2] * G[0];
  B[3] = A[0] * B[3] - A[3] * B[0];
  C[3] = A[0] * C[3] - A[3] * C[0];
  D[3] = A[0] * D[3] - A[3] * D[0];
  E[3] = A[0] * E[3] - A[3] * E[0];
  F[3] = A[0] * F[3] - A[3] * F[0];
  G[3] = A[0] * G[3] - A[3] * G[0];
  B[4] = A[0] * B[4] - A[4] * B[0];
  C[4] = A[0] * C[4] - A[4] * C[0];
  D[4] = A[0] * D[4] - A[4] * D[0];
  E[4] = A[0] * E[4] - A[4] * E[0];
  F[4] = A[0] * F[4] - A[4] * F[0];
  G[4] = A[0] * G[4] - A[4] * G[0];
  B[5] = A[0] * B[5] - A[5] * B[0];
  C[5] = A[0] * C[5] - A[5] * C[0];
  D[5] = A[0] * D[5] - A[5] * D[0];
  E[5] = A[0] * E[5] - A[5] * E[0];
  F[5] = A[0] * F[5] - A[5] * F[0];
  G[5] = A[0] * G[5] - A[5] * G[0];
  solve_5 (B + 1, C + 1, D + 1, E + 1, F + 1, G + 1);
  G[0] = (G[0] - B[0] * G[1] - C[0] * G[2] - D[0] * G[3] - E[0] * G[4]
          - F[0] * G[5]) / A[0];
#if EPSILON
  if (fabsl (G[0]) < LDBL_EPSILON)
    G[0] = 0.L;
#endif
}

#endif
