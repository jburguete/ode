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
#ifndef OPTIMIZE__H
#define OPTIMIZE__H 1

/**
 * \struct Optimize
 * \brief struct defining the data to perform a method optimization.
 */
typedef struct _Optimize Optimize;
struct _Optimize
{
  void (*print) (Optimize * optimize, FILE * file);
  ///< pointer to the function to print the results data.
  void (*print_maxima) (FILE * file, unsigned int nsteps, unsigned int order);
  ///< pointer to the function to print the maxima format data.
  void (*method) (Optimize * optimize);
  ///< pointer to the function to calculate the method variables.
  long double (*objective) (Optimize * optimize);
  ///< pointer to the function to calculate the objective function.
  gsl_rng *rng;                 ///< GSL pseudo-random numbers generator struct.
  ///< pointer to the array of GSL pseudo-random numbers generator structs.
  long double *coefficient;
  ///< array of method coefficientes.
  long double *random_data;
  ///< array of freedom degree values.
  long double *value_optimal;
  ///< pointer to the array of optimal values of the freedom degrees.
  long double *minimum;
  ///< pointer to the array of minimum values of the freedom degrees.
  long double *interval;
  ///< pointer to the array of intervals of the freedom degrees.
  const long double *minimum0;
  ///< pointer to the initial array of minimum values of the freedom degrees.
  const long double *interval0;
  ///< pointer to the initial array of intervals of the freedom degrees.
  long double *optimal;
  ///< pointer to the optimal objective function value.
  const unsigned int *random_type;
  ///< pointer the the array of random generation types for the freedom degrees.
  void *data;
  ///< pointer to additional method data.
  long double convergence_factor;       ///< convergence factor.
  long double search_factor;
  ///< factor to the coordinates search optimization algorithm.
  unsigned long long int nsimulations;
  ///< number of total simulations on Monte-Carlo optimization algorithm.
  unsigned long long int nrandom;
  ///< number of simulations per thread on Monte-Carlo optimization algorithm.
  unsigned int nsearch;
  ///< number of steps on coordinates search optimization algorithm.
  unsigned int niterations;     ///< iterations number.
  unsigned int nfree;           ///< number of freedom degrees.
  unsigned int size;            ///< total variables number.
  unsigned int type;            ///< method type.
};

typedef void (*OptimizeMethod) (Optimize * optimize);
typedef long double (*OptimizeObjective) (Optimize * optimize);
typedef void (*OptimizePrint) (Optimize * optimize, FILE * file);

#if PRINT_RANDOM
extern FILE *file_random;
extern FILE *file_random2;
#endif
extern int rank;
extern int nnodes;
extern unsigned nthreads;

void optimize_print_random (Optimize * optimize, FILE * file);
void optimize_step (Optimize * optimize);
void optimize_init (Optimize * optimize, gsl_rng * rng);
void optimize_delete (Optimize * optimize);
void optimize_bucle (Optimize * optimize);
void optimize_create (Optimize * optimize,
                      long double *optimal,
                      long double *value_optimal,
                      long double convergence_factor,
                      long double search_factor,
                      unsigned long long int nsimulations,
                      unsigned int nsearch, unsigned int niterations);

/**
 * Function to generate the freedom degree values.
 */
static inline void
optimize_generate_random (Optimize * optimize)  ///< Optimize struct.
{
  gsl_rng *rng;
  long double *data, *minimum, *interval;
  const unsigned int *type;
  unsigned int i, n;
  n = optimize->nfree;
  data = optimize->random_data;
  minimum = optimize->minimum;
  interval = optimize->interval;
  type = optimize->random_type;
  rng = optimize->rng;
  for (i = 0; i < n; ++i)
    switch (type[i])
      {
      case 0:
        data[i] = minimum[i] + interval[i] * gsl_rng_uniform (rng);
        break;
      case 1:
        data[i] = minimum[i] + interval[i] * random_zero (rng);
        break;
      default:
        data[i] = minimum[i] + interval[i] * random_one (rng);
      }
}

/**
 * Function to reduce the search variable intervals to increase convergence.
 */
static inline void
optimize_converge (Optimize * optimize) ///< Optimize struct.
{
  long double d, factor;
  unsigned int i, n;
  n = optimize->nfree;
  factor = optimize->convergence_factor;
  for (i = 0; i < n; ++i)
    {
      d = optimize->interval[i] *= factor;
      optimize->minimum[i] = fmaxl (0.L, optimize->value_optimal[i] - d * 0.5L);
    }
}

#endif
