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
  int (*method) (Optimize * optimize);
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
  long double *minimum0;
  ///< pointer to the initial array of minimum values of the freedom degrees.
  long double *interval0;
  ///< pointer to the initial array of intervals of the freedom degrees.
  long double *optimal;
  ///< pointer to the optimal objective function value.
  unsigned int *random_type;
  ///< pointer the the array of random generation types for the freedom degrees.
  void *data;
  ///< pointer to additional method data.
  long double convergence_factor;       ///< convergence factor.
  long double climbing_factor;
  ///< factor to the coordinates hill climbing optimization algorithm.
  unsigned long long int nsimulations;
  ///< number of total simulations on optimization algorithm.
  unsigned int thread;          ///< thread number.
  unsigned int nvariable;
  ///< number of total simulations per variable on optimization algorithm.
  unsigned int nclimbings;
  ///< number of steps on coordinates hill climbing optimization algorithm.
  unsigned int niterations;     ///< iterations number.
  unsigned int nfree;           ///< number of freedom degrees.
  unsigned int size;            ///< total variables number.
  unsigned int type;            ///< method type.
  unsigned int order;           ///< accuracy order.
  unsigned int nsteps;          ///< steps number.
};

typedef int (*OptimizeMethod) (Optimize * optimize);
typedef long double (*OptimizeObjective) (Optimize * optimize);
typedef void (*OptimizePrint) (Optimize * optimize, FILE * file);

extern FILE *file_variables;
extern int rank;
extern int nnodes;
extern unsigned nthreads;

void optimize_print_random (Optimize * optimize, FILE * file);
void optimize_step (Optimize * optimize);
void optimize_init (Optimize * optimize, gsl_rng * rng, unsigned int thread);
void optimize_delete (Optimize * optimize);
void optimize_bucle (Optimize * optimize);
void optimize_create (Optimize * optimize, long double *optimal,
                      long double *value_optimal);
int optimize_read (Optimize * optimize, xmlNode * node);

/**
 * Function to generate the freedom degree values.
 */
static inline void
optimize_generate_freedom (Optimize * optimize, ///< Optimize struct.
                           unsigned long long int ns)   ///< simulation number.
{
  gsl_rng *rng;
  long double *data, *minimum, *interval;
  unsigned long long int j;
  const unsigned int *type;
  unsigned int i, k, n;
  n = optimize->nfree;
  data = optimize->random_data;
  minimum = optimize->minimum;
  interval = optimize->interval;
  type = optimize->random_type;
  rng = optimize->rng;
  j = ns;
  for (i = 0; i < n; ++i)
    switch (type[i])
      {
      case RANDOM_TYPE_UNIFORM:
        data[i] = minimum[i] + interval[i] * gsl_rng_uniform (rng);
        break;
      case RANDOM_TYPE_BOTTOM:
        data[i] = minimum[i] + interval[i] * random_zero (rng);
        break;
      case RANDOM_TYPE_EXTREME:
        data[i] = minimum[i] + interval[i] * random_extreme (rng);
        break;
      case RANDOM_TYPE_TOP:
        data[i] = minimum[i] + interval[i] * random_one (rng);
        break;
      case RANDOM_TYPE_REGULAR:
        k = j % optimize->nvariable;
        j /= optimize->nvariable;
        data[i] = minimum[i] + interval[i] * k / (optimize->nvariable - 1);
        break;
      default:
        k = j % optimize->nvariable;
        j /= optimize->nvariable;
        data[i] = minimum[i]
          + interval[i] * (k + gsl_rng_uniform (rng)) / optimize->nvariable;
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
