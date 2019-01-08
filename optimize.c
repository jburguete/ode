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
 * \file optimize.c
 * \brief Source file with optimize functions.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <libxml/parser.h>
#include <glib.h>
#include <libintl.h>
#include <gsl/gsl_rng.h>
#if HAVE_MPI
#include <mpi.h>
#endif
#include "config.h"
#include "utils.h"
#include "optimize.h"

#define DEBUG_OPTIMIZE 0        ///< macro to debug.

GMutex mutex[1];                ///< GMutex struct.
FILE *file_variables = NULL;    ///< random variables file.
int rank;                       ///< MPI rank.
int nnodes;                     ///< MPI nodes number.
unsigned nthreads;              ///< threads number.

/**
 * Function to print the random variables on a file.
 */
void
optimize_print_random (Optimize * optimize,     ///< Optimize struct.
                       FILE * file)     ///< file.
{
  unsigned int i, n;
  n = optimize->nfree;
  for (i = 0; i < n; ++i)
    fprintf (file, "o%d:%.19Le;\n", i, optimize->value_optimal[i]);
  for (i = 0; i < n; ++i)
    fprintf (file, "m%d:%.19Le;\n", i, optimize->minimum[i]);
  for (i = 0; i < n; ++i)
    fprintf (file, "i%d:%.19Le;\n", i, optimize->interval[i]);
}

/**
 * Function to perform every optimization step.
 */
void
optimize_step (Optimize * optimize)     ///< Optimize struct.
{
  long double *is, *vo, *vo2, *random;
  long double o, o2, v, f;
  unsigned long long int ii, nrandom;
  unsigned int i, j, k, n, nfree;

#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_step: start\n");
#endif

  // save optimal values
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_step: save optimal values\n");
#endif
  nfree = optimize->nfree;
  o2 = *optimize->optimal;
  vo = (long double *) alloca (nfree * sizeof (long double));
  vo2 = (long double *) alloca (nfree * sizeof (long double));
  memcpy (vo, optimize->value_optimal, nfree * sizeof (long double));

  // optimization algorithm sampling
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_step: optimization algorithm sampling\n");
  fprintf (stderr, "optimize_step: nsimulations=%Lu\n", optimize->nsimulations);
#endif
  random = optimize->random_data;
  ii = optimize->nsimulations * (rank * nthreads + optimize->thread)
    / (nnodes * nthreads);
  nrandom = optimize->nsimulations * (rank * nthreads + optimize->thread + 1)
    / (nnodes * nthreads);
  for (; ii < nrandom; ++ii)
    {

      // random freedom degrees
#if DEBUG_OPTIMIZE
      fprintf (stderr, "optimize_step: random freedom degrees\n");
      fprintf (stderr, "optimize_step: simulation=%Lu\n", ii);
#endif
      optimize_generate_freedom (optimize, ii);

      // method coefficients
#if DEBUG_OPTIMIZE
      fprintf (stderr, "optimize_step: method coefficients\n");
#endif
      if (!optimize->method (optimize))
        o = INFINITY;
      else
        o = optimize->objective (optimize);
      if (o < o2)
        {
          o2 = o;
          memcpy (vo, random, nfree * sizeof (long double));
        }
      if (file_variables)
        {
          g_mutex_lock (mutex);
          print_variables (random, nfree, file_variables);
          fprintf (file_variables, "%.19Le\n", o);
          g_mutex_unlock (mutex);
        }
    }

  // array of intervals to climb around the optimal
#if DEBUG_OPTIMIZE
  fprintf (stderr,
           "optimize_step: array of intervals to climb around the optimal\n");
#endif
  is = (long double *) alloca (nfree * sizeof (long double));
  for (j = 0; j < nfree; ++j)
    is[j] = optimize->interval0[j] * optimize->climbing_factor;

  // hill climbing algorithm bucle
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_step: hill climbing algorithm bucle\n");
#endif
  memcpy (vo2, vo, nfree * sizeof (long double));
  n = optimize->nclimbings;
  for (i = 0; i < n; ++i)
    {
      memcpy (random, vo, nfree * sizeof (long double));
      for (j = k = 0; j < nfree; ++j)
        {
          v = vo[j];
          random[j] = v + is[j];
          if (!optimize->method (optimize))
            o = INFINITY;
          else
            o = optimize->objective (optimize);
          if (o < o2)
            {
              k = 1;
              o2 = o;
              memcpy (vo2, random, nfree * sizeof (long double));
            }
          if (file_variables)
            {
              g_mutex_lock (mutex);
              print_variables (random, nfree, file_variables);
              fprintf (file_variables, "%.19Le\n", o);
              g_mutex_unlock (mutex);
            }
          random[j] = fmaxl (0.L, v - is[j]);
          if (!optimize->method (optimize))
            o = INFINITY;
          else
            o = optimize->objective (optimize);
          if (o < o2)
            {
              k = 1;
              o2 = o;
              memcpy (vo2, random, nfree * sizeof (long double));
            }
          if (file_variables)
            {
              g_mutex_lock (mutex);
              print_variables (random, nfree, file_variables);
              fprintf (file_variables, "%.19Le\n", o);
              g_mutex_unlock (mutex);
            }
          random[j] = v;
        }


      // update optimal values and increase or reduce intervals if converging or
      // not
      if (!k)
        f = 0.5L;
      else
        {
          f = 1.2L;
          memcpy (vo, vo2, nfree * sizeof (long double));
        }
      for (j = 0; j < nfree; ++j)
        is[j] *= f;
    }

  // update optimal values
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_step: update optimal values\n");
#endif
  if (o2 < *optimize->optimal)
    {
      g_mutex_lock (mutex);
      *optimize->optimal = o2;
      memcpy (optimize->value_optimal, vo2, nfree * sizeof (long double));
      g_mutex_unlock (mutex);
    }

#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_step: end\n");
#endif
}

/**
 * Function to init required variables on an Optimize struct data.
 */
void
optimize_init (Optimize * optimize,     ///< Optimize struct.
               gsl_rng * rng,   ///< GSL pseudo-random number generator struct.
               unsigned int thread)     ///< thread number.
{
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_init: start\n");
#endif
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_init: nsimulations=%Lu nfree=%u size=%u\n",
           optimize->nsimulations, optimize->nfree, optimize->size);
#endif
  optimize->random_data
    = (long double *) g_slice_alloc (optimize->nfree * sizeof (long double));
  optimize->coefficient
    = (long double *) g_slice_alloc (optimize->size * sizeof (long double));
  optimize->minimum
    = (long double *) g_slice_alloc (optimize->nfree * sizeof (long double));
  optimize->interval
    = (long double *) g_slice_alloc (optimize->nfree * sizeof (long double));
  memcpy (optimize->minimum, optimize->minimum0,
          optimize->nfree * sizeof (long double));
  memcpy (optimize->interval, optimize->interval0,
          optimize->nfree * sizeof (long double));
  optimize->rng = rng;
  optimize->thread = thread;
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_init: end\n");
#endif
}

/**
 * Function to free the memory allocated by an Optimize struct.
 */
void
optimize_delete (Optimize * optimize)   ///< Optimize struct.
{
  g_slice_free1 (optimize->nfree * sizeof (long double), optimize->interval);
  g_slice_free1 (optimize->nfree * sizeof (long double), optimize->minimum);
  g_slice_free1 (optimize->size * sizeof (long double), optimize->coefficient);
  g_slice_free1 (optimize->nfree * sizeof (long double), optimize->random_data);
}

/**
 * Function to do the optimization bucle.
 */
void
optimize_bucle (Optimize * optimize)    ///< Optimize struct.
{
  GThread *thread[nthreads];
#if HAVE_MPI
  long double *vo;
  MPI_Status status;
#endif
  unsigned int i, j, nfree;

#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_bucle: start\n");
  fprintf (stderr, "optimize_bucle: nfree=%u\n", optimize->nfree);
#endif

  // Allocate local array of optimal values
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_bucle: allocate local array of optimal values\n");
#endif
  nfree = optimize->nfree;
#if HAVE_MPI
  vo = (long double *) alloca ((1 + nfree) * sizeof (long double));
  printf ("Rank=%d NNodes=%d\n", rank, nnodes);
#endif

  // Init some parameters
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_bucle: init some parameters\n");
#endif
  *optimize->optimal = INFINITY;
  for (i = 0; i < nfree; ++i)
    optimize->value_optimal[i]
      = optimize->minimum[i] + 0.5L * optimize->interval[i];

  // Iterate
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_bucle: iterate\n");
#endif
  for (i = 0; i < optimize->niterations; ++i)
    {

      // Optimization step parallelized for every node by GThreads
      if (nthreads > 1)
        {
          for (j = 0; j < nthreads; ++j)
            thread[j]
              = g_thread_new (NULL,
                              (GThreadFunc) (void (*)(void)) optimize_step,
                              (void *) (optimize + j));
          for (j = 0; j < nthreads; ++j)
            g_thread_join (thread[j]);
        }
      else
        optimize_step (optimize);

#if HAVE_MPI
      if (rank > 0)
        {

          // Secondary nodes send the optimal coefficients to the master node
          vo[0] = *optimize->optimal;
          memcpy (vo + 1, optimize->value_optimal,
                  nfree * sizeof (long double));
          MPI_Send (vo, 1 + nfree, MPI_LONG_DOUBLE, 0, 1, MPI_COMM_WORLD);

          // Secondary nodes receive the optimal coefficients
          MPI_Recv (vo, 1 + nfree, MPI_LONG_DOUBLE, 0, 1, MPI_COMM_WORLD,
                    &status);
          *optimize->optimal = vo[0];
          memcpy (optimize->value_optimal, vo + 1,
                  nfree * sizeof (long double));
        }
      else
        {
          printf ("rank=%d optimal=%.19Le\n", rank, *optimize->optimal);

          for (j = 1; j < nnodes; ++j)
            {

              // Master node receives the optimal coefficients obtained by
              // secondary nodes
              MPI_Recv (vo, 1 + nfree, MPI_LONG_DOUBLE, j, 1, MPI_COMM_WORLD,
                        &status);

              // Master node selects the optimal coefficients
              if (vo[0] < *optimize->optimal)
                {
                  *optimize->optimal = vo[0];
                  memcpy (optimize->value_optimal, vo + 1,
                          nfree * sizeof (long double));
                }
            }

          // Master node sends the optimal coefficients to secondary nodes
          vo[0] = *optimize->optimal;
          memcpy (vo + 1, optimize->value_optimal,
                  nfree * sizeof (long double));
          for (j = 1; j < nnodes; ++j)
            MPI_Send (vo, 1 + nfree, MPI_LONG_DOUBLE, j, 1, MPI_COMM_WORLD);
        }

#endif

      // Print the optimal coefficients
#if DEBUG_OPTIMIZE
      optimize_print_random (optimize, stderr);
      fprintf (stderr, "optimal=%.19Le\n", *optimize->optimal);
#endif

      // Updating coefficient intervals to converge
      optimize_converge (optimize);

      // Iterate
#if HAVE_MPI
      printf ("Rank %u\n", rank);
#endif
      printf ("Iteration %u Optimal %.19Le\n", i + 1, *optimize->optimal);
    }

#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_bucle: end\n");
#endif
}

/**
 * Function to create an Optimize struct data.
 */
void
optimize_create (Optimize * optimize,   ///< Optimize struct.
                 long double *optimal,
                 ///< pointer to the optimal objective function value.
                 long double *value_optimal)
                 ///< array of optimal freedom degree values.
{
  unsigned long long int nsimulations;
  unsigned int i, nfree;
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_create: start\n");
#endif
  optimize->optimal = optimal;
  optimize->value_optimal = value_optimal;
  nfree = optimize->nfree;
  optimize->nsimulations = nsimulations = optimize->nvariable;
  for (i = 1; i < nfree; ++i)
    optimize->nsimulations *= nsimulations;
  optimize->nclimbings *= nfree;
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_create nsimulations=%Lu nclimbings=%u nfree=%u\n",
           optimize->nsimulations, optimize->nclimbings, nfree);
  fprintf (stderr, "optimize_create: end\n");
#endif
}

/**
 * Function to read the Optimize struct data on a XML node.
 *
 * \return 1 on success, 0 on error.
 */
int
optimize_read (Optimize * optimize,     ///< Optimize struct.
               xmlNode * node)  ///< XML node.
{
  int code;
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_read: start\n");
#endif
  optimize->nvariable = xml_node_get_uint (node, XML_NSIMULATIONS, &code);
  if (code || !optimize->nvariable)
    {
      error_message = g_strdup (_("Bad simulations number"));
      goto exit_on_error;
    }
  optimize->nclimbings
    = xml_node_get_uint_with_default (node, XML_NCLIMBINGS, 0, &code);
  if (code)
    {
      error_message = g_strdup (_("Bad hill climbings number"));
      goto exit_on_error;
    }
  optimize->niterations = xml_node_get_uint (node, XML_NITERATIONS, &code);
  if (code || !optimize->niterations)
    {
      error_message = g_strdup (_("Bad iterations number"));
      goto exit_on_error;
    }
  optimize->convergence_factor
    = xml_node_get_float (node, XML_CONVERGENCE_FACTOR, &code);
  if (code || optimize->convergence_factor < LDBL_EPSILON)
    {
      error_message = g_strdup (_("Bad convergence factor"));
      goto exit_on_error;
    }
  optimize->climbing_factor
    = xml_node_get_float (node, XML_CLIMBING_FACTOR, &code);
  if (code || optimize->climbing_factor < LDBL_EPSILON)
    {
      error_message = g_strdup (_("Bad climging factor"));
      goto exit_on_error;
    }
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_read: end\n");
#endif
  return 1;

exit_on_error:
#if DEBUG_OPTIMIZE
  fprintf (stderr, "optimize_read: end\n");
#endif
  return 0;
}
