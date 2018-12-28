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
#include <math.h>
#include <libxml/parser.h>
#include <glib.h>
#include <libintl.h>
#include <gsl/gsl_rng.h>
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "steps.h"
#include "steps_3_2.h"
#include "steps_3_3.h"
#include "steps_4_2.h"
#include "steps_4_3.h"
#include "steps_5_2.h"
#include "steps_5_3.h"
#include "steps_5_4.h"
#include "steps_6_2.h"
#include "steps_6_3.h"
#include "steps_6_4.h"
#include "steps_6_5.h"
#include "steps_7_2.h"
#include "steps_7_3.h"
#include "steps_7_4.h"
#include "steps_7_5.h"
#include "steps_8_2.h"
#include "steps_8_3.h"
#include "steps_8_4.h"
#include "steps_8_5.h"
#include "steps_8_6.h"
#include "steps_9_2.h"
#include "steps_9_3.h"
#include "steps_9_4.h"
#include "steps_9_5.h"
#include "steps_9_6.h"
#include "steps_10_2.h"
#include "steps_10_3.h"
#include "steps_10_4.h"
#include "steps_10_5.h"
#include "steps_10_6.h"
#include "steps_11_2.h"
#include "steps_11_3.h"
#include "steps_11_4.h"
#include "steps_11_5.h"
#include "steps_11_6.h"

/**
 * Function to print a maxima format file to check the accuracy order of a
 * multi-steps method.
 */
void
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
void
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
void
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
void
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
void
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
void
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
void
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
void
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
void
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
void
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
void
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
void
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
 * Function to get the objective function of a 2 steps mult-steps method.
 * 
 * \return objective function value.
 */
long double
steps_objective_2 (Optimize * optimize) ///< Optimize struct.
{
  register long double *x;
  register long double k;
  x = optimize->coefficient;
  if (isnan (a0 (x)) || isnan (c0 (x)) || isnan (a1 (x)) || isnan (c1 (x)))
    return INFINITY;
  k = fminl (0.L, a0 (x));
  if (a1 (x) < 0.L)
    k += a1 (x);
  if (b0 (x) < 0.L)
    k += b0 (x);
  if (b1 (x) < 0.L)
    k += b1 (x);
  if (k < 0.L)
    return 20.L - k;
  return fmaxl (c0 (x), c1 (x));
}

/**
 * Function to get the objective function of a 3 steps mult-steps method.
 * 
 * \return objective function value.
 */
long double
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
long double
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
long double
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
long double
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
long double
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
long double
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
long double
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
long double
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
                                                                   c9
                                                                   (x))))))))));
}

/**
 * Function to get the objective function of a 11 steps mult-steps method.
 * 
 * \return objective function value.
 */
long double
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
                                                                   fmaxl (c9
                                                                          (x),
                                                                          c10
                                                                          (x)))))))))));
}

/**
 * Function to select the multi-steps method.
 *
 * \return 1 on success, 0 on error.
 */
int
steps_select (Optimize * optimize,      ///< Optimize struct.
              unsigned int nsteps,      ///< number of steps.
              unsigned int order)       ///< order of accuracy.
{
  optimize->size = 2 * nsteps;
  optimize->nfree = optimize->size - order - 1;
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
          optimize->minimum0 = steps_minimum_3_2;
          optimize->interval0 = steps_interval_3_2;
          optimize->random_type = steps_random_3_2;
          optimize->method = steps_3_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_3_3;
          optimize->interval0 = steps_interval_3_3;
          optimize->random_type = steps_random_3_3;
          optimize->method = steps_3_3;
          break;
        default:
          return 0;
        }
      break;
    case 4:
      optimize->print = steps_print_4;
      optimize->objective = steps_objective_4;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_4_2;
          optimize->interval0 = steps_interval_4_2;
          optimize->random_type = steps_random_4_2;
          optimize->method = steps_4_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_4_3;
          optimize->interval0 = steps_interval_4_3;
          optimize->random_type = steps_random_4_3;
          optimize->method = steps_4_3;
          break;
        default:
          return 0;
        }
      break;
    case 5:
      optimize->print = steps_print_5;
      optimize->objective = steps_objective_5;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_5_2;
          optimize->interval0 = steps_interval_5_2;
          optimize->random_type = steps_random_5_2;
          optimize->method = steps_5_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_5_3;
          optimize->interval0 = steps_interval_5_3;
          optimize->random_type = steps_random_5_3;
          optimize->method = steps_5_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_5_4;
          optimize->interval0 = steps_interval_5_4;
          optimize->random_type = steps_random_5_4;
          optimize->method = steps_5_4;
          break;
        default:
          return 0;
        }
      break;
    case 6:
      optimize->print = steps_print_6;
      optimize->objective = steps_objective_6;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_6_2;
          optimize->interval0 = steps_interval_6_2;
          optimize->random_type = steps_random_6_2;
          optimize->method = steps_6_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_6_3;
          optimize->interval0 = steps_interval_6_3;
          optimize->random_type = steps_random_6_3;
          optimize->method = steps_6_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_6_4;
          optimize->interval0 = steps_interval_6_4;
          optimize->random_type = steps_random_6_4;
          optimize->method = steps_6_4;
          break;
        case 5:
          optimize->minimum0 = steps_minimum_6_5;
          optimize->interval0 = steps_interval_6_5;
          optimize->random_type = steps_random_6_5;
          optimize->method = steps_6_5;
          break;
        default:
          return 0;
        }
      break;
    case 7:
      optimize->print = steps_print_7;
      optimize->objective = steps_objective_7;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_7_2;
          optimize->interval0 = steps_interval_7_2;
          optimize->random_type = steps_random_7_2;
          optimize->method = steps_7_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_7_3;
          optimize->interval0 = steps_interval_7_3;
          optimize->random_type = steps_random_7_3;
          optimize->method = steps_7_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_7_4;
          optimize->interval0 = steps_interval_7_4;
          optimize->random_type = steps_random_7_4;
          optimize->method = steps_7_4;
          break;
        case 5:
          optimize->minimum0 = steps_minimum_7_5;
          optimize->interval0 = steps_interval_7_5;
          optimize->random_type = steps_random_7_5;
          optimize->method = steps_7_5;
          break;
        default:
          return 0;
        }
      break;
    case 8:
      optimize->print = steps_print_8;
      optimize->objective = steps_objective_8;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_8_2;
          optimize->interval0 = steps_interval_8_2;
          optimize->random_type = steps_random_8_2;
          optimize->method = steps_8_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_8_3;
          optimize->interval0 = steps_interval_8_3;
          optimize->random_type = steps_random_8_3;
          optimize->method = steps_8_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_8_4;
          optimize->interval0 = steps_interval_8_4;
          optimize->random_type = steps_random_8_4;
          optimize->method = steps_8_4;
          break;
        case 5:
          optimize->minimum0 = steps_minimum_8_5;
          optimize->interval0 = steps_interval_8_5;
          optimize->random_type = steps_random_8_5;
          optimize->method = steps_8_5;
          break;
        case 6:
          optimize->minimum0 = steps_minimum_8_6;
          optimize->interval0 = steps_interval_8_6;
          optimize->random_type = steps_random_8_6;
          optimize->method = steps_8_6;
          break;
        default:
          return 0;
        }
      break;
    case 9:
      optimize->print = steps_print_9;
      optimize->objective = steps_objective_9;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_9_2;
          optimize->interval0 = steps_interval_9_2;
          optimize->random_type = steps_random_9_2;
          optimize->method = steps_9_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_9_3;
          optimize->interval0 = steps_interval_9_3;
          optimize->random_type = steps_random_9_3;
          optimize->method = steps_9_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_9_4;
          optimize->interval0 = steps_interval_9_4;
          optimize->random_type = steps_random_9_4;
          optimize->method = steps_9_4;
          break;
        case 5:
          optimize->minimum0 = steps_minimum_9_5;
          optimize->interval0 = steps_interval_9_5;
          optimize->random_type = steps_random_9_5;
          optimize->method = steps_9_5;
          break;
        case 6:
          optimize->minimum0 = steps_minimum_9_6;
          optimize->interval0 = steps_interval_9_6;
          optimize->random_type = steps_random_9_6;
          optimize->method = steps_9_6;
          break;
        default:
          return 0;
        }
      break;
    case 10:
      optimize->print = steps_print_10;
      optimize->objective = steps_objective_10;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_10_2;
          optimize->interval0 = steps_interval_10_2;
          optimize->random_type = steps_random_10_2;
          optimize->method = steps_10_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_10_3;
          optimize->interval0 = steps_interval_10_3;
          optimize->random_type = steps_random_10_3;
          optimize->method = steps_10_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_10_4;
          optimize->interval0 = steps_interval_10_4;
          optimize->random_type = steps_random_10_4;
          optimize->method = steps_10_4;
          break;
        case 5:
          optimize->minimum0 = steps_minimum_10_5;
          optimize->interval0 = steps_interval_10_5;
          optimize->random_type = steps_random_10_5;
          optimize->method = steps_10_5;
          break;
        case 6:
          optimize->minimum0 = steps_minimum_10_6;
          optimize->interval0 = steps_interval_10_6;
          optimize->random_type = steps_random_10_6;
          optimize->method = steps_10_6;
          break;
        default:
          return 0;
        }
      break;
    case 11:
      optimize->print = steps_print_11;
      optimize->objective = steps_objective_11;
      switch (order)
        {
        case 2:
          optimize->minimum0 = steps_minimum_11_2;
          optimize->interval0 = steps_interval_11_2;
          optimize->random_type = steps_random_11_2;
          optimize->method = steps_11_2;
          break;
        case 3:
          optimize->minimum0 = steps_minimum_11_3;
          optimize->interval0 = steps_interval_11_3;
          optimize->random_type = steps_random_11_3;
          optimize->method = steps_11_3;
          break;
        case 4:
          optimize->minimum0 = steps_minimum_11_4;
          optimize->interval0 = steps_interval_11_4;
          optimize->random_type = steps_random_11_4;
          optimize->method = steps_11_4;
          break;
        case 5:
          optimize->minimum0 = steps_minimum_11_5;
          optimize->interval0 = steps_interval_11_5;
          optimize->random_type = steps_random_11_5;
          optimize->method = steps_11_5;
          break;
        case 6:
#if OPTIMIZE_STEPS_11_6 == 1
          optimize->nfree = 10;
#elif OPTIMIZE_STEPS_11_6 == 2
          optimize->nfree = 4;
#endif
          optimize->minimum0 = steps_minimum_11_6;
          optimize->interval0 = steps_interval_11_6;
          optimize->random_type = steps_random_11_6;
          optimize->method = steps_11_6;
          break;
        default:
          return 0;
        }
      break;
    default:
      return 0;
    }
  return 1;
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
  s->method (s);
  snprintf (filename, 32, "steps-%u-%u.mc", nsteps, order);
  file = fopen (filename, "w");
  s->print (s, file);
  s->print_maxima (file, nsteps, order);
  fclose (file);

  // Free memory
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
