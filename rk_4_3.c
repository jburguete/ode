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
 * \file rk_4_3.c
 * \brief Source file to optimize Runge-Kutta 4 steps 3rd order methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "rk.h"
#include "rk_4_2.h"
#include "rk_4_3.h"

#define DEBUG_RK_4_3 0 ///< macro to debug.

///> array of minimum freedom degree values for the t-b coefficients of the 4
///> steps 3rd order Runge-Kutta method.
const long double minimum_tb_4_3[6] = { 0.L, 0.L, 0.L, 0.L, 0.L, 0.L};

///> array of minimum freedom degree intervals for the t-b coefficients of the 4
///> steps 3rd order Runge-Kutta method.
const long double interval_tb_4_3[6] = { 1.L, 1.L, 1.L, 1.L, 1.L, 1.L};

///> array of freedom degree random function types for the t-b coefficients of
///> the 4 steps 3rd order Runge-Kutta method.
const unsigned int random_tb_4_3[6] = { 2, 2, 2, 2, 2, 2};

/**
 * Function to print a maxima format file to check the accuracy order of a 4
 * steps 3rd order Runge-Kutta method.
 */
void
rk_print_maxima_4_3 (FILE * file,       ///< file.
		                 unsigned int nsteps, ///< steps number.
										 unsigned int order) ///< accuracy order.
{
  rk_print_maxima_4_2 (file, nsteps, order);
  fprintf (file, "b42*b21*t1+b43*(b31*t1+b32*t2)-1/6;\n");
  fprintf (file, "b41*t1^2+b42*t2^2+b43*t3^2-1/3;\n");
}

/**
 * Function to obtain the coefficients of a 4 steps 3rd order Runge-Kutta 
 * method.
 */
void
rk_tb_4_3 (Optimize * optimize)      ///< Optimize struct.
{
	long double *tb, *r;
#if DEBUG_RK_4_3
	fprintf (stderr, "rk_tb_4_3: start\n");
#endif
	tb = optimize->coefficient;
	r = optimize->random_data;
	t1 (tb) = r[0];
	t2 (tb) = r[1];
	b21 (tb) = r[2];
	t3 (tb) = r[3];
	b32 (tb) = r[4];
	b43 (tb) = r[5];
  t4 (tb) = 1.L;
	b42 (tb) = ((1.L / 3.L - b43 (tb) * sqr (t3 (tb))) 
			        - t1 (tb) * (0.5L - b43 (tb) * t3 (tb)))
						 / (t2 (tb) * (t2 (tb) - t1 (tb)));
	b41 (tb) = (0.5L - b42 (tb) * t2 (tb) - b43 (tb) * t3 (tb)) / t1 (tb);
	b31 (tb) = ((1.L / 6.L - b42 (tb) * b21 (tb) * t1 (tb)) / b43 (tb)
			        - b32 (tb) * t2 (tb)) / t1 (tb);
	rk_b_4 (tb);
#if DEBUG_RK_4_3
	fprintf (stderr, "rk_tb_4_3: end\n");
#endif
}

/**
 * Function to calculate the objective function of a 4 steps 3rd order 
 * Runge-Kutta method.
 *
 * \return objective function value.
 */
long double
rk_objective_tb_4_3 (RK * rk) ///< RK struct.
{
	long double *tb;
	long double o;
#if DEBUG_RK_4_3
	fprintf (stderr, "rk_objective_tb_4_3: start\n");
#endif
	tb = rk->tb->coefficient;
	if (isnan (b31 (tb)) || isnan (b41 (tb)) || isnan (b42 (tb)))
	  {
	    o = INFINITY;
		goto end;
	  }
	o = fminl (0.L, b20 (tb));
	if (b30 (tb) < 0.L)
		o += b30 (tb);
	if (b31 (tb) < 0.L)
		o += b31 (tb);
	if (b40 (tb) < 0.L)
		o += b40 (tb);
	if (b41 (tb) < 0.L)
		o += b41 (tb);
	if (b42 (tb) < 0.L)
		o += b42 (tb);
	if (o < 0.L)
	  {
		  o = 40.L - o;
			goto end;
		}
	o = 30.L + fmaxl (1.L, fmaxl (t1 (tb), fmaxl (t2 (tb), t3 (tb))));
	rk_bucle_ac (rk);
	o = fminl (o, *rk->ac0->optimal);
end:
#if DEBUG_RK_4_3
	fprintf (stderr, "rk_objective_tb_4_3: optimal=%Lg\n", o);
	fprintf (stderr, "rk_objective_tb_4_3: end\n");
#endif
	return o;
}
