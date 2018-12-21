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
 * \file ode.c
 * \brief Source file with the main function.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <alloca.h>
#include <time.h>
#include <unistd.h>
#include <glib.h>
#include <gsl/gsl_rng.h>
#if HAVE_MPI
#include <mpi.h>
#endif
#include "config.h"
#include "utils.h"
#include "optimize.h"
#include "steps.h"
#include "rk.h"


// OPTIMES
//
// RK_2_2 (cfl=1, cfl/steps=1/2)
// t1=1 t2=1
// a10=1 a20=a21=1/2
// b10=1 b20=b21=1/2
// c10=1 c20=0 c21=1
//
// RK_3_2 (cfl=2, cfl/steps=2/3)
// t1=1/2 t2=t3=1
// a10=1 a20=0 a21=1 a30=1/3 a31=0 a32=2/3
// b10=1/2 b20=b21=1/2 b30=b31=b32=1/3
// c10=1/2 c20=0 c21=1/2 c30=0 c31=0 c32=1/2
//
// RK_4_2 (cfl=3, cfl/steps=3/4)
// t1=1/3 t2=2/3 t3=t4=1
// a10=1 a20=0 a21=1 a30=a31=0 a32=1 a40=1/4 a41=a42=0 a43=3/4
// b10=1/3 b20=b21=1/2 b30=b31=b32=1/3
// c10=1/3 c20=0 c21=1/3 c30=c31=0 c32=1/3 c40=c41=c42=0 c43=1/3
//
// RK_5_2 (cfl=4, cfl/steps=4/5)
// t1=1/4 t2=1/2 t3=3/4 t4=t5=1
// a10=1 a20=0 a21=1 a30=a31=0 a32=1 a40=a41=a42=0 a43=1
// a50=1/5 a51=a52=a53=0 a54=4/5
// b10=b20=b21=b30=b31=b32=b40=b41=b42=b43=1/4 b50=b51=b52=b53=b54=1/5
// c10=1/4 c20=0 c21=1/4 c30=c31=0 c32=1/4 c40=c41=c42=0 c43=1/4
// c50=c51=c52=c53=0 c54=1/4
//
// RK_6_2 (cfl=5, cfl/steps=5/6)
// t1=1/5 t2=2/5 t3=3/5 t4=4/5 t5=t6=1
// a10=1 a20=0 a21=1 a30=a31=0 a32=1 a40=a41=a42=0 a43=1 a50=a51=a52=a53=0 a54=1
// a60=1/6 a61=a62=a63=a64=0 a65=5/6
// b10=b20=b21=b30=b31=b32=b40=b41=b42=b43b50=b51=b52=b53=b54=1/5
// b60=b61=b62=b63=b64=b65=1/6
// c10=1/5 c20=0 c21=1/5 c30=c31=0 c32=1/5 c40=c41=c42=0 c43=1/5
// c50=c51=c52=c53=0 c54=1/5 c60=c61=c62=c63=c64=0 c65=1/5
//
// RK_7_2 (cfl=6, cfl/steps=6/7)
// t1=1/6 t2=2/6 t3=3/6 t4=4/6 t5=5/6 t6=t7=1
// a10=1 a20=0 a21=1 a30=a31=0 a32=1 a40=a41=a42=0 a43=1 a50=a51=a52=a53=0 a54=1
// a60=a61=a62=a63=a64=0 a65=1 a70=1/7 a71=a72=a73=a74=a75=0 a76=6/7
// b10=b20=b21=b30=b31=b32=b40=b41=b42=b43b50=b51=b52=b53=b54=b60=b61=b62=b63
// =b64=b65=1/6 b70=b71=b72=b73=b74=b75=b76=1/7
// c10=1/6 c20=0 c21=1/6 c30=c31=0 c32=1/6 c40=c41=c42=0 c43=1/6
// c50=c51=c52=c53=0 c54=1/6 c60=c61=c62=c63=c64=0 c65=1/6
// c70=c71=c73=c74=c75=0 c76=1/6
//
// STEPS_3_2 (cfl=1/2)
// a0=3/4 a1=0 a2=1/4
// b0=3/2 b1=b2=0
// c0=2 c1=c2=0
//
// STEPS_4_2 (cfl=2/3)
// a0=8/9 a1=a2=0 a3=1/9
// b0=4/3 b1=b2=b3=0
// c0=3/2 c1=c2=c3=0
//
// STEPS_5_2 (cfl=3/4)
// a0=15/16 a1=a2=a3=0 a4=1/16
// b0=5/4 b1=b2=b3=b4=0
// c0=4/3 c1=c2=c3=c4=0
//
// STEPS_6_2 (cfl=4/5)
// a0=24/25 a1=a2=a3=a4=0 a5=1/25
// b0=6/5 b1=b2=b3=b4=b5=0
// c0=5/4 c1=c2=c3=c4=c5=0
//
// STEPS_7_2 (cfl=5/6)
// a0=35/36 a1=a2=a3=a4=0 a5=1/36
// b0=7/6 b1=b2=b3=b4=b5=0
// c0=6/5 c1=c2=c3=c4=c5=0
//
// STEPS_8_2 (cfl=6/7)
// a0=48/49 a1=a2=a3=a4=0 a5=1/49
// b0=8/7 b1=b2=b3=b4=b5=0
// c0=7/6 c1=c2=c3=c4=c5=0
//
// RK_3_3 (cfl=1, cfl/steps=1/3)
// t1=t3=1 t2=1/2
// a10=1 a20=3/4 a21=1/4 a30=1/3 a31=0 a32=2/3
// b10=1 b20=b21=1/4 b30=b31=1/6 b32=2/3
// c10=1 c20=0 c21=1 c30=0 c31=0 c32=1
//
// RK_4_3 (cfl=2, cfl/steps=1/2)
// t1=t3=1/2 t2=1
// a10=1 a20=0 a21=1 a30=2/3 a31=0 a32=1/3 a40=a41=a42=0 a43=1
// b10=b20=b21=1/2 b30=b31=b32=b40=b41=b42=1/6 b43=1/2
// c10=c20=c21=1/2 c30=0 c31=c32=c40=c41=c42=c43=1/2
//
// STEPS_4_3 (cfl=1/3)
// a0=16/27 a1=a2=0 a3=11/27
// b0=16/9 b1=b2=0 b3=4/9
// c0=3 c1=c2=0 c3=12/11
//
// STEPS_5_3 (cfl=1/2)
// a0=25/32 a1=a2=a3=0 a4=7/32
// b0=25/16 b1=b2=b3=0 b4=5/16
// c0=2 c1=c2=c3=0 c4=10/7
//
// RK_4_4
// b20=b30=b31=0
// t1=t2=b10=b21=1/2
// t3=t4=b32=1
// b40=b43=1/6
// b41=b42=1/3
//
// RK_5_4 (cfl=5.7383119999072648693e-01)
// a10=1
// a20=7.6291847536448934231e-01
// a21=2.3708152463551065768e-01
// a30=2.8277566215397998339e-03
// a31=9.4483076489543642366e-01
// a32=5.2341478483023776481e-02
// a40=6.6309349322657261432e-01
// a41=4.4608078655530437600e-04
// a42=3.3644655897571642926e-01
// a43=1.3867011155652078012e-05
// a50=6.1508837844244168818e-01
// a51=1.6215445203250317099e-01
// a52=1.5031540295088349761e-03
// a53=2.9995002701206228898e-03
// a54=2.1825451522542568299e-01
// b10=4.5307898578512145396e-01
// b20=1.0742850906542086858e-01
// b21=4.1289584245869493973e-01
// b30=4.3724780003845571699e-01
// b31=1.6681424138984873981e+00
// b32=9.0588475806896936944e-02
// b40=3.6664090929447934446e-02
// b41=1.3964776663947759666e-01
// b42=5.8495388577489235762e-01
// b43=1.9355255172522962676e-05
// b50=1.8734964237542205069e-01
// b51=3.0001844590932289148e-01
// b52=1.2903718344972795025e-01
// b53=3.4166232498545846666e-03
// b54=3.8017810501567252292e-01
// c10=4.5307898578512145396e-01
// c20=1.5535519922690283242e-05
// c21=1.7415774725317857109e+00
// c30=1.2525364871787093798e+00
// c31=1.7426727581493662919e+00
// c32=1.7307206145557779180e+00
// c40=4.7046995817421687883e-04
// c41=1.5854730818354572705e+00
// c42=1.7386197420604550316e+00
// c43=1.3957769958693602267e+00
// c50=1.6974106727520458494e-01
// c51=1.6275559956612365566e+00
// c52=7.2955698585694530597e-01
// c53=1.1376557995384390896e+00
// c54=1.7419025884665109550e+00
// 
// STEPS_6_4 (cfl=1.6475925238473621578e-01)
// a0=3.4246085571701207596e-01
// a1=a2=0
// a3=1.9179825943473607269e-01
// a4=9.3562124939009442453e-02
// a5=3.721787599092424089e-01
// b0=2.078553105578055306e+00
// b1=b2=b5=0
// b3=1.1641122222796929277e+00
// b4=5.6787174974870979874e-01
// c0=c3=c4=6.0694618695213438363e+00
// c1=c2=c5=0

/**
 * enum to define the method type.
 */
enum
{
  METHOD_TYPE_RUNGE_KUTTA = 1,  ///< Runge-Kutta method.
  METHOD_TYPE_STEPS = 2,        ///< multi-steps method.
  METHOD_TYPE_RUNGE_KUTTA_SIMPLE = 3,   ///< simple stable Runge-Kutta method.
  METHOD_TYPE_RUNGE_KUTTA_TIME = 4,
  ///< Runge-Kutta method with extended accuracy in time.
  METHOD_TYPE_RUNGE_KUTTA_TIME_SIMPLE = 5,
  ///< simple stable Runge-Kutta method with extended accuracy in time.
} MethodType;

unsigned int nsteps;            ///< steps number.
unsigned int order;             ///< accuracy order.

/**
 * Main function
 *
 * \return 0 on succes, error code on error.
 */
int
main (int argn,                 ///< arguments number.
      char **argc)              ///< argument chains array.
{
  char buffer[32];
  RK *rk;
  Optimize *s, *tb, *ac;
  gsl_rng *rng0, **rng;
  FILE *file;
  long double *value_optimal, *value_optimal2;
  long double optimal, optimal2, convergence_factor, convergence_factor2,
    search_factor, search_factor2;
  time_t d0;
  clock_t t0;
  unsigned long long int nsimulations, nsimulations2;
  unsigned long int seed;
  unsigned int i, j, k, type, niterations, niterations2, nsearch, nsearch2,
    nfree, nfree2;

#if HAVE_MPI
  // Init MPI
  MPI_Init (&argn, &argc);
#endif

  // Number of processors
  //nthreads = sysconf (_SC_NPROCESSORS_CONF);
  nthreads = 1;

  // Init the clock
  t0 = clock ();
  d0 = time (NULL);

#if HAVE_MPI

  // Nodes number
  MPI_Comm_size (MPI_COMM_WORLD, &nnodes);

  // Actual node
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

#else

  nnodes = 1;
  rank = 0;

#endif

  printf ("Rank=%d nnodes=%d nthreads=%u\n", rank, nnodes, nthreads);

  // Select the numerical model
  printf ("Selecting method\n");
  if ((argn != 10 && argn != 15)
      || sscanf (argc[1], "%u", &type) != 1
      || type < 1
      || type > 5
      || sscanf (argc[1], "%u", &type) != 1
      || sscanf (argc[2], "%u", &nsteps) != 1
      || sscanf (argc[3], "%u", &order) != 1
      || sscanf (argc[4], "%u", &niterations) != 1
      || sscanf (argc[5], "%Lu", &nsimulations) != 1
      || sscanf (argc[6], "%u", &nsearch) != 1
      || sscanf (argc[7], "%Lf", &convergence_factor) != 1
      || sscanf (argc[8], "%Lf", &search_factor) != 1
      || sscanf (argc[9], "%lu", &seed) != 1)
    {
      printf ("Usage is:\n"
              "./ode 2 steps order niterations nsimulations nsearch "
              "convergence_factor search_factor seed\n"
              "or\n"
              "./ode 1 steps order niterations nsimulations nsearch "
              "convergence_factor search_factor seed niterations2 "
              "nsimulations2 nsearch2 convergence_factor2 search_factor2\n");
      return 1;
    }
  if (type == 2 && argn != 10)
    {
      printf ("Bad steps parameters\n");
      return 2;
    }
  if (type == 3 && argn != 10)
    {
      printf ("Bad RK parameters\n");
      return 3;
    }
  if (type == 1 &&
      (argn != 15
       || sscanf (argc[10], "%u", &niterations2) != 1
       || sscanf (argc[11], "%Lu", &nsimulations2) != 1
       || sscanf (argc[12], "%u", &nsearch2) != 1
       || sscanf (argc[13], "%Lf", &convergence_factor2) != 1
       || sscanf (argc[14], "%Lf", &search_factor2) != 1))
    {
      printf ("Bad RK parameters\n");
      return 4;
    }


  printf ("type=%u steps=%u order=%u\n", type, nsteps, order);

  // Init a random numbers generator per node and thread
  printf ("Initing random numbers\n");
  rng = (gsl_rng **) g_slice_alloc (nnodes * nthreads * sizeof (gsl_rng *));
  rng0 = gsl_rng_alloc (gsl_rng_ranlxs2);
  gsl_rng_set (rng0, seed);
  for (i = k = 0; (int) i < nnodes; ++i)
    for (j = 0; j < nthreads; ++j, ++k)
      {
        rng[k] = gsl_rng_alloc (gsl_rng_taus2);
        gsl_rng_set (rng[k], gsl_rng_get (rng0));
      }

#if PRINT_RANDOM
  file_random = fopen ("random", "w");
#endif

  j = rank * nthreads;
  switch (type)
    {
    case METHOD_TYPE_RUNGE_KUTTA:
    case METHOD_TYPE_RUNGE_KUTTA_TIME:

#if PRINT_RANDOM
      file_random2 = fopen ("random2", "w");
#endif

      // Create optimize data
      rk = (RK *) alloca (nthreads * sizeof (RK));
      rk->strong = 1;
      if (type == METHOD_TYPE_RUNGE_KUTTA)
        rk->time_accuracy = 0;
      else
        rk->time_accuracy = 1;
      if (!rk_select (rk, nsteps, order))
        return 4;
      tb = rk->tb;
      ac = rk->ac0;
      nfree = tb->nfree;
      value_optimal
        = (long double *) g_slice_alloc (nfree * sizeof (long double));
      nfree2 = ac->nfree;
      value_optimal2
        = (long double *) g_slice_alloc (nfree2 * sizeof (long double));
      rk_create (rk, &optimal, value_optimal, &optimal2, value_optimal2,
                 convergence_factor, convergence_factor2, search_factor,
                 search_factor2, nsimulations, nsimulations2, nsearch, nsearch2,
                 niterations, niterations2);
      for (i = 1; i < nthreads; ++i)
        memcpy (rk + i, rk, sizeof (RK));
      for (i = 0; i < nthreads; ++i)
        rk_init (rk + i, rng[j + i]);

      // Method bucle
      printf ("Optimize bucle\n");
      rk_bucle_tb (rk);

      // Print the optimal coefficients
      printf ("Print the optimal coefficients\n");
      memcpy (tb->random_data, tb->value_optimal, nfree * sizeof (long double));
      tb->method (tb);
      memcpy (ac->random_data, ac->value_optimal,
              nfree2 * sizeof (long double));
      memcpy (rk->ac, ac, sizeof (Optimize));
      ac->method ((Optimize *) rk);
      snprintf (buffer, 32, "rk-%u-%u.mc", nsteps, order);
      file = fopen (buffer, "w");
      tb->print ((Optimize *) rk, file);
      tb->print_maxima (file, nsteps, order);
      fclose (file);

#if PRINT_RANDOM
      fclose (file_random2);
#endif

      // Free memory
      g_slice_free1 (nfree2 * sizeof (long double), value_optimal2);
      for (i = 0; i < nthreads; ++i)
        rk_delete (rk + i);
      break;

    case METHOD_TYPE_STEPS:

      // Create optimize data
      s = (Optimize *) alloca (nthreads * sizeof (Optimize));
      if (!steps_select (s, nsteps, order))
        return 2;
      nfree = s->nfree;
      value_optimal
        = (long double *) g_slice_alloc (nfree * sizeof (long double));
      optimize_create (s, &optimal, value_optimal, convergence_factor,
                       search_factor, nsimulations, nsearch, niterations);
      for (i = 1; i < nthreads; ++i)
        memcpy (s + i, s, sizeof (Optimize));
      for (i = 0; i < nthreads; ++i)
        optimize_init (s + i, rng[j + i]);

      // Method bucle
      fprintf (stderr, "Optimize bucle\n");
      optimize_bucle (s);

      // Print the optimal coefficients
      fprintf (stderr, "Print the optimal coefficients\n");
      memcpy (s->random_data, s->value_optimal, nfree * sizeof (long double));
      s->method (s);
      optimal = s->objective (s);
      snprintf (buffer, 32, "steps-%u-%u.mc", nsteps, order);
      file = fopen (buffer, "w");
      s->print (s, file);
      s->print_maxima (file, nsteps, order);
      fclose (file);

      // Free memory
      for (i = 0; i < nthreads; ++i)
        optimize_delete (s + i);
      break;

    case METHOD_TYPE_RUNGE_KUTTA_SIMPLE:
    case METHOD_TYPE_RUNGE_KUTTA_TIME_SIMPLE:

#if PRINT_RANDOM
      file_random2 = fopen ("random2", "w");
#endif
      // Create optimize data
      rk = (RK *) alloca (nthreads * sizeof (RK));
      rk->strong = 0;
      if (type == METHOD_TYPE_RUNGE_KUTTA_SIMPLE)
        rk->time_accuracy = 0;
      else
        rk->time_accuracy = 1;
      if (!rk_select (rk, nsteps, order))
        return 4;
      tb = rk->tb;
      nfree = tb->nfree;
      value_optimal
        = (long double *) g_slice_alloc (nfree * sizeof (long double));
      rk_create_simple (rk, &optimal, value_optimal, convergence_factor,
                        search_factor, nsimulations, nsearch, niterations);
      for (i = 1; i < nthreads; ++i)
        memcpy (rk + i, rk, sizeof (RK));
      for (i = 0; i < nthreads; ++i)
        rk_init (rk + i, rng[j + i]);

      // Method bucle
      printf ("Optimize bucle\n");
      rk_bucle_tb (rk);

      // Print the optimal coefficients
      printf ("Print the optimal coefficients\n");
      memcpy (tb->random_data, tb->value_optimal, nfree * sizeof (long double));
      tb->method (tb);
      snprintf (buffer, 32, "rk-%u-%u.mc", nsteps, order);
      file = fopen (buffer, "w");
      tb->print ((Optimize *) rk, file);
      tb->print_maxima (file, nsteps, order);
      fclose (file);

#if PRINT_RANDOM
      fclose (file_random2);
#endif

      // Free memory
      for (i = 0; i < nthreads; ++i)
        rk_delete (rk + i);

      break;

    default:
      printf ("Unknown method type\n");
      return 3;
    }

  printf ("optimal=%.19Le cpu time=%lg real time=%lu\n",
          optimal, (clock () - t0) / ((double) CLOCKS_PER_SEC),
          time (NULL) - d0);

#if PRINT_RANDOM
  fclose (file_random);
#endif

  // Free memory
  g_slice_free1 (nfree * sizeof (long double), value_optimal);
  j = nnodes * nthreads;
  for (i = 0; i < j; ++i)
    gsl_rng_free (rng[i]);
  g_slice_free1 (j * sizeof (gsl_rng *), rng);
  gsl_rng_free (rng0);

#if HAVE_MPI
  // Close MPI
  MPI_Finalize ();
#endif

  return 0;
}
