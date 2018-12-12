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
 * \file rk.h
 * \brief Header file with common variables and functions to optimize
 *   Runge-Kutta methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#ifndef RK__H
#define RK__H 1

/**
 * \struct RK
 * \brief struct to define Runge-Kutta data and methods.
 */
typedef struct
{
	Optimize tb[1];
	///< Optimize struct to define t-b coefficients data and methods.
	Optimize ac[1];
	///< Optimize struct to define a-c coefficients data and methods.
	Optimize ac0[1];
	///< Optimize struct to define initial a-c coefficients data and methods.
} RK;

#define RUNGE_KUTTA(o) ((Optimize *)((Optimize *)o)->data)
///< macro to get the Optimize external data.
#define t1(x) x[0]              ///< t1 Runge-Kutta coefficient.
#define t2(x) x[1]              ///< t2 Runge-Kutta coefficient.
#define b20(x) x[2]             ///< b20 Runge-Kutta coefficient.
#define b21(x) x[3]             ///< b21 Runge-Kutta coefficient.
#define t3(x) x[4]              ///< t3 Runge-Kutta coefficient.
#define b30(x) x[5]             ///< b30 Runge-Kutta coefficient.
#define b31(x) x[6]             ///< b31 Runge-Kutta coefficient.
#define b32(x) x[7]             ///< b32 Runge-Kutta coefficient.
#define t4(x) x[8]              ///< t4 Runge-Kutta coefficient.
#define b40(x) x[9]             ///< b40 Runge-Kutta coefficient.
#define b41(x) x[10]            ///< b41 Runge-Kutta coefficient.
#define b42(x) x[11]            ///< b42 Runge-Kutta coefficient.
#define b43(x) x[12]            ///< b43 Runge-Kutta coefficient.
#define t5(x) x[13]              ///< t4 Runge-Kutta coefficient.
#define b50(x) x[14]             ///< b40 Runge-Kutta coefficient.
#define b51(x) x[15]            ///< b41 Runge-Kutta coefficient.
#define b52(x) x[16]            ///< b42 Runge-Kutta coefficient.
#define b53(x) x[17]            ///< b43 Runge-Kutta coefficient.
#define b54(x) x[18]            ///< b43 Runge-Kutta coefficient.
#define a20(x) x[0]             ///< a20 Runge-Kutta coefficient.
#define a21(x) x[1]             ///< a21 Runge-Kutta coefficient.
#define c20(x) x[2]             ///< c20 Runge-Kutta coefficient.
#define c21(x) x[3]             ///< c21 Runge-Kutta coefficient.
#define a30(x) x[4]             ///< a30 Runge-Kutta coefficient.
#define a31(x) x[5]             ///< a31 Runge-Kutta coefficient.
#define a32(x) x[6]             ///< a32 Runge-Kutta coefficient.
#define c30(x) x[7]             ///< c30 Runge-Kutta coefficient.
#define c31(x) x[8]             ///< c31 Runge-Kutta coefficient.
#define c32(x) x[9]             ///< c32 Runge-Kutta coefficient.
#define a40(x) x[10]            ///< a40 Runge-Kutta coefficient.
#define a41(x) x[11]            ///< a41 Runge-Kutta coefficient.
#define a42(x) x[12]            ///< a42 Runge-Kutta coefficient.
#define a43(x) x[13]            ///< a43 Runge-Kutta coefficient.
#define c40(x) x[14]            ///< c40 Runge-Kutta coefficient.
#define c41(x) x[15]            ///< c41 Runge-Kutta coefficient.
#define c42(x) x[16]            ///< c42 Runge-Kutta coefficient.
#define c43(x) x[17]            ///< c43 Runge-Kutta coefficient.
#define a50(x) x[18]            ///< a50 Runge-Kutta coefficient.
#define a51(x) x[19]            ///< a51 Runge-Kutta coefficient.
#define a52(x) x[20]            ///< a52 Runge-Kutta coefficient.
#define a53(x) x[21]            ///< a53 Runge-Kutta coefficient.
#define a54(x) x[22]            ///< a54 Runge-Kutta coefficient.
#define c50(x) x[23]            ///< c50 Runge-Kutta coefficient.
#define c51(x) x[24]            ///< c51 Runge-Kutta coefficient.
#define c52(x) x[25]            ///< c52 Runge-Kutta coefficient.
#define c53(x) x[26]            ///< c53 Runge-Kutta coefficient.
#define c54(x) x[27]            ///< c54 Runge-Kutta coefficient.

extern const long double minimum_ac_rk_2[1];
extern const long double interval_ac_rk_2[1];
extern const unsigned int random_ac_rk_2[1];
extern const long double minimum_ac_rk_3[3];
extern const long double interval_ac_rk_3[3];
extern const unsigned int random_ac_rk_3[3];
extern const long double minimum_ac_rk_4[6];
extern const long double interval_ac_rk_4[6];
extern const unsigned int random_ac_rk_4[6];
extern const long double minimum_ac_rk_5[10];
extern const long double interval_ac_rk_5[10];
extern const unsigned int random_ac_rk_5[10];

void rk_print_tb_2 (long double * tb, char *label, FILE *file);
void rk_print_tb_3 (long double * tb, char *label, FILE *file);
void rk_print_tb_4 (long double * tb, char *label, FILE *file);
void rk_print_tb_5 (long double * tb, char *label, FILE *file);
void rk_print_2 (RK * rk, FILE * file);
void rk_print_3 (RK * rk, FILE * file);
void rk_print_4 (RK * rk, FILE * file);
void rk_print_5 (RK * rk, FILE * file);
void rk_print_maxima_2 (FILE * file);
void rk_print_maxima_3 (FILE * file);
void rk_print_maxima_4 (FILE * file);
void rk_print_maxima_5 (FILE * file);
void rk_ac_2 (RK * rk);
void rk_ac_3 (RK * rk);
void rk_ac_4 (RK * rk);
void rk_ac_5 (RK * rk);
long double rk_objective_ac_2 (RK * rk);
long double rk_objective_ac_3 (RK * rk);
long double rk_objective_ac_4 (RK * rk);
long double rk_objective_ac_5 (RK * rk);
void rk_init (RK * rk, gsl_rng * rng);
void rk_delete (RK * rk);
void rk_create (RK * rk,
                      long double *optimal,
                      long double *value_optimal,
                      long double *optimal2,
                      long double *value_optimal2,
                      long double convergence_factor,
                      long double convergence_factor2,
                      long double search_factor,
                      long double search_factor2,
                      unsigned long long int nsimulations,
                      unsigned long long int nsimulations2,
                      unsigned int nsearch, 
                      unsigned int nsearch2, 
											unsigned int niterations,
										 	unsigned int niterations2);
void rk_step_ac (RK * rk);
void rk_bucle_ac (RK * rk);
void rk_step_tb (RK * rk);
void rk_bucle_tb (RK * rk);
int rk_select (RK * rk, unsigned int nsteps, unsigned int order);

/**
 * Function to get \f$b_{i0}\f$ coefficients of the 2 steps Runge-Kutta methods.
 */
static inline void
rk_b_2 (long double *tb)        ///< array of Runge-Kutta t-b coefficients.
{
  b20 (tb) = t2 (tb) - b21 (tb);
}

/**
 * Function to get \f$b_{i0}\f$ coefficients of the 3 steps Runge-Kutta methods.
 */
static inline void
rk_b_3 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  rk_b_2 (tb);
  b30 (tb) = t3 (tb) - b31 (tb) - b32 (tb);
}

/**
 * Function to get \f$b_{i0}\f$ coefficients of the 4 steps Runge-Kutta methods.
 */
static inline void
rk_b_4 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  rk_b_3 (tb);
  b40 (tb) = t4 (tb) - b41 (tb) - b42 (tb) - b43 (tb);
}

/**
 * Function to get \f$b_{i0}\f$ coefficients of the 5 steps Runge-Kutta methods.
 */
static inline void
rk_b_5 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  rk_b_4 (tb);
  b50 (tb) = t5 (tb) - b51 (tb) - b52 (tb) - b53 (tb) - b54 (tb);
}

/**
 * Function to calculate the maximum CFL number in 2 steps Runge-Kutta.
 *
 * \return CFL value.
 */
static inline long double
rk_cfl_2 (long double *tb,      ///< array of Runge-Kutta t-b coefficients.
          long double *ac)      ///< array of Runge-Kutta a-c coefficients.
{
  return 1.L / fmaxl (t1 (tb), fmaxl (c20 (ac), c21 (ac)));
}

/**
 * Function to calculate the maximum CFL number in 3 steps Runge-Kutta.
 *
 * \return CFL value.
 */
static inline long double
rk_cfl_3 (long double *tb,      ///< array of Runge-Kutta t-b coefficients.
          long double *ac)      ///< array of Runge-Kutta a-c coefficients.
{
  return fminl (rk_cfl_2 (tb, ac),
                1.L / fmaxl (c30 (ac), fmaxl (c31 (ac), c32 (ac))));
}

/**
 * Function to calculate the maximum CFL number in 4 steps Runge-Kutta.
 *
 * \return CFL value.
 */
static inline long double
rk_cfl_4 (long double *tb,      ///< array of Runge-Kutta t-b coefficients.
          long double *ac)      ///< array of Runge-Kutta a-c coefficients.
{
  return fminl (rk_cfl_3 (tb, ac),
                1.L / fmaxl (c40 (ac),
									           fmaxl (c41 (ac), fmaxl (c42 (ac), c43 (ac)))));
}

/**
 * Function to calculate the maximum CFL number in 5 steps Runge-Kutta.
 *
 * \return CFL value.
 */
static inline long double
rk_cfl_5 (long double *tb,      ///< array of Runge-Kutta t-b coefficients.
          long double *ac)      ///< array of Runge-Kutta a-c coefficients.
{
  return 
		fminl (rk_cfl_4 (tb, ac),
           1.L / fmaxl (c50 (ac),
									      fmaxl (c51 (ac), 
															 fmaxl (c52 (ac), fmaxl (c53 (ac), c54 (ac))))));
}

#endif
