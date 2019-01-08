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
  unsigned int strong;          ///< boolean to cope with strong stability.
  unsigned int time_accuracy;
  ///< boolean to cope with extended time accuracy.
  unsigned int pair;            ///< boolean to build Runge-Kutta pairs.
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
#define e30(x) x[8]             ///< e30 Runge-Kutta pair coefficient.
#define e31(x) x[9]             ///< e31 Runge-Kutta pair coefficient.
#define t4(x) x[8]              ///< t4 Runge-Kutta coefficient.
#define b40(x) x[9]             ///< b40 Runge-Kutta coefficient.
#define b41(x) x[10]            ///< b41 Runge-Kutta coefficient.
#define b42(x) x[11]            ///< b42 Runge-Kutta coefficient.
#define b43(x) x[12]            ///< b43 Runge-Kutta coefficient.
#define e40(x) x[13]            ///< e40 Runge-Kutta pair coefficient.
#define e41(x) x[14]            ///< e41 Runge-Kutta pair coefficient.
#define e42(x) x[15]            ///< e42 Runge-Kutta pair coefficient.
#define t5(x) x[13]             ///< t5 Runge-Kutta coefficient.
#define b50(x) x[14]            ///< b50 Runge-Kutta coefficient.
#define b51(x) x[15]            ///< b51 Runge-Kutta coefficient.
#define b52(x) x[16]            ///< b52 Runge-Kutta coefficient.
#define b53(x) x[17]            ///< b53 Runge-Kutta coefficient.
#define b54(x) x[18]            ///< b53 Runge-Kutta coefficient.
#define e50(x) x[19]            ///< e50 Runge-Kutta pair coefficient.
#define e51(x) x[20]            ///< e51 Runge-Kutta pair coefficient.
#define e52(x) x[21]            ///< e52 Runge-Kutta pair coefficient.
#define e53(x) x[22]            ///< e53 Runge-Kutta pair coefficient.
#define t6(x) x[19]             ///< t6 Runge-Kutta coefficient.
#define b60(x) x[20]            ///< b60 Runge-Kutta coefficient.
#define b61(x) x[21]            ///< b61 Runge-Kutta coefficient.
#define b62(x) x[22]            ///< b62 Runge-Kutta coefficient.
#define b63(x) x[23]            ///< b63 Runge-Kutta coefficient.
#define b64(x) x[24]            ///< b64 Runge-Kutta coefficient.
#define b65(x) x[25]            ///< b65 Runge-Kutta coefficient.
#define e60(x) x[26]            ///< e60 Runge-Kutta pair coefficient.
#define e61(x) x[27]            ///< e61 Runge-Kutta pair coefficient.
#define e62(x) x[28]            ///< e62 Runge-Kutta pair coefficient.
#define e63(x) x[29]            ///< e63 Runge-Kutta pair coefficient.
#define e64(x) x[30]            ///< e64 Runge-Kutta pair coefficient.
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
#define a60(x) x[28]            ///< a60 Runge-Kutta coefficient.
#define a61(x) x[29]            ///< a61 Runge-Kutta coefficient.
#define a62(x) x[30]            ///< a62 Runge-Kutta coefficient.
#define a63(x) x[31]            ///< a63 Runge-Kutta coefficient.
#define a64(x) x[32]            ///< a64 Runge-Kutta coefficient.
#define a65(x) x[33]            ///< a65 Runge-Kutta coefficient.
#define c60(x) x[34]            ///< c60 Runge-Kutta coefficient.
#define c61(x) x[35]            ///< c61 Runge-Kutta coefficient.
#define c62(x) x[36]            ///< c62 Runge-Kutta coefficient.
#define c63(x) x[37]            ///< c63 Runge-Kutta coefficient.
#define c64(x) x[38]            ///< c64 Runge-Kutta coefficient.
#define c65(x) x[39]            ///< c65 Runge-Kutta coefficient.

void rk_print_tb (Optimize * tb, char *label, FILE * file);
void rk_print_maxima_2 (FILE * file);
void rk_print_maxima_3 (FILE * file);
void rk_print_maxima_4 (FILE * file);
void rk_print_maxima_5 (FILE * file);
void rk_print_maxima_6 (FILE * file);
void rk_bucle_ac (RK * rk);
int rk_run (xmlNode * node, gsl_rng ** rng);

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
 * Function to get \f$b_{i0}\f$ coefficients of the 6 steps Runge-Kutta methods.
 */
static inline void
rk_b_6 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  rk_b_5 (tb);
  b60 (tb) = t6 (tb) - b61 (tb) - b62 (tb) - b63 (tb) - b64 (tb) - b65 (tb);
}

/**
 * Function to get \f$\hat{b}_{40}\f$ coefficient of the 4 steps Runge-Kutta
 * pairs.
 */
static inline void
rk_e_4 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  e40 (tb) = 1.L - e41 (tb) - e42 (tb);
}

/**
 * Function to get \f$\hat{b}_{50}\f$ coefficient of the 5 steps Runge-Kutta
 * pairs.
 */
static inline void
rk_e_5 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  e50 (tb) = 1.L - e51 (tb) - e52 (tb) - e53 (tb);
}

/**
 * Function to get \f$\hat{b}_{60}\f$ coefficient of the 6 steps Runge-Kutta
 * pairs.
 */
static inline void
rk_e_6 (long double *tb)        ///< array of Runge-Kutta coefficients.
{
  e60 (tb) = 1.L - e61 (tb) - e62 (tb) - e63 (tb) - e64 (tb);
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

/**
 * Function to calculate the maximum CFL number in 5 steps Runge-Kutta.
 *
 * \return CFL value.
 */
static inline long double
rk_cfl_6 (long double *tb,      ///< array of Runge-Kutta t-b coefficients.
          long double *ac)      ///< array of Runge-Kutta a-c coefficients.
{
  return
    fminl (rk_cfl_5 (tb, ac),
           1.L / fmaxl (c60 (ac),
                        fmaxl (c61 (ac),
                               fmaxl (c62 (ac),
                                      fmaxl (c63 (ac),
                                             fmaxl (c64 (ac), c65 (ac)))))));
}

#endif
