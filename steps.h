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
 * \file steps.h
 * \brief Header file with common variables and functions to optimize
 *   multi-steps methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#ifndef STEPS__H
#define STEPS__H 1

#define c0(x) x[0]
///< c0 multi-steps coefficient.
#define a0(x) x[1]
///< a0 multi-steps coefficient.
#define c1(x) x[2]
///< c1 multi-steps coefficient.
#define a1(x) x[3]
///< a1 multi-steps coefficient.
#define c2(x) x[4]
///< c2 multi-steps coefficient.
#define a2(x) x[5]
///< a2 multi-steps coefficient.
#define c3(x) x[6]
///< c3 multi-steps coefficient.
#define a3(x) x[7]
///< a3 multi-steps coefficient.
#define c4(x) x[8]
///< c4 multi-steps coefficient.
#define a4(x) x[9]
///< a4 multi-steps coefficient.
#define c5(x) x[10]
///< c5 multi-steps coefficient.
#define a5(x) x[11]
///< a5 multi-steps coefficient.
#define c6(x) x[12]
///< c6 multi-steps coefficient.
#define a6(x) x[13]
///< a6 multi-steps coefficient.
#define c7(x) x[14]
///< c7 multi-steps coefficient.
#define a7(x) x[15]
///< a7 multi-steps coefficient.
#define c8(x) x[16]
///< c8 multi-steps coefficient.
#define a8(x) x[17]
///< a8 multi-steps coefficient.
#define c9(x) x[18]
///< c9 multi-steps coefficient.
#define a9(x) x[19]
///< a9 multi-steps coefficient.
#define c10(x) x[20]
///< c10 multi-steps coefficient.
#define a10(x) x[21]
///< a10 multi-steps coefficient.
#define b0(x) (a0(x)*c0(x))
///< b0 multi-steps coefficient.
#define b1(x) (a1(x)*c1(x))
///< b1 multi-steps coefficient.
#define b2(x) (a2(x)*c2(x))
///< b2 multi-steps coefficient.
#define b3(x) (a3(x)*c3(x))
///< b3 multi-steps coefficient.
#define b4(x) (a4(x)*c4(x))
///< b4 multi-steps coefficient.
#define b5(x) (a5(x)*c5(x))
///< b5 multi-steps coefficient.
#define b6(x) (a6(x)*c6(x))
///< b6 multi-steps coefficient.
#define b7(x) (a7(x)*c7(x))
///< b7 multi-steps coefficient.
#define b8(x) (a8(x)*c8(x))
///< b8 multi-steps coefficient.
#define b9(x) (a9(x)*c9(x))
///< b9 multi-steps coefficient.
#define b10(x) (a10(x)*c10(x))
///< b10 multi-steps coefficient.

void steps_print_maxima (FILE * file, unsigned int nsteps, unsigned int order);
void steps_print_1 (Optimize * optimize, FILE * file);
void steps_print_2 (Optimize * optimize, FILE * file);
void steps_print_3 (Optimize * optimize, FILE * file);
void steps_print_4 (Optimize * optimize, FILE * file);
void steps_print_5 (Optimize * optimize, FILE * file);
void steps_print_6 (Optimize * optimize, FILE * file);
void steps_print_7 (Optimize * optimize, FILE * file);
void steps_print_8 (Optimize * optimize, FILE * file);
void steps_print_9 (Optimize * optimize, FILE * file);
void steps_print_10 (Optimize * optimize, FILE * file);
void steps_print_11 (Optimize * optimize, FILE * file);
long double steps_objective_2 (Optimize * optimize);
long double steps_objective_3 (Optimize * optimize);
long double steps_objective_4 (Optimize * optimize);
long double steps_objective_5 (Optimize * optimize);
long double steps_objective_6 (Optimize * optimize);
long double steps_objective_7 (Optimize * optimize);
long double steps_objective_8 (Optimize * optimize);
long double steps_objective_9 (Optimize * optimize);
long double steps_objective_10 (Optimize * optimize);
long double steps_objective_11 (Optimize * optimize);
int steps_select (Optimize * optimize, unsigned int nsteps, unsigned int order);
int steps_run (xmlNode * node, gsl_rng **rng);

#endif
