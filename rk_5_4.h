/*
ODE: a program to get optime Runge-Kutta and multi-steps methods.

Copyright 2011-2019, Javier Burguete Tolosa.

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
 * \file rk_5_4.h
 * \brief Header file to optimize Runge-Kutta 5 steps 4th order methods.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2019.
 */
#ifndef RK_5_4__H
#define RK_5_4__H 1

int rk_tb_5_4 (Optimize * optimize);
int rk_tb_5_4t (Optimize * optimize);
int rk_tb_5_4p (Optimize * optimize);
int rk_tb_5_4tp (Optimize * optimize);
long double rk_objective_tb_5_4 (RK * rk);
long double rk_objective_tb_5_4t (RK * rk);
long double rk_objective_tb_5_4p (RK * rk);
long double rk_objective_tb_5_4tp (RK * rk);

#endif
