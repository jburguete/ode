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
 * \file config.h
 * \brief Configuration header file.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2018.
 */
#ifndef CONFIG__H
#define CONFIG__H 1

#define EPSILON 1
#define MAXIMA_PRECISION 20
#define OPTIMIZE_STEPS_11_6 0

#define XML_AC                 (const xmlChar *) "ac"
#define XML_BOTTOM             (const xmlChar *) "bottom"
#define XML_CLIMBING_FACTOR    (const xmlChar *) "climbing-factor"
#define XML_CONVERGENCE_FACTOR (const xmlChar *) "convergence-factor"
#define XML_EXTREME            (const xmlChar *) "extreme"
#define XML_INTERVAL           (const xmlChar *) "interval"
#define XML_MINIMUM            (const xmlChar *) "minimum"
#define XML_NCLIMBINGS         (const xmlChar *) "nclimbings"
#define XML_NITERATIONS        (const xmlChar *) "niterations"
#define XML_NO                 (const xmlChar *) "no"
#define XML_NSIMULATIONS       (const xmlChar *) "nsimulations"
#define XML_ORDER              (const xmlChar *) "order"
#define XML_ORTHOGONAL         (const xmlChar *) "orthogonal"
#define XML_PAIR               (const xmlChar *) "pair"
#define XML_RANDOM             (const xmlChar *) "random"
#define XML_REGULAR            (const xmlChar *) "regular"
#define XML_RUNGE_KUTTA        (const xmlChar *) "Runge-Kutta"
#define XML_STEPS              (const xmlChar *) "steps"
#define XML_STRONG             (const xmlChar *) "strong"
#define XML_TIME_ACCURACY      (const xmlChar *) "time-accuracy"
#define XML_TOP                (const xmlChar *) "top"
#define XML_TYPE               (const xmlChar *) "type"
#define XML_VARIABLE           (const xmlChar *) "variable"
#define XML_YES                (const xmlChar *) "yes"

#endif
