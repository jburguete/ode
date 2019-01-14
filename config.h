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
 * \file config.h
 * \brief Configuration header file.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2019.
 */
#ifndef CONFIG__H
#define CONFIG__H 1

#define EPSILON 1
///< using epsilon to avoid small numbers solving linear equations systems.
#define MAXIMA_PRECISION 20     ///< precision digits on maxima files.
#define OPTIMIZE_STEPS_11_6 0
///< special optimization for 11 steps 6th order multi-steps method.

#define XML_AC                 (const xmlChar *) "ac"
///< ac XML label.
#define XML_BOTTOM             (const xmlChar *) "bottom"
///< bottom XML label.
#define XML_CLIMBING_FACTOR    (const xmlChar *) "climbing-factor"
///< climbing-factor XML label.
#define XML_CONVERGENCE_FACTOR (const xmlChar *) "convergence-factor"
///< convergence-factor XML label.
#define XML_EXTREME            (const xmlChar *) "extreme"
///< extreme XML label.
#define XML_INTERVAL           (const xmlChar *) "interval"
///< interval XML label.
#define XML_MINIMUM            (const xmlChar *) "minimum"
///< minimum XML label.
#define XML_NCLIMBINGS         (const xmlChar *) "nclimbings"
///< nclimbings XML label.
#define XML_NITERATIONS        (const xmlChar *) "niterations"
///< niterations XML label.
#define XML_NO                 (const xmlChar *) "no"
///< no XML label.
#define XML_NSIMULATIONS       (const xmlChar *) "nsimulations"
///< nsimulations XML label.
#define XML_ORDER              (const xmlChar *) "order"
///< order XML label.
#define XML_ORTHOGONAL         (const xmlChar *) "orthogonal"
///< orthogonal XML label.
#define XML_PAIR               (const xmlChar *) "pair"
///< pair XML label.
#define XML_RANDOM             (const xmlChar *) "random"
///< random XML label.
#define XML_REGULAR            (const xmlChar *) "regular"
///< regular XML label.
#define XML_RUNGE_KUTTA        (const xmlChar *) "Runge-Kutta"
///< Runge-Kutta XML label.
#define XML_STEPS              (const xmlChar *) "steps"
///< steps XML label.
#define XML_STRONG             (const xmlChar *) "strong"
///< strong XML label.
#define XML_TIME_ACCURACY      (const xmlChar *) "time-accuracy"
///< time-accuracy XML label.
#define XML_TOP                (const xmlChar *) "top"
///< top XML label.
#define XML_TYPE               (const xmlChar *) "type"
///< type XML label.
#define XML_VARIABLE           (const xmlChar *) "variable"
///< variable XML label.
#define XML_YES                (const xmlChar *) "yes"
///< yes XML label.

#endif
