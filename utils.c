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
 * \file utils.c
 * \brief Source file with util functions.
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

gchar *error_message = NULL;    ///< error message string.

/**
 * Function to print an error message.
 */
void
show_error (const char *message)        ///< error message string.
{
  printf ("%s\n%s\n", _("ERROR!"), message);
}

/**
 * Function to get an integer number of a XML node property.
 *
 * \return Integer number value.
 */
int
xml_node_get_int (xmlNode * node,       ///< XML node.
                  const xmlChar * prop, ///< XML property.
                  int *error_code)      ///< Error code.
{
  int i = 0;
  xmlChar *buffer;
  buffer = xmlGetProp (node, prop);
  if (!buffer)
    *error_code = 1;
  else
    {
      if (sscanf ((char *) buffer, "%d", &i) != 1)
        *error_code = 2;
      else
        *error_code = 0;
      xmlFree (buffer);
    }
  return i;
}

/**
 * Function to get an unsigned integer number of a XML node property.
 *
 * \return Unsigned integer number value.
 */
unsigned int
xml_node_get_uint (xmlNode * node,      ///< XML node.
                   const xmlChar * prop,        ///< XML property.
                   int *error_code)     ///< Error code.
{
  unsigned int i = 0;
  xmlChar *buffer;
  buffer = xmlGetProp (node, prop);
  if (!buffer)
    *error_code = 1;
  else
    {
      if (sscanf ((char *) buffer, "%u", &i) != 1)
        *error_code = 2;
      else
        *error_code = 0;
      xmlFree (buffer);
    }
  return i;
}

/**
 * Function to get an unsigned integer number of a XML node property with a
 *   default value.
 *
 * \return Unsigned integer number value.
 */
unsigned int
xml_node_get_uint_with_default (xmlNode * node, ///< XML node.
                                const xmlChar * prop,   ///< XML property.
                                unsigned int default_value,
                                ///< default value.
                                int *error_code)        ///< Error code.
{
  unsigned int i;
  if (xmlHasProp (node, prop))
    i = xml_node_get_uint (node, prop, error_code);
  else
    {
      i = default_value;
      *error_code = 0;
    }
  return i;
}

/**
 * Function to get an long double number of a XML node property.
 *
 * \return Long double number value.
 */
long double
xml_node_get_float (xmlNode * node,     ///< XML node.
                    const xmlChar * prop,       ///< XML property.
                    int *error_code)    ///< Error code.
{
  long double x = 0.L;
  xmlChar *buffer;
  buffer = xmlGetProp (node, prop);
  if (!buffer)
    *error_code = 1;
  else
    {
      if (sscanf ((char *) buffer, "%Lg", &x) != 1)
        *error_code = 2;
      else
        *error_code = 0;
      xmlFree (buffer);
    }
  return x;
}

/**
 * Function to read the minimum, the interval and the random function type of a
 * variable on a XML node.
 *
 * \return 1 on success, 0 on error.
 */
int
read_variable (xmlNode * node,  ///< XML node.
               long double *minimum,    ///< minimum value.
               long double *interval,   ///< interval value.
               unsigned int *type,      ///< random function type.
               unsigned int n)  ///< variable number.
{
  char number[16];
  gchar *buffer;
  xmlChar *prop;
  int code;

  if (!node)
    {
      error_message = g_strdup (_("No XML node"));
      goto exit_on_error;
    }
  if (xmlStrcmp (node->name, XML_VARIABLE))
    {
      error_message = g_strdup (_("Bad XML node"));
      goto exit_on_error;
    }
  minimum[n] = xml_node_get_float (node, XML_MINIMUM, &code);
  if (code)
    {
      error_message = g_strdup (_("Bad minimum"));
      goto exit_on_error;
    }
  interval[n] = xml_node_get_float (node, XML_INTERVAL, &code);
  if (code)
    {
      error_message = g_strdup (_("Bad interval"));
      goto exit_on_error;
    }
  prop = xmlGetProp (node, XML_TYPE);
  if (!prop || !xmlStrcmp (prop, XML_RANDOM))
    type[n] = RANDOM_TYPE_UNIFORM;
  else if (!xmlStrcmp (prop, XML_BOTTOM))
    type[n] = RANDOM_TYPE_BOTTOM;
  else if (!xmlStrcmp (prop, XML_EXTREME))
    type[n] = RANDOM_TYPE_EXTREME;
  else if (!xmlStrcmp (prop, XML_TOP))
    type[n] = RANDOM_TYPE_TOP;
  else if (!xmlStrcmp (prop, XML_REGULAR))
    type[n] = RANDOM_TYPE_REGULAR;
  else if (!xmlStrcmp (prop, XML_ORTHOGONAL))
    type[n] = RANDOM_TYPE_ORTHOGONAL;
  else
    {
      xmlFree (prop);
      error_message = g_strdup (_("Bad random type function"));
      goto exit_on_error;
    }
  xmlFree (prop);
  return 1;

exit_on_error:
  snprintf (number, 16, "%u", n + 1);
  buffer = error_message;
  error_message
    = g_strconcat (_("Variable"), " ", number, ":\n", error_message, NULL);
  g_free (buffer);
  return 0;
}
