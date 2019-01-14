#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline void
rk_print_maxima (FILE * file,   ///< file.
                 unsigned int nsteps,   ///< steps number.
                 unsigned int ncoefficients,    ///< coefficients number.
                 unsigned int order,    ///< accuracy order.
                 char label)    ///< coefficient label.
{
  unsigned int i, j, k, l;
  // b_{ij}=1 (1st order)
  for (i = 0; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u+", label, nsteps, i);
  fprintf (file, "-1;\n");
  // b_{ij}t_j=1/2 (2nd order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u+", label, nsteps, i, i);
  fprintf (file, "-1/2;\n");
  if (order < 2)
    return;
  // b_{ij}t_j^2=1/3 (3rd order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^2+", label, nsteps, i, i);
  fprintf (file, "-1/3;\n");
  if (order < 3)
    return;
  // b_{ij}b_{jk}t_k=1/6 (3rd order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/6;\n");
  // b_{ij}t_j^3=1/4 (4th order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^3+", label, nsteps, i, i);
  fprintf (file, "-1/4;\n");
  if (order < 4)
    return;
  // b_{ij}b_{jk}b_{kl}t_l=1/24 (4th order)
  for (i = 3; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 2; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 1; k < j; ++k)
            fprintf (file, "b%u%u*t%u+", j, k, k);
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/24;\n");
  // b_{ij}b_{jk}t_k^2=1/12 (4th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u^2+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/12;\n");
  // b_{ij}t_jb_{jk}t_k=1/8 (4th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u*(", label, nsteps, i, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/8;\n");
  // b_{ij}t_j^4=1/5 (5th order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^4+", label, nsteps, i, i);
  fprintf (file, "-1/5;\n");
  if (order < 5)
    return;
  // b_{ij}t_j^2b_{jk}t_k^2=1/10 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u^2*(", label, nsteps, i, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/10;\n");
  // b_{ij}b_{jk}t_k^3=1/20 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u^3+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/20;\n");
  // b_{ij}(b_{jk}t_k)^2=1/20 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u+", i, j, j);
      fprintf (file, "0)^2+");
    }
  fprintf (file, "-1/20;\n");
  // b_{ij}t_jb_{jk}t_k^2=1/15 (5th order)
  for (i = 2; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u*(", label, nsteps, i, i);
      for (j = 1; j < i; ++j)
        fprintf (file, "b%u%u*t%u^2+", i, j, j);
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/8;\n");
  // b_{ij}t_jb_{jk}b_{kl}t_l=1/24 (5th order)
  for (i = 3; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*t%u*(", label, nsteps, i, i);
      for (j = 2; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 1; k < j; ++k)
            fprintf (file, "b%u%u*t%u+", j, k, k);
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-7/120;\n");
  // b_{ij}b_{jk}b_{kl}t_l^2=1/60 (5th order)
  for (i = 3; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 2; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 1; k < j; ++k)
            fprintf (file, "b%u%u*t%u^2+", j, k, k);
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/60;\n");
  // b_{ij}b_{jk}b_{kl}b_{lm}t_m=1/120 (5th order)
  for (i = 4; i < ncoefficients; ++i)
    {
      fprintf (file, "%c%u%u*(", label, nsteps, i);
      for (j = 3; j < i; ++j)
        {
          fprintf (file, "b%u%u*(", i, j);
          for (k = 2; k < j; ++k)
            {
              fprintf (file, "b%u%u*(", j, k);
              for (l = 1; l < k; ++l)
                fprintf (file, "b%u%u*t%u+", k, l, l);
              fprintf (file, "0)+");
            }
          fprintf (file, "0)+");
        }
      fprintf (file, "0)+");
    }
  fprintf (file, "-1/120;\n");
  // b_{ij}t_j^5=1/6 (6th order)
  for (i = 1; i < ncoefficients; ++i)
    fprintf (file, "%c%u%u*t%u^5+", label, nsteps, i, i);
  fprintf (file, "-1/6;\n");
}

static inline void
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
  fprintf (file, "-1b0;\n");

  // 1st order
  fprintf (file, "b0");
  for (i = 1; i < nsteps; ++i)
    fprintf (file, "+b%u", i);
  for (i = 1; i < nsteps; ++i)
    fprintf (file, "-%ub0*a%u", i, i);
  fprintf (file, "-1b0;\n");

  // high order
  for (j = 2, m = 1; j <= order; ++j, m = -m)
    {
      for (i = 1; i < nsteps; ++i)
        {
          for (k = 1, l = i; k < j; ++k)
            l *= i;
          fprintf (file, "-%ub0*a%u", l, i);
        }
      for (i = 1; i < nsteps; ++i)
        {
          for (k = 2, l = i * j; k < j; ++k)
            l *= i;
          fprintf (file, "+%ub0*b%u", l, i);
        }
      fprintf (file, "+%db0;\n", m);
    }
}

int
main (int argn, char **argc)
{
  int nsteps, order;
  if (argn != 4)
    return 1;
  nsteps = atoi (argc[2]);
  order = atoi (argc[3]);
  if (!strcmp (argc[1], "rk"))
    rk_print_maxima (stdout, nsteps, nsteps, order, 'b');
  else if (!strcmp (argc[1], "steps"))
    steps_print_maxima (stdout, nsteps, order);
  return 0;
}
