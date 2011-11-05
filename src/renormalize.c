#include <Rmath.h>

/* Renormalizes the vector w to sum to 1. If log \ne 0, the weights are assumed to be log weights */
void renormalize(double *w, int n, int log) {
  int i;

  if (log) {
    double max=w[0];
    for (i=1; i<n; i++) {
      if (w[i]>max) max = w[i];
    }
    for (i=0; i<n; i++) {
      w[i] = exp(w[i]-max);
    }
  }

  double sum=0;
  for (i=0; i<n; i++) {
    sum += w[i];
  }
  for (i=0; i<n; i++) {
    w[i] /= sum;
  }
}

void renormalize_wrap(double *w, int *n, int *log) {
  renormalize(w, *n, *log);
}

