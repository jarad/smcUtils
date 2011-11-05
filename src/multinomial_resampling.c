#include <R.h>
#include <Rmath.h>

/* Performs multinomial resampling on the weights w. Returns the indices for num_samples samples */
void multinomial_resampling(double *w, int n, int *indices, int num_samples) {
  int i, j;
  double cusum[n], random_uniform;

  /* Calculate the cumulative sum of the weight series */
  cusum[0] = w[0];
  for (i=1; i<n; i++) {
    cusum[i] = cusum[i-1]+w[i];
  }
  GetRNGstate();
  for (i=0; i<num_samples; i++) {
    random_uniform = runif(0,1);
    j=0;
    while (cusum[j]<random_uniform) { j++; }
    indices[i] = j;
  }
  PutRNGstate();
}

void multinomial_resampling_wrap(double *w, int *n, int *indices, int *num_samples) {
  multinomial_resampling(w, *n, indices, *num_samples);
}

