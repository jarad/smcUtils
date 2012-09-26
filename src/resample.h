
int compare_doubles(const void *, const void *);

void is_increasing_R(int *, const double *, int*);
int is_increasing(int , const double *);

void cumulative_sum_R(int *, double *);
int cumulative_sum(int , double *);

void rep2id_R(int *, int *, int *);
int rep2id(int *, int , int *);

void inverse_cdf_weights_R(int *, double *, int *, double *, int *);
int inverse_cdf_weights(int , double *, int , double *, int *);

void ess_R(int *, double *, double *);
double ess(int , double *);

void cov2_R(int *, double *, double *);
double cov2(int , double *);

void entropy_R(int *, double *, double *);
double entropy(int , double *);


void resample_R(int *, double *, int *, int *, int *);
int resample(int , double *, int , int *, int );

void multinomial_resample_R(int *, double *, int *, int *);
int multinomial_resample(int, double *, int, int *);

void stratified_resample_R(int *, double *, int *, int *);
int stratified_resample(int , double *, int , int *);

void systematic_resample_R(int *, double *, int *, int *);
int systematic_resample(int , double *, int , int *);

void residual_resample_R(int *, double *, int *, int *, int *);
int residual_resample(int , double *, int , int *, int );


