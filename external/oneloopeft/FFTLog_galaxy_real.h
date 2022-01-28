

void pgloops_nonpropag(struct fft_struct  *fft_input, double k, double z, int cleanup, double *pg_loops);
double pgloops_propag(struct fft_struct  *fft_input, double k, double z, int cleanup);


double complex M_cIdelta2cG2(double complex nu1, double complex nu2);
double complex M_cIcG2cG2(double complex nu1, double complex nu2);
double complex M_cIdelta2delta2(double complex nu1, double complex nu2);
double complex M_cIcG2(double complex nu1, double complex nu2);
double complex M_cIdelta2(double complex nu1, double complex nu2);
double complex M_cFcG2(double complex nu1);
