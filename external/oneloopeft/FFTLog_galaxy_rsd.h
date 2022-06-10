#include <complex.h>

void   rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *np_loops, double *p_loops);
void   rsd_1_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops, double *p_loops);
void   rsd_2_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops, double *p_loops);
void   rsd_3_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops, double *p_loops);
void   rsd_4_FFTLog(struct fft_struct *fft_input, double k, double mu, double *np_loops);

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup);
double P13_new(struct fft_struct *fft_input, double k, double z, int cleanup);

// Matrices

//0-th moment
double complex I2200(double complex n1, double complex n2, double mu);
double complex Idelta200(double complex n1, double complex n2, double mu);
double complex IG200(double complex n1, double complex n2, double mu);
double complex Idelta2delta200(double complex n1, double complex n2, double mu);
double complex IG2G200(double complex n1, double complex n2, double mu);
double complex Idelta2G200(double complex n1, double complex n2, double mu);

double complex I1300(double complex n1, double mu);
double complex FG200(double complex nu1, double mu);

//1-st moment
double complex I2201(double complex n1, double complex n2, double mu);
double complex Idelta201(double complex n1, double complex n2, double mu);
double complex IG201(double complex n1, double complex n2, double mu);
double complex FG201(double complex n1, double complex n2, double mu);
double complex J21101(double complex n1, double complex n2, double mu);
double complex Jdelta201(double complex n1, double complex n2, double mu);
double complex JG201(double complex n1, double complex n2, double mu);

double complex I1301(double complex n1, double mu);
double complex J12101(double complex n1, double mu);
double complex J11201(double complex n1, double mu);

//2-nd moment
double complex J21102(double complex n1, double complex n2, double mu);
double complex Jdelta202(double complex n1, double complex n2, double mu);
double complex JG202(double complex n1, double complex n2, double mu);
double complex I2211(double complex n1, double complex n2, double mu);
double complex J21111(double complex n1, double complex n2, double mu);
double complex N11a(double complex n1, double complex n2, double mu);
double complex N11b(double complex n1, double complex n2, double mu);
double complex N11c(double complex n1, double complex n2, double mu);

double complex J12102(double complex n1, double mu);
double complex I1311(double complex n1, double mu);
double complex J12111(double complex n1, double mu);
double complex J11211(double complex n1, double mu);

//3-rd moment
double complex J21112(double complex n1, double complex n2, double mu);
double complex N12a(double complex n1, double complex n2, double mu);
double complex N12b(double complex n1, double complex n2, double mu);

double complex J12112(double complex n1, double mu);

//4-th moment
double complex N22a(double complex n1, double complex n2, double mu);
double complex N22b(double complex n1, double complex n2, double mu);
double complex N22c(double complex n1, double complex n2, double mu);
