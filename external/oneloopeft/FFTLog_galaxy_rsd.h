#include <complex.h>

void   rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *loops);
void   rsd_1_FFTLog(struct fft_struct *fft_input, double k, double *loops);
void   rsd_2_FFTLog(struct fft_struct *fft_input, double k, double mu, double *loops);
void   rsd_3_FFTLog(struct fft_struct *fft_input, double k, double mu, double *loops);
void   rsd_4_FFTLog(struct fft_struct *fft_input, double k, double mu, double *loops);

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup);
double P13_new(struct fft_struct *fft_input, double k, double z, int cleanup);

// Matrices
double complex IF2(double complex n1, double complex n2);
double complex IF2S2(double complex n1, double complex n2);
double complex IPP(double complex n1, double complex n2);
double complex IS2(double complex n1, double complex n2);
double complex IS2S2(double complex n1, double complex n2);

double complex IF2G2(double complex n1, double complex n2);
double complex IF3G3(double complex n1);
double complex IG2(double complex n1, double complex n2);
double complex IS2G2(double complex n1, double complex n2);
double complex IF2p(double complex n1, double complex n2);
double complex IF2p2(double complex n1);
double complex IG2p(double complex n1);
double complex IPPp(double complex n1, double complex n2);
double complex IS2p(double complex n1, double complex n2);

double complex IF2pm3(double complex n1, double complex n2, double mu);
double complex IG2pm2(double complex n1, double complex n2, double mu);
double complex IPPpm(double complex n1, double complex n2, double mu);
double complex IS2pm(double complex n1, double complex n2, double mu);
double complex IG2G2(double complex n1, double complex n2);
double complex IG3(double complex n1);
double complex IG2pm(double complex n1, double complex n2);
double complex IG2pmm(double complex n1, double mu);
double complex IF2pm(double complex n1);
double complex IPPmm1(double complex n1, double complex n2, double mu);
double complex IPPmm2(double complex n1, double complex n2);
double complex IPPmm3(double complex n1, double complex n2, double mu);

double complex IG2pm3(double complex n1, double complex n2, double mu);
double complex IG2pmm3(double complex n1, double mu);
double complex IPPpm31(double complex n1, double complex n2, double mu);
double complex IPPpm32(double complex n1, double complex n2, double mu);

double complex IPPm41(double complex n1, double complex n2, double mu);
double complex IPPm42(double complex n1, double complex n2, double mu);
double complex IPPm43(double complex n1, double complex n2, double mu);