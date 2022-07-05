#include <complex.h>

void rsd_0_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double *np_loops, double *p_loops);
void rsd_1_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops, double *p_loops);
void rsd_2_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops, double *p_loops);
void rsd_3_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops, double *p_loops);
void rsd_4_FFTLog(struct oneloop_fftlog_workspace *fft_ws, double k, double mu, double *np_loops);

double P22_new(struct fft_struct *fft_input, double k, double z, int cleanup);
double P13_new(struct fft_struct *fft_input, double k, double z, int cleanup);

// Matrices

//0-th moment
double complex I2200(double complex n1, double complex n2);
double complex Idelta200(double complex n1, double complex n2);
double complex IG200(double complex n1, double complex n2);
double complex Idelta2delta200(double complex n1, double complex n2);
double complex IG2G200(double complex n1, double complex n2);
double complex Idelta2G200(double complex n1, double complex n2);

double complex I1300(double complex n1);
double complex FG200(double complex nu1);

//1-st moment
double complex I2201(double complex n1, double complex n2);
double complex Idelta201(double complex n1, double complex n2);
double complex IG201(double complex n1, double complex n2);
double complex FG201(double complex n1, double complex n2);
double complex J21101(double complex n1, double complex n2);
double complex Jdelta201(double complex n1, double complex n2);
double complex JG201(double complex n1, double complex n2);

double complex I1301(double complex n1);
double complex J12101(double complex n1);

//2-nd moment
double complex J21102x(double complex n1, double complex n2);
double complex J21102y(double complex n1, double complex n2);
double complex Jdelta202x(double complex n1, double complex n2);
double complex Jdelta202y(double complex n1, double complex n2);
double complex JG202x(double complex n1, double complex n2);
double complex JG202y(double complex n1, double complex n2);
double complex I2211(double complex n1, double complex n2);
double complex J21111(double complex n1, double complex n2);
double complex N11x(double complex n1, double complex n2);
double complex N11y(double complex n1, double complex n2);

double complex J12102x(double complex n1);
double complex J12102y(double complex n1);
double complex I1311(double complex n1);
double complex J12111(double complex n1);

//3-rd moment
double complex J21112x(double complex n1, double complex n2);
double complex J21112y(double complex n1, double complex n2);
double complex N12x(double complex n1, double complex n2);
double complex N12y(double complex n1, double complex n2);

double complex J12112x(double complex n1);
double complex J12112y(double complex n1);

//4-th moment
double complex N22x(double complex n1, double complex n2);
double complex N22y(double complex n1, double complex n2);
double complex N22z(double complex n1, double complex n2);
