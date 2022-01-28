

double complex Gamma(double complex z);
double complex Ifunc(double complex nu1, double complex nu2);
void FFT_compute_coeff(struct background * pba, struct primordial * ppm, struct fourier * pfo, double z, int nfft, double kmin_fft, double fft_bias, fftw_complex *biased_etam, fftw_complex *cmsym);

double FFT_kmax_Brent_solver(struct background * pba, struct primordial * ppm, struct fourier * pfo, double z, double kmin_fft, double fft_bias);
double FFT_kmax_Brent(double kmax, void *params);

double FFT_window(double k, double kmin, double kmax, double kleft, double kright);



