

double complex Gamma(double complex z);
double complex Ifunc(double complex nu1, double complex nu2);
void FFT_compute_coeff(struct Cosmology *cosmo, c_datablock * block,  void *config_in, struct interp_ptrs *ptrs, double z, 
                      int nfft, double kmin_fft, double fft_bias, fftw_complex *biased_etam, fftw_complex *cmsym);

double FFT_kmax_Brent_solver(Interpolator2D *PK, double z, double kmin_fft, double fft_bias);
double FFT_kmax_Brent(double kmax, void *params);

double FFT_window(double k, double kmin, double kmax, double kleft, double kright);



