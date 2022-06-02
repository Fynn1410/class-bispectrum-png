

void   rsd_0_FFTLog(struct fft_struct *fft_input, double k, double *loops);
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
double complex IF2p2(double complex n1, double complex n2);
double complex IG2p(double complex n1);
double complex IPPp(double complex n1, double complex n2);
double complex IS2p(double complex n1, double complex n2);
double complex IF2p21(double complex n1, double complex n2, double mu);
double complex IF2p22(double complex n1, double complex n2);
double complex IG2pm1(double complex n1, double complex n2);
double complex IG2pm2(double complex n1, double mu);
double complex IPPpm1(double complex n1, double complex n2);
double complex IPPpm2(double complex n1, double complex n2, double mu);
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
double complex IG2pmm1(double complex n1, double mu);
double complex IG2pmm2(double complex n1, double mu);
double complex IPPpm31(double complex n1, double complex n2, double mu);
double complex IPPpm32(double complex n1, double complex n2, double mu);
double complex IPPm41(double complex n1, double complex n2, double mu);
double complex IPPm42(double complex n1, double complex n2, double mu);