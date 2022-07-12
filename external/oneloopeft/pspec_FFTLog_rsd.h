
int rsd_oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 double k, double z, double f, double mu, long SPLIT, double * pk_nl);

void FFTLog_rsd_init(struct background *pba, struct primordial *ppm, struct fourier *pfo, double z);

void FFTLog_fill_bias_vector(struct fourier *pfo, double b1, double b2, double bG2, double btd, double cs2, double R2);