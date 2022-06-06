


int pm_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    double k,  double z, long SPLIT, double * pk_nl);

int pg_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    double k,  double z, long SPLIT, double * pk_nl);

int rsd_oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 double k,  double z, double mu, double f, long SPLIT, double * pk_nl);

