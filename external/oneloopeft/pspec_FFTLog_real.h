
int pm_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 int index_k, double z, long SPLIT, double *pk);

int pg_IR_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 int index_k, double z, double b1, double b2, double bG2, double btd, double R2, double cs2, long SPLIT, double *pk);