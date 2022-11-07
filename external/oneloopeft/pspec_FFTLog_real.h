int Real_Oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT);

int Real_Matter_IR_Resummed(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 int index_k, double z, long SPLIT, double *pk);

int Real_Galaxy_IR_Resummed_default(struct fourier *pfo, struct background *pba, struct primordial *ppm,
                    int index_k, double z, long SPLIT, double *pk);

int Real_Galaxy_IR_Resummed(struct fourier *pfo, struct background *pba, struct primordial *ppm,
                    int index_k, double z,
                    double b1, double b2, double bG2, double btd, double R2,
                    double cs2,
                    long SPLIT,
                    double *pk);