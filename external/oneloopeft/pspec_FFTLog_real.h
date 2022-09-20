int Real_Oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT);

int Real_Matter_IR_Resummed(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 int index_k, double z, long SPLIT, double *pk);

int Real_Galaxy_IR_Resummed(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                 int index_k, double z, long SPLIT, double *pk);