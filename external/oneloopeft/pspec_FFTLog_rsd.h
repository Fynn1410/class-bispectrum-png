
int rsd_oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT);

int RSD_IR_Ressummed(struct fourier *pfo, double f, int index_k, double z, double mu, double * result);

int RSD_Multipole(struct fourier *pfo, double f, int index_k, double z, int l, double * result);

double RSD_Multipole_integrand(double x, void *par);

double Legendre_Polynomial(int l, double mu);