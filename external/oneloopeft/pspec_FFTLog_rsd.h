
int rsd_oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT);

int RSD_IR_Ressummed(struct fourier *pfo, struct background *pba,
     int index_k, double z, double mu, 
     double b1, double b2, double bG2, double btd, 
     double c00, double c10, double c20, double c22, double c30, double c32, double c42,
     double * result);

int RSD_Multipole(struct fourier *pfo, struct background *pba, int index_k, double z, double b1, double b2, double bG2, double btd,
     double c00, double c10, double c20, double c22, double c30, double c32, double c42, int l, double * result);

double RSD_Multipole_integrand(double x, void *par);

double Legendre_Polynomial(int l, double mu);