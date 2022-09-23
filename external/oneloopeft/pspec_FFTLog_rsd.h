
int RSD_Oneloop_FFTLog(struct background *pba, struct primordial *ppm, struct fourier *pfo,
                    int index_k, double z, long SPLIT);

int RSD_IR_Ressummed_default(struct fourier *pfo, struct background *pba, 
                    int index_k, double z, double mu, 
                    double * result);

int RSD_IR_Ressummed(struct fourier *pfo, struct background *pba, 
                    int index_k, double z, double mu, 
                    double b1, double b2, double bG2, double btd,
                    double c00, double c10, double c20, double c22, double c30, double c32, double c42, 
                    double * result);

int RSD_Multipole_default(struct fourier *pfo, struct background *pba, 
                          int index_k, double z, int l, 
                          double * result);

int RSD_Multipole(struct fourier *pfo, struct background *pba, int index_k, double z, int l,
     double b1, double b2, double bG2, double btd,
     double c00, double c10, double c20, double c22, double c30, double c32, double c42, 
     double * result);

double RSD_Multipole_integrand(double x, void *par);

double Legendre_Polynomial(int l, double mu);