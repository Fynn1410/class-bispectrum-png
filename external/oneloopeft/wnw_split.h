int eft_ln_pk_nw_gfilter(struct background *pba, struct primordial *ppm, struct fourier *pfo, const int index_pk, const int index_k0, const int index_kmin, const int k_size, double *ln_pknw_array);
int eft_ln_pk_nw_gfilter_3d(struct background *pba, struct primordial *ppm, struct fourier *pfo, const int index_pk, const int index_k0, const int index_kmin, const int k_size, double *ln_pknw_array);
double eft_pk_nw_eisenstein_hu_factor(struct background *pba, struct primordial *ppm, struct fourier *pfo, const double k, const double k0);
double eft_pk_w_eisenstein_hu_factor(struct background *pba, struct primordial *ppm, struct fourier *pfo, const double k, const double k0);

double pk_nw_integrand(double x, void *par); 
double pk_Gfilter_nw( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double kf0, double z); 
double EH_PS_w( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double k0, double p0);
double EH_PS_nw( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k,double k0,double p0);

double T0( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k);
double T( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k);
double Tt0( struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double x1, double x2);




