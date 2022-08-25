
double pm_nowiggle_dst(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k,double z, int mode);
double pm_nowiggle_gfilter(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double z, int mode);
double pm_nowiggle_bspline(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double z, int mode);

double pm_nowiggle(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double z, double kf0, int cleanup, long SPLIT);
double pm_IR_LO(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double z,  long SPLIT);
double pm_IR_NLO(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double z, long SPLIT);

double IR_Sigma2_integrand(double x, void *par);
double IR_Sigma2(struct background * pba, struct primordial * ppm, struct fourier * pfo, double z,double kf0, long SPLIT);

double IR_del_Sigma2_integrand(double x, void *par);
double IR_del_Sigma2(struct background * pba, struct primordial * ppm, struct fourier * pfo,  double z, double kf0, long SPLIT);

double sigman2_integrand(double x, void *p);
double sigman(struct background * pba, struct primordial * ppm, struct fourier * pfo, double z, double k_min, double k_max, long n, long SPLIT);