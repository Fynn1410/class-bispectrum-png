
double eft_ir_sigma2(
              struct background * pba,
              struct primordial * ppm,
              struct fourier * pfo, 
              const double z, 
              const double k_split,
              const double k_bao
              );

double eft_ir_dsigma2(
              struct background * pba,
              struct primordial * ppm,
              struct fourier * pfo, 
              const double z, 
              const double k_split,
              const double k_bao
              );

int eft_ir_pk_lo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            double * ln_kvec,
            const int kvec_size,
            const double z,
            const double k_split,
            const double k_bao,
            double * out_pk
            );

int eft_ir_pk_nlo(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            enum linear_or_logarithmic mode,
            double * ln_kvec,
            const int kvec_size,
            const double z,
            const double k_split,
            const double k_bao,
            double * out_pk
            );

double eft_pk_moment(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo, 
            const int n,
            const double z, 
            enum eft_pk_type use_pk_type,
            const double k_split,
            const double k_bao
            );


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
