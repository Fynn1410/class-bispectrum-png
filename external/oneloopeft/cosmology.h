#ifndef COSMOLOGY_H
#define COSMOLOGY_H

 double Pk_dlnPk(struct background * pba, struct primordial * ppm, struct fourier * pfo, double k, double z , int mode);

double growth_D(struct background * pba, double z);
double growth_f(struct background * pba, double z);


#endif