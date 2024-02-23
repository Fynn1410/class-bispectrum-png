#ifndef __EFT_TOOLS__
#define __EFT_TOOLS__

extern struct fft_struct;

#define _ERR_ARG_OUT_OF_RANGE_ 1
#define _ERR_RES_OUT_OF_RANGE_ 2
extern int _errno_util;
extern double complex _err_last_arg;

double complex cGamma(const double complex z);
double rGamma(const double z);

#endif
