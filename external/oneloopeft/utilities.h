#ifndef __EFT_TOOLS__
#define __EFT_TOOLS__

#define _ERR_ARG_OUT_OF_RANGE_ 1
#define _ERR_RES_OUT_OF_RANGE_ 2

#include "carray.h"

extern int _errno_util;
extern class_complex _err_last_arg;

class_complex cGamma(class_complex z);
double rGamma(double z);

#endif
