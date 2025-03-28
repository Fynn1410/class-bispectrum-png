


/** utility function to calculate the triangle master integral
 * @param y                     Input: antiderivative variable
 * @param k12                   Input: wavenumber k1^2
 * @param k22                   Input: wavenumber k2^2
 * @param k32                   Input: wavenumber k3^2
 * @param M1                    Input: complex mass
 * @param M2                    Input: complex mass
 * @param M3                    Input: complex mass
 * @param out                   Output: pointer to output value
 * @return the error status
 */

 int util_tri_master(double y, double k12, double k22, double k32, class_complex M1, class_complex M2, class_complex M3, class_complex *out){
    class_complex M12 = M1*M1;      // M1^2
    class_complex M22 = M2*M2;      // M2^2
    class_complex M32 = M3*M3;      // M3^2
    double k14 = k12*k12;           // k1^4
    double k24 = k22*k22;           // k2^4
    double k34 = k32*k32;           // k3^4
    double CHOP_TOL = 1e-20;

	// antiderivative of the master integral integrand in y 
	double Num1 = 4.*k22*y+2.*k12-2.*k22-2.*k32;
	class_complex Num0 = -4.*k22*y+2.*M2-2.*M3+2.*k22;
	double DeltaR2 = -k12*y+k32*y-k32;
	class_complex DeltaR1 = -M2*y+M3*y+k12*y-k32*y+M1-M3+k32;
	class_complex DeltaR0 = M2*y-M3*y+M3;
	double DeltaS2 = -k14+2.*k12*k22+2.*k12*k32-k24+2.*k22*k32-k34;
	class_complex DeltaS1 =-4.*M1*k22-2.*M2*k12+2.*M2*k22+2.*M2*k32+2.*M3*k12+2.*M3*k22-2.*M3*k32-2.*k12*k22+2.*k24-2.*k22*k32;
	class_complex DeltaS0 =-M22+2.*M2*M3-2.*M2*k22-M32-2.*M3*k22-k24;

    class_complex DiakrS = sqrt(DeltaS1*DeltaS1 - 4.*DeltaS2*DeltaS0);

	class_complex solS1 = (-DeltaS1+DiakrS)/2./DeltaS2;
	class_complex solS2 = (-DeltaS1-DiakrS)/2./DeltaS2;

	class_complex cf2 = -(Num1*solS2+Num0)/DiakrS;
	class_complex cf1 = (Num1*solS1+Num0)/DiakrS;
		
    class_complex DiakrR = sqrt(DeltaR1*DeltaR1 - 4.*DeltaR2*DeltaR0);
				  
	class_complex solR1 = ((-DeltaR1+DiakrR)/2.)/DeltaR2;
	class_complex solR2 = ((-DeltaR1-DiakrR)/2.)/DeltaR2;

    class_complex temp1, temp2;
	if (cabs(cf1) < CHOP_TOL){
        // neglect cf1
        util_Fint(DeltaR2, solR1, solR2, solS2, &temp1);
        *out = cf2*temp1;
    } else if (cabs(cf2) < CHOP_TOL){
        // neglect cf2
        util_Fint(DeltaR2, solR1, solR2, solS1, &temp1);
        *out = cf1*temp1;
    } else{
        util_Fint(DeltaR2, solR1, solR2, solS2, &temp1);
        util_Fint(DeltaR2, solR1, solR2, solS1, &temp2);
        *out = cf2*temp1+cf1*temp2;
    }

    return _SUCCESS_;
}

int util_Fint(double aa, class_complex Y1, class_complex Y2, class_complex X0, class_complex *out){
    // calculates Fint. Most of the code of to find and count the branch cuts
    double CHOP_TOL = 1e-20;

	class_complex y1 = Y1;
	class_complex y2 = Y2;
	class_complex x0 = X0;

    // this is necessary because in gmpy2 sqrt(-1-0j) = -1j and not 1j -> TODO: not necessary in C?
	if (fabs(cimag(y2)) < CHOP_TOL){
        y2 = class_complex(creal(y2), 0.);
    }
	if (fabs(cimag(y1)) < CHOP_TOL){
        y1 = class_complex(creal(y1), 0.);
    }
	if (fabs(cimag(x0)) < CHOP_TOL){
        x0 = class_complex(creal(x0), 0.);
    }
		
	double rey1 = creal(y1);
	double imy1 = cimag(y1);
	double rey2 = creal(y2);
	double imy2 = cimag(y2);
	double rex0 = creal(x0);
	double imx0 = cimag(x0);
    
	double rey12 = rey1*rey1;
	double imy12 = imy1*imy1;
	double rey22 = rey2*rey2;
	double imy22 = imy2*imy2;
	double rex02 = rex0*rex0;
	double imx02 = imx0*imx0;

    int numbranchpoints = 0, signx0 = 0, sign = 0;

	double c = imy12*imy2*rex0 - imy1*imy22*rex0-imx02*imy2*rey1 + imx0*imy22*rey1-imy2*rex02*rey1 + imy2*rex0*rey12 + imx02*imy1*rey2-imx0*imy12*rey2+imy1*rex02*rey2-imx0*rey12*rey2-imy1*rex0*rey22+imx0*rey1*rey22;
	double a = imy1*rex0-imy2*rex0-imx0*rey1+imy2*rey1+imx0*rey2-imy1*rey2;
	double b = -imx02*imy1 + imx0*imy12+imx02*imy2-imy12*imy2-imx0*imy22+imy1*imy22-imy1*rex02+imy2*rex02+imx0*rey12-imy2*rey12-imx0*rey22+imy1*rey22;

    // if x0 is real there will always be a crossing through i or -i, which gives a cut of pi/2 instead of pi
    class_complex cutx0 = class_complex(0., 0.);

	if (0. < rex0 && rex0 < 1. && fabs(imx0) < CHOP_TOL){
		class_complex derivcritx0 = (y1 - y2)/2./sqrt(-(rex0-y1)*(rex0-y1))/(rex0-y2);	
		if (creal(derivcritx0) < 0.){
            signx0 = 1;
        } else {
            signx0 = -1;
        }
		cutx0 = signx0*_PI_/(sqrt(-rex0+y1+class_complex(0., 0.))*sqrt(rex0-y2+class_complex(0., 0.)));
    } else{
        cutx0 = class_complex(0., 0.);
    }

    // find remaining crossings of the imaginary axis
    int xsol_len;
    if (cabs(a) < CHOP_TOL){
        if (b != 0.){
            xsol_len = 1;
            double xsol[] = {- c / b};
        } else{
            xsol_len = 0;
        }		
    } else {
        if (b*b-4.*a*c > 0.){
            xsol_len = 1;
            double xsol[] = {(-b + sqrt(b*b-4.*a*c))/(2.*a),(-b - sqrt(b*b-4.*a*c))/(2.*a)};
        } else {
            //case where there is no intersection of the real axis (includes double zero)
            xsol_len = 0;
        }
    }

    int count=0;
    int i;
    if (xsol_len>0){
        for (i=0; i<xsol_len; i++){
            if (xsol[i] > CHOP_TOL && xsol[i] < 1. - CHOP_TOL && absc(class_complex(xsol[i], 0.)-x0)>CHOP_TOL){
                count+=1; 
            }
        }
        if (count>0){
            double xsol_new[count];
            for (i=0; i<xsol_len; i++){
                if (x > CHOP_TOL && x < 1. - CHOP_TOL && absc(class_complex(x, 0.)-x0)>CHOP_TOL){
                    xsol_new[i]=xsol[i];
                }
            }
        }
        class_complex atanarglist[count], abscrit[count], recrit[count];
        double temp;
        for (i=0; i<count; i++){
            temp = sqrt(class_complex(xsol_new[i], 0.)-y1)*sqrt(x0-y2)/(sqrt(-x0+y1)*sqrt(x-y2))
            atanarglist[i] = cabs(temp);
            abscrit[i] = creal(temp);
            if (abscrit[i] > 1 && cabs(recrit[i])<CHOP_TOL){
                numbranchpoints += 1;
                xbranch[i]=xsol_new[i];
            } else {
                xbranch[i]= class_complex(0., 0.);
            }
        } 
    }

    if (numbranchpoints==1){
        class_complex derivcrit_val;
        for (i=0; i<count; i++){
            if (xbranch[i]!=0){
                x = xbranch[i];
                derivcrit_val = sqrt(x0-y2)/sqrt(-x0+y1)*(1/(2.*sqrt(x-y1)*sqrt(x-y2)) -sqrt(x-y1)/(2.*(x-y2)*sqrt(x-y2)));
            }
        }

        if (creal(derivcrit_val) < 0){
            sign = 1;
        } else {
            sign = -1;
        }
        cut = sign*_PI_*2/(sqrt(-x0+y1)*sqrt(x0-y2));
    } else {
        cut = class_complex(0., 0.);
    }

    // TODO
	// class_complex prefac0 = class_complex(Prefactor(aa,y1,y2))
	// class_complex result = prefac0*(_SQRT_PI_/2.)*(cut + cutx0 + Antideriv(1,y1,y2,x0) - Antideriv(0,y1,y2,x0))


    return _SUCCESS_;
}

