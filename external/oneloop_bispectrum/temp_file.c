
/** Calculate arctan*prefactor, TODO: reference equation in thesis or paper
 * @param x                     Input: integration variable, evaluated at x=0 and x=1
 * @param y1                    Input: complex variable z_+
 * @param y2                    Input: complex variable z_-
 * @param x0                    Input: complex variable, evaluated for x_- or x_+
 * @param out                   Output: pointer to output value
 * @return the error status
 */

 int util_antideriv(
    struct gen_tri_integral *pti,
    double x_in,
    class_complex y1,
    class_complex y2,
    class_complex x0,
    class_complex *out
    ){
double CHOP_TOL = 1.e-14;
class_complex x = class_complex(x_in, 0.);

if (cabs(y1-x0)<CHOP_TOL || cabs(x0-y2)<CHOP_TOL){
printf("\nERROR in util_antideriv: division by zero. either x0_re_c=y1 or x0_re_c=y2.\n");
}

// case where x0 = y2 = 0 or 1
if (cabs(y2 - x0)<CHOP_TOL){
if (cabs(y2-x)<CHOP_TOL){
printf("\nspecial case: antideriv = 0");
*out = class_complex(0., 0.);
return _SUCCESS_;
}
*out = 2.*sqrt(x-y1)/(-x0+y1)/sqrt(x-y2);
return _SUCCESS_;
}

// case where x0 = y2 = 0 or 1
if (cabs(x0-y1)<CHOP_TOL){
printf("WARNING: switching var in Antideriv");
class_call(util_antideriv(pti,
                   x_in,
                   y2,
                   y1,
                   x0,
                   out),
     pti->error_message,
     pti->error_message);
return _SUCCESS_;
}

class_complex prefac = 2./(sqrt(-x0+y1)*sqrt(x0-y2));
class_complex temp = sqrt(x-y1)*sqrt(x0-y2)/sqrt(-x0+y1);
class_complex LimArcTan;

if (x_in==1. && cabs(y2-class_complex(1., 0.))<CHOP_TOL ){
LimArcTan = class_complex(0., 1.)*sqrt(-temp*temp) * _PI_/(2.*temp);
*out = prefac * LimArcTan;
// printf("\nInside util_antideriv: SPECIAL CASE 1");
// printf("\ntemp = %f + %fj", creal(temp), cimag(temp));
// printf("\nprefac = %f + %fj", creal(prefac), cimag(prefac));
// printf("\nLimArcTan = %f + %fj", creal(LimArcTan), cimag(LimArcTan));
// printf("\nout = %f + %fj\n\n", creal(*out), cimag(*out));
return _SUCCESS_;
}
if (x_in==0. && cabs(y2)<CHOP_TOL){
LimArcTan = sqrt(temp*temp) * _PI_/(2.*temp);
*out = prefac * LimArcTan;
// printf("\nInside util_antideriv: SPECIAL CASE 2");
// printf("\ntemp = %f + %fj", creal(temp), cimag(temp));
// printf("\nprefac = %f + %fj", creal(prefac), cimag(prefac));
// printf("\nLimArcTan = %f + %fj", creal(LimArcTan), cimag(LimArcTan));
// printf("\nout = %f + %fj\n\n", creal(*out), cimag(*out));
return _SUCCESS_;
}

class_complex arg_atan = temp/sqrt(x-y2);
class_complex atan_ = atan(arg_atan);
*out = prefac * atan_;

// printf("\nInside util_antideriv:");
// printf("\ny1 = %f + %fj", creal(y1), cimag(y1));
// printf("\ny2 = %f + %fj", creal(y2), cimag(y2));
// printf("\nx0 = %f + %fj", creal(x0), cimag(x0));
// printf("\nx_in = %f", x_in);
// printf("\nprefac = %f + %fj", creal(prefac), cimag(prefac));
// printf("\ntemp = %f + %fj", creal(temp), cimag(temp));
// printf("\nx-y2 = %f + %fj", creal(x-y2), cimag(x-y2));
// printf("\nsqrt(x-y2) = %f + %fj", creal(sqrt(x-y2)), cimag(sqrt(x-y2)));
// printf("\natan_ = %f + %fj", creal(atan_), cimag(atan_));
// printf("\narg atan_ = %f + %fj", creal(temp/sqrt(x-y2)), cimag(temp/sqrt(x-y2)));
// printf("\nout = %f + %fj\n\n", creal(*out), cimag(*out));

return _SUCCESS_;
}


/** Calculate TODO: not completely sure what this is exactly
* I thought this is s(z+, z-), but I am confused about the 1/sqrt(a) since a<0
* @param a                     Input: integer
* @param y1                    Input: complex variable z_+
* @param y2                    Input: complex variable z_-
* @param out                   Output: pointer to output value
* @return the error status
*/

int util_prefactor(
    struct gen_tri_integral *pti,
    double a_in,
    class_complex y1,
    class_complex y2,
    class_complex *out
    ){
class_complex y1_re, y1_im, y2_re, y2_im, m_y1_re, m_y2_re, a, temp_ay1y2, result;
double CHOP_TOL;

CHOP_TOL = 1.e-14;

a = class_complex(a_in, 0.);

y1_re = class_complex(creal(y1),    0.);
y1_im = class_complex(0,            cimag(y1));
y2_re = class_complex(creal(y2),    0.);
y2_im = class_complex(0.,           cimag(y2));

// the m_ is for "minus ..."
// sqrt(-real - I*0.) = - I*sqrt(real)          -> incorrect since imag part is exactly zero
// sqrt(-real + I*0.) = + I*sqrt(real)          -> correct
m_y1_re = class_complex(-creal(y1),    0.);
m_y2_re = class_complex(-creal(y2),    0.);


// compute s(y1, y2) = sqrt(-y1)*sqrt(-y2)/sqrt(a*y1*y2)
// special cases if the real or imaginary part of y1 and/or y2 is close to zero

// The first case requires a special treatment since both imaginary parts are zero
// It is important that we don't calculate sqrt(real_number - I*0.) as explained above
if (cabs(y2_im) < CHOP_TOL && cabs(y1_im) < CHOP_TOL){
if (cabs(y1_re) >= CHOP_TOL && cabs(y2_re) >= CHOP_TOL){
temp_ay1y2 = class_complex(a_in*creal(y1)*creal(y2), 0.);
result = sqrt(m_y1_re)*sqrt(m_y2_re)/sqrt(temp_ay1y2);
}
if (cabs(y1_re) < CHOP_TOL && cabs(y2_re) >= CHOP_TOL){
temp_ay1y2 = class_complex(-a_in*creal(y2), 0.);
result = sqrt(m_y2_re)/sqrt(temp_ay1y2);
}
if (cabs(y1_re) >= CHOP_TOL && cabs(y2_re) < CHOP_TOL){
temp_ay1y2 = class_complex(-a_in*creal(y1), 0.);
result = sqrt(m_y1_re)/sqrt(temp_ay1y2);
}
if (cabs(y1_re) < CHOP_TOL && cabs(y2_re) < CHOP_TOL){
temp_ay1y2 = class_complex(a_in, 0.);
result = 1./sqrt(temp_ay1y2);
}
} else if (cabs(y2_im) >= CHOP_TOL && cabs(y1_im) < CHOP_TOL){
if (cabs(y1_re) >= CHOP_TOL){
result = sqrt(m_y1_re)*sqrt(-y2)/sqrt(a*y1_re*y2);
}
if (cabs(y1_re) < CHOP_TOL){
result = sqrt(-y2)/sqrt(-a*y2);
}
} else if (cabs(y2_im) < CHOP_TOL && cabs(y1_im) >= CHOP_TOL){
if (cabs(y2_re) >= CHOP_TOL){
result = sqrt(-y1)*sqrt(m_y2_re)/sqrt(a*y1*y2_re);
}
if (cabs(y2_re) < CHOP_TOL){
result = sqrt(-y1)/sqrt(-a*y1);
}
} else {
// case where cabs(y2_im) >= CHOP_TOL && cabs(y1_im) >= CHOP_TOL
result = sqrt(-y1)*sqrt(-y2)/sqrt(a*y1*y2);
}

*out = result;

return _SUCCESS_;
}


/** Calculate F_int
* @param R2                    Input: negative double
* @param y1                    Input: complex double z_+
* @param y2                    Input: complex double z_-
* @param x0                    Input: complex double
* @param out                   Output: pointer to output value
* @return the error status
*/

int util_F_int(
struct gen_tri_integral *pti,
double R2,
class_complex y1,
class_complex y2,
class_complex x0,
class_complex *out
){
class_complex cutx0, cut, derivx0, deriv, x0_re_c, atan_arg, prefactor, antideriv0, antideriv1;
double CHOP_TOL, y1_re, y1_im, y2_re, y2_im, x0_re, x0_im, y1_re2, y1_im2, y2_re2, y2_im2, x0_re2, x0_im2, a, b, c;
double xcross_temp, xsol[2], xbranch[2]; // at most two branch cuts, so length of list is two
int idx, signx0, sign, n_branchpoints;

// just for debugging
double xsol_temp[2], scale=1.;


CHOP_TOL = 1.e-10;

// printf("\nCalling util_F_int");

// contribution of a branch cut crossing at the boundary <=> B=0 (only pi/2 gap)
cutx0   = class_complex(0., 0.);

// regular contribution from branch cut crossing (pi gap)
cut     = class_complex(0., 0.);

// initialize values for the lists that information on the branch cuts.
// branch cut crossings can only occur between 0 and 1, so the value -1 indicates no branch cut crossing
for (idx = 0; idx<2; idx++){
xsol_temp[idx] = 0;
xsol[idx] = -1.;
xbranch[idx] = -1.;
}

// initialize number of branch cut crossings
n_branchpoints = 0;

// split y1, y2, and x0 into real and imag parts
y1_re = creal(y1);
y1_im = cimag(y1);
y2_re = creal(y2);
y2_im = cimag(y2);
x0_re = creal(x0);
x0_im = cimag(x0);

// squared
y1_re2 = y1_re*y1_re;
y1_im2 = y1_im*y1_im;
y2_re2 = y2_re*y2_re;
y2_im2 = y2_im*y2_im;
x0_re2 = x0_re*x0_re;
x0_im2 = x0_im*x0_im;


// create special cases if the imaginary part is close to zero
if (abs(cimag(y1)) < CHOP_TOL*fmax(1.0, cabs(y1))){
y1 = class_complex(creal(y1), 0.);
}
if (abs(cimag(y2)) < CHOP_TOL*fmax(1.0, cabs(y2))){
y2 = class_complex(creal(y2), 0.);
}
if (abs(cimag(x0)) < CHOP_TOL*fmax(1.0, cabs(x0))){
x0 = class_complex(creal(x0), 0.);
}

a = y1_im*x0_re-y2_im*x0_re-x0_im*y1_re+y2_im*y1_re+x0_im*y2_re-y1_im*y2_re;
b = -x0_im2*y1_im + x0_im*y1_im2+x0_im2*y2_im-y1_im2*y2_im-x0_im*y2_im2+y1_im*y2_im2-y1_im*x0_re2+y2_im*x0_re2+x0_im*y1_re2-y2_im*y1_re2-x0_im*y2_re2+y1_im*y2_re2;
c = y1_im2*y2_im*x0_re - y1_im*y2_im2*x0_re-x0_im2*y2_im*y1_re + x0_im*y2_im2*y1_re-y2_im*x0_re2*y1_re + y2_im*x0_re*y1_re2 + x0_im2*y1_im*y2_re-x0_im*y1_im2*y2_re+y1_im*x0_re2*y2_re-x0_im*y1_re2*y2_re-y1_im*x0_re*y2_re2+x0_im*y1_re*y2_re2;


// case B=0 <=> 0 < x0_re < 1 and x0_im=0 <=> branch cut crossing at the boundary of the branch cut
// determine the sign of the extra contribution by taking the derivative of the argument of arctan
// w.r.t. x
if (CHOP_TOL < x0_re && x0_re < 1.-CHOP_TOL && abs(x0_im)<CHOP_TOL*fmax(1.0, cabs(x0))){
// printf("\nBranch cut crossing at B=0!!!");
x0_re_c = class_complex(x0_re, 0.);
// if (cabs(y1-x0_re_c)<CHOP_TOL || cabs(x0_re_c-y2)<CHOP_TOL){
//   printf("\nERROR in util_F_int: division by zero. either x0_re_c=y1 or x0_re_c=y2.\n");
//}
derivx0 = (y1 - y2)/2./sqrt(-(x0_re_c-y1)*(x0_re_c-y1))/(x0_re_c-y2);
if (creal(derivx0) < 0){
signx0 = 1;
} else {
signx0 = -1;
}
cutx0 = signx0*_PI_/(sqrt(-x0_re_c+y1)*sqrt(x0_re_c-y2));
}

// other possible branch cut crossings:
// x_1,2 = 1/(2a) * (-b +- sqrt( b^2 -4ac ) ) -> needs to be between 0 and 1

// special case if a = 0
if (abs(a)<CHOP_TOL){
if (abs(b)>CHOP_TOL){
if ( 0<-c/b && -c/b<1){
 xsol[0] = -c/b;
 xsol_temp[0] = -c/b;
}
}
} else {
// check that the sqrt is not imagniary
if (b*b-4*a*c > CHOP_TOL){
xcross_temp = (-b + sqrt(b*b-4*a*c))/(2.*a);           // first potential branch cut crossing
xsol_temp[0] = xcross_temp;
// check that the solutions are between 0 and 1
// and that x≠x0 since this corresponds to B=0
if ((CHOP_TOL < xcross_temp && xcross_temp < 1.-CHOP_TOL) && !(abs(xcross_temp-x0_re) < CHOP_TOL && x0_im < CHOP_TOL)){
 xsol[0] = xcross_temp;
}
xcross_temp = (-b - sqrt(b*b-4*a*c))/(2.*a);           // second potential branch cut crossing
xsol_temp[1] = xcross_temp;
if ((CHOP_TOL < xcross_temp && xcross_temp < 1.-CHOP_TOL) && !(abs(xcross_temp-x0_re) < CHOP_TOL && x0_im < CHOP_TOL)){
 xsol[1] = xcross_temp;
}
}
}
for (idx=0; idx<2; idx++){
if (xsol[idx] != -1.){
// if (cabs(y1-x0)<CHOP_TOL || cabs(x0-y2)<CHOP_TOL){
//     printf("\nERROR in util_F_int: division by zero. either x0=y1 or x0=y2\n you have y1-x0=%.10f + %.10fj, y2-x0=%.10f + %.10fj, cabs(y1-x0)=%.10f , cabs(y2-x0)=%.10f , CHOP_TOL = %.10f.\n", creal(x0-y1), cimag(x0-y1), creal(x0-y2), cimag(x0-y2), cabs(x0-y1), cabs(x0-y2), CHOP_TOL);
// }
atan_arg = sqrt(xsol[idx]-y1)*sqrt(x0-y2)/(sqrt(-x0+y1)*sqrt(xsol[idx]-y2));
// critereon for branch cut crossing: argument of atan lies on the imaginary axis either between -i*infty to -i or between i to i*infty
scale = fmax(1., cabs(y1)*cabs(y1));
scale = fmax(scale, cabs(y2)*cabs(y2));
scale = fmax(scale, cabs(x0)*cabs(x0));
if (cabs(atan_arg)>1. && abs(creal(atan_arg))<scale*CHOP_TOL){
 n_branchpoints += 1;
} else {
 // if the critereon is not met, set the entry of xsol back to -1 to indicate that there is not valid branch cut crossing
 xsol[idx] = -1;
}
}
}

// calculate the contributions from the branch cut crossings
// only relevant contribution if there is exactly one branch cut crossing
// printf("\nn_branchpoints=%d",n_branchpoints);
if (n_branchpoints==1){
// printf("\nbranch cut crossing!!!");
// printf("xsol[0] = %.10f, xsol[1] = %.10f", xsol[0], xsol[1]);
for (idx=0; idx<2; idx++){
if (xsol[idx] != -1.){
 //if (cabs(y1-xsol[idx])<CHOP_TOL || cabs(xsol[idx]-y2)<CHOP_TOL){
 //    printf("\nERROR in util_F_int: division by zero. either xsol[idx]=y1 or xsol[idx]=y2.\n");
 //}
 deriv = sqrt(x0-y2)/sqrt(-x0+y1)*(1./(2.*sqrt(xsol[idx]-y1)*sqrt(xsol[idx]-y2)) - sqrt(xsol[idx]-y1)/(2.*(xsol[idx]-y2)*sqrt(xsol[idx]-y2)));
}
}

// determine the sign of the extra contribution by taking the derivative of the argument of arctan
// w.r.t. x
if (creal(deriv)<0.){
sign = 1.;
} else {
sign = -1.;
}
//if (cabs(y1-x0)<CHOP_TOL || cabs(x0-y2)<CHOP_TOL){
//    printf("\nERROR in util_F_int: division by zero. either x0=y1 or x0=y2.\n");
//}
cut = sign*_PI_*2./(sqrt(-x0+y1)*sqrt(x0-y2));

// check for mistakes
if (xsol[0] == -1. && xsol[1] == -1.){
printf("\n ERROR: number of branch points is supposed to be equal to one, but both entries of xsol are -1");
printf("\nxsol[0] = %.10f, xsol[1] = %.10f", xsol[0], xsol[1]);
}
} else {
// I think unnecessary since we initialized it as zero
cut = class_complex(0., 0.);
}

class_call(util_prefactor(pti, 
               R2,
               y1,
               y2,
               &prefactor),
 pti->error_message,
 pti->error_message);

class_call(util_antideriv(pti,
               0,
               y1,
               y2,
               x0,
               &antideriv0),
 pti->error_message,
 pti->error_message);

class_call(util_antideriv(pti,
               1,
               y1,
               y2,
               x0,
               &antideriv1),
 pti->error_message,
 pti->error_message);


*out = prefactor * _PI_*_PI_/2. * (cut + cutx0 + antideriv1 - antideriv0);


// printf("\n\nAll relevant variables:\n");
//printf("\ncut = %f + %fj\ncutx0 = %f + %fj\nantideriv1 = %f + %fj\nantideriv0 = %f + %fj", creal(cut), cimag(cut), creal(cutx0), cimag(cutx0), creal(antideriv1), cimag(antideriv1), creal(antideriv0), cimag(antideriv0));
//for (idx=0; idx<2; idx++){
//    if (xsol[idx] != -1.){
//        printf("\nsqrt(x0-y2) = %f + %fj", creal(sqrt(x0-y2)), cimag(sqrt(x0-y2)));
//        printf("\nsqrt(-x0+y1) = %f + %fj", creal(sqrt(-x0+y1)), cimag(sqrt(-x0+y1)));
//        printf("\nsqrt(x-y1) = %f + %fj", creal(sqrt(xsol[idx]-y1)), cimag(sqrt(xsol[idx]-y1)));
//        printf("\nsqrt(x-y2) = %f + %fj\n", creal(sqrt(xsol[idx]-y2) ), cimag(sqrt(xsol[idx]-y2) ));
//    }
//}
//printf("\nderiv = %f + %fj", creal(deriv), cimag(deriv));
//printf("\nx0 = %f + %fj", creal(x0), cimag(x0));
//printf("\ny1 = %f + %fj", creal(y1), cimag(y1));
//printf("\ny2 = %f + %fj", creal(y2), cimag(y2));
//printf("\nx0_new = %f + %fj", creal(x0_new ), cimag(x0_new ));
//printf("\ny1_new = %f + %fj", creal(y1_new ), cimag(y1_new ));
//printf("\ny2_new = %f + %fj", creal(y2_new ), cimag(y2_new ));
//printf("\nxsol_new = %f + %fj", creal(xsol_new[0] ), cimag(xsol_new[0] ));
//printf("\nxsol_new = %f + %fj", creal(xsol_new[1] ), cimag(xsol_new[1] ));
//printf("\nR2 = %f", R2);
//printf("\na = %f", a);
//printf("\nb = %f", b);
//printf("\nc = %f", c);
//printf("\nb*b-4*a*c = %f", b*b-4*a*c);
//printf("\nxsol0 = %f", xsol[0]);
//printf("\nxsol1 = %f", xsol[1]);
//printf("\nxsol_temp = %f", xsol_temp[0]);
//printf("\nxsol_temp = %f", xsol_temp[1]);
//printf("\nn_branchpoints = %d", n_branchpoints);
//printf("\nn_branchpoints_new = %d", n_branchpoints_new);
//printf("\nangle_x0 = %f, angle_xsol = %f, angle_y2 = %f", angle_x0, angle_xsol, angle_y2);
//printf("scale = %f, scale*CHOP_TOL = %f", scale, scale*CHOP_TOL);
//printf("\ncabs(atan_arg) = %0.20f, creal(atan_arg) = %0.20f", cabs(atan_arg), creal(atan_arg));
//if (0. < x0_re && x0_re < 1. && abs(x0_im)<CHOP_TOL*fmax(1.0, cabs(x0))){
//    printf("\nderivx0 = %f", derivx0);
//}
//printf("\nprefactor = %f + %fj\n", creal(prefactor), cimag(prefactor));
//printf("\nout = %f + %fj\n", creal(*out), cimag(*out));

return _SUCCESS_;
}


/** Calculate the T_master contributions
* @param y                     Input: evaluation of the integral at 0 and 1
* @param k12                   Input: wavenumber k1^2
* @param k22                   Input: wavenumber k2^2
* @param k32                   Input: wavenumber k3^2
* @param M1                    Input: complex mass
* @param M1                    Input: complex mass
* @param M1                    Input: complex mass
* @param out                   Output: pointer to output value
* @return the error status
*/
int util_Tmaster_contr(
        struct gen_tri_integral *pti,
        int y,
        double k12,
        double k22,
        double k32,
        class_complex M1,
        class_complex M2,
        class_complex M3,
        class_complex *out
        ){
class_complex Num0, R1, R0, S1, S1_overk22, S0, DiscS_sqrt, DiscR_sqrt, solS1, solS2, solR1, solR2, c1, c2, Fint_temp1, Fint_temp2;
double Num1, R2, S2, k1, k2, k3, k14, k24, k34, CHOP_TOL;

// antiderivative of the master integral integrand in y 

CHOP_TOL = 1.e-14;

k1 = sqrt(k12);
k2 = sqrt(k22);
k3 = sqrt(k32);
k14 = k12*k12;
k24 = k22*k22;
k34 = k32*k32;

// ---------------- Eq. 5.79 ----------------
// implement these formulas separately to avoid minimize erros
if (y==0){
R2 = -k32;
R1 = M1-M3+k32;
R0 = M3;

Num1 = 2.*(k12-k22-k32);
Num0 = 2.*(M2-M3+k22);
} else if (y==1){
R2 = -k12;
R1 = M1-M2+k12;
R0 = M2;

Num1 = 2.*(k12+k22-k32);
Num0 = 2.*(M2-M3-k22);
} else {
printf("Error in util_Tmaster_contr. Case not considered. y=%d",y);
}


// S_2, S_1, S_0
S2 = -k14+2.*k12*k22+2.*k12*k32-k24+2.*k22*k32-k34;
S1 = -4.*M1*k22-2.*M2*k12+2.*M2*k22+2.*M2*k32+2.*M3*k12+2.*M3*k22-2.*M3*k32-2.*k12*k22+2.*k24-2.*k22*k32;
S0 = -M2*M2+2.*M2*M3-2.*M2*k22-M3*M3-2.*M3*k22-k24;

// Diakr(a,b,c)= b^2 - 4*a*c -> defined after ---------------- Eq. 5.94 ----------------
// Sqrt of Diakr
DiscS_sqrt = sqrt(S1*S1 - 4.*S2*S0);

if (cabs(DiscS_sqrt)<CHOP_TOL){
printf("\nERROR in util_Tmaster_contr: DiscS_sqrt=%.10f + %.10fj\n", creal(DiscS_sqrt), cimag(DiscS_sqrt));
}
if (abs(S2)<CHOP_TOL){
printf("\nERROR in util_Tmaster_contr: S2 = 4*k12*k22*(cos12^2 - 1)=%.10f\n", S2);
printf("\nk12 = %.10f, k22 = %.10f, cos12 = %.10f", k12, k22, (k32-k12-k12)/(2.*sqrt(k12)*sqrt(k22)));
}

//a*x^2 + b*x + c = 0 has the solutions
// x_1,2 = (-b +- sqrt(b^2 - 4*a*c))/(2*a)

// ---------------- Eq. 5.80 ----------------
// solS1 and solS2 are the solution of 
// DeltaS2*x^2 + DeltaS1*x + DeltaS0 = 0
solS1 = 0.5*(-S1 + DiscS_sqrt)/S2;              // x_+
solS2 = 0.5*(-S1 - DiscS_sqrt)/S2;              // x_-

c2 = -(Num1*solS2 + Num0)/DiscS_sqrt;           // c1
c1 =  (Num1*solS1 + Num0)/DiscS_sqrt;           // c2

DiscR_sqrt = sqrt(R1*R1 - 4.*R2*R0);            // sqrt( DeltaR1^2 - 4*DeltaR2*DeltaR0 )

solR1 = 0.5*(-R1+DiscR_sqrt)/R2;                // z_+
solR2 = 0.5*(-R1-DiscR_sqrt)/R2;                // z_-

if (cabs(DiscR_sqrt)<CHOP_TOL){
printf("\nERROR in util_Tmaster_contr: DiscR_sqrt=0\n");
}
if (abs(R2)<CHOP_TOL){
printf("\nERROR in util_Tmaster_contr: R2=0\n");
}



// ---------------- Eq. 5.82 ----------------
if (cabs(c1) < CHOP_TOL){
class_call(util_F_int(pti,
               R2,
               solR1,
               solR2,
               solS2,
               &Fint_temp2),
     pti->error_message,
     pti->error_message);

*out = c2*Fint_temp2;
return _SUCCESS_;

} else if (cabs(c2) < CHOP_TOL){
class_call(util_F_int(pti,
               R2,
               solR1,
               solR2,
               solS1,
               &Fint_temp1),
     pti->error_message,
     pti->error_message);

*out = c1*Fint_temp1;
return _SUCCESS_;

} else {
class_call(util_F_int(pti,
               R2,
               solR1,
               solR2,
               solS1,
               &Fint_temp1),
     pti->error_message,
     pti->error_message);

class_call(util_F_int(pti,
               R2,
               solR1,
               solR2,
               solS2,
               &Fint_temp2),
     pti->error_message,
     pti->error_message);
     

*out = c1*Fint_temp1 + c2*Fint_temp2;
return _SUCCESS_;

}
}