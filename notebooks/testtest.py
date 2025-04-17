import numpy as np
from gmpy2 import atan, log, sqrt, gamma


def almosteq(a, b, tol):
    return abs(a-b)<tol

def Antideriv(x, y1, y2, x0, CHOP_TOL):
    # calculates the antiderivative (arctan times a prefactor)
    
    # case where x0 = y2 = 0 or 1
    if almosteq(x0,y2,CHOP_TOL):
        if almosteq(x,y2,CHOP_TOL):
            # no idea why return zero
            return 0
        return 2.*sqrt(x-y1)/(-x0+y1)/sqrt(x-y2)
    
    if abs(x0-y1) < CHOP_TOL:
        print('WARNING: switching var in Antideriv')
        #x0 = y2 = 0 or 1
        return Antideriv(x,y2,y1,x0)

    prefac = 2/(sqrt(-x0+y1)*sqrt(x0-y2))
    temp = sqrt(x-y1)*sqrt(x0-y2)/sqrt(-x0+y1)
    LimArcTan = 0.+1j*0.

    # print('temp', temp)
    if x == 1 and almosteq(1.+1j*0., y2, CHOP_TOL):
        LimArcTan = 1j * sqrt(-temp**2) * np.pi/(2*temp)
        return  prefac * LimArcTan
    if x == 0 and almosteq(0.+1j*0., y2, CHOP_TOL):
        LimArcTan = sqrt(temp**2) * np.pi/(2*temp)
        return  prefac * LimArcTan

    return prefac*atan(temp/sqrt(x-y2))

def Prefactor(a_in, y1, y2, CHOP_TOL):
    #computes prefactor that shows up in Fint
    y2re = y2.real
    y1re = y1.real

    a = a_in


    if abs(y2.imag) < CHOP_TOL and abs(y1.imag) < CHOP_TOL:
        if abs(y1re) >= CHOP_TOL and abs(y2re) >= CHOP_TOL:
            print("\nTHIS IS THE CASE\n")
            print(f"-y1 = {-y1re + 1j*0.}")
            print(f"-y2 = {-y2re + 1j*0.}")
            print(f"sqrt(-y1re + 1j*0.) = {sqrt(-y1re + 1j*0.)}")
            print(f"sqrt(-y2re + 1j*0.) = {sqrt(-y2re + 1j*0.)}")
            print(f"sqrt(a*(y1re)*(y2re) + 1j*0.) = {sqrt(a*(y1re)*(y2re) + 1j*0.)}")
            return sqrt(-y1re + 1j*0.)*sqrt(-y2re + 1j*0.)/(sqrt(a*(y1re)*(y2re) + 1j*0.))
        if abs(y1re) < CHOP_TOL and abs(y2re) >= CHOP_TOL:
            return sqrt(-y2re + 1j*0.)/sqrt(-a*y2re + 1j*0.)
        if abs(y1re) >= CHOP_TOL and abs(y2re) < CHOP_TOL:
            return sqrt(-y1re + 1j*0.)/sqrt(-a*(y1re) + 1j*0.)
        if abs(y1re) < CHOP_TOL and abs(y2re) < CHOP_TOL:
            return 1/sqrt(a + 0.*1j)

    elif abs(y2.imag) >= CHOP_TOL and abs(y1.imag) < CHOP_TOL:
        if abs(y1re) >= CHOP_TOL: 
            return sqrt(-y1re + 1j*0.)*sqrt(-y2)/sqrt(a*y1re*y2)
        if abs(y1re) < CHOP_TOL:
            return sqrt(-y2)/sqrt(-a*y2)

    elif abs(y2.imag) < CHOP_TOL and abs(y1.imag) >= CHOP_TOL:
        if abs(y2re) > CHOP_TOL: 
            return sqrt(-y1)*sqrt(-y2re + 1j*0.)/(sqrt(a*y1*y2re))
        if abs(y2re) < CHOP_TOL:
            return sqrt(-y1)/sqrt(-a*y1)
    else:
        # case where abs(y2.imag) >= CHOP_TOL and abs(y1.imag) >= CHOP_TOL
        return sqrt(-y1)*sqrt(-y2)/sqrt(a*y1*y2)

# calculates Fint. Most of the code of to find and count the branch cuts
    
    # aa = DeltaR2
def Fint(aa, y1, y2, x0, CHOP_TOL):

    # this is necessary because in gmpy2 sqrt(-1-0j) = -1j and not 1j
    if abs(y2.imag) < CHOP_TOL:            # z_- is real
        y2 = y2.real + 1j*0.
    if abs(y1.imag) < CHOP_TOL:            # z_+ is real
        y1 = y1.real + 1j*0.
    if abs(x0.imag) < CHOP_TOL:            # x0 is real
        x0 = x0.real + 1j*0.


    # split y1=z_+, y2=z_1, x0 up into real an imaginary parts

    # z_+
    rey1 = y1.real
    imy1 = y1.imag

    # z_-
    rey2 = y2.real
    imy2 = y2.imag

    # x0
    rex0 = x0.real
    imx0 = x0.imag

    # initialize relevant qunatities
    # numbranchpoints: number of actual branch cut crossings
    # signx0:             
    # sign:             
    numbranchpoints = 0
    signx0 = 0
    sign = 0
    # xsol:             
    #    - description:     probably stores the x values for which there is a branch cut crossing, recall: x = Re(z), z: integration variable from 0 to 1
    #    - length:         2 (at most I think -> maye in C make it fixed size two and set one value to -1 if there isnt a branch cut crossing)
    # xbranch:    
    #     - description:     list of branch cut crossing points
    #    - length:         2 (at most two, could also be zero or one, but in C maybe fix it to two and set value to -1 if no branch cut crossing)        
    # atanarglist:
    #    - description:     argument of arc tan for all values in xsol
    #    - length:         len(xsol)
    # abscrit:
    #    - description:     abs value of arctanarg
    #    - length:         len(atanarglist) = len(xsol)
    # recrit: 
    #    - description:    real part of arctanarg
    #    - length:         len(atanarglist) = len(xsol)        
    # derivcrit:        
    xsol = []
    xbranch = []
    atanarglist = []
    abscrit = []
    recrit = []
    derivcrit = []


    # ---------------- Eq. 5.92 ----------------
    c = imy1**2*imy2*rex0 - imy1*imy2**2*rex0-imx0**2*imy2*rey1 + imx0*imy2**2*rey1-imy2*rex0**2*rey1 + imy2*rex0*rey1**2 + imx0**2*imy1*rey2-imx0*imy1**2*rey2+imy1*rex0**2*rey2-imx0*rey1**2*rey2-imy1*rex0*rey2**2+imx0*rey1*rey2**2
    a = imy1*rex0-imy2*rex0-imx0*rey1+imy2*rey1+imx0*rey2-imy1*rey2
    b = -imx0**2*imy1 + imx0*imy1**2+imx0**2*imy2-imy1**2*imy2-imx0*imy2**2+imy1*imy2**2-imy1*rex0**2+imy2*rex0**2+imx0*rey1**2-imy2*rey1**2-imx0*rey2**2+imy1*rey2**2
    
    cutx0 = 0. + 0.*1j

    # ---------------- Top of p. 60 ----------------
    # relevant for case B=0
    if 0 < rex0 < 1 and abs(imx0) < CHOP_TOL:                                # x0 real and 0<x0<1         <=>        B = 0
        derivcritx0 = (y1 - y2)/2/sqrt(-(rex0-y1)**2)/(rex0-y2)                # dA/dx
        if derivcritx0.real < 0:                                            # dA/dx determines the sign of branch cut crossig
            signx0 = 1                                                        # + pi/2
        else:
            signx0 = -1                                                        # - pi/2
        cutx0 = signx0*np.pi/(sqrt(-rex0+y1+0j)*sqrt(rex0-y2+0j))
    else:
        cutx0 = 0. + 0.*1j                                                        # no discontinuity
    
    if abs(a) < CHOP_TOL:                     # |a| = 0
        if b != 0:
            xsol = [- c / b]                # add the only solution to the list, solution of 0=bx+c
        else:
            xsol = []                        # no solution 
    else:
        #print(f"\nb*b-4*a*c = { b*b-4*a*c }")
        #print(f"\n(-b + sqrt(b*b-4*a*c))/(2.*a) = { (-b + sqrt(b*b-4*a*c))/(2.*a) }")
        #print(f"\n(-b - sqrt(b*b-4*a*c))/(2.*a) = { (-b - sqrt(b*b-4*a*c))/(2.*a) }")
        #print(f"\na = { a }, b = { b }, c = { c }")
        if b**2-4*a*c > 0:                    # full solution -> a≠0
            xsol = [(-b + sqrt(b**2-4*a*c))/(2*a),(-b - sqrt(b**2-4*a*c))/(2*a)]
        else:
            #case where there is no intersection of the real axis (includes double zero)
            xsol = []


    # check two things:
    #     - the entries of xsol are within the interval [0, 1]
    #     - z≠x0 since this correspond to B=0 which only has a pi/2 contribution
    xsol = [x for x in xsol if x > CHOP_TOL and x < 1. - CHOP_TOL and not(almosteq(x,x0,CHOP_TOL))]
    

    # Now, we still have to check if branch cut crossing belong to the correct arc
    # ---------------- below Eq. 5.94 ----------------
    # ---------------- below Eq. 5.90 ----------------

    # cut is the extra term from the branch cut crossing
    cut = 0. + 1j*0.
    if len(xsol) > 0:

        # ---------------- below Eq. 5.88 ----------------
        # A = arctanarg
        ##### In the paper there is sqrt(y1-x)/sqrt(x0-y1), here is sqrt(x-y1)/sqrt(-x0+y1) -> Is that problematic?
        atanarglist = [sqrt(x-y1)*sqrt(x0-y2)/(sqrt(-x0+y1)*sqrt(x-y2)) for x in xsol]

        # recall: 
        # branch cut for A from -i*infty to -i and from i to infty i
        # So at A.real = 0 and abs(A)≥1 (Note: abs(A)=1 is numerically impossible)


        # absolute values of the argument of arctan
        # abs(A)
        abscrit = [abs(atanarg) for atanarg in atanarglist]

        # real part of argument of arctan
        # A.real
        recrit = [atanarg.real for atanarg in atanarglist]
    
        for i in range(len(xsol)):
            # check conditions for branch cut crossing:
            # abs(A)>1 and A.real=0
            if abscrit[i] > 1 and abs(recrit[i])<CHOP_TOL:
                numbranchpoints += 1
                xbranch.append(xsol[i])

    # only  if there is one branch cut crossing we have to do something
    if numbranchpoints == 1:
        # determine the sign of the contribution through dA/dx
        # this implementation unnecessary since xbranch has length 1
        derivcrit = [sqrt(x0-y2)/sqrt(-x0+y1)*(1/(2*sqrt(x-y1)*sqrt(x-y2)) -sqrt(x-y1)/(2*(x-y2)*sqrt(x-y2))) for x in xbranch]
        if derivcrit[0].real < 0:
            sign = 1
        else:
            sign = -1
        cut = sign*np.pi*2/(sqrt(-x0+y1)*sqrt(x0-y2))
    # no additional contributions for zero or two branch cut crossings (aussuming B≠0)
    else:
        cut = 0. + 1j*0.

    prefac0 = Prefactor(aa,y1,y2, CHOP_TOL)


    print(f"Babis: prefac: aa={aa}, y1={y1}, y2={y2}, CHOP_TOL={CHOP_TOL}")

    result = prefac0*(np.pi**2/2.)*(cut + cutx0 + Antideriv(1,y1,y2,x0,CHOP_TOL) - Antideriv(0,y1,y2,x0,CHOP_TOL))

    # print(f"\nBabis: prefac = {prefac0}, cut = {cut}, cutx0 = {cutx0}, Antideriv1 = {Antideriv(1,y1,y2,x0,CHOP_TOL)}, Antideriv0 = {Antideriv(0,y1,y2,x0, CHOP_TOL)}")

    return result

def Diakr(a, b, c):
	return b**2-4*a*c

def TrMxy(y, k21, k22, k23, M1, M2, M3, CHOP_TOL):
	# antiderivative of the master integral integrand in y 


	# ---------------- Eq. 5.79 ----------------
	# N_1, N_0
	Num1 = 4*k22*y+2*k21-2*k22-2*k23
	Num0 = -4*k22*y+2*M2-2*M3+2*k22

	# R_2, R_1, R_0
	DeltaR2 = -k21*y+k23*y-k23
	DeltaR1 = -M2*y+M3*y+k21*y-k23*y+M1-M3+k23
	DeltaR0 = M2*y-M3*y+M3

	# S_2, S_1, S_0
	DeltaS2 = -k21**2+2*k21*k22+2*k21*k23-k22**2+2*k22*k23-k23**2
	DeltaS1 =-4*M1*k22-2*M2*k21+2*M2*k22+2*M2*k23+2*M3*k21+2*M3*k22-2*M3*k23-2*k21*k22+2*k22**2-2*k22*k23
	DeltaS0 =-M2**2+2*M2*M3-2*M2*k22-M3**2-2*M3*k22-k22**2

	# Diakr(a,b,c)= b^2 - 4*a*c -> defined after ---------------- Eq. 5.94 ----------------
	# Sqrt of Diakr

	DiakrS = sqrt(Diakr(DeltaS2, DeltaS1, DeltaS0)) 	# sqrt( DeltaS1^2 - 4*DeltaS2*DeltaS0 )

	#a*x^2 + b*x + c = 0 has the solutions
	# x_1,2 = (-b +- sqrt(b^2 - 4*a*c))/(2*a)

	# ---------------- Eq. 5.80 ----------------
	# solS1 and solS2 are the solution of 
	# DeltaS2*x^2 + DeltaS1*x + DeltaS0 = 0
	solS1 = (-DeltaS1+DiakrS)/2/DeltaS2				# x_+
	solS2 = (-DeltaS1-DiakrS)/2/DeltaS2  				# x_-

	cf2 = -(Num1*solS2+Num0)/DiakrS					# c1
	cf1 = (Num1*solS1+Num0)/DiakrS						# c2
		
	DiakrR = sqrt(Diakr(DeltaR2, DeltaR1, DeltaR0))	# sqrt( DeltaR1^2 - 4*DeltaR2*DeltaR0 )
				  
	solR1 = ((-DeltaR1+DiakrR)/2)/DeltaR2     			# z_+
	solR2 = ((-DeltaR1-DiakrR)/2)/DeltaR2 				# z_-

	# ---------------- Eq. 5.82 ----------------
	if abs(cf1) < CHOP_TOL:										# |c1| = 0 -> note that c1 and c2 are complex
		# neglect cf1
		return cf2*Fint(DeltaR2, solR1, solR2, solS2, CHOP_TOL)
	elif abs(cf2) < CHOP_TOL:									# |c2| = 0
		# neglect cf2
		return cf1*Fint(DeltaR2, solR1, solR2, solS1, CHOP_TOL)
	else:
		return cf2*Fint(DeltaR2, solR1, solR2, solS2, CHOP_TOL)+cf1*Fint(DeltaR2, solR1, solR2, solS1, CHOP_TOL)


def TriaMasterZeroMasses(k21, k22, k23):
	#case for triangle master integrals where all masses vanish
	return np.pi**3/sqrt(k21)/sqrt(k22)/sqrt(k23)


def TriaMaster(k21, k22, k23, M1, M2, M3, CHOP_TOL):
	#--- masses are squared as in the paper ---
	# calculates triangle master integral
	if M1 == 0.+1j*0. and M2 == 0.+1j*0. and M3 == 0.+1j*0.:
		return  TriaMasterZeroMasses(k21, k22, k23)
	
	return TrMxy(1, k21, k22, k23, M1, M2, M3, CHOP_TOL)-TrMxy(0, k21, k22, k23,M1, M2, M3, CHOP_TOL)

