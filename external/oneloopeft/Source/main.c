

/** @file main.c Documented main module, including functions to initilize and cleanup the cosmology structure
 *  and examples of calls to functions in other modules to compute the line clustering and shot power spectrum. 
 *  
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 *  In order to call any function from the package, the function calls should be placed in the marked section of main() function.
 */


#include "header.h"

struct globals gb;


int main(int argc,char *argv[])
{	
	extern struct globals gb;
 	
	clock_t t_start = clock();

 	initialize();

	long fish_dim;
	long int i,j,l, q, s,t;

	fish_dim = gb.Npars;	

	double pk_zmax = 14.;
	double pk_kmax = 200.;
	double M_min   = 1.e9;
	long mode_mf   = ST;

	/*
	* Define which lines to compute, and how many points to have for z-interpolation
	*/
	int nlines     = 0;
      int lines[7]   = {CO10, CO21, CO32, CO43, CO54, CO65, CII};
      int JJ[7]      = {1,2,3,4,5,6,0};
      size_t ninterp = 150;

	double gb_pars[] = {gb.logAs,gb.ns,gb.h,gb.Omega_b,gb.Omega_cdm,gb.sigFOG0};

	struct Cosmology Cx_ref;
	for(i=0;i<gb.Npars;i++){
		Cx_ref.cosmo_pars[i] = gb_pars[i];
	}

	/*
	* Initialize the cosmology structure, which includes CLASS cosmology and Line structures
	*/
	clock_t tic_r = clock();
	CL_Cosmology_initialize(&Cx_ref,pk_kmax,pk_zmax);//(&Cx_ref, pk_kmax, pk_zmax, nlines, lines, ninterp, M_min, mode_mf);
	printf("Reference Cosmology initialized\n");
	clock_t toc_r = clock();
	printf("Elapsed: %f seconds for ref cosmology\n", (double)(toc_r - tic_r)/CLOCKS_PER_SEC);

	
	/*
	* -----------------------------
	* Depending on what quantities ou want to compute, you need to modify this part of the main() function.
	* -----------------------------
	*/

	/* 
	 * Set the k and Pk arrays and setting the redshift
	 */
 	double z = 0.0;
    int nk = 200;
	double *k  = loginit_1Darray(nk, 1.e-3, 8.);
	double *PK = make_1Darray(nk);

	for(i=0; i<nk; i++){
		PK[i] = PS_hh_G(&Cx_ref, k[i], z, z, LOOP, NOIR, DST, 0.0);
	}

	free(k);
 	free(PK);

	/* 
	* -----------------------------
	* Modify up to here 
	* -----------------------------
	*/

	cleanup(&Cx_ref);
	
	clock_t t_finish = clock();
	printf("	Elapsed: %f seconds for the calculation\n", (double)(t_finish - t_start)/CLOCKS_PER_SEC);
	return 0;
}


/**
 * Initialize the path to the required directories, set the values of cosmological parmaeters, and initialize the interpolator of the SFR(M,z) 
 * from tabulated data provided in gb.SFR_filename.  
 * 
 * The global structure "gb" have several elements to hold the paths to project source directory, input, and output folders, 
 * and values of cosmological parmaeters. 
 * 
 * @return void 
 */
void initialize()
{
	extern struct globals gb;

	getcwd(gb.project_home , sizeof(gb.project_home));  //This gives the path to the source directory
	chdir("..");  ///Change the path to the parent directory
	getcwd(gb.project_home , sizeof(gb.project_home));

	sprintf(gb.output_dir,"%s/Output", gb.project_home);
	sprintf(gb.data_dir,"%s/Input", gb.project_home);
	sprintf(gb.data_priors,"%s/Priors", gb.data_dir);

	sprintf(gb.SFR_filename,"%s/release-sfh_z0_z8_052913/sfr/sfr_release.dat",gb.data_dir);  

	//logSFR_alloc_init();

	///// cosmological parameters corresponding to initial conditions of HV simulations
	gb.c 		 = 2.99792458e5;  /// In units of km/s
	gb.h 		 = 0.677;
	gb.Omega_cdm = 0.11923/pow(gb.h,2.);  ////omega_cdm = Omega_cdm h^2 ;
	gb.Omega_b 	 = 0.02247/pow(gb.h,2.);   ///omega_b = Omega_b h^2;
	gb.Omega_r 	 = 0.0000910958;  ////radiation = photons + neutrinos
	gb.Omega_g 	 = 5.3956421715871286e-05;   ////photons, input for Class
	gb.Omega_nu  = 0.00;    ////neutrinos
	gb.ns        = 0.96824;
	gb.As 	 	 = 2.1085e-9;
	gb.logAs     = log(gb.As*pow(10.,10.));  ///3.0665
	gb.kp 		 = 0.05;  //in unit of Mpc^-1

	gb.sigFOG0 	 = 250.; 	

	gb.PS_kmin 	 = 1.e-4;
	gb.PS_kmax 	 = 200.;

	gb.Npars   	 = 6;
	return;
}


/**
 * Free the memory allocated to cosmology structure and SFR interpolator
 * 
 * @return void 
 */
void cleanup(struct Cosmology *Cx)
{
	extern struct globals gb;	

	Cosmology_free(Cx);

	return;
}
	
 	
	
