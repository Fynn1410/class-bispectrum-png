

/** @file ps_halo_1loop.c Documented real-space, direct integration computation of 1loop contributions of the halo/galaxy power spectrum
 * See arXiv:2010.14523 for explicit expressions
 *
 * Azadeh Moradinezhad Dizgah, November 4th 2021
 *
 * This module computes the 1loop halo/galaxy power spectrum in real-space via direct numerical integration.
 * IR-resummation and EFT counter terms are included. In addition to loops due to gravitational loops, terms arising only in the presence of local PNG are
 * also included. The explicit expressions of all the loops are given in 2010.14523.
 *
 * In summary, the following functions can be called from other modules:
 * -# PS_hh_G()
 * -# Compute_Gloops()
 * -# F2_s()
 * -# F3_s()
 * -# S2_s()
 * -# F2()
 * -# S2()
 */

#include "header.h"
#include <time.h>
//struct globals gb;

/**
 * Compute the contributions up to 1loop to halo power spectrum for Gaussian initial conditions
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber
 * @param z            Input: redshift of interest
 * @param M            Input: halo mass, used in computing the halo bias
 * @param mode_pt      Input: switch to decide whether to compute tree-level halo power spectrum or the 1loop
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param mode_mf      Input: switch to set the theoretical model of the mass function used to compute the halo biases
 * @return G loop contributions of P_h
 */

int PS_hh_G(
               struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               short has_loop,
               short has_ir,
               long SPLIT, 
               double * pk_nl)
{
      int cleanup = 0;

      double pm_lin = 0., pm_lin_IR = 0., pm_1loop_IR = 0., pm_22 = 0., pm_13 = 0., pm_1loop =0., pm_ct = 0., ph_tot = 0.;

      //DL
      //double *bias_arr = make_1Darray(4);
      //halo_bias(Cx, M, z, mode_mf, bias_arr);
      // double b1  = bias_arr[0];
      // double b2  = bias_arr[1];
      // double bG2 = bias_arr[2];
      // double btd = bias_arr[3];

      //DL CLASS-PT values page 30
      double b1  = 2.0;
      double b2  = -1.0;
      double bG2 = 0.1;
      double btd = -0.1;
      double R2  =  5.0;
      double cs2 = 0.2;

      //printf("Bias Vals: %12.6e %12.6e %12.6e %12.6e %12.6e  %12.6e \n",k,M,b1,b2,bG2,btd);

      clock_t start, matter, halo;
      start = clock();

      if (has_loop == _TRUE_) {

          double *ps_hloops, *ps_mloops;
          ps_hloops = make_1Darray(6);
          // ps_mloops = make_1Darray(2);

          if (has_ir == _TRUE_) {
            // Compute_G_loops(pba, ppm, pfo, k, z, _TRUE_, MATTER, SPLIT, ps_mloops);
            // matter = clock();
            // fprintf(stderr,"Matter Loops done! %f sec. needed.\n", (double)(matter-start) / CLOCKS_PER_SEC);
            Compute_G_loops(pba, ppm, pfo, k, z, _TRUE_, HALO, SPLIT, ps_hloops);
            halo = clock();
            // fprintf(stderr,"Matter Loops done! %f sec. needed.\n", (double)(halo-matter) / CLOCKS_PER_SEC);
          }
          else {
            Compute_G_loops(pba, ppm, pfo, k, z, NOIR, HALO, SPLIT, ps_hloops);
            // Compute_G_loops(pba, ppm, pfo, k, z, NOIR, MATTER, SPLIT, ps_mloops);
          }
          pm_lin_IR       = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
          double pb1b2    = b1 * b2  * ps_hloops[0];
          double pb1bg2   = 2. * b1 * bG2 * ps_hloops[1];
          double pb22     = 0.25 * pow(b2, 2.) * ps_hloops[2];
          double pbg22    = pow(bG2, 2.)  * ps_hloops[3];
          double pb2bg2   = b2 * bG2 * ps_hloops[4];
          double pb1b3nl  = 2. * b1 * (bG2 + 2./5. * btd) * ps_hloops[5];
          double pR2      = - 2. * b1 * R2 * pow(k, 2.) * pm_lin_IR;
          double ph_loops =  pb1b2 + pb1bg2 + pb22+ pbg22 + pb2bg2 + pb1b3nl + pR2;

          double khat     = 1. * pba->h;

          double p_mm;
          PS_mm_G(pba,ppm,pfo,k,z,has_loop,has_ir,SPLIT,&p_mm);

          // if (has_ir == _FALSE_) {
          //   pm_lin   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
          //   pm_22    = ps_mloops[0];
          //   pm_13    = ps_mloops[1];
          //   pm_1loop = pm_lin + pm_22 + pm_13;
          //   pm_ct    = - 2. * b1 * (R2 + cs2 * b1) * pow(k, 2.) * pm_lin;
          //   ph_tot   = pow(b1, 2.) * pm_1loop + pm_ct + ph_loops;
          // }
          // else {
          //   pm_lin      = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
          //   pm_22       = ps_mloops[0];
          //   pm_13       = ps_mloops[1];
          //   pm_lin_IR   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
          //   pm_1loop_IR = pm_IR_NLO(pba, ppm, pfo, k, z, SPLIT);
          //   pm_ct       = - 2. * b1 * (R2 + cs2 * b1) * pow(k, 2.) *  pm_lin_IR; //pow(k, 2.)/(1.+pow(k/khat,2.))*
          //   ph_tot      = pow(b1, 2.) * pm_1loop_IR + pm_ct + ph_loops;
          // }

          // fprintf(stderr,"Loops done in %f sec.\n", (double)(halo-start) / CLOCKS_PER_SEC);

          ph_tot = pow(b1, 2.) * p_mm  + ph_loops;

          FILE *fpa;
          char file_name[50];
          sprintf(file_name, "data/pg_DI.txt");
          fpa = fopen(file_name, "a");
          fprintf(fpa, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",\
              k, pm_lin_IR, pow(b1, 2.) * p_mm, pb1b2, pb1bg2, pb22, pbg22, pb2bg2, pb1b3nl, ph_loops, ph_tot);
          fclose(fpa);

          free(ps_hloops);
          free(ps_mloops);
      }
      else {
          ph_tot  = pow(b1, 2.) * Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
      }

      //free(bias_arr);
      *pk_nl = ph_tot;

      return _SUCCESS_;
}

/**
 * Compute one loop matter power spectrum for Gaussian initial conditions
 *
 * @param ppr          Input: pointer to precision structure
 * @param pba          Input: pointer to background structure
 * @param ppt          Input: pointer to perturbation structure
 * @param ppm          Input: pointer to primordial structure
 * @param pfo          Input: pointer to fourier structure
 * @param index_pk     Input: index for type of power spectrum (total matter m, baryon+CDM bc)
 * @param k            Input: wavenumber
 * @param z            Input: redshift
 * @param mpt          Input: flag to decide whether to compute tree-level halo power spectrum or the 1loop
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param pk_nl        Output: pointer to output non-linear P(k)
 * @return the error status
 */

int PS_mm_G(
            struct background * pba,
            struct primordial * ppm,
            struct fourier * pfo,
            double k,
            double z,
            short has_loop,
            short has_ir,
            long SPLIT,
            double * pk_nl) {

  double pm_lin = 0.;
  double pm_lin_IR = 0.;
  double pm_22 = 0.;
  double pm_13 = 0.;
  double ph_tot = 0.;
  double pm_1loop_IR = 0.;
  double pm_ct = 0.;
  double *ps_mloops;
  double khat;
  double cs2;

  if (has_loop == _TRUE_) {

    ps_mloops = make_1Darray(2);
    Compute_G_loops(pba, ppm, pfo, k, z, has_ir, MATTER, SPLIT, ps_mloops);

    khat = 1. * pba->h;
    cs2 = 0.2;

    pm_lin   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
    pm_22    = ps_mloops[0];
    pm_13    = ps_mloops[1];

    if (has_ir == _FALSE_){
      pm_ct    = - 2. * cs2 * pow(k, 2.) * pm_lin;
      ph_tot   = pm_lin + pm_22 + pm_13 + pm_ct;
      // fprintf(stderr,"%e %e %e %e %e %e\n",k, pm_lin, pm_13, pm_22, pm_ct, ph_tot);
    }
    else {
      double p_nowiggle = pm_nowiggle(pba, ppm, pfo, k, z, 0., 0, SPLIT);
      double p_wiggle   = pm_lin - p_nowiggle;
      double sigma2     = IR_Sigma2(pba, ppm, pfo, z, 0., SPLIT);
      double sup        = exp(-k * k * sigma2);

      pm_lin_IR   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      pm_1loop_IR = p_nowiggle + sup * p_wiggle * (1. + k * k * sigma2) + pm_22 + pm_13; 
      pm_ct       = - 2. * cs2 * pow(k, 2.) *  pm_lin_IR; // pow(k, 2.)/(1.+pow(k/khat,2.))
      ph_tot      = pm_1loop_IR + pm_ct;

      FILE *fpa;
      char file_name[50];
      sprintf(file_name, "data/pm_DI.txt");
      fpa = fopen(file_name, "a");
      fprintf(fpa,"%e %e %e %e %e %e\n",k , pm_lin_IR, pm_22, pm_13, pm_ct, ph_tot);
      fclose(fpa);
    }

    free(ps_mloops);
  }
  else {
    ph_tot  =  Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
  }

  //fprintf(stderr,"%e => %e\n",k, ph_tot);
  *pk_nl = ph_tot;

  return _SUCCESS_;
}

int PS_hh_1(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT)
{
      int cleanup = 0;

      double pm_lin = 0., pm_lin_IR = 0., pm_1loop_IR = 0., pm_22 = 0., pm_13 = 0., pm_1loop =0., pm_ct = 0., ph_tot = 0.;

      //DL CLASS-PT values page 30
      double b1  = 2.0;
      double b2  = -1.0;
      double bG2 = 0.1;
      double btd = -0.1;
      double R2  =  5.0;
      double cs2 = 0.2;

      clock_t start, matter, halo;
      start = clock();

      double *ps_hloops;
      ps_hloops = make_1Darray(10);

      Compute_1_loops(pba, ppm, pfo, k, z, mu, SPLIT, ps_hloops);
      halo = clock();
      // fprintf(stderr,"Halo Loops for the 1-st moment done! %f sec. needed.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      double I2201     = ps_hloops[0];
      double I1301p3101     = ps_hloops[1];
      double Idelta201 = ps_hloops[2];
      double IG201     = ps_hloops[3];
      double FG201     = ps_hloops[4];

      double J12101    = ps_hloops[5];
      double J11201    = ps_hloops[6];
      double J21101    = ps_hloops[7];
      double Jdelta201 = ps_hloops[8];
      double JG201     = ps_hloops[9];

      double ph_loops  = I2201 + I1301p3101 + Idelta201 + IG201 + FG201;
      double plos_loops  = J12101 + J11201 + J21101 + Jdelta201 + JG201;
      // pm_lin_IR =  pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      pm_lin_IR =  Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
      double khat     = 1. * pba->h;

      // fprintf(stderr,"Loops done in %f sec.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      FILE *fpa;
      char file_name[50];
      sprintf(file_name, "data/1_moment_DI_%g.txt",mu);
      fpa = fopen(file_name, "a");
      fprintf(fpa, "%12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e %12.12e\n",\
          k, pm_lin_IR, I2201, I1301p3101, Idelta201, IG201, FG201, J12101, J11201, J21101, Jdelta201, JG201, ph_loops, plos_loops);
      fclose(fpa);

      free(ps_hloops);
      //free(ps_mloops);
      //*pk_nl = ph_tot;

      return _SUCCESS_;
}

int PS_hh_2(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT)
{
      int cleanup = 0;

      double pm_lin = 0., pm_lin_IR = 0., pm_1loop_IR = 0., pm_22 = 0., pm_13 = 0., pm_1loop =0., pm_ct = 0., ph_tot = 0.;

      //DL CLASS-PT values page 30
      double b1  = 2.0;
      double b2  = -1.0;
      double bG2 = 0.1;
      double btd = -0.1;
      double R2  =  5.0;
      double cs2 = 0.2;

      clock_t start, matter, halo;
      start = clock();

      double *ps_hloops;
      ps_hloops = make_1Darray(10);

      Compute_2_loops(pba, ppm, pfo, k, z, mu, SPLIT, ps_hloops);
      halo = clock();
      // fprintf(stderr,"Halo Loops for the 1-st moment done! %f sec. needed.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      double J12102    = ps_hloops[0];
      double J21102    = ps_hloops[1];
      double Jdelta202 = ps_hloops[2];
      double JG202     = ps_hloops[3];

      double I2211  = ps_hloops[4];
      double I1311  = ps_hloops[5];
      double J12111 = ps_hloops[6];
      double J11211 = ps_hloops[7];
      double J21111 = ps_hloops[8];
      double N11    = ps_hloops[9];

      double ph_loops  = J12102 + J21102 + Jdelta202 + JG202;
      double plos_loops  = I2211 + I1311 + J12111 + J11211 + J21111 + N11;
      // pm_lin_IR =  pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      pm_lin_IR =  Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
      double khat     = 1. * pba->h;

      // fprintf(stderr,"Loops done in %f sec.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      FILE *fpa;
      char file_name[50];
      sprintf(file_name, "data/2_moment_DI_%g.txt",mu);
      fpa = fopen(file_name, "a");
      fprintf(fpa, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",\
          k, pm_lin_IR, J12102, J21102, Jdelta202, JG202, I2211, I1311, J12111, J11211, J21111, N11, ph_loops, plos_loops);
      fclose(fpa);

      free(ps_hloops);
      //free(ps_mloops);
      //*pk_nl = ph_tot;

      return _SUCCESS_;
}

int PS_hh_3(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT)
{
      int cleanup = 0;

      double pm_lin = 0., pm_lin_IR = 0., pm_1loop_IR = 0., pm_22 = 0., pm_13 = 0., pm_1loop =0., pm_ct = 0., ph_tot = 0.;

      //DL CLASS-PT values page 30
      double b1  = 2.0;
      double b2  = -1.0;
      double bG2 = 0.1;
      double btd = -0.1;
      double R2  =  5.0;
      double cs2 = 0.2;

      clock_t start, matter, halo;
      start = clock();

      double *ps_hloops;
      ps_hloops = make_1Darray(3);

      Compute_3_loops(pba, ppm, pfo, k, z, mu, SPLIT, ps_hloops);
      halo = clock();
      // fprintf(stderr,"Halo Loops for the 1-st moment done! %f sec. needed.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      double J21112 = ps_hloops[0];
      double J12112 = ps_hloops[1];

      double N12 = ps_hloops[2];

      double ph_loops  = J21112 + J12112 + N12;
      // pm_lin_IR =  pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      pm_lin_IR =  Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
      double khat     = 1. * pba->h;

      // fprintf(stderr,"Loops done in %f sec.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      FILE *fpa;
      char file_name[50];
      sprintf(file_name, "data/3_moment_DI_%g.txt",mu);
      fpa = fopen(file_name, "a");
      fprintf(fpa, "%12.6e %12.6e %12.6e %12.6e %12.6e %12.6e\n",\
          k, pm_lin_IR, J21112, J12112, N12, ph_loops);
      fclose(fpa);

      free(ps_hloops);
      //free(ps_mloops);
      //*pk_nl = ph_tot;

      return _SUCCESS_;
}

int PS_hh_4(   struct background * pba,
               struct primordial * ppm,
               struct fourier * pfo,
               double k,
               double z,
               double mu,
               long SPLIT)
{
      int cleanup = 0;

      clock_t start, matter, halo;
      start = clock();

      double *ps_hloops;
      ps_hloops = make_1Darray(1);

      Compute_4_loops(pba, ppm, pfo, k, z, mu, SPLIT, ps_hloops);
      halo = clock();
      // fprintf(stderr,"Halo Loops for the 1-st moment done! %f sec. needed.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      double N22 = ps_hloops[0];

      // double pm_lin_IR =  pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      double pm_lin_IR =  Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
      double khat     = 1. * pba->h;

      // fprintf(stderr,"Loops done in %f sec.\n", (double)(halo-start) / CLOCKS_PER_SEC);

      FILE *fpa;
      char file_name[50];
      sprintf(file_name, "data/4_moment_DI_%g.txt",mu);
      fpa = fopen(file_name, "a");
      fprintf(fpa, "%12.6e %12.6e %12.6e\n",\
          k, pm_lin_IR, N22);
      fclose(fpa);

      free(ps_hloops);
      //free(ps_mloops);
      //*pk_nl = ph_tot;

      return _SUCCESS_;
}

/**
 * Compute the loop contributions dure to nonlinear evolution of matter fluctuations and nonlinear halo bias, present for Gaussian initial conditions
 * The function G_loop_integrands() defines the integrand and Compute_G_loops() computes the integrals
 *
 * @param Cx           Input: pointer to cosmology structure
 * @param k            Input: wavenumber
 * @param z            Input: redshift of interest
 * @param M            Input: halo mass, used in computing the halo bias
 * @param IR_switch    Input: switch to decide whether to perform IR resummation or no
 * @param hm_switch    Input: switch to decide whether to compute the 1loop terms due to matter or bias
 * @param SPLIT        Input: switch to set the method to perform the wiggle-nowiggle split of matter power spectrum
 * @param result       Output: an output array containing the results of the 1loop terms,
 *                             has 2 elements for hm_switch=MATTER, and 6 elements for hm_switch=HALO
 * @return void
 */

int Compute_G_loops(
                    //struct Cosmology *Cx,
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    short has_ir,
                    long hm_switch,
                    long SPLIT,
                    double *result)
{
      struct integrand_parameters2 par;

      int ncomp = 0;

      if(hm_switch == HALO)
        ncomp = 6;
      else if(hm_switch == MATTER)
        ncomp = 2;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
      double plin_IR_k;

      if(has_ir == _TRUE_) {
        plin_IR_k   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      }
      else{plin_IR_k = 0;}


      // Parsing information to the integrand
      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = k;
      par.p5  = z;
      par.p6  = log(1.e-5);
      par.p7  = log(1000.);
      if (has_ir == _TRUE_)
        par.p13 = WIR;
      else
        par.p13 = NOIR;
      par.p14 = hm_switch;
      par.p15 = SPLIT;
      par.p8  = plin_IR_k;  /* Note: since cuhre integrator is parallelized, if evaluating all the
                             * IR resummed power spectra in the integrand, it will build the interpolators
                             * For each thread and creats a mess. Here we evaluate the p(k) (which builds the interpolator)
                             * The first time it is called, and then call the p(q) and p(kmq) and p(kpq) inside the integrand
                             */

      // int ndim = 2, nvec = 1, verbose = 0, last = 4, seed = 0,
      //     mineval = 0, maxeval = 2e6, key1 =50, key2 = 50, key3 = 1, maxpass = 200,
      //       border =0, maxchisq = 10, mindeviation=0.25,
      //       ngiven = 0, ldxgiven = ndim, nextra = 0;
      // int nregions, neval;

      // Divonne(ndim, ncomp, G_loop_integrands, &par, nvec,
      //           RelErr, AbsErr, verbose, seed,
      //          mineval, maxeval,  key1, key2, key3, maxpass,
      //         border, maxchisq, mindeviation,
      //         ngiven, ldxgiven, NULL, nextra, NULL,
      //        NULL, NULL, &nregions, &neval, fail, result, error, prob);

      int key;
      if(hm_switch == MATTER)
          key = 7;
      else if(hm_switch == HALO)
          key = 13;

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e8;
      int nregions, neval;

      Cuhre(ndim,ncomp, G_loop_integrands, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, fail, result, error, prob);

      // for(int i =0; i<ncomp; i++)
      //     printf("Gloops integral : %d %12.6e %12.6e %12.6e %12.6e %d \n", i, k, result[i], error[i], prob[i], fail[i]);

      return _SUCCESS_;
}

static int G_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);

      struct background *pba = pij.pba;
      struct primordial *ppm = pij.ppm;
      struct fourier *pfo    = pij.pfo;
      double k               = pij.p4;
      double z               = pij.p5;
      double logqmin         = pij.p6;
      double logqmax         = pij.p7;
      long IR_switch         = pij.p13;
      long hm_switch         = pij.p14;
      long SPLIT             = pij.p15;
      double plin_IR_k       = pij.p8;

      //fprintf(stderr,"%e, %e, %e,%e \n",k,z,logqmin,logqmax);
      double logq  = (x[0] * (logqmax - logqmin) + logqmin);
      double cos   = (2. * x[1] - 1.);
      double q     = exp(logq);
      double kmq   = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double kpq   = pow(fabs(pow(q, 2.) + pow(k, 2.) + 2. * q * k * cos), 1./2.);
      // fprintf(stderr, "%ld", hm_switch);
      static int cleanup = 0;

      /// Model used in 1907.06666, the integrals are given in the appendix, Eq. A1, note that my S2_s = sigma^2(q,k-1) and F2_s = F2(q,k-q) in their notation.
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle.
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to
      /// devide each 0 power spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).

      /// The integrands below correspond to the follwing bias combinaions:
      /*  ff[0] : b1b2
          ff[1] : b1bG2
          ff[2] : b2b2
          ff[3] : bG2bG2
          ff[4] : b2bG2
          ff[5] : b1b3nl, b3nl = 2*(bG2+2/5*btd)

          ff[0] : P22_m
          ff[1] : P13_m
      */

      ////Defining regularized version of p22 integrand as in Senatore et al. We need to first define heavisde step function

      double theta_kmq = 0, theta_kpq = 0.;
      double plin_k;
      double plin_q;
      double plin_kmq;
      double plin_kpq;
      double plin_IR_q;
      double plin_IR_kmq;

      if(kmq>=q)
        theta_kmq = 1.;
      else
        theta_kmq = 0.;

      if(kpq>=q)
        theta_kpq = 1.;
      else
        theta_kpq = 0.;

      int nn =0;
      if(hm_switch == HALO)
         nn = 6;
      else if(hm_switch == MATTER)
         nn = 2;

      if(kmq <= exp(logqmax) && kmq>=exp(logqmin) && kpq <= exp(logqmax) && kpq>=exp(logqmin)){
        if(IR_switch == NOIR){
          plin_k   =  Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);
          plin_q   =  Pk_dlnPk(pba, ppm, pfo, q, z, LPOWER);
          plin_kmq =  Pk_dlnPk(pba, ppm, pfo, kmq, z, LPOWER);
          plin_kpq =  Pk_dlnPk(pba, ppm, pfo, kpq, z, LPOWER);

            if(hm_switch == HALO)
            {
                ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * F2_s(q, k, cos);
                ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * F2_s(q, k, cos) * S2_s(q, k, cos);
                ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * (plin_kmq - plin_q);
                ff[3] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * pow(S2_s(q, k, cos), 2.);
                ff[4] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * S2_s(q,k,cos);
                ff[5] = 1./pow(2.*M_PI,3.) * 4. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_k * plin_q   * S2_s(q, k, cos) * F2(q, k, -cos);
            }
            else if(hm_switch == MATTER)
            {
                ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_q * plin_kmq * pow(F2_s(q, k, cos), 2.);
                ff[1] = 1./pow(2.*M_PI,3.) * 6. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_k * plin_q * F3_s(q, k, cos);
            }
            for(int i =0;i<nn;i++){
              if (isnan(ff[i]))
                 printf("%d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",i, k, q, kmq, plin_q, plin_kmq, ff[i]);
            }

        }
        else if(IR_switch == WIR){
          plin_IR_q   = pm_IR_LO(pba, ppm, pfo, q, z, SPLIT);
          plin_IR_kmq = pm_IR_LO(pba, ppm, pfo, kmq, z, SPLIT);
          double plin_IR_kpq = pm_IR_LO(pba, ppm, pfo, kpq, z, SPLIT);

          if(hm_switch == HALO)
          {
              ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * F2_s(q, k, cos);
              ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * F2_s(q, k, cos) * S2_s(q, k, cos);
              // ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * (plin_IR_kmq - plin_IR_q);
              ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * (plin_IR_kmq - plin_IR_q);
              ff[3] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(S2_s(q, k, cos), 2.);
              ff[4] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * S2_s(q, k, cos);
              ff[5] = 1./pow(2.*M_PI,3.) * 4. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q * S2_s(q, k, cos) * F2(q, k, -cos);
          }
          else if(hm_switch == MATTER)
          {
              // ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(F2_s(q, k, cos), 2.) * theta_kmq\
              //       + 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kpq * pow(F2_s(q, k, -cos), 2.) * theta_kpq;
              ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(F2_s(q, k, cos), 2.);
              ff[1] = 1./pow(2.*M_PI,3.) * 6. * 2. * M_PI * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q * F3_s(q, k, cos);
          }


          for(int i =0;i<nn;i++){
            if (isnan(ff[i]))
                printf("%d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",i, k, q, kmq, plin_IR_q, plin_IR_kmq, ff[i]);

          }
        }
    }
    else{
      for(int i =0;i <nn; i++)
        ff[i] = 0.;
    }

      return 0;
}

int Compute_1_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result)
{
      struct integrand_parameters2 par;

      int ncomp = 10;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
      double plin_IR_k;

      // plin_IR_k   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      plin_IR_k   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);

      // Parsing information to the integrand
      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = k;
      par.p5  = z;
      par.p6  = log(1.e-4);
      par.p7  = log(1000.);
      par.p13 = WIR;
      par.p15 = SPLIT;
      par.p9  = mu;
      par.p8  = plin_IR_k;  /* Note: since cuhre integrator is parallelized, if evaluating all the
                             * IR resummed pow spectra in the integrand, it will build the interpolators
                             * For each thread and creats a mess. Here we evaluate the p(k) (which builds the interpolator)
                             * The first time it is called, and then call the p(q) and p(kmq) and p(kpq) inside the integrand
                             */

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e8;
      int nregions, neval;

      int key = 7;

      Cuhre(ndim,ncomp, G1_loop_integrands, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, fail, result, error, prob);

      // for(int i =0; i<ncomp; i++)
      //     printf("Gloops integral : %d %12.6e %12.6e %12.6e %12.6e %d \n", i, k, result[i], error[i], prob[i], fail[i]);

      return _SUCCESS_;
}

static int G1_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);

      struct background *pba = pij.pba;
      struct primordial *ppm = pij.ppm;
      struct fourier *pfo    = pij.pfo;
      double k               = pij.p4;
      double z               = pij.p5;
      double logqmin         = pij.p6;
      double logqmax         = pij.p7;
      long SPLIT             = pij.p15;
      double plin_IR_k       = pij.p8;
      double mu              = pij.p9;

      //fprintf(stderr,"%e, %e, %e,%e \n",k,z,logqmin,logqmax);
      double logq  = (x[0] * (logqmax - logqmin) + logqmin);
      double cos   = (2. * x[1] - 1.);
      double q     = exp(logq);
      double kmq   = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double kpq   = pow(fabs(pow(q, 2.) + pow(k, 2.) + 2. * q * k * cos), 1./2.);
      // fprintf(stderr, "%ld", hm_switch);
      static int cleanup = 0;

      /// Model used in 1907.06666, the integrals are given in the appendix, Eq. A1, note that my S2_s = sigma^2(q,k-1) and F2_s = F2(q,k-q) in their notation.
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle.
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to
      /// devide each 0 pow spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).

      /// The integrands below correspond to the follwing bias combinaions:
      /*  ff[0] : b1b2
          ff[1] : b1bG2
          ff[2] : b2b2
          ff[3] : bG2bG2
          ff[4] : b2bG2
          ff[5] : b1b3nl, b3nl = 2*(bG2+2/5*btd)

          ff[0] : P22_m
          ff[1] : P13_m
      */

      ////Defining regularized version of p22 integrand as in Senatore et al. We need to first define heavisde step function

      double theta_kmq = 0, theta_kpq = 0.;
      double plin_k;
      double plin_q;
      double plin_kmq;
      double plin_kpq;
      double plin_IR_q;
      double plin_IR_kmq;

      if(kmq>=q)
        theta_kmq = 1.;
      else
        theta_kmq = 0.;

      if(kpq>=q)
        theta_kpq = 1.;
      else
        theta_kpq = 0.;

      int nn =10;

      if(kmq <= exp(logqmax) && kmq>=exp(logqmin) && kpq <= exp(logqmax) && kpq>=exp(logqmin)){

        // plin_IR_q   = pm_IR_LO(pba, ppm, pfo, q, z, SPLIT);
        plin_IR_q   = Pk_dlnPk(pba, ppm, pfo, q, z, LPOWER);
        // plin_IR_kmq = pm_IR_LO(pba, ppm, pfo, kmq, z, SPLIT);
        plin_IR_kmq   = Pk_dlnPk(pba, ppm, pfo, kmq, z, LPOWER);

        // double plin_IR_kpq = pm_IR_LO(pba, ppm, pfo, kpq, z, SPLIT);
        double plin_IR_kpq   = Pk_dlnPk(pba, ppm, pfo, kpq, z, LPOWER);


        double kmq2 = pow(k,2.)+pow(q,2.)-2.*cos*k*q;

        ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * F2_s(q, k, cos) * G2_s(q, k, cos);
        ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *  (F3_s(q, k, cos) + G3_s(q, k, cos));
        ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * G2_s(q, k, cos);
        ff[3] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * S2_s(q, k, cos) * G2_s(q, k, cos);
        ff[4] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_k   * S2_s(q, k, cos) * F2(q, k, -cos);
      
        ff[5] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   LoS1(q,k,cos,mu) * G2(q,k,-cos);
        ff[6] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   LoS2(q,k,cos,mu) * F2(q,k,-cos);
        ff[7] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS2(q,k,cos,mu) * F2_s(q, k, cos);
        ff[8] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS2(q,k,cos,mu);
        ff[9] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS2(q,k,cos,mu) * S2_s(q, k, cos);

        // ff[5] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   -1/7*(M_PI*(k - q*cos)*(-6*k*q + 7*(pow(k,2) + pow(q,2))*cos - 8*k*q*pow(cos,2))*mu)/(k*q*(kmq2));
        // ff[6] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (M_PI*cos*(10*k*q - 7*(pow(k,2) + pow(q,2))*cos + 4*k*q*pow(cos,2))*mu)/(7.*k*pow(q,2));
        // ff[7] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (pow(k,2)*M_PI*cos*(7*k*cos + q*(3 - 10*pow(cos,2)))*mu)/(7.*pow(q,2)*(kmq2));
        // ff[8] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (2*M_PI*cos*mu)/q;
        // ff[9] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (2*pow(k,2)*M_PI*cos*(-1 + pow(cos,2))*mu)/(q*(kmq2));


        // ff[5] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (k-q*cos)/(kmq2) *    G2(q, k, -cos);
        // ff[6] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (q*cos)/(pow(q,2.)) * F2(q, k, -cos);
        // ff[7] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (q*cos)/(pow(q,2.)) * F2_s(q, k, cos);
        // ff[8] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (q*cos)/(pow(q,2.));
        // ff[9] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (q*cos)/(pow(q,2.)) * S2_s(q, k, cos);
      }
      else{
        for(int i =0;i <nn; i++)
          ff[i] = 0.;
      }

      return 0;
}

int Compute_2_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result)
{
      struct integrand_parameters2 par;

      int ncomp = 10;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
      double plin_IR_k;

      // plin_IR_k   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      plin_IR_k   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);

      // Parsing information to the integrand
      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = k;
      par.p5  = z;
      par.p6  = log(1.e-5);
      par.p7  = log(1000.);
      par.p13 = WIR;
      par.p15 = SPLIT;
      par.p9  = mu;
      par.p8  = plin_IR_k;  /* Note: since cuhre integrator is parallelized, if evaluating all the
                             * IR resummed pow spectra in the integrand, it will build the interpolators
                             * For each thread and creats a mess. Here we evaluate the p(k) (which builds the interpolator)
                             * The first time it is called, and then call the p(q) and p(kmq) and p(kpq) inside the integrand
                             */

      // int ndim = 2, nvec = 1, verbose = 0, last = 4, seed = 0,
      //     mineval = 0, maxeval = 2e6, key1 =50, key2 = 50, key3 = 1, maxpass = 200,
      //       border =0, maxchisq = 10, mindeviation=0.25,
      //       ngiven = 0, ldxgiven = ndim, nextra = 0;
      // int nregions, neval;

      // Divonne(ndim, ncomp, G2_loop_integrands, &par, nvec,
      //           RelErr, AbsErr, verbose, seed,
      //          mineval, maxeval,  key1, key2, key3, maxpass,
      //         border, maxchisq, mindeviation,
      //         ngiven, ldxgiven, NULL, nextra, NULL,
      //        NULL, NULL, &nregions, &neval, fail, result, error, prob);
      
      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e6;
      int nregions, neval;

      int key = 7;

      Cuhre(ndim,ncomp, G2_loop_integrands, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, fail, result, error, prob);

      // for(int i =0; i<ncomp; i++)
          // printf("Gloops integral : %d %12.6e %12.6e %12.6e %12.6e %d \n", i, k, result[i], error[i], prob[i], fail[i]);

      return _SUCCESS_;
}

static int G2_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);

      struct background *pba = pij.pba;
      struct primordial *ppm = pij.ppm;
      struct fourier *pfo    = pij.pfo;
      double k               = pij.p4;
      double z               = pij.p5;
      double logqmin         = pij.p6;
      double logqmax         = pij.p7;
      long SPLIT             = pij.p15;
      double plin_IR_k       = pij.p8;
      double mu              = pij.p9;

      //fprintf(stderr,"%e, %e, %e,%e \n",k,z,logqmin,logqmax);
      double logq  = (x[0] * (logqmax - logqmin) + logqmin);
      double cos   = (2. * x[1] - 1.);
      double q     = exp(logq);
      double kmq   = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double kpq   = pow(fabs(pow(q, 2.) + pow(k, 2.) + 2. * q * k * cos), 1./2.);
      // fprintf(stderr, "%ld", hm_switch);
      static int cleanup = 0;

      /// Model used in 1907.06666, the integrals are given in the appendix, Eq. A1, note that my S2_s = sigma^2(q,k-1) and F2_s = F2(q,k-q) in their notation.
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle.
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to
      /// devide each 0 pow spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).

      /// The integrands below correspond to the follwing bias combinaions:
      /*  ff[0] : b1b2
          ff[1] : b1bG2
          ff[2] : b2b2
          ff[3] : bG2bG2
          ff[4] : b2bG2
          ff[5] : b1b3nl, b3nl = 2*(bG2+2/5*btd)

          ff[0] : P22_m
          ff[1] : P13_m
      */

      ////Defining regularized version of p22 integrand as in Senatore et al. We need to first define heavisde step function

      double theta_kmq = 0, theta_kpq = 0.;
      double plin_k;
      double plin_q;
      double plin_kmq;
      double plin_kpq;
      double plin_IR_q;
      double plin_IR_kmq;

      if(kmq>=q)
        theta_kmq = 1.;
      else
        theta_kmq = 0.;

      if(kpq>=q)
        theta_kpq = 1.;
      else
        theta_kpq = 0.;

      int nn =10;

      if(kmq <= exp(logqmax) && kmq>=exp(logqmin) && kpq <= exp(logqmax) && kpq>=exp(logqmin)){

        // plin_IR_q   = pm_IR_LO(pba, ppm, pfo, q, z, SPLIT);
        plin_IR_q   = Pk_dlnPk(pba, ppm, pfo, q, z, LPOWER);
        // plin_IR_kmq = pm_IR_LO(pba, ppm, pfo, kmq, z, SPLIT);
        plin_IR_kmq   = Pk_dlnPk(pba, ppm, pfo, kmq, z, LPOWER);

        // double plin_IR_kpq = pm_IR_LO(pba, ppm, pfo, kpq, z, SPLIT);
        double plin_IR_kpq   = Pk_dlnPk(pba, ppm, pfo, kpq, z, LPOWER);

        double kmq2 = pow(k,2.)+pow(q,2.)-2.*cos*k*q;

        // ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.)) * G2(q, k, -cos);
        // ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.)) * F2_s(q, k, cos);
        // ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.));
        // ff[3] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.)) * S2_s(q, k, cos);

        // ff[0] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   -1/14*(M_PI*(-6*k*q + 7*(pow(k,2) + pow(q,2))*cos - 8*k*q*pow(cos,2))*(q*(-1 + pow(cos,2)) + (q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2)))/(k*pow(q,2)*(kmq2));
        // ff[1] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (pow(k,2)*M_PI*(7*k*cos + q*(3 - 10*pow(cos,2)))*(q*(-1 + pow(cos,2)) + (q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2)))/(14.*pow(q,2)*pow(kmq2,2));
        // ff[2] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (M_PI*q*(-1 + pow(cos,2)) + M_PI*(q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2))/(q*(kmq2));
        // ff[3] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (pow(k,2)*M_PI*(-1 + pow(cos,2))*(q*(-1 + pow(cos,2)) + (q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2)))/(q*pow(kmq2,2));

        ff[0] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   LoS3(q,k,cos,mu) * G2(q,k,-cos);
        ff[1] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS3(q,k,cos,mu) * F2_s(q,k,cos);
        ff[2] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS3(q,k,cos,mu);
        ff[3] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS3(q,k,cos,mu) * S2_s(q,k,cos);

        ff[4] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow(G2_s(q, k, cos),2.);
        ff[5] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   G3_s(q, k, cos);

        ff[6] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   LoS1(q,k,cos,mu) * G2(q,k,-cos);
        ff[7] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   LoS2(q,k,cos,mu) * F2(q,k,-cos);
        ff[8] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS2(q,k,cos,mu) * G2_s(q, k, cos);
        ff[9] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS4(q,k,cos,mu);

        // ff[6] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   -1/7*(M_PI*(k - q*cos)*(-6*k*q + 7*(pow(k,2) + pow(q,2))*cos - 8*k*q*pow(cos,2))*mu)/(k*q*(kmq2));
        // ff[7] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (M_PI*cos*(10*k*q - 7*(pow(k,2) + pow(q,2))*cos + 4*k*q*pow(cos,2))*mu)/(7.*k*pow(q,2));
        // ff[8] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (pow(k,2)*M_PI*cos*(7*k*cos - q*(1 + 6*pow(cos,2)))*mu)/(7.*pow(q,2)*(kmq2));
        // ff[9] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (k*M_PI*(-((k - 2*q*cos)*(-1 + pow(cos,2))) + (-k + 4*q*cos + 3*k*pow(cos,2) - 6*q*pow(cos,3))*pow(mu,2)))/(pow(q,2)*(kmq2));

        // ff[6] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (k-q*cos)/(kmq2) * G2(q, k, -cos);
        // ff[7] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q *   (q*cos)/(pow(q,2.)) * F2(q, k, -cos);
        // ff[8] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (q*cos)/(pow(q,2.)) * G2_s(q, k, cos);
        // ff[9] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (q*cos)/(pow(q,2.)) * ((q*cos)/(pow(q,2.)) + (k-q*cos)/(kmq2));
      }
      else{
        for(int i =0;i <nn; i++)
          ff[i] = 0.;
      }

      return 0;
}

int Compute_3_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result)
{
      struct integrand_parameters2 par;

      int ncomp = 3;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
      double plin_IR_k;

      // plin_IR_k   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      plin_IR_k   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);

      // Parsing information to the integrand
      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = k;
      par.p5  = z;
      par.p6  = log(1.e-5);
      par.p7  = log(1000.);
      par.p13 = WIR;
      par.p15 = SPLIT;
      par.p9  = mu;
      par.p8  = plin_IR_k;  /* Note: since cuhre integrator is parallelized, if evaluating all the
                             * IR resummed pow spectra in the integrand, it will build the interpolators
                             * For each thread and creats a mess. Here we evaluate the p(k) (which builds the interpolator)
                             * The first time it is called, and then call the p(q) and p(kmq) and p(kpq) inside the integrand
                             */

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e6;
      int nregions, neval;

      int key = 13;

      Cuhre(ndim,ncomp, G3_loop_integrands, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, fail, result, error, prob);

      // for(int i =0; i<ncomp; i++)
      //     printf("Gloops integral : %d %12.6e %12.6e %12.6e %12.6e %d \n", i, k, result[i], error[i], prob[i], fail[i]);

      return _SUCCESS_;
}

static int G3_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);

      struct background *pba = pij.pba;
      struct primordial *ppm = pij.ppm;
      struct fourier *pfo    = pij.pfo;
      double k               = pij.p4;
      double z               = pij.p5;
      double logqmin         = pij.p6;
      double logqmax         = pij.p7;
      long SPLIT             = pij.p15;
      double plin_IR_k       = pij.p8;
      double mu              = pij.p9;

      //fprintf(stderr,"%e, %e, %e,%e \n",k,z,logqmin,logqmax);
      double logq  = (x[0] * (logqmax - logqmin) + logqmin);
      double cos   = (2. * x[1] - 1.);
      double q     = exp(logq);
      double kmq   = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double kpq   = pow(fabs(pow(q, 2.) + pow(k, 2.) + 2. * q * k * cos), 1./2.);
      // fprintf(stderr, "%ld", hm_switch);
      static int cleanup = 0;

      /// Model used in 1907.06666, the integrals are given in the appendix, Eq. A1, note that my S2_s = sigma^2(q,k-1) and F2_s = F2(q,k-q) in their notation.
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle.
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to
      /// devide each 0 pow spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).

      /// The integrands below correspond to the follwing bias combinaions:
      /*  ff[0] : b1b2
          ff[1] : b1bG2
          ff[2] : b2b2
          ff[3] : bG2bG2
          ff[4] : b2bG2
          ff[5] : b1b3nl, b3nl = 2*(bG2+2/5*btd)

          ff[0] : P22_m
          ff[1] : P13_m
      */

      ////Defining regularized version of p22 integrand as in Senatore et al. We need to first define heavisde step function

      double theta_kmq = 0, theta_kpq = 0.;
      double plin_k;
      double plin_q;
      double plin_kmq;
      double plin_kpq;
      double plin_IR_q;
      double plin_IR_kmq;

      if(kmq>=q)
        theta_kmq = 1.;
      else
        theta_kmq = 0.;

      if(kpq>=q)
        theta_kpq = 1.;
      else
        theta_kpq = 0.;

      int nn =3;

      if(kmq <= exp(logqmax) && kmq>=exp(logqmin) && kpq <= exp(logqmax) && kpq>=exp(logqmin)){

        // plin_IR_q   = pm_IR_LO(pba, ppm, pfo, q, z, SPLIT);
        plin_IR_q   = Pk_dlnPk(pba, ppm, pfo, q, z, LPOWER);
        // plin_IR_kmq = pm_IR_LO(pba, ppm, pfo, kmq, z, SPLIT);
        plin_IR_kmq   = Pk_dlnPk(pba, ppm, pfo, kmq, z, LPOWER);

        // double plin_IR_kpq = pm_IR_LO(pba, ppm, pfo, kpq, z, SPLIT);
        double plin_IR_kpq   = Pk_dlnPk(pba, ppm, pfo, kpq, z, LPOWER);

        double kmq2 = pow(k,2.)+pow(q,2.)-2.*cos*k*q;

        // ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.)) * G2_s(q, k, cos);
        // ff[1] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q   * (k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.)) * G2(q, k, -cos);
        // ff[2] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (k-q*cos)/(kmq2) * pow((q*cos)/(pow(q,2.)),2.);
      
        // ff[0] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (pow(k,2)*M_PI*(7*k*cos - q*(1 + 6*pow(cos,2)))*(q*(-1 + pow(cos,2)) + (q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2)))/(14.*pow(q,2)*pow(kmq2,2));
        // ff[1] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q   * -1/14*(M_PI*(-6*k*q + 7*(pow(k,2) + pow(q,2))*cos - 8*k*q*pow(cos,2))*(q*(-1 + pow(cos,2)) + (q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2)))/(k*pow(q,2)*(kmq2));
        // ff[2] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (-(M_PI*(k - 3*q*cos)*(-1 + pow(cos,2))*mu) + M_PI*(q*cos*(3 - 5*pow(cos,2)) + k*(-1 + 3*pow(cos,2)))*pow(mu,3))/(pow(q,2)*(kmq2));

        ff[0] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS3(q,k,cos,mu) * G2_s(q,k,cos);
        ff[1] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_k * plin_IR_q   * LoS3(q,k,cos,mu) * G2(q,k,-cos);
        ff[2] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS5(q,k,cos,mu);
      }
      else{
        for(int i =0;i <nn; i++)
          ff[i] = 0.;
      }

      return 0;
}

int Compute_4_loops(
                    struct background * pba,
                    struct primordial * ppm,
                    struct fourier * pfo,
                    double k,
                    double z,
                    double mu,
                    long SPLIT,
                    double *result)
{
      struct integrand_parameters2 par;

      int ncomp = 1;

      double AbsErr = 0.0;                            // Required absolute error (0.0 for none)
      double RelErr = 1.e-3;                          // Required relative error (1.0e-2)
      int    fail[ncomp];
      double error[ncomp];
      double prob[ncomp];
      double plin_IR_k;

      // plin_IR_k   = pm_IR_LO(pba, ppm, pfo, k, z, SPLIT);
      plin_IR_k   = Pk_dlnPk(pba, ppm, pfo, k, z, LPOWER);

      // Parsing information to the integrand
      par.ppm = ppm;
      par.pba = pba;
      par.pfo = pfo;
      par.p4  = k;
      par.p5  = z;
      par.p6  = log(1.e-5);
      par.p7  = log(1000.);
      par.p13 = WIR;
      par.p15 = SPLIT;
      par.p9  = mu;
      par.p8  = plin_IR_k;  /* Note: since cuhre integrator is parallelized, if evaluating all the
                             * IR resummed pow spectra in the integrand, it will build the interpolators
                             * For each thread and creats a mess. Here we evaluate the p(k) (which builds the interpolator)
                             * The first time it is called, and then call the p(q) and p(kmq) and p(kpq) inside the integrand
                             */

      int ndim = 2,  nvec = 1, verbose = 0, last = 4, mineval = 0, maxeval = 1e8;
      int nregions, neval;

      int key = 13;

      Cuhre(ndim,ncomp, G4_loop_integrands, &par, nvec,
              RelErr, AbsErr, verbose | last,
               mineval, maxeval, key,
               NULL, NULL, &nregions, &neval, fail, result, error, prob);

      // for(int i =0; i<ncomp; i++)
      //     printf("Gloops integral : %d %12.6e %12.6e %12.6e %12.6e %d \n", i, k, result[i], error[i], prob[i], fail[i]);

      return _SUCCESS_;
}

static int G4_loop_integrands(
                             const int *ndim,
                             const cubareal x[],
                             const int *ncomp,
                             cubareal ff[],
                             void *p)
{

      struct integrand_parameters2 pij;
      pij = *((struct integrand_parameters2 *)p);

      struct background *pba = pij.pba;
      struct primordial *ppm = pij.ppm;
      struct fourier *pfo    = pij.pfo;
      double k               = pij.p4;
      double z               = pij.p5;
      double logqmin         = pij.p6;
      double logqmax         = pij.p7;
      long SPLIT             = pij.p15;
      double plin_IR_k       = pij.p8;
      double mu              = pij.p9;

      //fprintf(stderr,"%e, %e, %e,%e \n",k,z,logqmin,logqmax);
      double logq  = (x[0] * (logqmax - logqmin) + logqmin);
      double cos   = (2. * x[1] - 1.);
      double q     = exp(logq);
      double kmq   = pow(fabs(pow(q, 2.) + pow(k, 2.) - 2. * q * k * cos), 1./2.);
      double kpq   = pow(fabs(pow(q, 2.) + pow(k, 2.) + 2. * q * k * cos), 1./2.);
      // fprintf(stderr, "%ld", hm_switch);
      static int cleanup = 0;

      /// Model used in 1907.06666, the integrals are given in the appendix, Eq. A1, note that my S2_s = sigma^2(q,k-1) and F2_s = F2(q,k-q) in their notation.
      /// Factor of 2. * (logqmax - logqmin) is due to change of variable from 0 to logarithmic k, and a factor of 2*PI is due to integration over azimuthal angle.
      /// Note that to compare the theoretical predictions against Emiliano's measurement, since he is using a different notation for Fourier transform, I need to
      /// devide each 0 pow spectrum by a factor of 1/pow(2.*M_PI,3.), which I do in my pk_lin() function. If using another notation for Fourier transform
      /// (the one that I usually use, which has a factor of 1/pow(2*M_PI,3) in the definition), you need to multiply these integrands by a factor of 1/pow(2*M_PI,3).

      /// The integrands below correspond to the follwing bias combinaions:
      /*  ff[0] : b1b2
          ff[1] : b1bG2
          ff[2] : b2b2
          ff[3] : bG2bG2
          ff[4] : b2bG2
          ff[5] : b1b3nl, b3nl = 2*(bG2+2/5*btd)

          ff[0] : P22_m
          ff[1] : P13_m
      */

      ////Defining regularized version of p22 integrand as in Senatore et al. We need to first define heavisde step function

      double theta_kmq = 0, theta_kpq = 0.;
      double plin_k;
      double plin_q;
      double plin_kmq;
      double plin_kpq;
      double plin_IR_q;
      double plin_IR_kmq;

      if(kmq>=q)
        theta_kmq = 1.;
      else
        theta_kmq = 0.;

      if(kpq>=q)
        theta_kpq = 1.;
      else
        theta_kpq = 0.;

      int nn =1;

      if(kmq <= exp(logqmax) && kmq>=exp(logqmin) && kpq <= exp(logqmax) && kpq>=exp(logqmin)){

        // plin_IR_q   = pm_IR_LO(pba, ppm, pfo, q, z, SPLIT);
        plin_IR_q   = Pk_dlnPk(pba, ppm, pfo, q, z, LPOWER);
        // plin_IR_kmq = pm_IR_LO(pba, ppm, pfo, kmq, z, SPLIT);
        plin_IR_kmq   = Pk_dlnPk(pba, ppm, pfo, kmq, z, LPOWER);

        // double plin_IR_kpq = pm_IR_LO(pba, ppm, pfo, kpq, z, SPLIT);
        double plin_IR_kpq   = Pk_dlnPk(pba, ppm, pfo, kpq, z, LPOWER);

        double kmq2 = pow(k,2.)+pow(q,2.)-2.*cos*k*q;

        // ff[0] = 1./pow(2.*M_PI,3.) * 2. * 2. * M_PI  * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * pow((k-q*cos)/(kmq2) * (q*cos)/(pow(q,2.)),2.);
        
        // ff[0] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * (M_PI*(3*pow(q,2)*pow(-1 + pow(cos,2),2) - 2*(-1 + pow(cos,2))*(2*pow(k,2) - 12*k*q*cos + 3*pow(q,2)*(-1 + 5*pow(cos,2)))*pow(mu,2) + (8*k*q*cos*(3 - 5*pow(cos,2)) + 4*pow(k,2)*(-1 + 3*pow(cos,2)) + pow(q,2)*(3 - 30*pow(cos,2) + 35*pow(cos,4)))*pow(mu,4)))/(4.*pow(kmq2,2));
        
        ff[0] = 1./pow(2.*M_PI,3.) * 2. * (logqmax - logqmin) * pow(q, 3.) * plin_IR_q * plin_IR_kmq * LoS6(q,k,cos,mu);
      
      }
      else{
        for(int i =0;i <nn; i++)
          ff[i] = 0.;
      }

      return 0;
}

//*********************************************************************************////
//****** Second-order kernels *******////
//*********************************************************************************////

////This is F2(q,k-q) since only this configuration apears in the loops ////
double F2_s(double k1,double k2,double mu)
{
    // double kmq2 = pow(q,2.)+pow(k,2.)-2.*q*k*mu;
    double f = 0.;
    // if(k2<0.01)
    //     f = 3.*pow(k2/k1,2.) - 5./7. * pow(k2/k1,2.)* pow(mu,2.) + 13./14.* pow(k2/k1,3.) * mu - 10./7. * pow(k2/k1,3.)*pow(mu,3.);
    // else
    f = pow(k2,2.)*(3.*k1 + 7.*k2*mu - 10.*k1*pow(mu,2.))/(14.*k1*(pow(k1,2.) + pow(k2,2.) - 2.*k1*k2*mu));
    // f = (5*(q + k*mu) + ((2*pow(k,2) + 5*pow(q,2))*(-q + k*mu))/kmq2)/(14.*q);

    return f;
}

////This is G2(q,k-q) since only this configuration apears in the loops ////
double G2_s(double q,double k,double mu)
{
  double kmq2   = pow(q,2.)+pow(k,2.)-2.*q*k*mu;
  // double f = 3./14. * (1. + k*mu) - (1./14.*q)*(3.*pow(q,2.)+4.*pow(k,2.))*(q-k*mu)/(pow(q,2.)+pow(k,2.)-2.*q*k*mu);
  double f = (3*(q + k*mu) + ((4*pow(k,2) + 3*pow(q,2))*(-q + k*mu))/kmq2)/(14.*q);
  return f;
}

////This is S2(k1,k2-k1) since only this configuration apears in the loops ////
double S2_s(double k1,double k2,double mu)
{
      double f = 0.;
      f = pow(-k1 + k2*mu,2.)/(pow(k1,2.) + pow(k2,2.) - 2.*k1*k2*mu) - 1.;

      return f;
}

////This is the symmetrized F3(q,-q,kmq) given in 1603.04405
double F3_s(double q, double k, double mu)
{

    double f = 0.;
    double kmq2  = pow(q,2.)+pow(k,2.)-2.*q*k*mu;
    double kpq2  = pow(q,2.)+pow(k,2.)+2.*q*k*mu;
    double kdotq = k*q*mu;

    // f = 1./kmq2*(5.*pow(k,2.)/63. - (11.*kdotq)/54. - pow(k,2.)*pow(kdotq,2.)/(6.*pow(q,4.))\
    //    +  19.*pow(kdotq,3.)/(63.*pow(q,4.))  - 23.*pow(k,2.)*kdotq/(378.*pow(q,2.))\
    //    - 23.*pow(kdotq,2.)/(378.*pow(q,2.)) + pow(kdotq,3.)/(9.*pow(k*q,2.)));

    // f = -1./(756.*k*pow(q,2.))*(98.*pow(k,3.)*pow(mu,2.)-14.*k*pow(q,2.)*(6+pow(mu,2.)+(2*pow(k,2.)+7.*pow(q,2.))*(-((q-k*mu)*(-6.*k*q+7.*(pow(k,2.)+pow(q,2.))*mu-8.*k*q*pow(mu,2.))/kmq2 \
    // + (q+k*mu)*(6.*k*q+7.*(pow(k,2.)+pow(q,2.))*mu+8.*k*q*pow(mu,2.))/kpq2))));

    // f = -1/756*  (98*pow(k,3)*pow(mu,2) - 14*k*pow(q,2)*(6 + pow(mu,2)) + (2*pow(k,2) \
    //     + 7*pow(q,2)) *(-(((q - k*mu)*(-6*k*q + 7*(pow(k,2) + pow(q,2))*mu - 8*k*q*pow(mu,2)))/kmq2) \
    //     + ((q + k*mu)*(6*k*q + 7*(pow(k,2) + pow(q,2))*mu + 8*k*q*pow(mu,2)))/kpq2))/(k*pow(q,2));

    f = 1./kmq2*((5.*pow(k,2.))/63. - kdotq/54. - (pow(k,2.)*pow(kdotq,2.))/(6.*pow(q,4.)) \
      + (71.*pow(kdotq,3.))/(189.*pow(q,4.)) - (4.*pow(kdotq,4.))/(27.*pow(k,2.)*pow(q,4.))\
      + (47.*pow(k,2.)*kdotq)/(378.*pow(q,2.)) - (163.*pow(kdotq,2.))/(378.*pow(q,2.))\
      + (5.*pow(kdotq,3.))/(27.*pow(k,2.)*pow(q,2.)));

    return f;
}

////This is the symmetrized G3(q,-q,k) given in 1603.04405
double G3_s(double q, double k, double mu)
{

    double f = 0.;
    double kmq2  = pow(q,2.)+pow(k,2.)-2.*q*k*mu;
    double kpq2  = pow(q,2.)+pow(k,2.)+2.*q*k*mu;
    double kdotq = k*q*mu;

    // f = 1./(252.*pow(q,2.)*k) * (-30.*pow(k,3.)*pow(mu,2.) + 2.*pow(q,2.)*k*(6.+pow(mu,2.)) \
    //     + ((q-k*mu) * (-12.*q*pow(k,3.) + 7.*pow(q,4.)*mu + 13.*pow(q*k,2.)*mu) + 6.*pow(k,4.)*mu - 2.*pow(q,3.)*k*(3.+4.*pow(mu,2.)))/kmq2 \
    //     - (q+k*mu)*(12.*q*pow(k,3.) + 7.*pow(q,4)*mu + 13.*pow(q*k,2.) + 6.*pow(k,4.)*mu + 2.*pow(q,3.)*k*(3.+4*pow(mu,2.)))/kpq2);

    // f = (-30*pow(k,3)*pow(mu,2) + 2*k*pow(q,2)*(6 + pow(mu,2)) + ((q - k*mu)*(-12*pow(k,3)*q + 6*pow(k,4)*mu \
    //       + 13*pow(k,2)*pow(q,2)*mu + 7*pow(q,4)*mu - 2*k*pow(q,3)*(3 + 4*pow(mu,2))))/pow(kmq2,2) \
    //       - ((q + k*mu)*(12*pow(k,3)*q + 6*pow(k,4)*mu + 13*pow(k,2)*pow(q,2)*mu + 7*pow(q,4)*mu + 2*k*pow(q,3)*(3 + 4*pow(mu,2))))\
    //       /kpq2)/(252.*k*pow(q,2));

    f = 1./kmq2* (-pow(k,2.)/21. + (11.*kdotq)/126. - (pow(k,2.)*pow(kdotq,2.))/(6.*pow(q,4.))\
      + (17.*pow(kdotq,3.))/(63.*pow(q,4.)) - (4.*pow(kdotq,4.))/(63.*pow(k,2)*pow(q,4.))\
      + (29.*pow(k,2)*kdotq)/(126.*pow(q,2.)) - (7.*pow(kdotq,2.))/(18.*pow(q,2.))\
      + (5.*pow(kdotq,3.))/(63.*pow(k,2.)*pow(q,2.)));

    return f;
}

////This is S2(k1,k2) ////
double S2(double mu)
{
      double f = 0.;
      f = pow(mu,2.) - 1.;

      return f;
}

double F2(double k1,double k2,double mu)
{
      double f=0.;

      f = 5./7. + 1./2. * (k1/k2 + k2/k1) * mu + 2./7. * pow(mu,2.);

      return f;
}

double G2(double k1,double k2,double mu)
{
  // double f = (1./14.)*(6 + mu * (7.*(k1/k2 + k2/k1) + 8.*mu));
  double f = 3./7. + 1./2. * (k1/k2 + k2/k1) * mu + 4./7. * pow(mu,2.);
  return f;
}

// cos -> angle between q and k (gets integrated over), mu -> angle between k and los (fix for all loops)

//(k-q)_p /(|k-q|²)
double LoS1(double q, double k, double cos, double mu)
{
  double f = (2*M_PI*(k - q*cos)*mu)/(pow(k,2) + pow(q,2) - 2*k*q*cos);
  return f;
}

//(q)_p /(|q|²)
double LoS2(double q, double k, double cos, double mu)
{
  double f = (2*M_PI*cos*mu)/q;
  return f;
}

//(k-q)_p q_p /(|k-q|² |q|²)
double LoS3(double q, double k, double cos, double mu)
{
  double f = (M_PI*q*(-1 + pow(cos,2)) + M_PI*(q + 2*k*cos - 3*q*pow(cos,2))*pow(mu,2))/(q*(pow(k,2) + pow(q,2) - 2*k*q*cos));
  return f;
}

//q_p /(|q|²) * ((k-q)_p /(|k-q|²) + q_p /(|q|²))
double LoS4(double q, double k, double cos, double mu)
{
  double f = (k*M_PI*(-((k - 2*q*cos)*(-1 + pow(cos,2))) \
            + (-k + 4*q*cos + 3*k*pow(cos,2) - 6*q*pow(cos,3))*pow(mu,2)))/(pow(q,2)*(pow(k,2) + pow(q,2) - 2*k*q*cos));
  return f;
}

//(k-q)_p /(|k-q|²) (q_p /(|q|²))²
double LoS5(double q, double k, double cos, double mu)
{
  double f = (-(M_PI*(k - 3*q*cos)*(-1 + pow(cos,2))*mu) \
              + M_PI*(q*cos*(3 - 5*pow(cos,2)) + k*(-1 + 3*pow(cos,2)))*pow(mu,3))/(pow(q,2)*(pow(k,2) + pow(q,2) - 2*k*q*cos));
  return f;
}

//((k-q)_p q_p /(|k-q|² |q|²))²
double LoS6(double q, double k, double cos, double mu)
{
  double f = (M_PI*(3*pow(q,2)*pow(-1 + pow(cos,2),2) - 2*(-1 + pow(cos,2))*(2*pow(k,2) - 12*k*q*cos + 3*pow(q,2)*(-1 + 5*pow(cos,2)))*pow(mu,2) \
            + (8*k*q*cos*(3 - 5*pow(cos,2)) + 4*pow(k,2)*(-1 + 3*pow(cos,2)) + pow(q,2)*(3 - 30*pow(cos,2) + 35*pow(cos,4)))*pow(mu,4)))/(4.*pow(q,2)*pow(pow(k,2) + pow(q,2) - 2*k*q*cos,2));
  return f;
}