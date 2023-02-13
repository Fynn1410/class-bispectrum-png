#include "header.h"

int eft_init(struct background * pba,
            struct fourier * pfo,
            struct eft * peft,
            double z,
            short compute_rsd_spectrum,
            double k_min,
            int fft_k_points) {

  peft->z0 = z;
  peft->has_rsd = compute_rsd_spectrum;
  peft->k_size_fft = fft_k_points;

  class_call(eft_indices(peft), peft->error_message, peft->error_message);

  return _SUCCESS_;
}


int eft_indices(struct eft * peft) {

  int i = 0;
  /** - define all indices for the different spectra contributions (guaranteed to be consecutive) */

  /** Real space moments */
  class_define_index(peft->index_I2200,           _TRUE_, i, 1);
  class_define_index(peft->index_I1300,           _TRUE_, i, 1);
  class_define_index(peft->index_Idelta200,       _TRUE_, i, 1);
  class_define_index(peft->index_IG200,           _TRUE_, i, 1);
  class_define_index(peft->index_Idelta2delta200, _TRUE_, i, 1);
  class_define_index(peft->index_IG2G200,         _TRUE_, i, 1);
  class_define_index(peft->index_Idelta2G200,     _TRUE_, i, 1);
  class_define_index(peft->index_FG200,           _TRUE_, i, 1);
  
  if (peft->has_rsd) {
    /** - 1-st moment of RSD expansion */
    class_define_index(peft->index_I2201,     _TRUE_, i, 1);
    class_define_index(peft->index_Idelta201, _TRUE_, i, 1);
    class_define_index(peft->index_IG201,     _TRUE_, i, 1);
    class_define_index(peft->index_J21101,    _TRUE_, i, 1);
    class_define_index(peft->index_Jdelta201, _TRUE_, i, 1);
    class_define_index(peft->index_JG201,     _TRUE_, i, 1);
    class_define_index(peft->index_FG201,     _TRUE_, i, 1);
    class_define_index(peft->index_I1301p3101,_TRUE_, i, 1);
    class_define_index(peft->index_J12101,    _TRUE_, i, 1);

    /** - 2-nd moment of RSD expansion */
    class_define_index(peft->index_J21102x,   _TRUE_, i, 1);
    class_define_index(peft->index_J21102y,   _TRUE_, i, 1);
    class_define_index(peft->index_Jdelta202x,_TRUE_, i, 1);
    class_define_index(peft->index_Jdelta202y,_TRUE_, i, 1);
    class_define_index(peft->index_JG202x,    _TRUE_, i, 1);
    class_define_index(peft->index_JG202y,    _TRUE_, i, 1);
    class_define_index(peft->index_I2211,     _TRUE_, i, 1);
    class_define_index(peft->index_J21111,    _TRUE_, i, 1);
    class_define_index(peft->index_N11x,      _TRUE_, i, 1);
    class_define_index(peft->index_N11y,      _TRUE_, i, 1);
    class_define_index(peft->index_J12102x,   _TRUE_, i, 1);
    class_define_index(peft->index_J12102y,   _TRUE_, i, 1);
    class_define_index(peft->index_I1311,     _TRUE_, i, 1);
    class_define_index(peft->index_J12111p11211,_TRUE_, i, 1);
    //class_define_index(peft->index_J11211,    _TRUE_, i, 1);

    /** - 3-rd moment of RSD expansion */
    class_define_index(peft->index_J21112x,   _TRUE_, i, 1);
    class_define_index(peft->index_J21112y,   _TRUE_, i, 1);
    class_define_index(peft->index_N12x,      _TRUE_, i, 1);
    class_define_index(peft->index_N12y,      _TRUE_, i, 1);
    class_define_index(peft->index_J12112x,   _TRUE_, i, 1);
    class_define_index(peft->index_J12112y,   _TRUE_, i, 1);

    /** - 4-th moment of RSD expansion */
    class_define_index(peft->index_N22x,      _TRUE_, i, 1);
    class_define_index(peft->index_N22y,      _TRUE_, i, 1);
    class_define_index(peft->index_N22z,      _TRUE_, i, 1);
  }

  peft->index_num = i;

  /** - store symmetry information */
  class_alloc(peft->symmetry, peft->index_num*sizeof(short), peft->error_message);

  peft->symmetry[peft->index_P22] = mat_symmetric;
  peft->symmetry[peft->index_P13] = vec;

  peft->symmetry[peft->index_Idelta200]       = mat_symmetric;
  peft->symmetry[peft->index_IG200]           = mat_symmetric;
  peft->symmetry[peft->index_Idelta2delta200] = mat_symmetric;
  peft->symmetry[peft->index_IG2G200]         = mat_symmetric;
  peft->symmetry[peft->index_Idelta2G200]     = mat_symmetric;
  peft->symmetry[peft->index_FG200]           = vec;

  peft->symmetry[peft->index_I2201]     = mat_symmetric;
  peft->symmetry[peft->index_I1301p3101]= mat_none;
  peft->symmetry[peft->index_Idelta201] = mat_symmetric;
  peft->symmetry[peft->index_IG201]     = mat_symmetric;
  peft->symmetry[peft->index_FG201]     = vec;
  peft->symmetry[peft->index_J12101]    = vec;
  peft->symmetry[peft->index_J21101]    = mat_none;
  peft->symmetry[peft->index_Jdelta201] = mat_none;
  peft->symmetry[peft->index_JG201]     = mat_none;

  peft->symmetry[peft->index_J12102x]     = vec;
  peft->symmetry[peft->index_J12102y]     = vec;
  peft->symmetry[peft->index_J21102x]     = mat_symmetric;
  peft->symmetry[peft->index_J21102y]     = mat_symmetric;
  peft->symmetry[peft->index_Jdelta202x]  = mat_symmetric;
  peft->symmetry[peft->index_Jdelta202y]  = mat_symmetric;
  peft->symmetry[peft->index_JG202x]      = mat_symmetric;
  peft->symmetry[peft->index_JG202y]      = mat_symmetric;
  peft->symmetry[peft->index_I2211]       = mat_symmetric;
  peft->symmetry[peft->index_I1311]       = vec;
  peft->symmetry[peft->index_J12111p11211]= vec;
  peft->symmetry[peft->index_J21111]      = mat_none;
  peft->symmetry[peft->index_N11x]        = mat_none;
  peft->symmetry[peft->index_N11y]        = mat_none;
  //peft->symmetry[peft->index_J11211]      = vec;

  peft->symmetry[peft->index_J21112x] = mat_symmetric;
  peft->symmetry[peft->index_J21112y] = mat_symmetric;
  peft->symmetry[peft->index_J12112x] = vec;
  peft->symmetry[peft->index_J12112y] = vec;
  peft->symmetry[peft->index_N12x]    = mat_none;
  peft->symmetry[peft->index_N12y]    = mat_none;

  peft->symmetry[peft->index_N22x] = mat_none;
  peft->symmetry[peft->index_N22y] = mat_none;
  peft->symmetry[peft->index_N22z] = mat_none;


  /** - allocate TODO */
  class_alloc(peft->spectra_contributions, index_num*sizeof(double *), peft->error_message);
  class_alloc(peft->fft_matrices, index_num*sizeof(double complex *), peft->error_message);
  class_alloc(peft->fft_matrices_size, index_num*sizeof(int), peft->error_message);


  for (i = 0; i < peft->index_num; i++) {
    class_alloc(peft->spectra_contributions[i], peft->k_size, peft->error_message);

    switch (peft->symmetry[i]) {
    case vec:
      peft->fft_matrices_size[i] = peft->k_size_fft; break;
    case mat_none:
      peft->fft_matrices_size[i] = peft->k_size_fft * peft->k_size_fft; break;
    case mat_symmetric:
      peft->fft_matrices_size[i] = peft->k_size_fft * (peft->k_size_fft + 1) / 2; break;
    default:
      peft->fft_matrices_size[i] = 0; break;
    }

    class_alloc(peft->fft_matrices[i], peft->fft_matrices_size[i]*sizeof(double complex), peft->error_message);
  }

  /** - allocate arrays for different P_lin types, pk_index_num is the last entry in enum eft_pk_type */
  class_alloc(peft->pk_l, pk_index_num*sizeof(double *), peft->error_message);

  for (i = 0; i < pk_index_num; i++) {
    class_alloc(peft->pk_l[i], peft->k_size*sizeof(double), peft->error_message);
  }


  return _SUCCESS_;
}