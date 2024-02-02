/** #if DIRECT_INTEGRATION : make will only compile this file if DIRECT_INTEGRATION=yes is given */

#include "direct_integration.h"

// static int indexed_double_cmp(const void * a, const void * b) {
//   struct indexed_double * a_id = (struct indexed_double *)a;
//   struct indexed_double * b_id = (struct indexed_double *)b;
//   if ((*a_id).value < (*b_id).value)
//       return -1;
//   else if ((*a_id).value > (*b_id).value)
//       return 1;
//   else
//       return 0;
// }

// /** non-tensor argument: (P_lin^IR,RSD)_i = P_lin^IR,RSD(k_i, mu_i) with mu_size = k_size */
// int eft_linear_spectrum_rsd_sort(
//                   struct background * pba,
//                   struct primordial * ppm,
//                   struct fourier * pfo,
//                   struct eft * peft,
//                   enum linear_or_logarithmic mode,
//                   const double * const ln_kvec,
//                   const double * const muvec,
//                   const int vec_size,
//                   const int index_pk_type,
//                   double * out_pk) {

//   int index_k;
//   double ln_kvec_sorted[vec_size];
//   struct indexed_double ln_kvec_indexed[vec_size];

//   /** - sort the input k-values */
//   // class_alloc(ln_kvec_indexed, kvec_size*sizeof(struct indexed_double), peft->error_message);
//   for (index_k = 0; index_k < kvec_size; index_k++) {
//     ln_kvec_indexed[index_k].value = ln_kvec[index_k];
//     ln_kvec_indexed[index_k].index = index_k;
//   }
//   qsort(ln_kvec_indexed, kvec_size, sizeof(struct indexed_double), indexed_double_cmp);
//   for (index_k = 0; index_k < kvec_size; index_k++) {
//     ln_kvec_sorted[index_k] = ln_kvec_indexed[index_k].value;
//   }
//   /** - mu does not require sorting since it is not interpolated (dependence is analytical) */
  

//   switch (index_pk_type)
//   {
//   case pkmu_rsd_ir_resummed_lo:
//     class_call(eft_ir_pk_rsd_lo(pba, ppm, pfo, mode,
//                                 ln_kvec_sorted,
//                                 kvec_size,
//                                 peft->z0,
//                                 peft->f_z0,
//                                 muvec,
//                                 muvec_size,
//                                 peft->Sigma2_ir,
//                                 peft->dSigma2_ir,
//                                 out_pk),
//                         pfo->error_message,
//                         peft->error_message);
//     break;

//   case pkmu_rsd_ir_resummed_nlo:
//     class_call(eft_ir_pk_rsd_nlo(pba, ppm, pfo, mode,
//                                  ln_kvec,
//                                  kvec_size,
//                                  peft->z0,
//                                  peft->f_z0,
//                                  muvec,
//                                  muvec_size,
//                                  peft->Sigma2_ir,
//                                  peft->dSigma2_ir,
//                                  out_pk),
//                         pfo->error_message,
//                         peft->error_message);
//     break;
  
//   default:
//     class_stop(pfo->error_message, "index_pk_type = %d not recognized.", index_pk_type);
//     break;
//   }

//   return _SUCCESS_;
// }


// integrands should return -999 for immediate cancellation