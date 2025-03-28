import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import pi

                              # get reduced Hubble for conversions to 1/Mpc

# galaxy bias:
b1, b2, bG2 = 1., 0., 0.
# b1, b2, bG2 = 2., -1., 0.1

# noise bias:
d1, d2, P_eps = 1., 1., 1e3
#d1, d2, P_eps = 1., 0.1, 1e3


k_size      = 30
cos_size    = 20
mu_size     = 10
z_size      = 1

k_min = 0.001
k_max = 0.1

mu_min = -0.5
mu_max = 0.5

cos_min = 0
cos_max = 1

z = 0.

print("Full cos loop")
print(f"Iterations: {k_size**2*cos_size*mu_size**2*z_size*1e-6:.2f} * 10^6")
print(f"Memory: {k_size**2*cos_size*mu_size**2*z_size*8/(1024)**2:.2f} MB")
print("\nFull k loop")
print(f"Iterations: {k_size**3*mu_size**2*z_size*1e-6:.2f} * 10^6")
print(f"Memory: {k_size*(k_size+1)*(k_size+2)*mu_size*(mu_size+1)*z_size*8/(1024)**2:.2f} MB")
print("\nOptimized k loop")
print(f"Iterations: {k_size*(k_size+1)*(k_size+2)*mu_size*(mu_size+1)*z_size/12*1e-6:.2f} * 10^6")
print(f"Memory: {k_size*(k_size+1)*(k_size+2)*mu_size*(mu_size+1)*z_size/12*8/(1024)**2:.2f} MB")
print("\nTypes loop")
print(f"Iterations: {k_size*mu_size**2*z_size*1e-4:.2f} * 10^4")
print(f"Memory: {k_size*mu_size**2*z_size/12*8/(1024)**2:.3f} MB")


k = np.logspace(-3, 0, k_size)
# np.linspace(0.001, 0.05, k_size)
cos = np.linspace(-0.8, 0.9, cos_size)
mu = np.linspace(0., 0.9, mu_size)
z = np.linspace(0., 1., mu_size)



dcos = (cos_max-cos_min) / cos_size
dk = (k_max-k_min)/k_size
dlog_k = np.log10(k_max/k_min) / k_size
dmu = (mu_max-mu_min) / mu_size



counter = 0
k1 = []
k2 = []
k3 = []
k1_full = []
k2_full = []
k3_full = []
k1_angle = []
k2_angle = []
k3_angle = []
k1_new = []
k2_new = []
k3_new = []
for index_k1 in range(k_size):
    k1_temp = k_min*10**(dlog_k * index_k1)
    for index_k2 in range(index_k1+1):
        k2_temp = k_min*10**(dlog_k * index_k2)
        for index_k3 in range(index_k2+1):
            k3_temp = k_min*10**(dlog_k * index_k3)

            if k3_temp + k2_temp > k1_temp:
                k1.append(k1_temp)
                k2.append(k2_temp)
                k3.append(k3_temp)
            else:
                k1_full.append(k1_temp)
                k2_full.append(k2_temp)
                k3_full.append(k3_temp)

for index_k1 in range(k_size):
    k1_temp = k_min*10**(dlog_k * index_k1)
    for index_k2 in range(index_k1+1):
        k2_temp = k_min*10**(dlog_k * index_k2)

        for cos_idx in range(cos_size):
            cos12 = cos_min + cos_idx*dcos
            k3_temp = np.sqrt(k1_temp**2 + k2_temp**2 - 2*cos12*k1_temp*k2_temp)
            k1_angle.append(k1_temp)
            k2_angle.append(k2_temp)
            k3_angle.append(k3_temp)

for index_diag in range(k_size):
    k1_temp = k_min*10**(dlog_k * index_diag)
    #k2_temp = k_min + index_diag*dk
    for index_k3 in range(index_diag+1):
        k3_temp = k_min*10**(dlog_k * index_k3)
        # adjust this loop range
        k2_limit = -np.abs(k3_temp - k1_temp/2) + k1_temp/2
        index_k2_max = int((k2_limit - k_min)/dk)
        for index_k2 in range(index_k2_max+1):
            k2_temp = k1_temp - index_k2*dk
            k1_new.append(k1_temp)
            k2_new.append(k2_temp)
            k3_new.append(k3_temp)

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.scatter(k1, k2, k3, c='b', marker='o', alpha=0.7)  # Scatter plot with blue points
#ax.scatter(k1_full, k2_full, k3_full, c='r', marker='.', alpha=0.7)  # Scatter plot with blue points
#ax.scatter(k1_angle, k2_angle, k3_angle, c='black', marker='o', alpha=0.7)  # Scatter plot with blue points

ax.scatter(k1_new, k2_new, k3_new, c='black', marker='o', alpha=0.7)  # Scatter plot with blue points


# Customize the plot
#ax.set_xlim(k_min, k_max)
#ax.set_ylim(k_min, k_max)
#ax.set_zlim(k_min, k_max)
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')

plt.show()








# Old stuff from fourier.c

"""

int fourier_bispectrum_multipoles_treelevel_at_k_and_z(
                                            struct background *pba,
                                            struct primordial * ppm,
                                            struct fourier * pfo,
                                            double b1,
                                            double b2,
                                            double bG2,
                                            double d1,
                                            double d2,
                                            double P_eps,
                                            double k1,
                                            double k2,
                                            double k3,
                                            double cos12,
                                            int l,
                                            int angle_size,
                                            double z,
                                            int use_IR_resum,
                                            int index_pk,
                                            double * B_l
                                            ) {

  int last_index = angle_size-1;
  int idx, idx1, idx2;
  double mu1, mu2, mu3, cos_phi, sin12, d_mu, d_phi, B_l_temp, ln_bispectrum_treelevel;
  double *P_l;
  double **B_ggg;

  class_alloc(P_l, angle_size*sizeof(double), pfo->error_message);
  class_alloc(B_ggg, angle_size*sizeof(double*), pfo->error_message);
  for (idx=0; idx<angle_size; idx++){
    class_alloc(B_ggg[idx], angle_size*sizeof(double), pfo->error_message);
  }

  // the function contains k3 and cos12 as arguments. if k3 is zero, then the angle between k1 and k2 is used (cos12)
  // Note: cos12 is not the angle between the vectors k1 and k2, but the angle within the triangle
  if (k3 == 0.) {
    k3 = sqrt(k1*k1 + k2*k2 - 2.*k1*k2*cos12);
  } else {
    cos12 = 0.5*(k1*k1 + k2*k2 - k3*k3);
  }
  sin12 = sqrt(1-cos12*cos12); // note that sin12 > 0 since the angle must be between 0 and pi

  d_mu = 2./(angle_size-1);
  d_phi = 2.*_PI_/(angle_size-1);


  // loop over both angles and fill the bispectrum for those and fixed k values
  for (idx1=0; idx1<angle_size; idx1++){
    mu1 = -1. + d_mu * idx1;    // since we integrate over mu1, this sampling is fine
    

    // define Legendre polynomials for different (even) multipoles up to l=6
    if (l==0) {
      P_l[idx1] = 1.;
    } else if (l==2) {
      P_l[idx1] = 0.5*(3.*mu1*mu1 - 1.);
    } else if (l==4) {
      P_l[idx1] = 0.125*(35.*mu1*mu1*mu1*mu1 - 30.*mu1*mu1 + 3.);
    } else if (l==6) {
      P_l[idx1] = 1./16.*(231.*mu1*mu1*mu1*mu1*mu1*mu1 - 315.*mu1*mu1*mu1*mu1 + 105.*mu1*mu1 - 5.);
    } else {
        class_test(l!=0 && l!=2 && l!=4 && l!=6, pfo->error_message, "invalid multipole, you have l=%d",l);
    }

    for (idx2=0; idx2<angle_size; idx2++) {
      cos_phi = cos(idx2*d_phi);  // since we integrate over dphi, we need to sample the cosine like this
      mu2 = - mu1*cos12 - sqrt(1-mu1*mu1)*sin12*cos_phi;

      // Note: here we set cos12=0. because we already calculated k3. If we would provide both to the function, we would obtain an error
      class_call(fourier_bispectrum_treelevel_at_k_and_z(pba, ppm, pfo, b1, b2, bG2, d1, d2, P_eps, k1, k2, k3, 0., mu1, mu2, z, use_IR_resum, index_pk, &ln_bispectrum_treelevel), 
                                                        pfo->error_message,
                                                        pfo->error_message);
      // TODO: we do not need log B, but B -> issues with the sign
      B_ggg[idx1][idx2] = ln_bispectrum_treelevel; // exp(ln_bispectrum_treelevel);
    }
  }

  // Legendre Gauss Quadrature:
  int n=5;
  double x_i[] = {-0.906179845938664, -0.538469310105683, 0., 0.538469310105683, 0.906179845938664};        // evaluation points
  double w_i[] = {0.236926885055189, 0.478628670499366, 128./225., 0.478628670499366, 0.236926885055189};   // weights

  B_l_temp = 0;
  for (idx1=0; idx1<angle_size-1; idx1++){
    phi = _PI_*(x_i[idx1]+1);
    for (idx2=0; idx2<angle_size-1; idx2++){
      mu1 = x_i[idx2];
      mu2 = - mu1*cos12 - sqrt(1-mu1*mu1)*sin12*cos(phi);
      if (l==0) {
        P_l = 1.;
      } else if (l==2) {
        P_l = 0.5*(3.*mu1*mu1 - 1.);
      } else if (l==4) {
        P_l = 0.125*(35.*mu1*mu1*mu1*mu1 - 30.*mu1*mu1 + 3.);
      } else if (l==6) {
        P_l = 1./16.*(231.*mu1*mu1*mu1*mu1*mu1*mu1 - 315.*mu1*mu1*mu1*mu1 + 105.*mu1*mu1 - 5.);
      } else {
          class_test(l!=0 && l!=2 && l!=4 && l!=6, pfo->error_message, "invalid multipole, you have l=%d",l);
      }

      class_call(fourier_bispectrum_treelevel_at_k_and_z(pba, ppm, pfo, b1, b2, bG2, d1, d2, P_eps, k1, k2, k3, 0., mu1, mu2, z, use_IR_resum, index_pk, &ln_bispectrum_treelevel), 
                                                        pfo->error_message,
                                                        pfo->error_message);
      B_l_temp += w_i[idx1]*w_i[idx2]*P_l*ln_bispectrum_treelevel;
    }
  }
  B_l = B_l_temp*0.25*(2*l+1);











  // 2D Trapezoid Rule // TODO: Check if this is correct // TODO: add the Legendre Polynomial

  // Corner terms:
  B_l_temp = (B_ggg[0][0]+B_ggg[0][last_index])*P_l[0]  + (B_ggg[last_index][0]+B_ggg[last_index][last_index])*P_l[last_index];

  // Edge Terms:
  for (idx=0; idx<angle_size-1; idx++){
    B_l_temp += 2.*(B_ggg[idx][0] + B_ggg[idx][last_index])*P_l[idx] + 2.*B_ggg[0][idx]*P_l[0] + 2.*B_ggg[last_index][idx]*P_l[last_index];
  }

  // Interior Terms:
  for (idx1=0; idx1<angle_size-1; idx1++){
    for (idx2=0; idx2<angle_size-1; idx2++){
      B_l_temp += 4.*B_ggg[idx1][idx2]*P_l[idx1];
    }
  }

  *B_l = 0.25*d_phi*d_mu * B_l_temp/(2.*_PI_); // Note: the 1/2pi comes from the integration measure from d_phi/2pi

  // free all the memory
  free(P_l);
  for (idx=0; idx<angle_size; idx++){
    free(B_ggg[idx]);
  }
  free(B_ggg);

  return _SUCCESS_;
}

"""