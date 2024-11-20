from classy import Class
import numpy as np
import matplotlib.pyplot as plt

settings = {
    'output':'mPk',
    'non_linear':'oneloopPT',
    'has_rsd': 'yes',
    'fourier_verbose': 3,
    'ic':'ad',
    'accurate_lensing':1,
    'T_cmb':2.7255,
    'Omega_b':0.04886698938,
    'N_eff':3.046,
    'Omega_cdm':0.25929733618, # 0.11923 / h^2
    'YHe':0.24,
    'tau_reio':0.0568,
    'k_pivot':0.002,
    'A_s':2.08e-9,
    'n_s':0.97,
    'P_k_max_1/Mpc': 200.,
    'z_max_pk':2.,
    'z_pk': 0.,
    'h':0.67810,
    'write_warnings':'yes'}

settings.update({'eft_verbose': 3,
                 'eft_compute_loop_matrices':'yes',
                 'eft_write_loop_matrices':'no',
                 'eft_direct_integration':'no'})

cosmo = Class()
cosmo.set(settings)
cosmo.compute()


id_z = 0
z = np.array([0.])

bias = np.array([[2., -0.7, 0.4, 0.1]]) # , [1.5, -0.4, 0.2, 0.05]
counter = np.array([[0., 0., 0., 0., 0., 0., 0.]])  # , [0., 0., 0., 0., 0., 0., 0.]
b1, b2, bG2, btd = bias[id_z]
c00,c10,c22,c32,c20,c30,c42 = counter[id_z]

zout, mu_grid_int, k_grid_int = cosmo.eft_get_output_sampling(z[id_z])
zint = np.array([zout])

k = k_grid_int[0]
ln_k = np.log(k)
mu = np.linspace(0., 1., num=5, endpoint=True)
unit_fin = np.array([1., 0., 0., 0.], dtype='float64')[:, None, None]
unit_UV = np.array([0., 1., 0., 0.], dtype='float64')[:, None, None]
mu_c = mu[None, :, None]
k_c = k[None, None, :]
mu_of_z = np.broadcast_to(mu, (len(z), len(mu)))
k_of_z_mu = np.broadcast_to(k, (len(z), len(mu), len(k)))

pkmu_rsd       = cosmo.eft_pkmu_rsd_grid(mu_of_z, k_of_z_mu, z, bias, counter, 'Pdd_hh_rsd')
pkmu_rsd_no_IR = cosmo.eft_pkmu_rsd_grid(mu_of_z, k_of_z_mu, z, bias, counter, 'Pdd_hh_rsd_no_IR_resummation')
pk_lin_int     = cosmo.eft_pk_linear_real_grid(k_of_z_mu, zint, 'Pk_linear')[id_z, 0]
pk_lin = pk_lin_int[None, None, :]
sigma_v2 = cosmo.eft_sigma_sq(zout, -1, 'Pk_linear')/3.

_, _, _, I2200          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 0)
_, _, _, I1300          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 1)
I1300 *= pk_lin
_, _, _, Idelta200      = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 2)
_, _, _, IG200          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 3)
_, _, _, Idelta2delta200= cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 4)
_, _, _, IG2G200        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 5)
_, _, _, Idelta2G200    = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 6)
_, _, _, FG200          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 7)
FG200 *= pk_lin

_, _, _, I2201          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 8)
_, _, _, Idelta201      = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 9)
_, _, _, IG201          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 10)
_, _, _, J21101         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 11)
J21101 = J21101 * mu_c
_, _, _, Jdelta201      = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 12)
Jdelta201 = Jdelta201 * mu_c
_, _, _, JG201          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 13)
JG201 = JG201 * mu_c
_, _, _, FG201          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 14)
FG201 *= pk_lin
_, _, _, I1301p3101     = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 15)
I1301p3101 *= pk_lin
_, _, _, J12101         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 16)
J12101 = J12101 * pk_lin * mu_c
_, _, _, J11201         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 17)
J11201 = J11201 * pk_lin * mu_c

_, _, _, J21102x        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 18)
_, _, _, J21102y        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 19)
J21102 = J21102x + mu_c**2 * J21102y
_, _, _, Jdelta202x     = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 20)
_, _, _, Jdelta202y     = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 21)
Jdelta202 = Jdelta202x + mu_c**2 * Jdelta202y
_, _, _, JG202x         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 22)
_, _, _, JG202y         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 23)
JG202 = JG202x + mu_c**2 * JG202y
_, _, _, I2211          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 24)
_, _, _, J21111         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 25)
J21111 = J21111 * mu_c
_, _, _, N11x           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 26)
_, _, _, N11y           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 27)
N11 = N11x + mu_c**2 * N11y
_, _, _, J12102x        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 28)
_, _, _, J12102y        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 29)
J12102 = pk_lin * (J12102x + mu_c**2 * J12102y)
_, _, _, I1311          = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 30)
I1311 *= pk_lin
_, _, _, J12111         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 31)
J12111 = J12111 * pk_lin * mu_c
_, _, _, J11211         = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 32)
J11211 = J11211 * pk_lin * mu_c

_, _, _, J21112x        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 33)
_, _, _, J21112y        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 34)
J21112 = J21112x + mu_c**2 * J21112y
_, _, _, N12x           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 35)
_, _, _, N12y           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 36)
N12 = mu_c * (N12x + mu_c**2 * N12y)
_, _, _, J12112x        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 37)
_, _, _, J12112y        = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 38)
J12112 = pk_lin * (J12112x + mu_c**2 * J12112y)

_, _, _, N22x           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 39)
_, _, _, N22y           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 40)
_, _, _, N22z           = cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_linear', 41)
N22 = N22x + mu_c**2 * N22y + mu_c**4 * N22z


D2  = (cosmo.scale_independent_growth_factor(z[id_z]) / cosmo.scale_independent_growth_factor(zout))**2
D4  = D2**2
f_z = cosmo.scale_independent_growth_factor_f(z[id_z])


Prsd0_loop = D4 * (2.*b1*b1 * (I2200 + 3.*I1300) + 2.*b1*b2 * Idelta200 + 2.*bG2*bG2 * IG2G200 + 4.*b1*bG2 * IG200        \
                    + 0.5*b2*b2 * Idelta2delta200 + 2.*b2*bG2 * Idelta2G200 + 8.*b1*(bG2 + 0.4*btd) * FG200)
Prsd1_loop = D4*f_z * (2.*mu_c*mu_c * (2.*b1 * I2201 + 3.*b1 * I1301p3101 + b2 * Idelta201 + 2.*bG2 * IG201 + 4.*(bG2 + 0.4*btd) * FG201)   \
                        + 4.*mu_c*k_c * (b1*b1 * (J12101 + J11201 + J21101) + 0.5*b1*b2 * Jdelta201 + b1*bG2 * JG201))
Prsd2_loop = D4*f_z*f_z * (mu_c*mu_c*k_c*k_c * (2.*b1 * (2.*J12102 + J21102) + b2 * Jdelta202 + 2.*bG2 * JG202)
                            + 2.*mu_c*mu_c*mu_c*mu_c * (I2211 + 3.*I1311) + 4.*mu_c*mu_c*mu_c*k_c * b1 * (J12111 + J11211 + J21111) + mu_c*mu_c*k_c*k_c * b1*b1 * N11);
Prsd3_loop = D4*f_z*f_z*f_z * (2.*mu_c*mu_c*mu_c*mu_c*k_c*k_c * (J21112 + 2.*J12112) + 2.*mu_c*mu_c*mu_c*k_c*k_c*k_c * b1 * N12)
Prsd4_loop = D4*f_z*f_z*f_z*f_z * 0.5*mu_c*mu_c*mu_c*mu_c*k_c*k_c*k_c*k_c * N22

# add counterterm contributions to UV-divergent parts
Prsd0_loop += unit_UV * D2 * c00 * k_c*k_c * pk_lin
Prsd1_loop += unit_UV * D2*f_z * c10 * mu_c*mu_c*k_c*k_c * pk_lin
Prsd2_loop += unit_UV * D2*f_z*f_z * (c20 + c22 * mu_c*mu_c) * mu_c*mu_c*k_c*k_c * pk_lin
Prsd3_loop += unit_UV * D2*f_z*f_z*f_z * (c30 + c32 * mu_c*mu_c) * mu_c*mu_c*mu_c*mu_c*k_c*k_c * pk_lin
Prsd4_loop += unit_UV * D2*f_z*f_z*f_z*f_z * c42 * mu_c*mu_c*mu_c*mu_c*mu_c*mu_c*k_c*k_c * pk_lin

# add sigma_v to finite parts
Prsd2_loop += unit_fin * -D4*f_z*f_z * b1*b1 * sigma_v2 * mu_c*mu_c*k_c*k_c * pk_lin
Prsd3_loop += unit_fin * -2.*D4*f_z*f_z*f_z * b1 * sigma_v2 * mu_c*mu_c*mu_c*mu_c*k_c*k_c * pk_lin
Prsd4_loop += unit_fin * -D4*f_z*f_z*f_z*f_z * sigma_v2 * mu_c*mu_c*mu_c*mu_c*mu_c*mu_c*k_c*k_c * pk_lin

pkmu_loop = Prsd0_loop + Prsd1_loop + Prsd2_loop + Prsd3_loop + Prsd4_loop
pkmu_loop_sum = np.zeros((len(mu), len(k)), dtype='float64')
for index_part in range(3, -1, -1):
    pkmu_loop_sum += pkmu_loop[index_part]

pkmu = D2 * ((b1 + f_z*mu**2)**2)[:, None] * pk_lin_int[None, :] + pkmu_loop_sum

plt.figure(1)
plt.loglog(k, (bias[id_z, 0]**2.)*pk_lin_int*D2, label='$b_1^2 \cdot$linear')
plt.loglog(k, pkmu_rsd[id_z, 0], label='1-loop real')
plt.loglog(k, pkmu_rsd[id_z, -1], label='1-loop $\mu = {0:.2f}$'.format(mu_of_z[0, -1]))
plt.loglog(k, pkmu[0], label='1-loop real pieces')
plt.loglog(k, pkmu[-1], label='1-loop $\mu = {0:.2f}$ pieces'.format(mu[-1]))
#plt.xlim(1.e-3, 1.)
#plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_rsd_spectra_contributions.pdf')

id_mu = -1
plt.figure(2)
plt.loglog(k, pkmu_loop[0, id_mu], label='finite', color='b')
plt.loglog(k, -pkmu_loop[0, id_mu], ls='--', label='finite', color='b')
plt.loglog(k, pkmu_loop[1, id_mu], label='UV divergence', color='r')
plt.loglog(k, -pkmu_loop[1, id_mu], ls='--', label='UV divergence', color='r')
plt.loglog(k, pkmu_loop[2, id_mu], label='IR divergence', color='g')
plt.loglog(k, -pkmu_loop[2, id_mu], ls='--', label='IR divergence', color='g')
plt.loglog(k, pkmu_loop[3, id_mu], label='Pole divergence', color='k')
plt.loglog(k, -pkmu_loop[3, id_mu], ls='--', label='Pole divergence', color='k')
#plt.xlim(1.e-3, 1.)
#plt.ylim(10., 50.e3)
plt.title("$\mu = {0:.2f}$".format(mu[id_mu]))
plt.legend()
plt.savefig('plots/oneloop_rsd_spectra_contributions_parts.pdf')