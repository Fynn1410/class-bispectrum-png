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


id_z = 1
z = np.array([1., 2.])
k_grid = np.geomspace(1.e-3, 1., num=100)
ln_k_grid = np.log(k_grid)
mu_grid = np.linspace(0., 1., num=5, endpoint=True)
mu_of_z = np.broadcast_to(mu_grid, (len(z), len(mu_grid)))
k_of_z_mu = np.broadcast_to(k_grid, (len(z), len(mu_grid), len(k_grid)))

bias = np.array([[2., -0.7, 0.4, 0.1], [1.5, -0.4, 0.2, 0.05]])
counter = np.array([[0., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0., 0.]])

pk_lin = np.empty((len(z), len(k_grid)), dtype='float64')
for index_z, zz in enumerate(z):
  for index_k, kk in enumerate(k_grid):
    pk_lin[index_z, index_k] = cosmo.pk_lin(kk, zz)

pkmu_rsd = cosmo.eft_pkmu_rsd_grid(mu_of_z, k_of_z_mu, z, bias, counter, 'Pdd_hm_rsd')

pkmu_rsd_spline = cosmo.eft_pkmu_rsd_spline(mu_of_z, z, bias, counter, 'Pdd_hm_rsd')
pkmu_rsd_interp = np.moveaxis(pkmu_rsd_spline(ln_k_grid), 0, 2)

ln_k_spline = pkmu_rsd_spline.x
pkmu_rsd_interp_spline_grid = np.moveaxis(pkmu_rsd_spline(ln_k_spline), 0, 2)
k_spline = np.exp(ln_k_spline)

k_spl_of_z_mu = np.broadcast_to(k_spline, (len(z), len(mu_grid), len(k_spline)))
pkmu_rsd_spline_grid = cosmo.eft_pkmu_rsd_grid(mu_of_z, k_spl_of_z_mu, z, bias, counter, 'Pdd_hm_rsd')

#print(pkmu_rsd_spline_grid[id_z, 0, 0])
#print(pkmu_rsd_spline.c[:, 0, id_z, 0])

plt.figure(1)
plt.loglog(k_grid, bias[id_z, 0]*pk_lin[id_z], label='$b_1 \cdot$linear')
plt.loglog(k_grid, pkmu_rsd[id_z, 0], label='1-loop real')
plt.loglog(k_grid, pkmu_rsd[id_z, -1], label='1-loop $\mu = {0:.2f}$'.format(mu_of_z[0, -1]))
plt.loglog(k_grid, pkmu_rsd_interp[id_z, 0], ls='--', label='1-loop real spline')
plt.loglog(k_spline, pkmu_rsd_interp_spline_grid[id_z, 0], ls='', marker='.', ms=2, label='1-loop real spline breakpoints')
plt.loglog(k_grid, pkmu_rsd_interp[id_z, -1], ls='--', label='1-loop $\mu = {0:.2f}$'.format(mu_of_z[0, -1]))
plt.xlim(1.e-3, 1.)
#plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_rsd_spectra_hm.pdf')