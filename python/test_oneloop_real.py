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
k_of_z = np.broadcast_to(k_grid, (len(z), len(k_grid)))

bias = np.array([[2., -0.7, 0.4, 0.1], [1.5, -0.4, 0.2, 0.05]])
counter = np.array([[0., 0.], [0., 0.]])

pk_real = cosmo.eft_pk_real_grid(k_of_z, z, bias, counter, 'Pdd_mm_real_no_IR_resummation')
I2200 =   cosmo.eft_pk_real_grid(k_of_z, z, bias, counter, 'Pdd_mm_22_no_IR_resummation')
I1300 =   cosmo.eft_pk_real_grid(k_of_z, z, bias, counter, 'Pdd_mm_13_no_IR_resummation')

pk_lin = np.empty_like(pk_real)
for index_z, zz in enumerate(z):
  for index_k, kk in enumerate(k_grid):
    pk_lin[index_z, index_k] = cosmo.pk_lin(kk, zz)

pk_real_pieces = pk_lin + 2.*(I2200 + 3.*I1300)
print(np.sqrt( np.sum((pk_real_pieces - pk_real)**2.) / float(len(z)*len(k_grid)) ) )

pk_nl, k_nl, z_nl = cosmo.get_pk_and_k_and_z(nonlinear=True)
id_z_nl = np.argmin(np.abs(z_nl - z[id_z]))

plt.figure(1)
plt.loglog(k_grid, pk_lin[id_z], label='linear')
plt.loglog(k_grid, pk_real[id_z], label='1-loop')
plt.loglog(k_grid, pk_real_pieces[id_z], ls='--', label='1-loop pieces')
plt.loglog(k_nl, pk_nl[:, id_z_nl], label='CLASS NL z={0:.3f}'.format(z_nl[id_z_nl]))
plt.xlim(1.e-3, 1.)
plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_real_spectra.pdf')