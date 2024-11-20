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

pk_real = cosmo.eft_pk_real_grid(k_of_z, z, bias, counter, 'Pdd_mm_real')
actual_z1, muvec1, kvec1, bare_I2200 =   cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_IR_resummed_LO', 0)
actual_z2, muvec2, kvec2, bare_I1300 =   cosmo.eft_pkmu_rsd_contribution(z[id_z], 'Pk_IR_resummed_LO', 1)
bare_I2200 = np.squeeze(bare_I2200)
bare_I1300 = np.squeeze(bare_I1300)
pk_lin_lo_int = np.squeeze( cosmo.eft_pk_linear_real_grid(np.broadcast_to(kvec1, (len(z),) + kvec1.shape), z, 'Pk_IR_resummed_LO') )
pk_lin_nlo_int = np.squeeze( cosmo.eft_pk_linear_real_grid(np.broadcast_to(kvec1, (len(z),) + kvec1.shape), z, 'Pk_IR_resummed_NLO') )

k_grid_int = np.squeeze(kvec1)
print("Actual z = {0:.3f}".format(actual_z1))

pk_lin = np.empty_like(pk_real)
growth_D = np.empty_like(z)
for index_z, zz in enumerate(z):
  growth_D[index_z] = cosmo.scale_independent_growth_factor(zz)
  for index_k, kk in enumerate(k_grid):
    pk_lin[index_z, index_k] = cosmo.pk_lin(kk, zz)
growth_D0 = cosmo.scale_independent_growth_factor(actual_z1)

I2200 = ((growth_D/growth_D0)**4)[None, :, None] * bare_I2200[:, None, :]
I1300 = ((growth_D/growth_D0)**2)[None, :, None] * pk_lin_lo_int[None,...] * bare_I1300[:, None, :]

pk_loop_pieces = 2.*(I2200 + 3.*I1300)
pk_real_pieces = pk_lin_nlo_int + np.sum(pk_loop_pieces, axis=0)
#print(np.sqrt( np.sum((pk_real_pieces - pk_real)**2.) / float(len(z)*len(k_grid)) ) )

pk_nl, k_nl, z_nl = cosmo.get_pk_and_k_and_z(nonlinear=True)
id_z_nl = np.argmin(np.abs(z_nl - z[id_z]))

plt.figure(1)
plt.loglog(k_grid, pk_lin[id_z], label='linear')
plt.loglog(k_grid, pk_real[id_z], label='1-loop')
plt.loglog(k_grid_int, pk_lin_lo_int[id_z], ls='--', label='linear internal')
plt.loglog(k_grid_int, pk_real_pieces[id_z], ls='--', label='1-loop pieces')
plt.loglog(k_nl, pk_nl[:, id_z_nl], label='CLASS NL z={0:.3f}'.format(z_nl[id_z_nl]))
plt.xlim(1.e-3, 1.)
plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_real_spectra.pdf')

plt.figure(2)
plt.loglog(k_grid_int, I2200[0, id_z], label='finite')
plt.loglog(k_grid_int, -I2200[0, id_z], ls='--', label='finite')
plt.loglog(k_grid_int, I2200[1, id_z], label='UV divergence')
plt.loglog(k_grid_int, -I2200[1, id_z], ls='--', label='UV divergence')
plt.loglog(k_grid_int, I2200[2, id_z], label='IR divergence')
plt.loglog(k_grid_int, -I2200[2, id_z], ls='--', label='IR divergence')
plt.loglog(k_grid_int, I2200[3, id_z], label='pole divergence')
plt.loglog(k_grid_int, -I2200[3, id_z], ls='--', label='pole divergence')
#plt.xlim(1.e-3, 1.)
#plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_pieces_I2200.pdf')

plt.figure(3)
plt.loglog(k_grid_int, I1300[0, id_z], label='finite')
plt.loglog(k_grid_int, -I1300[0, id_z], ls='--', label='finite')
plt.loglog(k_grid_int, I1300[1, id_z], label='UV divergence')
plt.loglog(k_grid_int, -I1300[1, id_z], ls='--', label='UV divergence')
plt.loglog(k_grid_int, I1300[2, id_z], label='IR divergence')
plt.loglog(k_grid_int, -I1300[2, id_z], ls='--', label='IR divergence')
plt.loglog(k_grid_int, I1300[3, id_z], label='pole divergence')
plt.loglog(k_grid_int, -I1300[3, id_z], ls='--', label='pole divergence')
#plt.xlim(1.e-3, 1.)
#plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_pieces_I1300.pdf')

plt.figure(4)
plt.loglog(k_grid_int, pk_loop_pieces[0, id_z], label='finite')
plt.loglog(k_grid_int, -pk_loop_pieces[0, id_z], ls='--', label='finite')
plt.loglog(k_grid_int, pk_loop_pieces[1, id_z], label='UV divergence')
plt.loglog(k_grid_int, -pk_loop_pieces[1, id_z], ls='--', label='UV divergence')
plt.loglog(k_grid_int, pk_loop_pieces[2, id_z], label='IR divergence')
plt.loglog(k_grid_int, -pk_loop_pieces[2, id_z], ls='--', label='IR divergence')
plt.loglog(k_grid_int, pk_loop_pieces[3, id_z], label='pole divergence')
plt.loglog(k_grid_int, -pk_loop_pieces[3, id_z], ls='--', label='pole divergence')
plt.loglog(k_grid_int, np.sum(pk_loop_pieces[:, id_z], axis=0), label='sum')
plt.loglog(k_grid_int, -np.sum(pk_loop_pieces[:, id_z], axis=0), ls='--', label='sum')
#plt.xlim(1.e-3, 1.)
#plt.ylim(10., 50.e3)
plt.legend()
plt.savefig('plots/oneloop_pieces_pk_loop_pieces.pdf')