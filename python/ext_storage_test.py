# import classy module
from classy import Class
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

# create instance of the class "Class"
LambdaCDM = Class()
# pass input parameters
LambdaCDM.set({'output':'mPk','P_k_max_1/Mpc':200,'z_max_pk':3.5})
LambdaCDM.set({'non_linear':'oneloopPT'})
LambdaCDM.set({'has_rsd':'yes'})
LambdaCDM.set({'eft_compute_loop_matrices':'yes'})
LambdaCDM.set({'eft_enable_mu_approximation':'yes'})
LambdaCDM.set({'eft_write_loop_matrices':'no'})
LambdaCDM.set({'eft_direct_integration':'no'})
LambdaCDM.set({'fourier_verbose': 3})
LambdaCDM.set({'eft_verbose': 3})
# run class
t1_start = perf_counter()
LambdaCDM.compute()
t1 = perf_counter() - t1_start

# vector of z values
z = np.array([0.,1.],'float64')
z_size = z.shape[0]
# number of mu values from 0 to 1
mu_size = 3
# number of k values, k_min and k_max in 1/Mpc
k_size = 100
kmin=0.001071519
kmax=0.9332543
# parameters of oneloop model
biases = np.tile(np.array([1,-0.5,0.3,0.8], dtype='float64'), (z_size, 1))
counterterms = np.tile(np.array([-10,20,20,20,0,0,0], dtype='float64'), (z_size, 1))

mu = np.tile(np.linspace(0., 1., mu_size, endpoint=True, dtype='float64'), (z_size, 1))
k = np.tile(np.geomspace(kmin, kmax, k_size, endpoint=True, dtype='float64'), (z_size, mu_size, 1))
print(mu.shape)
print(k.shape)

tpk1_start = perf_counter()
pkmuz = LambdaCDM.eft_pkmu_rsd_grid(mu, k, z, biases, counterterms, 'Pdd_hh_rsd')
tpk1 = perf_counter() - tpk1_start

print(pkmuz[0,0])

# recompute using ext_storage
LambdaCDM.set({'h': 0.67})
t2_start = perf_counter()
LambdaCDM.compute()
t2 = perf_counter() - t2_start

tpk2_start = perf_counter()
pkmuz = LambdaCDM.eft_pkmu_rsd_grid(mu, k, z, biases, counterterms, 'Pdd_hh_rsd')
tpk2 = perf_counter() - tpk2_start

tpkinter_start = perf_counter()
for i in range(100):
  pkmuz = LambdaCDM.eft_pkmu_rsd_grid(mu, k, z, biases, counterterms, 'Pdd_hh_rsd')
tpkinter = perf_counter() - tpkinter_start

print(pkmuz[0,0])

print("Compute(): Run1: {0:.3f}, Run2: {1:.3f}".format(t1, t2))
print("eft_pkmu_rsd_grid() [{0:d}x{1:d}x{2:d} points]: Run1: {3:.3f}, Run2: {4:.3f}".format(z_size, mu_size, k_size, tpk1, tpk2))
print("eft_pkmu_rsd_grid() bias/counterterm change [{0:d}x{1:d}x{2:d} points]: {3:.6f}".format(z_size, mu_size, k_size, tpkinter/100.))