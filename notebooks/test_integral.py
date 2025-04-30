# import classy module
from classy import Class
CHOP_TOL = 1.e-12
from testtest import TriaMaster

# create instance of the class "Class"
LambdaCDM = Class()
# pass input parameters
LambdaCDM.set({'omega_b':0.0223828,'omega_cdm':0.1201075,'h':0.67810,'A_s':2.100549e-09,'n_s':0.9660499,'tau_reio':0.05430842})
LambdaCDM.set({'omega_b':0.02237,'omega_cdm':0.1200,'h':0.6736,'A_s':2.0830e-9,'n_s':0.9649,'tau_reio':0.0568, 'k_pivot': 0.05, 'N_ur':3.046, 'T_cmb': 2.7255})
LambdaCDM.set({'output':'mPk','lensing':'no','P_k_max_1/Mpc':3.0, 'z_max_pk': 1.0})

# run class
LambdaCDM.compute()

# get reduced Hubble for conversions to 1/Mpc
h = LambdaCDM.h()
#print(h)
#h = 2/3
# uncomment to get plots displayed in notebook

import matplotlib.pyplot as plt
from math import pi
from matplotlib import rc
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
import numpy as np

import pandas as pd


# Mathematica directory
dir = "/Users/fynn/Desktop/Uni/Master/Masters-thesis/04 Mathematica"

# Enable LaTeX rendering
rcParams['figure.dpi'] = 120
rc('text', usetex=False)
rc('font', family='serif', size=15)

rcParams.update({
    'font.size': 15,
    'axes.titlesize': 15,
    'axes.labelsize': 15,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 15,
})

def compare_integrals(x1, x2, y1real, y1imag, y2real, y2imag, ylabel):
    factor_real = np.mean(y1real/y2real)
    factor_imag = np.mean(y1imag/y2imag)
    print("factor=", factor_real, factor_imag)


    residuals_real = 100 * (y1real-y2real)/y1real
    residuals_imag = 100 * (y1imag-y2imag)/y1imag

    fig = plt.figure(figsize=(8, 6))
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.1)

    # Main plot (data + fit)
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(x1, np.abs(y1real), label="real analytical")
    ax1.plot(x2, np.abs(y2real), ls="--", label="real mathematica numerical")
    ax1.plot(x1, np.abs(y1imag), label="imag analytical")
    ax1.plot(x2, np.abs(y2imag), ls="--", label="imag mathematica numerical")
    ax1.set_yscale("log")
    ax1.set_title(rf"{ylabel}")
    ax1.legend()
    ax1.grid()
    ax1.tick_params(labelbottom=False)  # Hide x-axis labels for top plot

    # Residual plot (same x-axis)
    ax2 = fig.add_subplot(gs[1], sharex=ax1)  # Shares x-axis with ax1
    ax2.axhline(0, color='black', linestyle='--', linewidth=1)
    ax2.plot(x1, (residuals_real), color='black', label="Residuals Real")
    ax2.plot(x1, (residuals_imag), color='black', label="Residuals Imag")
    ax2.set_xlabel(rf"k $\left[h/\mathrm{{Mpc}}\right]$")
    ax2.set_ylabel(r"Res. $\%$")
    #ax2.set_xscale("log")
    ax2.legend()
    ax2.grid()
    plt.show()

k2_peak_n = np.array([0., -3.4e-2, -1e-3, -7.6e-5, -1.56e-5, -3.4e-2, -1e-3, -7.6e-5, -1.56e-5, 0.])*h**2
k2_UV_n = np.array([1e-4, 6.9e-2, 8.2e-3, 1.3e-3, 1.35e-5, -6.9e-2, -8.2e-3, -1.3e-3, -1.35e-5, 0.])*h**2 # original 

M_n = -k2_peak_n + 1j*k2_UV_n
M_n[0] = 1e-4 + 1j*0.

cos12 = -0.5


# read Mathematica file
T_master_mathematica = np.array(pd.read_csv(dir+"/T_master_integral/T_master_mathematica_singleindex.csv", names=["M1", "M2", "M3", "k1", "k2", "Real", "Imag"]))

M1_arr = T_master_mathematica[:,0]
M2_arr = T_master_mathematica[:,1]
M3_arr = T_master_mathematica[:,2]
k1_math = T_master_mathematica[:,3]
k2_math = T_master_mathematica[:,4]
T_math_real = T_master_mathematica[:,5]
T_math_imag = T_master_mathematica[:,6]
if isinstance(T_math_imag[0], str):
    for idx, temp in enumerate(T_math_imag):
        T_math_imag[idx] =  eval(temp.replace("*^", "**"))


k_len = len(k1_math)
M1_len = len(M1_arr)
M2_len = len(M2_arr)
M3_len = len(M3_arr)



T_master_arr_CLASS = np.zeros(k_len, dtype=complex)
T_master_arr_Babis = np.zeros(k_len, dtype=complex)
x = np.arange(k_len)

for idx in range(k_len):

    k1_in = k1_math[idx]
    k2_in = k2_math[idx]
    k3_in_sq = k1_in**2 + k2_in**2 + 2*k1_in*k2_in*cos12
    M1_idx = int(M1_arr[idx])
    M2_idx = int(M2_arr[idx])
    M3_idx = int(M3_arr[idx])
    T_master_arr_CLASS[idx] = LambdaCDM.get_T_master(k1_in**2, k2_in**2, k3_in_sq, M_n[M1_idx], M_n[M2_idx], M_n[M3_idx])
    T_master_arr_Babis[idx] = TriaMaster(k1_in**2, k2_in**2, k3_in_sq, M_n[M1_idx], M_n[M2_idx], M_n[M3_idx], CHOP_TOL)




x2 = k1_math/h
y0real = T_master_arr_Babis.real + 1.e-8
y0imag = T_master_arr_Babis.imag + 1.e-8
y1real = T_master_arr_CLASS.real + 1.e-8
y1imag = T_master_arr_CLASS.imag + 1.e-8
y2real = T_math_real + 1.e-8
y2imag = T_math_imag + 1.e-8
ylabel = "T_master"

y1real = np.array(y1real)
y1imag = np.array(y1imag)
y2real = np.array(y2real)
y2imag = np.array(y2imag)


cut = 1.

mask = np.abs(y1imag/y2imag-1) < cut
#mask = (M1_arr == 0) & (M2_arr == 0) & (M3_arr == 0)


for idx in range(k_len):
    print("")
    print("")
    print(f"Real: M1 = {M1_arr[idx]}, M2 = {M2_arr[idx]}, M3 = {M3_arr[idx]}, k1 = {k1_math[idx]}, k2 = {k2_math[idx]}, CLASS = {y1real[idx]}, Babis = {y0real[idx]}, Mathematica = {y2real[idx]}, chi [%] = {np.abs(y1real[idx]/y2real[idx]-1)*100}, chi [%] = {np.abs(y0real[idx]/y2real[idx]-1)*100}")
    print(f"Imag: M1 = {M1_arr[idx]}, M2 = {M2_arr[idx]}, M3 = {M3_arr[idx]}, k1 = {k1_math[idx]}, k2 = {k2_math[idx]}, CLASS = {y1imag[idx]}, Babis = {y0imag[idx]}, Mathematica = {y2imag[idx]}, chi [%] = = {np.abs(y1imag[idx]/y2imag[idx]-1)*100}, chi [%] = = {np.abs(y0imag[idx]/y2imag[idx]-1)*100}")
    




M_len = 9*9*9
k_len = 3*3*3
tot_len = M_len*k_len
T_master_arr_CLASS = np.zeros(tot_len, dtype=complex)
T_master_arr_Babis = np.zeros(tot_len, dtype=complex)

k_arr = np.array([0.1])*h

idx = 0
for k1_in in k_arr:
    for k2_in in k_arr: 
        for k3_in in k_arr: 
            for M1_idx in range(9):
                for M2_idx in range(9):
                    for M3_idx in range(9):
                        T_master_arr_CLASS[idx] = LambdaCDM.get_T_master(k1_in**2, k2_in**2, k3_in**2, M_n[M1_idx], M_n[M2_idx], M_n[M3_idx])
                        T_master_arr_Babis[idx] = TriaMaster(k1_in**2, k2_in**2, k3_in**2, M_n[M1_idx], M_n[M2_idx], M_n[M3_idx], CHOP_TOL)
                        if idx==338:
                            print(k1_in, k2_in, k3_in)
                            print(M1_idx, M2_idx, M3_idx)
                            print(f"Real: M1 = {M1_idx}, M2 = {M2_idx}, M3 = {M3_idx}, k1 = {k1_in}, k2 = {k2_in}, cos12={(k3_in**2 - k1_in**2 - k2_in**2)/(2*k1_in*k2_in)}, CLASS = {T_master_arr_CLASS[idx].real}, Babis = {T_master_arr_Babis[idx].real}, chi [%] = {np.abs(T_master_arr_Babis[idx].real/T_master_arr_CLASS[idx].real-1)*100}")
                            print(f"Real: M1 = {M1_idx}, M2 = {M2_idx}, M3 = {M3_idx}, k1 = {k1_in}, k2 = {k2_in}, cos12={(k3_in**2 - k1_in**2 - k2_in**2)/(2*k1_in*k2_in)}, CLASS = {T_master_arr_CLASS[idx].imag}, Babis = {T_master_arr_Babis[idx].imag}, chi [%] = = {np.abs(T_master_arr_Babis[idx].imag/T_master_arr_CLASS[idx].imag-1)*100}")
                            

                            print("")
                        idx += 1

#x1 = k/h
x2 = k1_math/h
y0real = T_master_arr_Babis.real + 1.e-8
y0imag = T_master_arr_Babis.imag + 1.e-8
y1real = T_master_arr_CLASS.real + 1.e-8
y1imag = T_master_arr_CLASS.imag + 1.e-8
y2real = T_math_real + 1.e-8
y2imag = T_math_imag + 1.e-8
ylabel = "T_master"

y1real = np.array(y1real)
y1imag = np.array(y1imag)
y2real = np.array(y2real)
y2imag = np.array(y2imag)


cut = 1.

mask = np.abs(y1imag/y2imag-1) < cut
#mask = (M1_arr == 0) & (M2_arr == 0) & (M3_arr == 0)


x = np.arange(tot_len)
compare_integrals(x, x, y1real, y1imag, y0real, y0imag, ylabel)

#for idx in range(k_len):
#    print("")
#    print("")
#    print(f"Real: M1 = {M1_arr[idx]}, M2 = {M2_arr[idx]}, M3 = {M3_arr[idx]}, k1 = {k1_math[idx]}, k2 = {k2_math[idx]}, CLASS = {y1real[idx]}, Babis = {y0real[idx]}, Mathematica = {y2real[idx]}, chi [%] = {np.abs(y1real[idx]/y2real[idx]-1)*100}, chi [%] = {np.abs(y0real[idx]/y2real[idx]-1)*100}")
#    print(f"Imag: M1 = {M1_arr[idx]}, M2 = {M2_arr[idx]}, M3 = {M3_arr[idx]}, k1 = {k1_math[idx]}, k2 = {k2_math[idx]}, CLASS = {y1imag[idx]}, Babis = {y0imag[idx]}, Mathematica = {y2imag[idx]}, chi [%] = = {np.abs(y1imag[idx]/y2imag[idx]-1)*100}, chi [%] = = {np.abs(y0imag[idx]/y2imag[idx]-1)*100}")
    

##

# Real: M1 = 7.0, M2 = 3.0, M3 = 4.0, k1 = 0.03, k2 = 0.015, CLASS = 1149129.326369968, Mathematica = 349707.5734816661, chi [%] = 228.59720907079887
# Imag: M1 = 7.0, M2 = 3.0, M3 = 4.0, k1 = 0.03, k2 = 0.015, CLASS = 2510405.4934871565, Mathematica = -128629.92741454949, chi [%] = = 2051.6496230279313

