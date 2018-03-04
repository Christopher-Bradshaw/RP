#!/usr/bin/env python3
import numpy as np

# Collisional transition rates.
# Table 4.10, page 76. cm^3 s^-1
q_p_2s_2p = 1e-4 * (2.51 + 2.23) / 2
q_e_2s_2p = 1e-4 * (0.22 + 0.35) / 2

# Two photon transition probability
# Page 73, s^-1
A_2s_1s = 8.23

# Effective recombination coefficient
# Table 4.11, page 77. cm^3 s^-1
alpha_eff_2s = 0.838e-13

# Spectral distribution
# Table 4.11, page 77. erg Hz^-1
g_v = 1 # this is given in very broad bins. I think we want to find this elsewhere


# Number densities of protons and electrons
N_p = 1
N_e = 1

def two_photon_gamma_nu():

    bottom = 1 + (N_p * q_p_2s_2p + N_e * q_e_2s_2p) / A_2s_1s

    return alpha_eff_2s * g_v / bottom


# I'm not sure what H8 means? I'm guessing the 8th line?
lyman_intensities = np.array([32.7])
balmer_intensities = np.array([2.87, 1, 0.466, 0.256, 0.158])
paschen_intensities = np.array([0.352, 0.354*0.466, 0.354*0.256])

# Table 4.2, page 66.
j_Hb_const=1.24e-25

def j_Hb():
    return j_Hb_const * N_p * N_e / 4 * np.pi
