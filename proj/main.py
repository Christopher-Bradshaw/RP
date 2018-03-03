

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

def two_photon_gamma_nu():

    bottom = 1 + (N_p * q_p_2s_2p + N_e * q_e_2s_2p) / A_2s_1s

    return alpha_eff_2s * g_v / bottom



    pass
