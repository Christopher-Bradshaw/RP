#!/usr/bin/env python3
import numpy as np


# Spectral distribution
# Table 4.11, page 77. erg Hz^-1
g_v = 1 # this is given in very broad bins. I think we want to find this elsewhere




# I'm not sure what H8 means? I'm guessing the 8th line?
lyman_intensities = np.array([32.7])
balmer_intensities = np.array([2.87, 1, 0.466, 0.256, 0.158])
paschen_intensities = np.array([0.352, 0.354*0.466, 0.354*0.256])

# Table 4.2, page 66.
j_Hb_const=1.24e-25

def j_Hb():
    return j_Hb_const * N_p * N_e / 4 * np.pi
