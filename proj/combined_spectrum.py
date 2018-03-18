#!/usr/bin/env python3
"""
Calculate I_v from a bunch of J_v.
Use the diff eq that Tyler solved

I_v = S_v - 1/a_v e ^ (-a_v d_gal)
get to Flux by integrating over solid angle (size/distance)
"""

import two_photon_spectra
import hydrogen_lines
import matplotlib.pyplot as plt
import numpy as np

from tyler import free_free

n_p, n_e, T = 1, 1, 10000

nus_j, j_nus = [], [] # the nus used for the js and the js at those nus
nus_a, a_nus = [], []

# Two photon
two_photon_spectra.use_lambda = False
t_n, t_j = two_photon_spectra.get_two_photon_emission_coefficient(n_p, n_e)
nus_j.append(t_n[1:])
j_nus.append(t_j[1:])

# Hydrogen Lines
h_nus, h_j = hydrogen_lines.get_line_emission_coeff()
nus_j.append(h_nus)
j_nus.append(h_j)


# If other things have specific NUs they need on the x axis, add them in here




# Interpolated everything into same nu grid
min_nu = np.min(np.hstack(nus_j))#, np.min(np.hstack(nus_a)))
max_nu = np.max(np.hstack(nus_j))#, np.max(np.hstack(nus_a)))
output_nus = np.linspace(min_nu, max_nu, num=1000000)

# If they can take an array of nus, do so here.
j_nus.append(free_free.j_nu(output_nus,T, n_p, n_e))
nus_j.append(output_nus)
a_nus.append(free_free.alpha_nu(output_nus,T, n_p, n_e))
nus_a.append(output_nus)


output_js = np.sum([np.interp(output_nus, nus_j[i], j_nus[i]) for i in range(len(nus_j))], axis=0)

fig, ax = plt.subplots()
ax.set(yscale="log")
ax.plot(output_nus, output_js)
plt.show()
