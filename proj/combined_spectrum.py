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

nus_for_j_arr, j_nus = [], [] # the nus used for the js and the js at those nus
nus_for_a_arr, a_nus = [], []

# Two photon
two_photon_spectra.use_lambda = False
t_n, t_j = two_photon_spectra.get_two_photon_jv(n_p, n_e)
nus_for_j_arr.append(t_n[1:])
j_nus.append(t_j[1:])

# Hydrogen Lines
h_nus, h_j = hydrogen_lines.get_line_emission_jv(n_p, n_e)
nus_for_j_arr.append(h_nus)
j_nus.append(h_j)

h_nus, j_a = hydrogen_lines.get_line_emission_av()
nus_for_a_arr.append(h_nus)
a_nus.append(j_a)


# If other things have specific NUs they need on the x axis, add them in here


# Interpolated everything into same nu grid
min_nu = min(np.min(np.hstack(nus_for_j_arr)), np.min(np.hstack(nus_for_a_arr)))
max_nu = max(np.max(np.hstack(nus_for_j_arr)), np.max(np.hstack(nus_for_a_arr)))
output_nus = np.linspace(min_nu, max_nu, num=1000000)

# If they can take an array of nus, do so here.
j_nus.append(free_free.j_nu(output_nus,T, n_p, n_e))
nus_for_j_arr.append(output_nus)
a_nus.append(free_free.alpha_nu(output_nus,T, n_p, n_e))
nus_for_a_arr.append(output_nus)

output_js = np.sum([np.interp(output_nus, nus_for_j_arr[i], j_nus[i]) for i in range(len(nus_for_j_arr))], axis=0)
output_as = np.sum([np.interp(output_nus, nus_for_a_arr[i], a_nus[i]) for i in range(len(nus_for_a_arr))], axis=0)

output_js *= 1e53 # I think we need to multiple by the stellar lum?
# I_v = S_v - 1/a_v e ^ (-a_v d_gal)
d_gal = 3e21 # 1kpc
I_v = output_js / output_as - 1/output_as * np.exp(-output_as * d_gal)

fig, ax = plt.subplots()
ax.set(yscale="log", title="Intensity", xlabel="Nu", ylabel="ergs")
ax.plot(output_nus, I_v)
plt.show()
