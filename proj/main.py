#!/usr/bin/env python3

import matplotlib.pyplot as plt

import constants as c
from two_photon_spectra import get_two_photon_emission_coefficient
import boundfree

T_neb = 10000. #[Kelvin] Nebula temperature, maybe make it an argument?

def main():
    np, ne = 1e-1, 1e-1 # emailed brant, we can update this

    fig, ax = plt.subplots()

    lmbda, coeff = get_two_photon_emission_coefficient(np, ne)
    ax.plot(lmbda*1e4, coeff*c.c/lmbda)

    coeff_bf, lmbda_bf, nu_bf = boundfree.get_bound_free_emission_coefs(T_neb, n_levels = 6, print_progress = False)
    ax.plot(lmbda_bf*1e4, coeff_bf*c.c/lmbda_bf)

    #Increasing the number of quantum levels to calculate the bound free spectrum slows the code down a lot. Try to keep
    # at or below 8 for testing
    
    #Same example just all the arguments listed
    #boundfree.get_bound_free_emission_coefs(T, lams = np.linspace(.1,10,1000) * 1e-4, n_levels = 8, print_progress = True)

    # Load and plot your stuff here. NB the units of microns on the x axis and the weird units on the y axis

    ax.set(
            yscale="log", xscale="log", xlim=(1e-1, 1e1), ylim=(1e-26, 2e-23),
            xlabel="Wavelength (microns)", ylabel=r"$\nu \gamma_{\nu}\ (erg\ cm^3\ s^{-1})$",
    )
    ax.annotate(s="T = 10,000K", xy=(0.75, 0.80), xycoords="axes fraction")
    plt.savefig("4.1.png")


if __name__ == "__main__":
    main()
