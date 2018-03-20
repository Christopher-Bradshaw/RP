#!/usr/bin/env python3

import matplotlib.pyplot as plt

import constants as c
from two_photon_spectra import get_two_photon_emission_coefficient
import boundfree

from tyler import free_free

T_neb = 10000. #[Kelvin] Nebula temperature, maybe make it an argument?

def main():
    np, ne = 1e-1, 1e-1 # emailed brant, we can update this

    fig, ax = plt.subplots()

    lmbda, coeff = get_two_photon_emission_coefficient(np, ne)
    ax.plot(lmbda*1e4, coeff*c.c/lmbda)

    coeff_bf, lmbda_bf, nu_bf = boundfree.get_bound_free_emission_coefs(T_neb, n_levels = 6, print_progress = False)
    ax.plot(lmbda_bf*1e4, coeff_bf*c.c/lmbda_bf)


    coeff_ff = free_free.gamma_nu(c.c/lmbda_bf, T_neb)
    ax.plot(lmbda_bf*1e4, coeff_ff*c.c / lmbda_bf)

    ax.set(
            yscale="log", xscale="log", xlim=(1e-1, 1e1), ylim=(1e-26, 2e-23),
            xlabel="Wavelength (microns)", ylabel=r"$\nu \gamma_{\nu}\ (erg\ cm^3\ s^{-1})$",
    )
    ax.annotate(s="T = 10,000K", xy=(0.75, 0.80), xycoords="axes fraction")
    plt.savefig("4.1.png")
    plt.show()


if __name__ == "__main__":
    main()
