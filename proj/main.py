#!/usr/bin/env python3

import matplotlib.pyplot as plt

import constants as c
from two_photon_spectra import get_two_photon_emission_coefficient

def main():
    np, ne = 1e-1, 1e-1 # emailed brant, we can update this

    fig, ax = plt.subplots()

    lmbda, coeff = get_two_photon_emission_coefficient(np, ne)
    ax.plot(lmbda*1e4, coeff*c.c/lmbda)

    # Load and plot your stuff here. NB the units of microns on the x axis and the weird units on the y axis

    ax.set(yscale="log", xscale="log", xlim=(1e-1, 1e1), ylim=(1e-26, 2e-23), xlabel="Wavelength (microns)", ylabel=r"$\nu \gamma_{\nu}\ (erg\ cm^3\ s^{-1})$")
    ax.annotate(s="T = 10,000K", xy=(0.75, 0.80), xycoords="axes fraction")
    plt.show()


if __name__ == "__main__":
    main()
