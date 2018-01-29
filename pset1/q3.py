#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

k = 1.381 * 1e-16 # Boltzmann's const in ergs
c = 2.998 * 1e10 # Speed of light in cm/s
h = 6.626 * 1e-27 # Planck const in ergs*s

def planck_function(nu, T):
    thingy = np.exp(h * nu / (k * T))
    return 2 * h * (nu**3 / c**2) / (thingy - 1)


if __name__ == "__main__":
    print(planck_function(1e11, 1e4))

    nus = np.logspace(4, 20, num=500)
    temps = np.power(10, np.arange(5))

    fig, ax = plt.subplots()
    for temp in temps:
        ax.plot(nus, planck_function(nus, temp))
    ax.set(xscale="log", yscale="log")
    ax.set_ylim(bottom=1e-20)
    plt.show()

