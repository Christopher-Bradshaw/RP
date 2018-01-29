#!/usr/bin/env python3
import numpy as np
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt

k = 1.381 * 1e-16 # Boltzmann's const in ergs
c = 2.998 * 1e10 # Speed of light in cm/s
h = 6.626 * 1e-27 # Planck const in ergs*s

# Very safe margins so we don't run into numerical issues
n_min, n_max = 1e-50, 1e50

# Calculate the maximum nu value we can have given float precision
def nu_limits(T):
    max_limit = np.log(n_max) * k * T / h
    min_limit = n_min * k * T / h # need the safety margin here as the numerator will be small
    return min_limit, max_limit

def numerically_integrated_planck_function(T):
    min_nu, max_nu = nu_limits(T)
    integral, _ = scipy.integrate.quad(planck_function, min_nu, max_nu, args=(T))
    return integral

# 1.58a in Rybicki
def analytically_integrated_planck_function(T):
    return 2 * (np.pi * k * T)**4 / (15 * c**2 * h**3)

def planck_function(nu, T):
    return 2 * h * (nu**3 / c**2) / (np.exp(h * nu / (k * T)) - 1)

def analytical_derivative_planck_function(nu, T):
    factor = np.exp(h * nu / (k * T))
    return 2 * h * nu**2 / (c**2 * (factor - 1)) * (
            3 - (nu * h * factor) / (k * T * (factor - 1)))

def a():
    # assert our numerical integration is good
    for temp in np.logspace(0, 20):
        assert np.isclose(
                analytically_integrated_planck_function(temp),
                numerically_integrated_planck_function(temp))

    # And plot the planck function for fun
    _, ax = plt.subplots()
    for temp in np.power(10, np.arange(5)):
        nus = np.logspace(4, np.log10(nu_limits(temp)[1]), num=500)
        ax.plot(nus, planck_function(nus, temp), label=temp)
    ax.set(xscale="log", yscale="log")
    ax.set_ylim(bottom=1e-20)
    ax.legend()
    plt.show(block=False)

def b():
    T = 1e4
    min_nu, max_nu = nu_limits(T)
    min_nu = 1 # being reasonable...
    x = scipy.optimize.brentq(
            analytical_derivative_planck_function, min_nu, max_nu, args=(T,),
    )

    print(x)

    # And plot the derivative for fun
    _, ax = plt.subplots()
    for temp in np.power(10, np.arange(5)):
        nus = np.logspace(4, np.log10(nu_limits(temp)[1]), num=500)
        ax.plot(nus, analytical_derivative_planck_function(nus, temp), label=temp)
    ax.set(xscale="log")#, yscale="log")
    ax.legend()
    plt.show(block=False)


if __name__ == "__main__":
    a()
    b()
    input()
