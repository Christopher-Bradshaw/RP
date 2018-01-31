#!/usr/bin/env python3
import numpy as np
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt

k = 1.381 * 1e-16 # Boltzmann's const in ergs
c = 2.998 * 1e10 # Speed of light in cm/s
h = 6.626 * 1e-27 # Planck const in ergs*s

# Very safe margins so we don't run into numerical issues
n_min, n_max = 1e-150, 1e150

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

def wavelength_planck_function(lmdba, T):
    return 2 * h * (c**2 / lmdba**5) / (np.exp(h * c / (lmdba * k * T)) - 1)

def planck_function(nu, T):
    return 2 * h * (nu**3 / c**2) / (np.exp(h * nu / (k * T)) - 1)

def analytical_derivative_planck_function(nu, T):
    factor = np.exp(h * nu / (k * T))
    return 2 * h * nu**2 / (c**2 * (factor - 1)) * (
            3 - (nu * h * factor) / (k * T * (factor - 1)))

def analytical_temp_derivative_planck_function(nu, T):
    factor = np.exp(h * nu / (k * T))
    return (2 * h**2 * nu**4 / (c**2 * k * T**2)) * factor / (factor - 1)**2

# returns the nu at which the function peaks
def find_planck_peaks(T, f):
    # Working in log space
    # Start low, increase until derivative goes negative
    # Backtrack one step to make sure the previous step didn't straddle the peak
    # Halve the step size and keep going until we have constrained it enough
    lower, step = 1, 0.5
    while True:
        while f(10**(lower + step), T) > f(10**lower, T):
            lower = lower + step
        lower -= step
        step /= 2
        if step < 1e-10:
            return 10**lower

def part_a():
    # assert our numerical integration is good
    for temp in np.logspace(0, 5):
        assert np.isclose(
                analytically_integrated_planck_function(temp),
                numerically_integrated_planck_function(temp),
               rtol=1e-5,
        )

    # And plot the planck function for fun
    _, ax = plt.subplots()
    for temp in np.power(10, np.arange(5)):
        nus = np.logspace(4, np.log10(nu_limits(temp)[1]), num=500)
        ax.plot(nus, planck_function(nus, temp), label=temp)
    ax.set(xscale="log", yscale="log")
    ax.set_ylim(bottom=1e-20)
    ax.legend()
    plt.show(block=False)

def part_b():
    # I start running into numerical issues in the brent func at higher temps...
    # My code runs fine, I just have nothing to compare it to.
    nu_peaks, lambda_peaks, temps = [], [], np.logspace(0, 5)
    for temp in temps:
        min_nu, max_nu = nu_limits(temp)
        min_nu = 1 # being reasonable...

        d_nu_zero = scipy.optimize.brentq(
                analytical_derivative_planck_function, min_nu, max_nu, args=(temp,),
        )
        manual_nu_peak = find_planck_peaks(temp, planck_function)
        assert np.isclose(
                d_nu_zero,
                manual_nu_peak,
               rtol=1e-5,
        )
        nu_peaks.append(planck_function(manual_nu_peak, temp) * (k*temp)/(h*manual_nu_peak))

        # Won't do the checking for lambda...
        manual_lambda_peak = find_planck_peaks(temp, wavelength_planck_function)
        lambda_peaks.append(wavelength_planck_function(manual_lambda_peak, temp) * (manual_lambda_peak*k*temp) / (h*c))

    # And plot the derivative for fun
    _, ax = plt.subplots()
    for temp in np.power(10, np.arange(5)):
        nus = np.logspace(4, np.log10(nu_limits(temp)[1]), num=500)
        ax.plot(nus, analytical_derivative_planck_function(nus, temp), label=temp)
    ax.set(xscale="log")#, yscale="log")
    ax.legend()
    plt.show(block=False)

    # And plot the nu peaks
    _, ax = plt.subplots()
    ax.plot(temps, nu_peaks)
    ax.set(xscale="log", yscale="log")
    plt.show(block=False)

    # And lambda peaks
    _, ax = plt.subplots()
    ax.plot(temps, nu_peaks)
    ax.set(xscale="log", yscale="log")
    plt.show(block=False)

def part_c():
    _, ax = plt.subplots()
    for temp in np.power(10, np.arange(5)):
        nus = np.logspace(4, np.log10(nu_limits(temp)[1]), num=500)
        ax.plot(nus, analytical_temp_derivative_planck_function(nus, temp), label=temp)
    ax.set(xscale="log")#, yscale="log")
    plt.show(block=False)


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    input()
