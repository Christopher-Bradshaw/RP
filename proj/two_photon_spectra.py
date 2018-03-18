import numpy as np
import matplotlib.pyplot as plt
import constants as c


erg_to_ev = 6.242e+11
use_lambda = True

# Returns a list of nus and the ergs/Hz at those nus. This
# fully (we can linearly interpolate) describes the two photon
# spectrum
def get_two_photon_spectrum_vs_nu():
    eps = 1e-99 # We need this to save us from divide by 0 later
    nus_up = 1e14 * np.array([eps + 1.233 * i for i in range(10)])
    nu_midpoint = 12.34e14
    ergs_hz_up = 1e-27 * np.array([
        0.00, 0.303, 0.978, 1.836, 2.78,
        3.78, 4.80, 5.80, 6.78, 7.74
    ]) # for the increasing nus, the ergs/hz

    nus_down = 2*nu_midpoint - np.copy(nus_up)
    # More energy as their nus are greater
    ergs_hz_down = np.copy(ergs_hz_up) * (nus_down/nus_up)

    nus = np.concatenate((nus_up, [nu_midpoint], nus_down[::-1]))
    ergs_hz = np.concatenate((ergs_hz_up, [8.62e-27], ergs_hz_down[::-1]))

    return nus, ergs_hz

def get_two_photon_spectrum_vs_lambda():
    nus, ergs_hz = get_two_photon_spectrum_vs_nu()
    lambdas = c.c / nus

    # ergs_m = ergs_hz * c.c * np.power(lambdas, -2)
    return lambdas[::-1], ergs_hz[::-1]

def sanity_check_lambda():
    lambdas, ergs_m = get_two_photon_spectrum_vs_lambda()

    tot_ergs = np.trapz(ergs_m[:-1], lambdas[:-1])
    # Ignore the last point because that is at lambda = inf
    print("Total ergs:", tot_ergs)
    print("Total eV: ", tot_ergs * erg_to_ev)

    plt.plot(lambdas[:-1], ergs_m[:-1])
    plt.show()

def sanity_check_nu():
    nus, ergs_hz = get_two_photon_spectrum_vs_nu()

    tot_ergs = np.trapz(ergs_hz, nus)
    print("Total ergs:", tot_ergs)
    print("Total eV: ", tot_ergs * erg_to_ev)

    plt.plot(nus, ergs_hz)
    plt.show()

# Equation 4.29 in osterbrock
def get_two_photon_emission_coefficient(Np, Ne):
    if use_lambda:
        x, g_v = get_two_photon_spectrum_vs_lambda()
    else:
        x, g_v = get_two_photon_spectrum_vs_nu()
    denominator = 1 + (Np * c.q_p_2s_2p + Ne * c.q_e_2s_2p) / c.A_2s_1s

    gamma_2p = c.alpha_eff_2s * g_v / denominator

    return x, gamma_2p

def sanity_check_emission_coeff():
    x, coeff = get_two_photon_emission_coefficient(1, 1) # TODO find Np, Ne
    _, ax = plt.subplots()
    if use_lambda:
        coeff *= c.c/x
        x *= 1e4 # lambda in units of microns
        ax.set(yscale="log", xscale="log", xlim=(1e-1, 1e1))
    else:
        ax.set(yscale="log", xlim=(3e14, 10e14), ylim=(1e-40, 1e-38))
    ax.plot(x, coeff)
    plt.show()

# 4.28 in osterbrock
def get_two_photon_jv(Np, Ne, plot=False):
    x, gamma = get_two_photon_emission_coefficient(Np, Ne)
    jv = Np * Ne * gamma / (4 * np.pi)
    if plot:
        _, ax = plt.subplots()
        ax.plot(x, jv)
        ax.set(yscale="log")
        plt.show()
    return x, jv

if __name__ == "__main__":
    # sanity_check_lambda()
    # sanity_check_nu()
    sanity_check_emission_coeff()
    # get_two_photon_jv(1, 1)
