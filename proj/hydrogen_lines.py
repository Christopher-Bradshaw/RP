import math
import numpy as np
import constants as c
from scipy.stats import norm
import matplotlib.pyplot as plt

def get_base_array():
    return np.linspace(1e14, 1e15, num=4000000)

def ev_to_nu(ev):
    return ev / 4.135667662e-15 # plancks const in ev

# We are assuming thermal line broading
# I think I need more broadnign from somewhere
def get_line_width(nu):
    return math.sqrt(c.k * c.T / c.mH) * nu/c.c

def get_line(nu, height):
    width = get_line_width(nu)
    num_sds = 5
    min_nu, max_nu = nu - num_sds*width, nu + num_sds*width
    nus = np.linspace(min_nu, max_nu, num=100)
    rv = norm(loc = nu, scale = width)
    ergs = rv.pdf(nus)*height

    ergs[0], ergs[-1] = 0, 0 # terminate at tails
    return nus, ergs

def get_lines(lines):
    all_ergs, all_nus = [], []
    for line in lines:
        nus, ergs = get_line(ev_to_nu(line["energies"]), line["intensities"])
        all_ergs.append(ergs)
        all_nus.append(nus)
    return np.array(all_nus), np.array(all_ergs)


def get_line_emission_coeff():
    lines = choose_line_energies()
    all_nus, all_emission_hbeta = get_lines(lines)
    all_ergs = 1.24e-25 * all_emission_hbeta # magic number is hBeta

    for nus in all_nus:
        start, end = nus[0], nus[1]
        start_lt_starts = start < all_nus[:,0]
        end_lt_starts = end < all_nus[:,0]
        # If the start is less than some other start, then so must the end be
        # Therefore all xors must be false
        # If this is true, we cen just sort by the nus and plot
        assert np.all(np.logical_not(np.logical_xor(end_lt_starts, start_lt_starts)))

    order = np.argsort(all_nus.flatten())
    nus = all_nus.flatten()[order]
    ergs = all_ergs.flatten()[order]

    # Separate
    # _, ax = plt.subplots()
    # for i, ergs in enumerate(all_ergs):
    #     ax.plot(all_nus[i], ergs)
    # ax.set(yscale="log")#, xlim=(3e14, 10e14), ylim=(1e-40, 1e-38))
    # tot_ergs = np.trapz(all_ergs[1], all_nus[1])
    # print("Integration of H beta:", tot_ergs)
    # plt.show(block=False)

    # Together
    # _, ax = plt.subplots()
    # ax.set(yscale="log")#, xlim=(3e14, 10e14), ylim=(1e-40, 1e-38))
    # ax.plot(nus, ergs)
    # plt.show()

    return nus, ergs

def choose_line_energies():
    highest_level = 21
    # In case B, all Lyman series photons above lyman alpha are converted
    # In case B, we ignore everything that goes to level 1 (because it will be absorbed immediately)
    # Note that at the moment I am using the table with no collisional transitions
    # It doesn't really change it to be honest
    energy_levels = np.array([-13.6 / n**2 for n in range(1, highest_level)])
    balmer_energies = energy_levels[2:] - energy_levels[1]
    # 3, 4, 5 \n, 6 - 10\n 11 - 15 \n 16 - 20
    balmer_intensities = np.array([2.87, 1, 0.466,
        0.256, 0.158, 0.105, 0.073, 0.0529,
        *np.linspace(0.0529, 0.0154, num=6)[1:],
        *np.linspace(0.0154, 0.0064, num=6 )[1:]])
    paschen_energies = energy_levels[3:] - energy_levels[2]
    # 4, 5 \n, 6 - 10\n 11 - 15 \n 16 - 20
    paschen_intensities = np.array([0.352, 0.354,
        0.354, 0.352, 0.350, 0.350, 0.350,
        *np.linspace(0.350, 0.344, num=6)[1:],
        *np.linspace(0.344, 0.344, num=6)[1:]] * balmer_intensities[1:])

    energies = np.concatenate((balmer_energies, paschen_energies))
    intensities = np.concatenate((balmer_intensities, paschen_intensities))
    lines = np.zeros(len(energies), dtype=[("energies", np.float), ("intensities", np.float)])
    lines["energies"] = energies
    lines["intensities"] = intensities
    return(lines)

def main():
    get_line_emission_coeff()

if __name__ == "__main__":
    main()
