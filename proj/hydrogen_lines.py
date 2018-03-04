import math
import numpy as np
import constants as c
from scipy.stats import norm
import matplotlib.pyplot as plt

def get_base_array():
    # we want to go from 1e14 to 1e15. 901 points gives us slices of 1e12
    return np.linspace(1e14, 1e15, num=901)

def ev_to_nu(ev):
    return ev / 4.135667662e-15 # plancks const in ev

# We are assuming thermal line broading
def get_line_width(nu):
    return math.sqrt(c.k * c.T / (c.mH * c.c**2)) * nu

def get_line(nu, height):
    width = get_line_width(nu)
    num_sds = 5
    start, stop = nu - num_sds*width, nu + num_sds*width
    nus = np.linspace(start, stop)

    rv = norm(loc = nu, scale = width)
    ergs = rv.pdf(nus)*height
    ergs[0], ergs[-1] = 0, 0 # just zero it out
    # plt.plot(nus, rv.pdf(nus))
    # plt.show()

    return nus, ergs

def get_lines(lines):
    lines = np.sort(lines, order="energies")
    all_nus, all_ergs = [], []
    for line in lines:
        nus, ergs = get_line(ev_to_nu(line["energies"]), line["intensities"])
        all_nus.append(nus)
        all_ergs.append(ergs)
    nus = np.concatenate(all_nus)
    ergs = np.concatenate(all_ergs)
    return nus, ergs

def get_line_emission_coeff(lines):
    nus, ergs = get_lines(lines)
    ergs = ergs * 1.24e-25
    _, ax = plt.subplots()
    ax.plot(nus, ergs)
    ax.set(yscale="log", xlim=(3e14, 10e14), ylim=(1e-40, 1e-38))
    plt.show()
    return nus, ergs


def main():
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
        *np.linspace(0.0529, 0.0154, num=5, endpoint=False),
        *np.linspace(0.0154, 0.0064, num=5, endpoint=False)])
    paschen_energies = energy_levels[3:] - energy_levels[2]
    # 4, 5 \n, 6 - 10\n 11 - 15 \n 16 - 20
    paschen_intensities = np.array([0.352, 0.354,
        0.354, 0.352, 0.350, 0.350, 0.350,
        *np.linspace(0.350, 0.344, num=5, endpoint=False),
        *np.linspace(0.344, 0.344, num=5, endpoint=False)] * balmer_intensities[1:])

    energies = np.concatenate((balmer_energies, paschen_energies))
    intensities = np.concatenate((balmer_intensities, paschen_intensities))
    lines = np.zeros(len(energies), dtype=[("energies", np.float), ("intensities", np.float)])
    lines["energies"] = energies
    lines["intensities"] = intensities

    get_line_emission_coeff(lines)

if __name__ == "__main__":
    main()
