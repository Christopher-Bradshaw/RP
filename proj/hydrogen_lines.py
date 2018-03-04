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
    nus = get_base_array()
    width = get_line_width(nu)
    assert width / (nus[1] - nus[0]) > 10 # enough detail to integrate
    num_sds = 5
    start, stop = nu - num_sds*width, nu + num_sds*width
    assert start > nus[0] and stop < nus[-1]

    rv = norm(loc = nu, scale = width)
    ergs = rv.pdf(nus)*height

    return ergs

def get_lines(lines):
    output_ergs = np.zeros_like(get_base_array())
    all_ergs = []
    for line in lines:
        ergs = get_line(ev_to_nu(line["energies"]), line["intensities"])
        all_ergs.append(ergs)
        output_ergs += ergs
    return output_ergs, all_ergs

def get_line_emission_coeff(lines):
    nus = get_base_array()
    ergs, all_ergs = get_lines(lines)
    ergs = ergs * 1.24e-25
    all_ergs = [i*1.24e-25 for i in all_ergs]

    # Summed
    _, ax = plt.subplots()
    ax.plot(nus, ergs)
    ax.set(yscale="log", xlim=(3e14, 10e14), ylim=(1e-40, 1e-38))
    # plt.show()

    # Separate
    _, ax = plt.subplots()
    for i in all_ergs:
        ax.plot(nus, i)
    ax.set(yscale="log", xlim=(3e14, 10e14), ylim=(1e-40, 1e-38))
    print(all_ergs[1])
    tot_ergs = np.trapz(all_ergs[1], nus)
    print(tot_ergs)

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

    get_line_emission_coeff(lines)

if __name__ == "__main__":
    main()
