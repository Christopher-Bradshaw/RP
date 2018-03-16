import numpy as np
import math as math
import matplotlib.pyplot as plt

## Constants 
Z = 1                    ## units: unitless      (atomic number)
e = 4.80320425*1e-10     ## units: StatCoulombs  (charge of electron)
me = 9.10938356*1e-28    ## units: g             (mass of electron)           
c = 2.99792458*1e10      ## units: cm / s        (speed of light)
h = 6.62607004*1e-27     ## units: erg * s       (Planck constant)
Ry = 2.1798741e-11       ## units: erg           (Rydberg constant (Ionization energy)
kb = 1.380658e-16        ## units: erg/K         (Bolztmann constant)


def b(m, eta, rho, s, ll):
    """
    eq. C8 from Karzas and Latter 1961
    b_s values to calculate G1
    """
    
    if s == 0:
        return 1.
    
    if s == 1:
        return  2. * m * eta / ll
    
    #if s >= 2:
    bb_1 = b(m, eta, rho, s-1, ll) #(b_(s-1)
    #print bb_1
    bb_2 = b(m, eta, rho, s-2, ll) #(b_(s-2)

    coef = -1./ (s * (s + (2*ll) - 1))
    bbA = (4. * eta) * (s - 1 - m) * bb_1
    bbB = (2*m + 2 - s) * (2*m + 2*ll + 1 - s) * bb_2

    bb = coef * (bbA + bbB)

    return bb

def g(mm,eta,rho,ll):
    """
    ## eq. C5 Karzas and Latter 1961
    """
    m = -mm
    gg = 0.
    if m>=0:
        for s in range(2*m+1):
            gg = gg + b(m,eta,rho,s,ll)*rho**s 
    if m<0:
        print('m < 0 ?')
    return gg


def sigma_minus(ni, n, l, eta, rho):
    """
    ## eq. 36 Karzas and Latter 1961 for transitions to a lower l state
    """

    if l == 0:
        ## if l is already zero it can't go lower
        sigma = 0

        return sigma

    else:
        
        sigA = (np.pi * e**2) / (me * c * ni)
        sigB = (2. ** (4.*l)) / 3.
        sigC = (l**2.) * np.math.factorial(n + l)

        sigD = 1.
        for l_val in np.arange(1,l):
            #end not inclusive (l -1)
            sigD *= ( (l_val**2) + (eta**2))

        sigE = np.math.factorial( (2. * l) + 1) * np.math.factorial( (2. * l) - 1) * np.math.factorial( n - l - 1.)
        sigE = float(sigE)

        invcotrho = np.pi/2. - np.arctan(rho)
        sigF = np.exp( -4. * eta * invcotrho)  * (rho ** ((2. * l) + 2.))
        sigG = (1 - np.exp(- 2. * np.pi * eta)) * ((1. + rho**2) ** (2*n - 2.))

        #Gl portion
        sigH = g( l + 1 - n,  eta, rho, l) - (((1 + rho**2)**-2.) *  g( l - 1 - n, eta, rho, l))
        sigH = sigH**2
        
        sigma = sigA * sigB * sigC * (sigD / sigE) * (sigF / sigG) * sigH
    
        return sigma

def sigma_plus(ni, n, l, eta, rho):
    """
     ## eq. 36 Karzas and Latter 1961 for transitions to a upper l state
    """
    sigA = (np.pi * e**2) / (me * c * ni)
    sigB = (2. ** (4.*l + 6.)) / 3.
    sigC = ((l + 1.)**2.) * np.math.factorial(n + l)

    sigD = 1.
    for l_val in np.arange(0,l+1):
        #end not inclusive (l)
        sigD *= ( ((l_val+1 )**2.) + (eta**2.))

    sigEA = (2.*l + 1.) * np.math.factorial(2.*l + 1) * np.math.factorial(2.*l + 2.) *  np.math.factorial(n-l-1.)
    sigEB = (( (l + 1.)**2.) + (eta**2.) )**2.
    sigE = sigEA * sigEB

    invcotrho = np.pi/2. - np.arctan(rho)
    sigF = np.exp( -4. * eta * invcotrho)  * (rho ** ((2. * l) + 4.)) * (eta**2.)
    sigG = (1 - np.exp(- 2. * np.pi * eta)) * ((1. + rho**2.) ** (2.*n))

    #Gl portion
    sigHA = (l + 1. - n) * g( l + 1 - n,  eta, rho, l+1)
    sigHB = (((l + 1. + n)/ (1. + rho**2)) ) * g( l - n,  eta, rho, l+1)

    sigH = (sigHA + sigHB)**2.

    sigma = sigA * sigB * sigC * (sigD / sigE) * (sigF / sigG) * sigH

    return sigma


def gauntf(ni, n, l, eta, rho):

    """
     ## eq. 40 Karzas and Latter 1961 for bound free gaunt factor which
    is the ratio of the total bound free cross section( equation 38) and
    the Karmer's cross section (equation 39).

    The cross section is the sum of the up and down l transitions within
    a given energy level.
    """

    
    sigma_bf = sigma_plus(ni, n, l, eta, rho) + sigma_minus(ni, n, l, eta, rho)

    sigma_KA = (2.**4) / (3. * (3.**.5))
    sigma_KB = (e**2) / (me * c * ni) / n
    sigma_KC = ((rho**2.) / (1 + rho**2))**2.

    sigma_K = sigma_KA * sigma_KB * sigma_KC

    g_bf = sigma_bf/sigma_K

    return g_bf


def gauntlavg(n, nu, eta, rho):

    """
    ## The average bound free gaunt factor for a given energy level n.

    The average is weighted by the number of possible states in a given l
    of a energy level. Osterbrock and Finland talks about this, but I am 
    confused - Brittany M.
    """

    gtot = 0. 
    for l in range(n):
        gtot = gtot + gauntf(nu ,n, l, eta, rho)*((2*l+1))
    return gtot/n**2



def bf_emission(n, lam, T):
    """
    Calculating the bound free emission coefficients following equation 1 from
    Brown and Matthews 1970 using, the photoionization cross section a(v) 
    from Seaton 1960.
    """

    nu = c/lam
    E = h * nu
    eta = np.sqrt( (Z**2) * Ry / E )
    rho = eta/n

    I_n = Ry * (Z**2) / (n**2)  #Ionization potential of energy level n

    if E < I_n:
        # the photon of an energy below the ionization energy can't actually
        # ionize the gas.
        sigma_bf = 0
        
    else:
        # Calculating the photoionization cross section from the first
        # equation of Seaton 1960. The gaunt factor used is the averge
        # bound free gaunt factor

        epsilon = (nu * h / (Z**2) / Ry) - (1./n**2)
        sigma_A = 7.907e-18 
        sigma_B =  (n / (Z**2)) * ((1 + (n**2) * epsilon))**-3.

        g_avgbf =  gauntlavg(n, nu, eta, rho)
        #print g_avgbf

        sigma_bf = sigma_A * sigma_B * g_avgbf

    #Calculating the Emission Coefficient
    gamma_A = ((2./np.pi)**.5) * np.exp(I_n / kb /T  ) / ( (c**2) * ( (me * kb) * T)**(3./2.) )

    gamma_B = 2 * (n**2) * h * ((h * nu)**3) * sigma_bf  *  np.exp( - (h  / kb) * nu  /T)

    gamma = gamma_A * gamma_B #[erg/cm^3] or [erg/cm^3 s Hz]

    return gamma



def get_bound_free_emission_coefs(T, lams = np.linspace(.1,10,1000) * 1e-4, n_levels = 8, print_progress = False):
    """
    Inputs 
    T - Temperature of the gas                                              #[Kelvin]
    lams - grid of wavelengths to calculate your spectrum                   #[cm]
    n_levels = Number of levels to consider transitions to/from. 
               ~10 is fine but if you have time do 15. Increasing n
               allows more contribution from lower energy states/photons.

    print_progress - if you want updates, my code is very slow

    Outputs
    gammas [erg cm^3] or [erg cm^3 Hz * s-1]
    lams [cm] - if you input  a grid it will still pop this out
    nus  [Hz] - associated frequencies of the wavelengths
    
    you may need to np.argsort for a vs. frequency addition with another spectrum.
    """
    
    print 'Evaluating the bound free emission coefficients with ' + str(n_levels) + ' quantum levels'
    
    nus = c/lams
    
    gammas = np.zeros(len(lams))

    for lam_i, lam in enumerate(lams):

        if print_progress == True:
            print str(lam_i + 1) + ' of ' +  str(len(lams)) + ' wavelength elements'
            
        gamma = 0
        for n_i in np.arange(n_levels)+1:
            #print n_i
            gamma += bf_emission(n_i, lam, T)

            gammas[lam_i] = gamma

    return gammas, lams, nus
