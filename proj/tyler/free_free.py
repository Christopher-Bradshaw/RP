import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d
#

# data_2 = np.loadtxt('table_2.dat')
#
# d2_1 = data_2[:,:5]
# d2_2 = data_2[:,5:]
#
# data_2 = np.concatenate((d2_1,d2_2),axis=0)
# np.savetxt('table_2_ordered.dat', data_2)

data_2 = np.loadtxt('table_2_ordered.dat')

logGamma_data = data_2[:,0]
gaunt_data = data_2[:,1]

get_gaunt_gamma = interp1d( logGamma_data, gaunt_data )

nPoints = 100
logGamma_range = np.linspace(-4, 4, nPoints )
gaunt_interp = get_gaunt_gamma( logGamma_range )

# plt.figure(0)
# plt.clf()
# plt.plot(logGamma_range, gaunt_interp )
# plt.show()



data_1 = np.loadtxt('table_1.dat')
logGamma_data = data_1[1:,0]
logU_data = data_1[0,1:]
gaunt_data = data_1[1:,1:]

get_gaunt_u_gamma = interp2d(  logU_data, logGamma_data, gaunt_data, kind='cubic' )

nPoints = 100
logGamma_range = np.linspace(-4, 4, nPoints )
logU_range = np.linspace(-4, 4, nPoints )

gaunt_interp = get_gaunt_u_gamma( logU_range, logGamma_range )

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

for i in range( 0,nPoints, 2 ):
    ax1.plot( logU_range, gaunt_interp[:,i])
    ax2.plot( logGamma_range, gaunt_interp[i,:])

# plt.show()
fig1.savefig( 'ganunt_1.png')
fig2.savefig( 'ganunt_2.png')

# Calculating Emission Coefficient
def gamma_nu(nu, T):
    h = 6.62606885e-27 # erg * s
    kB = 1.3806485279e-16 # erg/K
    Ry = 2.1798723e-11 # erg
    e = 4.8032047e-10 # statC

    e4 = e**4
    # e4 = 5.322607e-38 # e^4
    m = 9.109383e-28 # g
    # m =1.6726219e-24 #g
    c = 2.99792458e10 # cm/s
    c3 = c*c*c
    # c3 = 2.69440024173740e25 # c^3
    #nu_0 = Ry/h # Brant said this is probably the ionization frequency
    Z = 2 # Hydrogen
    logU = np.log10((h * nu)/(kB * T))
    logGammaSquared = np.log10((Z*Z * Ry)/(kB * T))
    # Note that this j_nu is per N_+ and N_e
    # return  (32 * Z*Z * e4 * h)/(3 * m*m * c3) * np.sqrt((np.pi * Ry)/(3 * kB * T)) * np.exp((-h * nu)/(kB * T)) * get_gaunt_u_gamma(logU, logGammaSquared)
    return  6.8e-38 /np.sqrt(T) * np.exp((-h * nu)/(kB * T)) * get_gaunt_u_gamma(  logU, logGammaSquared )

def j_nu(nu, T, N_i, N_e):
    return 1/(4 * np.pi) * N_i * N_e * gamma_nu(nu, T)

def alpha_nu(nu, T, N_i, N_e):
    h = 6.62606885e-27 # erg * s
    kB = 1.3806485279e-16 # erg/K
    logU = np.log10((h * nu)/(kB * T))
    logGammaSquared = np.log10((Z*Z * Ry)/(kB * T))
    return 3.7e8/np.sqrt(T) * N_i * N_e * (1 - np.exp((-h * nu)/(kB * T)))/np.power(nu, 3) * get_gaunt_u_gamma(logU, logGammaSquared)


nu_array = np.logspace(13, 16, 1000)
gamma_array = np.zeros(nu_array.shape)
for i, nu in enumerate(nu_array):
    gamma_array[i] = gamma_nu(nu, 10000)

plt.plot(3e14/nu_array, nu* gamma_array)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-1, 10)
plt.xlabel('Wavelength (microns)')
plt.ylim(1e-26, 1e-22)
plt.ylabel('Nu * Gamma_nu  (free-free)')
plt.savefig('free_free.png')
