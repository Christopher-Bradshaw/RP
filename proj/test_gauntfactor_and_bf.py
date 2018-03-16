import boundfree
import matplotlib.pyplot as plt
import numpy as np

#constants
#not the best way sorry
Z = boundfree.Z
Ry = boundfree.Ry
c = boundfree.c
h = boundfree.h

logEZ2Ry = np.linspace(-6,13, 200)
EZ2Ry = 10**logEZ2Ry  
etas = np.sqrt(1./EZ2Ry)
ns = np.arange(1, 7)  #energy levels to cover s (l = 0)
ns2 = np.arange(2, 7)  #energy levels to cover p (l = 1)
g_bfs = np.zeros((len(etas), len(ns)))
g_bfs_p = np.zeros((len(etas), len(ns)))

#getting associated wavelengths with the energies
Es = (Z**2) * Ry / (etas**2)   #[erg]
nus = Es/ h          #[Hz]
lams = (c/nus) * 1e4 #[microns]


print 'calculating the gaunt factors up to n level ' + str(np.max(ns))
for n_i, n_val in enumerate(ns):
    print n_i+1
    for eta_i, eta_val in enumerate(etas):

        rho_val = eta_val/n_val
        E_val = (Z**2) * Ry / (eta_val**2)
        nu_val = E_val/h
        #lam_val = gauntfunc.c/nu_val
        #print str(lam_val * 1e4) + ' microns'
        
        g_bf = boundfree.gauntf(nu_val, n_val, 0, eta_val, rho_val) #1s

        g_bfs[eta_i,n_i] = g_bf


# Reproducing Figures from Karzas and Latter 1961 for the gaunt factor at individual n and l states.
# more figures can be reproduced but showing that the p orbitals can be reproduced we know the up and down work. s only tests l+ transitions
plt.figure('gbfs')
plt.subplot(141)
plt.title('Figure 7')

for n_i, n_val in enumerate(ns):
    plt.plot(logEZ2Ry , g_bfs[:, n_i], label = str(n_val) + 's')

plt.legend()
plt.yscale('log')
plt.ylim([1e-1, 1e2])
plt.xlim([-4,3])
plt.ylabel('g_bf')
plt.xlabel('log10(E/Z^2/Ry)')


plt.subplot(142)
plt.title('Figure 8')

for n_i, n_val in enumerate(ns):
    plt.plot(logEZ2Ry , g_bfs[:, n_i], label = str(n_val) + 's')

plt.legend()
plt.yscale('log')
plt.ylim([1e-5, 1e2])
plt.xlim([1,13])
plt.xlabel('log10(E/Z^2/Ry)')


print 'calculating the average gaunt factors up to n level ' + str(np.max(ns2))
for n_i, n_val in enumerate(ns2):
    print n_i
    for eta_i, eta_val in enumerate(etas):

        rho_val = eta_val/n_val
        E_val = (Z**2) * Ry / (eta_val**2)
        nu_val = E_val/h
        #lam_val = gauntfunc.c/nu_val
        #print str(lam_val * 1e4) + ' microns'
        
        g_bf = boundfree.gauntf(nu_val, n_val, 1, eta_val, rho_val) #1s

        g_bfs_p[eta_i,n_i] = g_bf


plt.subplot(143)
plt.title('Figure 9')

for n_i, n_val in enumerate(ns2):
    plt.plot(logEZ2Ry , g_bfs_p[:, n_i], label = str(n_val) + 'p')

plt.legend()
plt.yscale('log')
plt.ylim([1e-2, 1e1])
plt.xlim([-4,3])
plt.xlabel('log10(E/Z^2/Ry)')


plt.subplot(144)
plt.title('Figure 10')

for n_i, n_val in enumerate(ns2):
    plt.plot(logEZ2Ry , g_bfs_p[:, n_i], label = str(n_val) + 'p')

plt.legend()
plt.yscale('log')
plt.ylim([1e-7, 1e0])
plt.xlim([1,13])
plt.xlabel('log10(E/Z^2/Ry)')
plt.show()



#reproducing the average gaunt factor plots

#for figure 19
ns3 = np.arange(1, 7) #15 takes goddamn forever
g_avgbfs = np.zeros((len(etas), len(ns3)))
for n_i, n_val in enumerate(ns3):
    print n_i+1
    for eta_i, eta_val in enumerate(etas):

        rho_val = eta_val/n_val
        E_val = (Z**2) * Ry / (eta_val**2)
        nu_val = E_val/ h

        g_avgbf = boundfree.gauntlavg(n_val, nu_val, eta_val, rho_val)

        g_avgbfs[eta_i,n_i] = g_avgbf


# Figures 9 - 18
plt.figure('gbfs avg')
plt.subplot(121)
plt.title(' Figure 19 ')

for n_i, n_val in enumerate(ns3):
    plt.plot(logEZ2Ry , g_avgbfs[:, n_i], label = 'n = '+ str(n_val))

plt.legend()
#plt.yscale('log')
plt.ylim([.91, 1.12])
plt.xlim([-6,1])
plt.ylabel('g_bf average')
plt.xlabel('log10(E/Z^2/Ry)')

plt.subplot(122)
plt.title('average g_bf over relevant wavelengths')

for n_i, n_val in enumerate(ns3):
    plt.plot(lams , g_avgbfs[:, n_i], label = 'n = '+ str(n_val))

plt.legend()
plt.xscale('log')
plt.ylim([.91, 1.12])
plt.xlim([1e-1,1e1])
plt.xlabel('microns')
plt.show()


#Calculate Emission Coefficients
T = 10000.
gammas, lams, nus = boundfree.get_bound_free_emission_coefs(T, lams = np.linspace(.1,10,1000) * 1e-4, n_levels = 8, print_progress = True)

plt.figure('bound free emission ceofficients')
plt.subplot(211)
plt.title("T = 10,000 K nebula")
plt.plot(lams * 1e4, gammas)
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\gamma_{\nu}$ (erg * cm$^3$)')
plt.xlabel('microns')

plt.subplot(212)
plt.plot(lams * 1e4, gammas * nus)
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\gamma_{\nu}$ * $\nu$ (erg  cm$^3$  Hz)')
plt.xlabel('microns')
plt.ylim([1e-26, 1e-23])
plt.show()
