#!/usr/bin/python3.7

import matplotlib
import matplotlib.pyplot as plt

class td_plot():
    def wf_plot(self, x, psi, axis, fig):
        plt.xlabel('x ($\\mathrm{\AA}$)')
        #plt.ylabel('probability density ($\\mathrm{\AA}^{-1}$)')
        plt.ylabel("$\psi(q)$")
        #axis.plot(x, f*scaling + yoffset, color=COLOUR1)
        axis.plot(x, f, color='blue')
        #axis.fill_between(x, f*scaling + yoffset, yoffset, f > 0., color=COLOUR1, alpha=0.5)
        #axis.fill_between(x, f*scaling + yoffset, yoffset, f < 0., color=COLOUR2, alpha=0.5)

        plt.legend()
        plt.show()
        return

    def dens_plot(r, densities, energies):
        plt.xlabel('x ($\\mathrm{\AA}$)')
        plt.ylabel('probability density ($\\mathrm{\AA}^{-1}$)')

        energies = ['E = {: >5.2f} eV'.format(eigenvalues[i].real / e) for i in range(3)]
        plt.plot(r * 1e+10, densities[0], color='blue',  label=energies[0])
        plt.plot(r * 1e+10, densities[1], color='green', label=energies[1])
        plt.plot(r * 1e+10, densities[2], color='red',   label=energies[2])

        plt.legend()
        plt.show()
        return
