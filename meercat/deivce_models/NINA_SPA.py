import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

F0 = 4.97e9
L_UNIT = 3.29e-8
C_UNIT = 3.9e-15
E_UNIT = 3.29e-24

class NINA_SPA():

    def __init__(self):
        pass

    def estimate_kerrnull_flux(self):

        Npoints = len(self.spa.nlind.flux)
        phis = self.spa.nlind.flux

        omega_rough = self.spa.omega * F0 * 2 * np.pi
        g3_rough = self.spa.gn[0] * F0 * 2 * np.pi
        g4_rough = self.spa.gn[1] * F0 * 2 * np.pi

        K_rough = 12*(g4_rough-5*g3_rough**2/omega_rough)

        phi_rough = phis[np.where(np.abs(K_rough[Npoints//2:3*Npoints//4]) == np.min(np.abs(K_rough[Npoints//2:3*Npoints//4])))[0]+Npoints//2]

        phi_fine = np.linspace(-0.03, 0.03, 1001) + phi_rough

        omega_fine = np.interp(phi_fine, phis, self.spa.omega * F0 * 2 * np.pi)
        g3_fine = np.interp(phi_fine, phis, self.spa.gn[0] * F0 * 2 * np.pi)
        g4_fine = np.interp(phi_fine, phis, self.spa.gn[1] * F0 * 2 * np.pi)

        K_fine = 12 * (g4_fine - 5 * g3_fine ** 2 / omega_fine)

        return phi_fine[np.where(np.abs(K_fine)==np.min(np.abs(K_fine)))[0]]/(2*np.pi)

    def get_gns_at_flux(self, phi, plot=False):

        phis = self.spa.nlind.flux

        omega = np.interp(phi*2*np.pi, phis, self.spa.omega * F0 * 2 * np.pi)
        g3 = np.interp(phi*2*np.pi, phis, self.spa.gn[0] * F0 * 2 * np.pi)
        g4 = np.interp(phi*2*np.pi, phis, self.spa.gn[1] * F0 * 2 * np.pi)
        g5 = np.interp(phi*2*np.pi, phis, self.spa.gn[2] * F0 * 2 * np.pi)
        g6 = np.interp(phi*2*np.pi, phis, self.spa.gn[3] * F0 * 2 * np.pi)
        g7 = np.interp(phi*2*np.pi, phis, self.spa.gn[4] * F0 * 2 * np.pi)
        g8 = np.interp(phi*2*np.pi, phis, self.spa.gn[5] * F0 * 2 * np.pi)
        g9 = np.interp(phi * 2 * np.pi, phis, self.spa.gn[6] * F0 * 2 * np.pi)
        g10 = np.interp(phi * 2 * np.pi, phis, self.spa.gn[7] * F0 * 2 * np.pi)

        if plot:
            fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))
            axs[0].plot(self.spa.nlind.flux / 2 / np.pi, self.spa.omega * F0)
            axs[0].set_xlabel(r'$\bar{\varphi}_\mathrm{e}/2\pi$')
            axs[0].set_ylabel(r'$\omega / 2 \pi$', rotation=0, labelpad=10)
            axs[0].grid()

            Npoints = len(self.spa.nlind.flux)

            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(self.spa.gn[0][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g3')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(self.spa.gn[1][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g4')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(self.spa.gn[2][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g5')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(self.spa.gn[3][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g6')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi,
                            np.abs(self.spa.gn[4][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g7')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi,
                            np.abs(self.spa.gn[5][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g8')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi,
                            np.abs(self.spa.gn[6][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g9')
            axs[1].semilogy(self.spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi,
                            np.abs(self.spa.gn[7][Npoints // 2:3 * Npoints // 4]) * F0,
                            label='g10')
            axs[1].legend()
            axs[1].set_ylim([1e-5, 1e9])
            axs[1].set_xlabel('Relative flux')
            axs[1].set_ylabel('Nonlinearity (Hz)')

            plt.plot([phi, phi],
                     [1e-2, 3e10], 'k--')
            plt.grid()
            plt.show()

        return omega, g3, g4, g5, g6, g7, g8, g9, g10

    def get_Lj_at_flux(self, phi):

        omega, g3, g4, g5, g6, g7, g8, g9, g10 = self.get_gns_at_flux(phi, plot=False)

        return 1/(omega**2*self.C) - self.L


    def get_transitions_at_flux(self, phi, Ns, plot=False):

        omega, g3, g4, g5, g6, g7, g8, g9, g10 = self.get_gns_at_flux(phi)

        s = qt.destroy(int(Ns))
        H0 = omega * s.dag() * s + g3 * (s.dag() + s) ** 3 + g4 * (s.dag() + s) ** 4 + g5 * (s.dag() + s) ** 5 + g6 * (
                    s.dag() + s) ** 6 + g7 * (s.dag() + s) ** 7 + g8 * (s.dag() + s) ** 8 # + g9 * (
                         #s.dag() + s) ** 9 + g10 * (s.dag() + s) ** 10

        transitions = np.diff(H0.eigenstates()[0][0:Ns // 2])

        if plot:

            fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 2.5), dpi=200)
            axs[0].plot(np.diff(H0.eigenstates()[0][0:Ns // 2]) / (2 * np.pi * 1e9), 'o-')
            axs[0].set_ylabel('Transition Energies (GHz)')
            axs[0].grid()
            axs[0].set_xlabel('Photon number')
            # axs[0].set_ylim([0, omega[i]*2])

            axs[1].plot(np.diff(np.diff(H0.eigenstates()[0][0:Ns // 2])) / (2 * np.pi * 1e6), 'o-')
            axs[1].set_ylabel('Anaharmonicity (MHz)')
            axs[1].grid()
            axs[1].set_xlabel('Photon number')
            # axs[1].set_ylim([-12*g4[i]*2, 12*g4[i]*2])
            fig.tight_layout()
            plt.show()

        return transitions
