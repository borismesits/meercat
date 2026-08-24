from ninatool.internal.structures import loop, Nlosc
from ninatool.internal.elements import L, J
import numpy as np
import matplotlib.pyplot as plt
import warnings
from meercat.physics.JJ import *

class SNAILresonator():
    '''
    For now, just restricting to SNAILs with two large-area juctions (Manhattan style fab)
    '''

    def __init__(self, L=1e-9, C=1e-12, array_num=1, order=9):
        self.Ic_small_area = None
        self.Ic_large_area = None
        self.alpha = None
        self.order = order
        self.array_num = array_num
        self.L = L
        self.C = C

    def define_large_area_Lj(self, Lj):
        self.Ic_large_area = Lj_to_Ic(Lj)

    def define_large_area_Ic(self, Ic):
        self.Ic_large_area = Ic

    def define_large_area_Ej(self, Ej):
        self.Ic_large_area = Ej_to_Ic(Ej)

    def define_small_area_Lj(self, Lj):
        self.Ic_small_area = Lj_to_Ic(Lj)

    def define_small_area_Ic(self, Ic):
        self.Ic_small_area = Ic

    def define_small_area_Ej(self, Ej):
        self.Ic_small_area = Ej_to_Ic(Ej)

    def define_alpha(self, alpha):

        if alpha < 0:
            warnings.warn('alpha cannot be negative')
            raise Exception

        elif alpha > 1:
            self.alpha = 1/alpha
            warnings.warn(
                'alpha is defined to be less than 1 in our convention! setting alpha to reciprocal of provided value')

        else:
            self.alpha = alpha

    def try_to_build_model(self):

        if self.alpha != None and self.Ic_large_area != None and self.Ic_small_area != None:
            warnings.warn('Cannot define all three of alpha, large junction Ej/Ic/Lj, and small junction Ej/Ic/Lj')
            raise Exception

        else:
            if self.alpha == None:
                self.alpha = self.Ic_small_area/self.Ic_large_area
            elif self.Ic_small_area == None:
                self.Ic_small_area = self.Ic_large_area * self.alpha
            elif self.Ic_large_area == None:
                self.Ic_large_area = self.Ic_small_area / self.alpha

            f0 = 4.97e9
            Lunit = 3.29e-8
            Cunit = 3.9e-15
            Eunit = 3.29e-24

            # LEFT BRANCH ELEMENTS
            J0 = J(ic=Ic_to_Ej(self.Ic_small_area) / Eunit, order=self.order, name='J0') # NINA Ic is really Ej. dimensionless units are for nerds
            left_elements = [J0]

            # RIGHT BRANCH ELEMENTS
            J1 = J(ic=Ic_to_Ej(self.Ic_large_area) / Eunit, order=self.order, name='J1')
            J2 = J(ic=Ic_to_Ej(self.Ic_large_area), order=self.order, name='J2')
            right_elements = [J1, J2]

            # DEFINES THE SNAIL AS A LOOP INSTANCE
            snail = loop(left_branch=left_elements,
                         right_branch=right_elements,
                         stray_inductance=True,
                         name='SNAIL')

            snail.Nloops = self.array_num
            snail.Lstray.L0 = self.L / Lunit

            # DEFINES A SNAIL PARAMETRIC AMPLIFIER BY SHUNTING THE SNAIL WITH A CAPACITOR
            self.spa = Nlosc(snail)
            self.spa.C = self.C / Cunit
            self.Npoints = len(self.spa.nlind.flux)

    def get_H0_at_flux(self, phi):
        self.try_to_build_model()



    def get_gs_at_flux(self, phi):
        pass







def evaluate_design(L, C, alpha, Ic, M, i, nbar0=3, gain=10, g_Delta=0.1, plot=True, Ns=100):
    # default unit values
    f0 = 4.97e9
    Lunit = 3.29e-8
    Cunit = 3.9e-15
    Eunit = 3.29e-24

    hbar = 1.0545718e-34
    e = 1.60218e-19
    phi0 = np.pi * hbar / e  # this is the unreduced flux quantum, needed for next expression
    EJ = Ic * phi0 / (2 * np.pi)

    order = 9
    # LEFT BRANCH ELEMENTS
    J0 = J(ic=EJ / Eunit, order=order, name='J0')
    left_elements = [J0]

    # RIGHT BRANCH ELEMENTS
    J1 = J(ic=EJ / alpha / Eunit, order=order, name='J1')
    J2 = J(ic=EJ / alpha / Eunit, order=order, name='J2')
    right_elements = [J1, J2]

    # DEFINES THE SNAIL AS A LOOP INSTANCE
    snail = loop(left_branch=left_elements,
                 right_branch=right_elements,
                 stray_inductance=True,
                 name='SNAIL')

    snail.Nloops = M
    snail.Lstray.L0 = L_HFSS / Lunit

    # DEFINES A SNAIL PARAMETRIC AMPLIFIER BY SHUNTING THE SNAIL WITH A CAPACITOR
    spa = Nlosc(snail)
    spa.C = C_HFSS / Cunit

    Npoints = len(spa.nlind.flux)

    omega = spa.omega[Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g3 = spa.gn[0][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g4 = spa.gn[1][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g5 = spa.gn[2][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g6 = spa.gn[3][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g7 = spa.gn[4][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g8 = spa.gn[5][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g9 = spa.gn[6][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi
    g10 = spa.gn[7][Npoints // 2:3 * Npoints // 4] * f0 * 2 * np.pi

    if plot:
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))
        axs[0].plot(spa.nlind.flux / 2 / np.pi, spa.omega * f0)
        axs[0].set_xlabel(r'$\bar{\varphi}_\mathrm{e}/2\pi$')
        axs[0].set_ylabel(r'$\omega$', rotation=0, labelpad=10)
        axs[0].grid()

        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g3) / (2 * np.pi), label='g3')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g4) / (2 * np.pi), label='g4')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g5) / (2 * np.pi), label='g5')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g6) / (2 * np.pi), label='g6')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g7) / (2 * np.pi), label='g7')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g8) / (2 * np.pi), label='g8')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g9) / (2 * np.pi), label='g9')
        axs[1].semilogy(spa.nlind.flux[Npoints // 2:3 * Npoints // 4] / 2 / np.pi, np.abs(g10) / (2 * np.pi),
                        label='g10')
        axs[1].legend()
        axs[1].set_ylim([1e-5, 1e9])
        axs[1].set_xlabel('Relative flux')
        axs[1].set_ylabel('Nonlinearity (Hz)')

        plt.plot([spa.nlind.flux[Npoints // 2 + i] / 2 / np.pi, spa.nlind.flux[Npoints // 2 + i] / 2 / np.pi],
                 [1e-2, 3e10], 'k--')
        plt.grid()
        plt.show()

        print(spa.nlind.flux[i] / 2 / np.pi)

        print(snail.Lstray.L0 * Lunit)
        print(spa.C * Cunit)

        print(phi0 / Ic)

    s = destroy(Ns)

    s_classical = np.linspace(-100, 100, 1001)

    H0 = omega[i] * s.dag() * s + g3[i] * (s.dag() + s) ** 3 + g4[i] * (s.dag() + s) ** 4 + g5[i] * (s.dag() + s) ** 5 + \
         g6[i] * (s.dag() + s) ** 6 + g7[i] * (s.dag() + s) ** 7 + g8[i] * (s.dag() + s) ** 8 + g9[i] * (
                     s.dag() + s) ** 9 + g10[i] * (s.dag() + s) ** 10

    H_classical = omega[i] * s_classical ** 2 + g3[i] * s_classical ** 3 + g4[i] * s_classical ** 4 + g5[
        i] * s_classical ** 5 + g6[i] * s_classical ** 6 + g7[i] * s_classical ** 7 + g8[i] * s_classical ** 8 + g9[
                      i] * s_classical ** 9 + g10[i] * s_classical ** 10

    # H0 = omega[i]*s.dag()*s + g3[i]*(s.dag() + s)**3 + g4[i]*(s.dag()*s.dag()*s*s)

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

    return spa, omega[i], g3[i], g4[i], g5[i], g6[i]
