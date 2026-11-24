import scipy.optimize
from ninatool.internal.structures import loop, Nlosc
from ninatool.internal.elements import L, J
import numpy as np
import matplotlib.pyplot as plt
import warnings
from meercat.physics.JJ import *
import qutip as qt
import matplotlib as mpl
mpl.use('Qt5Agg')
from meercat.deivce_models.NINA_SPA import NINA_SPA

F0 = 4.97e9
L_UNIT = 3.29e-8
C_UNIT = 3.9e-15
E_UNIT = 3.29e-24

class SNAILresonator(NINA_SPA):
    '''
    For now, just restricting to SNAILs with two large-area juctions (Manhattan style fab)
    '''

    def __init__(self, L=1e-9, C=1e-12, array_num=1, order=9, num_large_area_JJs=2):

        super().__init__()

        self.Ic_small_area = None
        self.Ic_large_area = None
        self.alpha = None
        self.order = order
        self.array_num = int(array_num)
        self.num_large_area_JJs = num_large_area_JJs
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

    def change_C(self, C):
        self.C = C
        self.spa.C = C / C_UNIT

    def change_L(self, L):
        try:
            self.L = L[0]
        except:
            self.L = L

        self.snail.Lstray.L0 = self.L / L_UNIT

        self.spa = Nlosc(self.snail)
        self.spa.C = self.C / C_UNIT

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

            # LEFT BRANCH ELEMENTS
            J0 = J(ic=Ic_to_Ej(self.Ic_small_area) / E_UNIT, order=self.order, name='J0') # NINA Ic is really Ej. dimensionless units are for nerds
            left_elements = [J0]

            # RIGHT BRANCH ELEMENTS
            if self.num_large_area_JJs == 2:
                J1 = J(ic=Ic_to_Ej(self.Ic_large_area) / E_UNIT, order=self.order, name='J1')
                J2 = J(ic=Ic_to_Ej(self.Ic_large_area) / E_UNIT, order=self.order, name='J2')
                right_elements = [J1, J2]
            if self.num_large_area_JJs == 3:
                J1 = J(ic=Ic_to_Ej(self.Ic_large_area) / E_UNIT, order=self.order, name='J1')
                J2 = J(ic=Ic_to_Ej(self.Ic_large_area) / E_UNIT, order=self.order, name='J2')
                J3 = J(ic=Ic_to_Ej(self.Ic_large_area) / E_UNIT, order=self.order, name='J3')
                right_elements = [J1, J2, J3]

            # DEFINES THE SNAIL AS A LOOP INSTANCE
            self.snail = loop(left_branch=left_elements,
                         right_branch=right_elements,
                         stray_inductance=True,
                         name='SNAIL')

            self.snail.Nloops = self.array_num
            self.snail.Lstray.L0 = self.L / L_UNIT

            # DEFINES A SNAIL PARAMETRIC AMPLIFIER BY SHUNTING THE SNAIL WITH A CAPACITOR
            self.spa = Nlosc(self.snail)
            self.spa.C = self.C / C_UNIT
            self.Npoints = len(self.spa.nlind.flux)

            Lj_small_area = Ic_to_Lj(self.Ic_small_area)
            Lj_large_area = Ic_to_Lj(self.Ic_large_area)
            self.total_Lj_at_nullflux = self.array_num/(1/(2*Lj_large_area)+1/Lj_small_area)

    def copy(self):

        copy = SNAILresonator(L=self.L, C=self.C, array_num=self.array_num, order=self.order)
        copy.define_large_area_Ic(self.Ic_large_area)
        copy.define_small_area_Ic(self.Ic_small_area)

        return copy

    def tune_freq_with_C(self, target_freq):

        def cost_function(C, args):
            target_freq = args[0]
            snailres_copy = args[1]
            snailres_copy.change_C(C)

            omega, g3, g4, g5, g6, g7, g8, g9, g10 = snailres_copy.get_gns_at_flux(snailres_copy.estimate_kerrnull_flux())
            cost = np.abs(omega / (2 * np.pi) - target_freq)
            return cost

        copy = self.copy()
        copy.try_to_build_model()
        args = [target_freq, copy]
        x0 = [copy.C]

        result = scipy.optimize.minimize(cost_function, x0, args=args, method='Nelder-Mead', bounds=[(1e-15, 1e-6)],
                                         tol=1e-17)
        self.change_C(result.x)
        return result.x[0]

    def tune_freq_with_L(self, target_freq):

        def cost_function(L, args):
            target_freq = args[0]
            snailres_copy = args[1]
            snailres_copy.change_L(L)

            omega, g3, g4, g5, g6, g7, g8, g9, g10 = snailres_copy.get_gns_at_flux(snailres_copy.estimate_kerrnull_flux())
            cost = np.abs(omega / (2 * np.pi) - target_freq)
            return cost

        copy = self.copy()
        copy.try_to_build_model()
        args = [target_freq, copy]
        x0 = [copy.L]

        result = scipy.optimize.minimize(cost_function, x0, args=args, method='Nelder-Mead', bounds=[(1e-13,1e-3)], tol=1e-2)
        self.change_L(result.x[0])
        return result.x[0]

if __name__ == '__main__':
    plt.close('all')
    snailres = SNAILresonator(L=1.98e-9, C=0.093e-12, array_num=3)
    snailres.define_large_area_Ic(0.44e-6)
    snailres.define_alpha(0.18)
    snailres.try_to_build_model()
    omega, g3, g4, g5, g6, g7, g8, g9, g10 = snailres.get_gns_at_flux(snailres.estimate_kerrnull_flux(), plot=True)
    transitions = snailres.get_transitions_at_flux(snailres.estimate_kerrnull_flux(), 50, plot=True)
    Lj_HFSS = snailres.get_Lj_at_flux(snailres.estimate_kerrnull_flux())
    transitions_spread = np.max(transitions) - np.min(transitions)
    metric = np.abs(transitions_spread / g3)
    print(metric)
