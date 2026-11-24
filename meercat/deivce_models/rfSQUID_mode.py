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

class rfSQUIDresonator(NINA_SPA):

    def __init__(self, L=1e-9, C=1e-12, array_num=1, order=9):

        super().__init__()

        self.Ic = None
        self.Lshunt = None
        self.order = order
        self.array_num = int(array_num)
        self.L = L
        self.C = C

    def define_Lj(self, Lj):
        self.Ic = Lj_to_Ic(Lj)

    def define_Ic(self, Ic):
        self.Ic = Ic

    def define_Ej(self, Ej):
        self.Ic = Ej_to_Ic(Ej)

    def define_Lshunt(self, Lshunt):
        self.Lshunt = Lshunt

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

        # LEFT BRANCH ELEMENTS
        J0 = J(ic=Ic_to_Ej(self.Ic) / E_UNIT, order=self.order, name='J0') # NINA Ic is really Ej. dimensionless units are for nerds
        left_elements = [J0]

        # RIGHT BRANCH ELEMENTS
        L0 = L(L=self.Lshunt/L_UNIT, order=self.order, name='L0')
        right_elements = [L0]

        # DEFINES THE SNAIL AS A LOOP INSTANCE
        self.squid = loop(left_branch=left_elements,
                     right_branch=right_elements,
                     stray_inductance=True,
                     name='rfSQUID')

        self.squid.Nloops = self.array_num
        self.squid.Lstray.L0 = self.L / L_UNIT

        # DEFINES A SNAIL PARAMETRIC AMPLIFIER BY SHUNTING THE SNAIL WITH A CAPACITOR
        self.spa = Nlosc(self.squid)
        self.spa.C = self.C / C_UNIT
        self.Npoints = len(self.spa.nlind.flux)

    def copy(self):

        copy = rfSQUIDresonator(L=self.L, C=self.C, array_num=self.array_num, order=self.order)
        copy.define_Ic(self.Ic)
        copy.define_Lshunt(self.Lshunt)

        return copy

    def tune_freq_with_C(self, target_freq):

        def cost_function(C, args):
            target_freq = args[0]
            acopy = args[1]
            acopy.change_C(C)

            omega, g3, g4, g5, g6, g7, g8, g9, g10 = acopy.get_gns_at_flux(acopy.estimate_kerrnull_flux())
            cost = np.abs(omega / (2 * np.pi) - target_freq)
            return cost

        copy = self.copy()
        copy.try_to_build_model()
        args = [target_freq, copy]
        x0 = [copy.C]

        result = scipy.optimize.minimize(cost_function, x0, args=args, method='Nelder-Mead', bounds=[(1e-15, 1e-6)],
                                         tol=1e-17)
        return result.x[0]

    def tune_freq_with_L(self, target_freq):

        def cost_function(L, args):
            target_freq = args[0]
            acopy = args[1]
            acopy.change_L(L)

            omega, g3, g4, g5, g6, g7, g8, g9, g10 = acopy.get_gns_at_flux(acopy.estimate_kerrnull_flux())
            cost = np.abs(omega / (2 * np.pi) - target_freq)
            return cost

        copy = self.copy()
        copy.try_to_build_model()
        args = [target_freq, copy]
        x0 = [copy.L]

        result = scipy.optimize.minimize(cost_function, x0, args=args, method='Nelder-Mead', bounds=[(1e-13,1e-3)], tol=1e-2)
        return result.x

if __name__ == '__main__':
    plt.close('all')
    squidres = rfSQUIDresonator(L=0.01e-9, C=1.5e-12, array_num=20)
    squidres.define_Lj(2e-9)
    squidres.define_Lshunt(1e-9)
    squidres.try_to_build_model()
    result = squidres.tune_freq_with_C(5e9)
    squidres.change_C(result)
    omega, g3, g4, g5, g6, g7, g8, g9, g10 = squidres.get_gns_at_flux(squidres.estimate_kerrnull_flux(), plot=True)
    transitions = squidres.get_transitions_at_flux(squidres.estimate_kerrnull_flux()-0.005, 500, plot=True)
    transitions_spread = np.max(transitions) - np.min(transitions)
    metric = np.abs(transitions_spread / g3)
    print(metric)

