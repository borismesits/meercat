from meercat.simulators.meq import MasterEquationSim
from meercat.simulators.meq_program import Program
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from meercat.util.plot import plot_expect_with_variance
import matplotlib as mpl
mpl.use('Qt5Agg')
from tqdm.contrib import itertools

class Sqz_Disp_Meas(Program):

    def __init__(self):
        super().__init__()

    def body(self):
        self.delay(time=0.1)
        self.play_arb_pulse('DISP_O')
        self.delay()
        self.play_arb_pulse('CONV_O_S')

K4 = 1 * 2 * np.pi
qb_anharm = 275*2*np.pi
chi = 1.7*2*np.pi

sim = MasterEquationSim()
sim.add_mode('Qubit', 3,  K4=qb_anharm)
sim.add_mode('Readout', 20)
sim.add_mode('SNAIL', 20, K4=K4)
sim.add_mode('Output', 20)
q, r, s, o = sim.build_H0()
sim.add_cross_Kerr('Qubit','Readout', chi=chi)

DISP_O_eps = 20
g4 = 0.01 * 2 * np.pi
DISP_SUBH_eta = 4

g3 = 10 * 2 * np.pi
O_S_gdelta = 0.1
S_R_gdelta = 0.1
SQZ_eta = 0.1

CONV_O_S_eta = 0.8
CONV_S_R_eta = 0.8

# demo sequence
if True:
    prog = Sqz_Disp_Meas()
    prog.add_pulse_to_library('SUBH_S_G4', 6 * DISP_SUBH_eta ** 3 * g4 * s, 0.05, 0.001, 0, 0)
    prog.add_pulse_to_library('DISP_O', DISP_O_eps * o, 0.1, 0.01, 0, 0)
    prog.add_pulse_to_library('SQZ', 3*g3*s*s*SQZ_eta, 0.025, 0.001, 0, 0)
    prog.add_pulse_to_library('CONV_O_S', 3*g3*O_S_gdelta*s.dag()*o*CONV_O_S_eta, 0.1, 0.01, 0, 0, stark_op=2*CONV_O_S_eta**2*s.dag()*s*K4)
    prog.add_pulse_to_library('CONV_S_R', 3*g3*S_R_gdelta*s.dag()*r*CONV_S_R_eta, 0.1, 0.01, 0, 0, stark_op=2*CONV_S_R_eta**2*s.dag()*s*K4)

    sim.set_coherent_IC('Readout', 0)  # change qubit initial condition here
    sim.set_fock_IC('Qubit', 0) # change qubit initial condition here

    t = np.linspace(0, 0.5, 1001)
    qt.settings.auto_tidyup = False
    result = sim.simulate(t, prog)

    plot_expect_with_variance(3, t, [q+q.dag(), 1j*q-1j*q.dag(), q.dag()*q], ['xq', 'pq', 'nq'], result.states)
    plot_expect_with_variance(20, t, [r + r.dag(), 1j * r - 1j * r.dag(), r.dag() * r], ['xr', 'pr', 'nr'],
                              result.states)
    plot_expect_with_variance(20, t, [s + s.dag(), 1j * s - 1j * s.dag(), s.dag() * s], ['xs', 'ps', 'ns'],
                              result.states)
    plot_expect_with_variance(29, t, [o+o.dag(), 1j*o-1j*o.dag(), o.dag()*o], ['xo', 'po', 'no'], result.states)

