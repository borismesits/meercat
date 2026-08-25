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
        self.play_arb_pulse('SUBH_S_G4')
        self.delay()
        self.play_arb_pulse('CONV_S_R')
        self.delay()
        self.play_arb_pulse('CONV_S_R')
        self.delay()
        # self.delay(time=0.1)
        # self.play_arb_pulse('SUBH_S_G4')
        # self.delay()
        # self.play_arb_pulse('SQZ')
        # self.delay()
        # self.play_arb_pulse('CONV_O_S')


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

SUBH_eta = 4

g3 = 10 * 2 * np.pi
O_S_gdelta = 0.1
S_R_gdelta = 0.1
SQZ_eta = 0.1

CONV_O_S_eta = 0.8
CONV_S_R_eta = 0.8

# demo sequence
if True:
    prog = Sqz_Disp_Meas()
    prog.add_pulse_to_library('SUBH_S_G4', 6 * SUBH_eta ** 3 * g4 * s, 0.05, 0.001, 0, 0)
    prog.add_pulse_to_library('DISP_O', DISP_O_eps * o, 0.1, 0.01, 0, 0)
    prog.add_pulse_to_library('SQZ', 3*g3*s*s*SQZ_eta, 0.025, 0.001, 0, 0)
    prog.add_pulse_to_library('CONV_O_S', 3*g3*O_S_gdelta*s.dag()*o*CONV_O_S_eta, 0.1, 0.01, 0, 0, stark_op=2*CONV_O_S_eta**2*s.dag()*s*K4)
    prog.add_pulse_to_library('CONV_S_R', 3*g3*S_R_gdelta*s.dag()*r*CONV_S_R_eta, 0.1, 0.01, 0, 0, stark_op=2*CONV_S_R_eta**2*s.dag()*s*K4)

    sim.set_coherent_IC('Readout', 0)  # change qubit initial condition here
    sim.set_fock_IC('Qubit', 0) # change qubit initial condition here

    t = np.linspace(0, 0.5, 1001)
    qt.settings.auto_tidyup = False
    result = sim.simulate(t, prog)

    fig, axs = plt.subplots(nrows=4, ncols=3)
    axs[0, 0].plot(t, qt.expect(q + q.dag(), result.states))
    axs[0, 1].plot(t, qt.expect(1j * q - 1j * q.dag(), result.states))
    axs[0, 2].plot(t, qt.expect(q.dag()*q, result.states))
    axs[1, 0].plot(t, qt.expect(r + r.dag(), result.states))
    axs[1, 1].plot(t, qt.expect(1j * r - 1j * r.dag(), result.states))
    axs[1, 2].plot(t, qt.expect(r.dag()*r, result.states))
    axs[2, 0].plot(t, qt.expect(s + s.dag(), result.states))
    axs[2, 1].plot(t, qt.expect(1j * s - 1j * s.dag(), result.states))
    axs[2, 2].plot(t, qt.expect(s.dag()*s, result.states))
    axs[3, 0].plot(t, qt.expect(o + o.dag(), result.states))
    axs[3, 1].plot(t, qt.expect(1j * o - 1j * o.dag(), result.states))
    axs[3, 2].plot(t, qt.expect(o.dag()*o, result.states))

