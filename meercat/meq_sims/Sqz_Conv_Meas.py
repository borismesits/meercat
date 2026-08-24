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

sim = MasterEquationSim()
sim.add_mode('Qubit', 3)
sim.add_mode('Readout', 20)
sim.add_mode('SNAIL', 20, K4=K4)
sim.add_mode('Output', 20)
q, r, s, o = sim.build_H0()

DISP_O_eps = 20
g4 = 0.01 * 2 * np.pi
DISP_SUBH_eta = 4

g3 = 1 * 2 * np.pi
SQZ_eta = 1

CONV_O_S_eta = 1.5

# demo sequence
if True:
    prog = Sqz_Disp_Meas()
    prog.add_pulse_to_library('SUBH_S_G4', 6 * DISP_SUBH_eta ** 3 * g4 * s, 0.05, 0.001, 0, 0)
    prog.add_pulse_to_library('DISP_O', DISP_O_eps * o, 0.1, 0.01, 0, 0)
    prog.add_pulse_to_library('SQZ', 3 * g3 * s * s * SQZ_eta, 0.025, 0.001, 0, 0)
    prog.add_pulse_to_library('CONV_O_S', 3*g3*s.dag()*o*CONV_O_S_eta, 0.3, 0.01, 0, 0, stark_op=2*CONV_O_S_eta**2*s.dag()*s*K4)
    prog.add_pulse_to_library('CONV_O_S', 3 * g3 * s.dag() * o * CONV_O_S_eta, 0.3, 0.01, 0, 0, stark_op=2 * CONV_O_S_eta ** 2 * s.dag() * s * K4)

    sim.set_fock_IC('Qubit', 1)

    t = np.linspace(0, 0.5, 1001)
    qt.settings.auto_tidyup = False
    result = sim.simulate(t, prog)

    plot_expect_with_variance(3, t, [q+q.dag(), 1j*q-1j*q.dag(), q.dag()*q], ['xq', 'pq', 'nq'], result.states)
    plot_expect_with_variance(20, t, [r + r.dag(), 1j * r - 1j * r.dag(), r.dag() * r], ['xr', 'pr', 'nr'],
                              result.states)
    plot_expect_with_variance(20, t, [s + s.dag(), 1j * s - 1j * s.dag(), s.dag() * s], ['xs', 'ps', 'ns'],
                              result.states)
    plot_expect_with_variance(29, t, [o+o.dag(), 1j*o-1j*o.dag(), o.dag()*o], ['xo', 'po', 'no'], result.states)

# amplitude chevron
if False:
    CONV_O_S_etas = np.linspace(0, 5, 31)
    CONV_O_S_freqs = np.linspace(-200, 200, 31)

    final_o_means = np.zeros([len(CONV_O_S_etas), len(CONV_O_S_freqs)])
    final_o_rms = np.zeros([len(CONV_O_S_etas), len(CONV_O_S_freqs)])
    final_o_i = np.zeros([len(CONV_O_S_etas), len(CONV_O_S_freqs)])
    final_o_q = np.zeros([len(CONV_O_S_etas), len(CONV_O_S_freqs)])

    for i, j in itertools.product(range(0, len(CONV_O_S_etas)), range(0, len(CONV_O_S_freqs)), colour='blue'):
        prog = Sqz_Disp_Meas()

        CONV_O_S_eta = CONV_O_S_etas[i]
        CONV_O_S_freq = CONV_O_S_freqs[j]
        prog.add_pulse_to_library('SUBH_S_G4', 6 * DISP_SUBH_eta ** 3 * g4 * s, 0.05, 0.001, 0, 0, stark_op=2*DISP_SUBH_eta**2*s.dag()*s*K4)
        prog.add_pulse_to_library('DISP_O', DISP_O_eps * o, 0.1, 0.01, 0, 0)
        prog.add_pulse_to_library('SQZ', 3 * g3 * s * s * SQZ_eta, 0.025, 0.001, 0, 0, stark_op=2*SQZ_eta**2*s.dag()*s*K4)
        prog.add_pulse_to_library('CONV_O_S',3 * g3 * s.dag() * o * CONV_O_S_eta, 0.05, 0.001, CONV_O_S_freq, 0, stark_op=2*CONV_O_S_eta**2*s.dag()*s*K4)

        t = np.linspace(0, 0.5, 1001)
        qt.settings.auto_tidyup = False
        result = sim.simulate(t, prog)

        final_o_i[i, j] = qt.expect(o+o.dag(), result.states[-1])
        final_o_q[i, j] = qt.expect(1j*o-1j*o.dag(), result.states[-1])
        final_o_means[i, j] = np.abs(qt.expect(o, result.states[-1]))
        final_o_rms[i, j] = np.sqrt(qt.expect(o.dag()*o, result.states[-1]))

    fig, axs = plt.subplots(2,2)
    axs[0, 0].pcolor(CONV_O_S_freqs, CONV_O_S_etas, final_o_i)
    axs[1, 0].pcolor(CONV_O_S_freqs, CONV_O_S_etas, final_o_q)
    axs[0, 1].pcolor(CONV_O_S_freqs, CONV_O_S_etas, final_o_means)
    axs[1, 1].pcolor(CONV_O_S_freqs, CONV_O_S_etas, final_o_rms)
