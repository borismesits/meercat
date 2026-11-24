from meercat.simulators.meq import MasterEquationSim
from meercat.simulators.waveform import *
from meercat.simulators.meq_program import Program
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from meercat.util.plot import plot_expect_with_variance
import matplotlib as mpl
mpl.use('Qt5Agg')
from tqdm.contrib import itertools

class Disp_Meas(Program):

    def __init__(self, config):
        super().__init__()
        self.cfg = config

    def body(self):
        super().body()
        self.delay(0.5)
        self.play_arb_pulse('DISP_R')
        if self.cfg['reverse'] == True:
            self.delay(0.5)
            self.play_arb_pulse('DISP_R_inv')
        self.delay(0.5)


config = {}
config['reverse'] = False

sim = MasterEquationSim()
sim.add_mode('Readout', 100, omega0=0 * 2 * np.pi, K4=0.2 * 2 * np.pi)

r = sim.build_H0()[0]

prog = Disp_Meas(config)
prog.add_pulse_to_library('DISP_R', 6 * r, TanhWaveform(0.3, 0.1),  ampl=1)
prog.add_pulse_to_library('DISP_R_inv', 6 * r, TanhWaveform(0.3, 0.1), ampl=1, phase=1*np.pi)

t = np.linspace(0, prog.get_final_time(), 101)
result = sim.simulate(t, prog, [])

from meercat.util.plot import plot_IQnbar_vs_time, make_Qfunc_movie
plot_IQnbar_vs_time(t, r, result.states, name='readout')

i=-1
fig, axes = plt.subplots(1, 2, figsize=(5, 2.5))
axes[0].cla()
qt.plot_fock_distribution(result.states[i], fig=fig, ax=axes[0])
qt.plot_wigner(result.states[i], fig=fig, ax=axes[1])
plt.show()
