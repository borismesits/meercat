import qutip as qt
import numpy as np
import matplotlib.pyplot as plt

class MasterEquationSim():

    def __init__(self):

        self.modes = {}
        self.c_ops = []

    def add_mode(self, name, num_levels, omega0=0, K4=0, K6=0, K8=0, kappa=0, gamma_phi=0, fock0=0):

        self.modes[name] = {'num_levels': num_levels,
                            'omega0': omega0,
                            'K4': K4,
                            'K6': K6,
                            'K8': K8,
                            'kappa': kappa,
                            'gamma_phi': gamma_phi,
                            'psi0': qt.basis(num_levels, fock0)}

    def add_cross_Kerr(self, mode1_name, mode2_name, chi=1):

        op1 = self.modes[mode1_name]['op']
        op2 = self.modes[mode2_name]['op']

        self.H0 += op1.dag()*op1*op2.dag()*op2*chi

    def build_H0(self):

        self.num_modes = len(self.modes.keys())
        ops = []

        for mode_name in self.modes: #this builds up a hilbert space with arbitrary number of subspaces
            sub_ops = []
            for mode_name_2 in self.modes:
                if mode_name == mode_name_2:
                    sub_ops.append(qt.destroy(self.modes[mode_name_2]['num_levels']))
                else:
                    sub_ops.append(qt.qeye(self.modes[mode_name_2]['num_levels']))

            self.modes[mode_name]['op'] = qt.tensor(sub_ops)
            ops.append(self.modes[mode_name]['op'])

        for mode_name in self.modes:

            op = self.modes[mode_name]['op']
            omega0 = self.modes[mode_name]['omega0']
            K4 = self.modes[mode_name]['K4']
            K6 = self.modes[mode_name]['K6']
            K8 = self.modes[mode_name]['K8']

            try:
                self.H0 += omega0*op.dag() * op + K4 * op.dag()**2 * op**2 + K6 * op.dag()**3 * op**3 + K8 * op.dag()**4 * op**4
            except:
                self.H0 = omega0*op.dag() * op + K4 * op.dag()**2 * op**2 + K6 * op.dag()**3 * op**3 + K8 * op.dag()**4 * op**4

        return ops

    def set_fock_IC(self, mode_name, fock_num):
        self.modes[mode_name]['psi0'] = qt.basis(self.modes[mode_name]['num_levels'], fock_num)

    def set_coherent_IC(self, mode_name, alpha):
        self.modes[mode_name]['psi0'] = qt.coherent(self.modes[mode_name]['num_levels'], alpha)

    def simulate(self, t, program, c_ops):

        self.program = program
        self.program.body()
        self.build_TD_H(t)

        sub_psi0s = [] # build up the initial condition
        for mode_name in self.modes:
            sub_psi0 = self.modes[mode_name]['psi0']
            sub_psi0s.append(sub_psi0)
        self.psi0 = qt.tensor(sub_psi0s)

        # for mode_name in self.modes:
        #     op = self.modes[mode_name]['kappa']
        #     if self.modes[mode_name]['kappa'] > 0:
        #         se

        qt.settings.auto_tidyup = False
        self.result = qt.mesolve(self.H, self.psi0, t, c_ops)
        return self.result

    def build_TD_H(self, t):
        '''
        TODO: there's a detail that is wrong and I don't know how important it is
        For a subharmonic drive, the drive rate is not proportional to eta (as well as other parametric processes)
        This means that, for a ramped pulse, the ramp of the subh rate is not the same as the ramp of the pulse
        But I assume that they are the same here
        I'm ignoring because this is intended to be a rough simulation at the moment
        '''
        self.H = [self.H0]

        for pulse_num in self.program.sequence_dict:
            pulse = self.program.sequence_dict[pulse_num]
            pulse_waveform = pulse['waveform']
            envelope = pulse['ampl']*pulse_waveform.return_waveform(t, pulse['start'], pulse['waveform'].length) * np.exp(-1j*(pulse['freq']*t+pulse['phase']))
            op = pulse['op']
            stark_op = None

            if np.max(np.abs(op.full())) > 0: # something breaks when adding zero operators to TD solver, so just check here
                H_td = [op.dag(), envelope]
                H_td_c = [op, np.conjugate(envelope)]

                self.H.append(H_td)
                self.H.append(H_td_c)

                if stark_op != None:
                    H_td_stark = [stark_op, np.abs(envelope)**2]
                    self.H.append(H_td_stark)

    def plot_transitions_for_modes(self):
        pass # TODO

if __name__ == '__main__':

    sim = MasterEquationSim()
    sim.add_mode('Amplifier', 3)
    sim.add_mode('Output', 3)
    a, o = sim.build_H0()
    sim.add_cross_Kerr('Amplifier', 'Output', chi=1)

