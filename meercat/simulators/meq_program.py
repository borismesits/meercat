import numpy as np

class Program():

    def __init__(self):
        self.pulse_dict = {}
        self.sequence_dict = {}
        self.waveform_dict = {}

        self.time_counter = 0
        self.last_pulse_counter = 0
        self.pulse_num_counter = 0

    def play_arb_pulse(self, pulse_name, t=0):
        pulse = self.pulse_dict[pulse_name]
        # TODO, add options to change non-waveform parameters here. useful for something like X vs Y gates, just change pulse phase
        self.sequence_dict[self.pulse_num_counter] = pulse.copy()
        # if phase != None:
        #     self.sequence_dict[self.pulse_num_counter]['phase'] = phase
        self.sequence_dict[self.pulse_num_counter]['start'] = self.time_counter + t
        self.pulse_num_counter += 1
        self.last_pulse_counter = np.max([self.last_pulse_counter, pulse['waveform'].length+t])

    def delay(self, t=0):
        self.time_counter += self.last_pulse_counter
        self.last_pulse_counter = 0
        self.time_counter += t
        pass

    def add_pulse_to_library(self, name, op, waveform, freq=0, phase=0, ampl=1, stark_op2=None, stark_op4=None, eta_scaling=1):

        self.pulse_dict[name] = {'op': op, # the main operator of the parametric process
                                 'waveform': waveform,
                                 'freq': freq,
                                 'phase': phase,
                                 'ampl': ampl, # equivalent to eta
                                 'stark_op2': stark_op2, # this operator gets activated with the mod square of pump amplitude
                                 'stark_op4': stark_op4, # this operator gets activated with the mod^4 of pump amplitude
                                 'eta_scaling': eta_scaling} # some parametric processes scale quadratically or otherwise pump amplitude

    def plot_program(self):
        pass # TODO

    def get_final_time(self):
        self.body()
        return self.time_counter

    def body(self):
        # resets the counter and sequence dict, in case body() is called multiple times
        self.sequence_dict = {}
        self.time_counter = 0
