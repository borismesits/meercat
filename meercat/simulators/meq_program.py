import numpy as np

class Program():

    def __init__(self):
        self.pulse_dict = {}
        self.sequence_dict = {}

        self.time_counter = 0
        self.last_pulse_counter = 0
        self.pulse_num_counter = 0

    def play_arb_pulse(self, pulse_name, t=0):
        self.sequence_dict[self.pulse_num_counter] = self.pulse_dict[pulse_name].copy()
        self.sequence_dict[self.pulse_num_counter]['start'] = self.time_counter + t
        self.pulse_num_counter += 1
        self.last_pulse_counter = np.max([self.last_pulse_counter, self.pulse_dict[pulse_name]['length']+t])
        print(self.time_counter)

    def delay(self, time=0):
        self.time_counter += self.last_pulse_counter
        self.last_pulse_counter = 0
        self.time_counter += time
        print(self.time_counter)
        pass

    def add_pulse_to_library(self, name, op, length, ramp, freq, phase, stark_op=None):
        self.pulse_dict[name] = {'op': op,
                                 'length': length,
                                 'ramp': ramp,
                                 'freq': freq,
                                 'phase': phase,
                                 'stark_op': stark_op}

    def body(self):
        pass
