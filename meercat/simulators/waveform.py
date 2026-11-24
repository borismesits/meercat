import numpy as np

# all waveforms normalized to max amplitude of 1

# these are generic
def const_waveform(t, start, length):

    wave = t*0
    wave[np.where(start < t and t < (start + length))] = 1
    return wave


def tanh_waveform(t, start, length, ramp):
    return (np.tanh((t - start) / ramp) - np.tanh((t - (start + length)) / ramp)) / 2


def sechtanh_waveform(t, start, length, mu):
    beta = 1 / length
    t0 = start + length/2
    phi = mu * np.log(np.cosh((t - t0) * beta))
    return 1 / np.cosh((t - t0) * beta) * np.exp(-1j * phi)


# these classes are used since one might like several different types of, say, tanh pulse, in one program
class Waveform:
    def __init__(self):
        pass

    def return_waveform(self):
        pass


class ConstWaveform(Waveform):
    def __init__(self, length):
        super().__init__()
        self.length = length

    def return_waveform(self, t, start, length):
        return const_waveform(t, start, length)


class TanhWaveform(Waveform):
    def __init__(self, length, ramp):
        super().__init__()
        self.length = length
        self.ramp = ramp

    def return_waveform(self, t, start, length):
        return tanh_waveform(t, start, length, self.ramp)


class SechTanhWaveform(Waveform):
    def __init__(self, length, mu):
        super().__init__()
        self.length = length
        self.mu = mu

    def return_waveform(self, t, start, length):
        return sechtanh_waveform(t, start, length, self.mu)



