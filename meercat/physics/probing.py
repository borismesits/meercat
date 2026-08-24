import numpy as np

def RN2Lj(RN):
    e = 1.6e-19
    Delta = 180e-6 * e
    hbar = 6.63e-34 / (2 * np.pi)

    Lj = hbar * RN / (np.pi * Delta)

    return Lj


def RN2Ic(RN):
    Phi0 = 2.067e-15
    phi0 = Phi0 / (2 * np.pi)
    Lj = RN2Lj(RN)
    return phi0 / Lj


print(RN2Ic(3000-180))