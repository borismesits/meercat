import numpy as np
hbar = 1.0545718e-34
e = 1.60218e-19
Phi0 = np.pi*hbar/e
phi0 = hbar/(2*e)
'''
beware distinction between unreduced (Phi0) and reduced flux quantum (phi0)
'''

def Lj_to_Ej(Lj):
    return phi0**2/Lj

def Ej_to_Lj(Ej):
    return phi0**2/Ej

def Lj_to_Ic(Lj):
    return phi0/Lj

def Ic_to_Lj(Ic):
    return phi0/Ic

def Ic_to_Ej(Ic):
    return phi0*Ic

def Ej_to_Ic(Ic):
    return Ic/phi0