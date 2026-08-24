import matplotlib.pyplot as plt
import numpy as np

def get_L_C_from_HFSS(LJ0, LJ1, f0, f1):
    '''
    LJ - the lumped L
    f0 - frequency when LJ is 0
    f1 - frequency when LJ is nominal
    '''

    w0 = f0*2*np.pi
    w1 = f1*2*np.pi
    L = (w1**2 *LJ1 - w0**2*LJ0)/np.abs(w0**2 - w1**2)
    C = 1/(w0**2*(L+LJ0))

    return L, C

def get_L_C_from_Z_f0(Z, f0):

    omega = f0*2*np.pi
    C = 1/(omega*Z)
    L = 1/(C*omega**2)

    return L, C

def get_L_C_from_LJ_Pm_f0(LJ, Pm, f0):
    omega = f0*2*np.pi
    L = (1-Pm)/Pm*LJ
    C = 1/((L+LJ)*omega**2)
    return L, C