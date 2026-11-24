import matplotlib.pyplot as plt
import numpy as np
from qutip import (about, basis, coherent, destroy, displace, expect, mesolve,
                   num, plot_fock_distribution, plot_wigner, squeeze, variance, Options)

N = 30

# operators
a = destroy(N)
n = num(N)
x = (a + a.dag()) / np.sqrt(2)
p = -1j * (a - a.dag()) / np.sqrt(2)

disp = 0
phi = 0.0
Delta = 50*2*np.pi
epsilon_R = 0.0
epsilon_I = 0.3
omega_d = 50*2*np.pi
K = 5*2*np.pi

tlist = np.linspace(0, 1, 1001)

varepsilon = 5*np.exp(1j*(omega_d*tlist + phi))

H0 = Delta*a.dag()*a

# H = [H0, [a, np.conjugate(varepsilon)], [a.dag(), varepsilon]]
H = [H0, [a, np.exp(-1j*(omega_d*tlist + phi))]]

c_ops = [0.0*a]

N = 30

H0 = Delta*a.dag()*a + K*a.dag()*a.dag()*a*a

H1 = a.dag()-a

H = [H0, [H1, varepsilon-np.conjugate(varepsilon)]]

c_ops = [0.0*a]

psi0 = coherent(N, disp)

result = mesolve(H, psi0, tlist, c_ops, [])

plt.figure(figsize=(10,7))
plt.plot(tlist,expect(a.dag()*a, result.states),'k-',label='<n>',linewidth=1)

plt.legend()
plt.grid()
plt.xlabel('Time')
plt.show()
