import matplotlib.pyplot as plt
import numpy as np

r = np.linspace(0,3,100)
X = lambda zeta: 1/np.sqrt((1-r**2)**2+(2*zeta*r)**2)

zetas = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0]
lines = ['-', '--', '-.', 'dotted', (0, (5,5)), (0, (1,1))]

plt.figure("AmplitudeVZeta")
plt.xlabel(r'$r = \frac{\omega}{\omega_n}$')
plt.ylabel(r'$\~ X$')
for i, zeta in enumerate(zetas):
    plt.plot(r, X(zeta), label = r'$\zeta = $'+str(zeta), color = 'black', linestyle = lines[i])
plt.ylim([0, 5])
plt.xlim([0,3])
plt.legend()
plt.show()
    