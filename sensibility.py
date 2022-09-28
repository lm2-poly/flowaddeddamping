#Performing a sensitivity analysis
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-09-05

import matplotlib.pyplot as plt
import pandas as pd

filepath = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\sensibility\results\NACA00'
suffix = ['02-100','03-100','04-100','02-90','03-90','04-90','02-80','03-80','04-80']
lines = ['solid', 'dashed', 'dashdot', 'dotted', (0, (1,10)), (0, (5,10)), (0,(1,1)), (0,(5,1)), (0, (3, 5, 1, 5, 1, 5))]
thick = [2,3,4]
chord = [100,90,80]
markers = ['o','s','v','^','*','D','X','+','1']

lengths = []
slope = []

plt.figure("Sensitivity analysis")
plt.ylabel(r'$\zeta_{added}$')
plt.xlabel(r'$V$ [m/s]')
plt.title("Adimensional added damping according to velocity of different airfoils")
plt.grid(True)
for i in range(len(suffix)):
    thickness = thick[i%3]
    length = chord[i//3]
    filename = [filepath,suffix[i],'.csv']
    filename = "".join(filename)
    data = pd.read_csv(filename, sep=',', header=None).to_numpy()
    data[:,0] = data[:,0]/1000
    plt.plot(data[:,0], data[:,3], linestyle = lines[i], label = 'Relative thickness = '+str(thickness)+'%, effective length = '+str(length)+'%')
    lengths.append(length)
    slope.append((data[20,3]-data[5,3])/(data[20,0]-data[5,0]))
plt.legend(loc='best')

plt.figure("slope vs ERT")
plt.xlabel("Effective length [%]")
plt.ylabel(r"$\zeta / V$ [s/m]")
plt.title("Added damping to velocity slope according to effective length and relative thickness")
plt.grid(True)
for i in range(len(lengths)):
    thickness = thick[i % 3]
    length = chord[i // 3]
    plt.scatter(lengths[i],slope[i], marker = markers[i%3], c="black")
for i in range(len(thick)):
    thickness = thick[i]
    lengthis = [lengths[i],lengths[3+i],lengths[6+i]]
    slopes = [slope[i], slope[3+i], slope[6+i]]
    plt.plot(lengthis,slopes,label = 'Relative thickness = '+str(thickness)+'%', linestyle = lines[i])
plt.legend(loc='best')
plt.show()
