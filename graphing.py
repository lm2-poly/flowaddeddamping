#Plotting the important data for the article
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-09-05

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

nmodes = 3
lines = ['-','--','-.']
widths = [1.0, 1.0, 1.0]

#Structure typique d'une ligne: [U, UR, omega, zeta, Re(lambda), Im(lambda)]

plt.figure("Damping")
NACA_exp = r'TestFiles\Results\NACA0003_exp.csv'
NACA_exp = pd.read_csv(NACA_exp, sep=';', header = None).to_numpy()
NACA_exp[:,1] = NACA_exp[:,1]*95/3.4
NACA_shifted = NACA_exp.copy()
NACA_shifted[:,3] -= NACA_shifted[0,3]*np.ones(len(NACA_shifted[:,3]))
NACA_nas = []
for i in range(0,nmodes):
    NACA_nas.append(r'TestFiles\Results\NACA0003_Curves_hydroelastic_aeroelastic'+str(i)+'.csv')
    NACA_nas[i] = pd.read_csv(NACA_nas[i], sep=',', header = None).to_numpy()
    NACA_nas[i][:,1] = NACA_nas[i][:,1]*95/3.4
    
F0_exp = r'TestFiles\Results\F0_exp.csv'
F0_exp = pd.read_csv(F0_exp, sep=';', header = None).to_numpy()
F0_exp[:,1] = F0_exp[:,1]*250/12
F0_shifted = F0_exp.copy()
F0_shifted[:,3] -= F0_shifted[1,3]*np.ones(len(F0_shifted[:,3]))
F0_nas = []
for i in range(0,nmodes):
    F0_nas.append(r'TestFiles\Results\F0_Curves_hydroelastic_aeroelastic'+str(i)+'.csv')
    F0_nas[i] = pd.read_csv(F0_nas[i], sep=',', header = None).to_numpy()    
    F0_nas[i][:,1] = F0_nas[i][:,1]*250/12
  

F1_exp = r'TestFiles\Results\F1_exp.csv'
F1_exp = pd.read_csv(F1_exp, sep=';', header = None).to_numpy()
F1_exp[:,1] = F1_exp[:,1]*250/12
F1_shifted = F1_exp.copy()
F1_shifted[:,3] -= F1_shifted[1,3]*np.ones(len(F1_shifted[:,3]))
F1_nas = []
for i in range(0,nmodes):
    F1_nas.append(r'TestFiles\Results\F1_Curves_hydroelastic_aeroelastic'+str(i)+'.csv')
    F1_nas[i] = pd.read_csv(F1_nas[i], sep=',', header = None).to_numpy()
    F1_nas[i][:,1] = F1_nas[i][:,1]*250/12


ax = plt.subplot(2, 3, 1)
plt.xlim(0,max(NACA_nas[0][:,1]))
plt.scatter(NACA_exp[:,1], NACA_exp[:,2]/NACA_exp[0,2], marker="x", c='black', label="Numerical")
plt.arrow(NACA_exp[0,1]+2, NACA_exp[0,2]/NACA_exp[0,2]-0.1,3,-0.1)
plt.text(NACA_exp[0,1]+5, NACA_exp[0,2]/NACA_exp[0,2]-0.3,"Numerical")
plt.scatter([None],[None], marker = "v", c="black", label = "Numerical - shifted")
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][:,2]),NACA_nas[i][:,2]/NACA_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(NACA_nas[i][0,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][0,2]),NACA_nas[i][0,2]/NACA_nas[0][0,2]+0.05,"Hydroelastic mode "+str(i+1))
plt.ylabel(r'$\Omega_i$')
plt.ylim([0.6,4])
ax.set_xticklabels([])
#plt.legend()

ax = plt.subplot(2, 3, 2)
plt.scatter(F0_exp[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_exp[:,2]/F0_exp[0,2], marker="o", c='black', label="Experimental")
plt.arrow(F0_exp[0,1]+0.05, F0_exp[0,2]/F0_exp[0,2]+0.1,0.2,0.8)
plt.text(F0_exp[0,1]+0.3, 1.8,"Experimental")
plt.scatter([None],[None], marker = "^", c="black", label = "Experimental - shifted")
fviv = 11/(0.012*F0_exp[0,2])
plt.plot([0,fviv,2*fviv,3*fviv],[0,1,2,3], c='black', linestyle = ':', linewidth = 1.0)
plt.text(0.65*fviv, 2.1, "Vortex shedding frequency")
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1]*1/(F0_nas[0][0,2]/F0_nas[i][:,2]),F0_nas[i][:,2]/F0_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    if i != 0:
        plt.text(F0_nas[i][0,1]*1/(F0_nas[0][0,2]/F0_nas[i][0,2]),F0_nas[i][0,2]/F0_nas[0][0,2]+0.1,"Hydroelastic mode "+str(i+1))
    else:
        plt.text(1.75,0.7,'Hydroelastic mode 1')
plt.xlim(0,max(F0_nas[0][:,1]))
plt.ylim([0.4,4])
ax.set_xticklabels([])
#plt.legend()

ax = plt.subplot(2, 3, 3)
plt.xlim(0,max(F1_nas[0][:,1]))
plt.scatter(F1_exp[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_exp[:,2]/F1_exp[0,2], marker="s", c='black', label="Experimental")
plt.arrow(F1_exp[0,1]+0.1, F1_exp[0,2]/F1_exp[0,2]+0.1,0.175,0.4)
plt.text(F1_exp[0,1]+0.3, 1.5,"Experimental")
plt.scatter([None],[None], marker = "P", c="black", label = "Experimental - shifted")
fviv = 8/(0.012*F1_exp[0,2])
plt.plot([0,fviv,2*fviv,3*fviv,4*fviv],[0,1,2,3,4], c='black', linestyle = ':', linewidth = 1.0)
plt.text(2, 3.7, "Vortex shedding frequency")
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1]*1/(F1_nas[0][0,2]/F1_nas[i][:,2]),F1_nas[i][:,2]/F1_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    if i != 0:
        plt.text(F1_nas[i][0,1]*1/(F1_nas[0][0,2]/F1_nas[i][0,2]),F1_nas[i][0,2]/F1_nas[0][0,2]+0.1,"Hydroelastic mode "+str(i+1))
    else:
        plt.text(2, 0.7, "Hydroelastic mode 1")
ax.set_xticklabels([])
plt.ylim([0.4,4])
#plt.legend()




plt.subplot(2, 3, 4)
plt.title("a)",y=-0.32)
plt.xlim(0,max(NACA_nas[0][:,1]))
plt.scatter(NACA_exp[:,1], NACA_exp[:,3], marker="x", c='black', label="Numerical")
plt.scatter(NACA_shifted[:,1], NACA_shifted[:,3], marker="v", c='black', label="Numerical - shifted")
plt.arrow(NACA_shifted[0,1]+1, NACA_shifted[0,3]-0.01,3,-0.03)
plt.text(NACA_shifted[0,1]+5, NACA_shifted[0,3]-0.05,"Numerical - shifted")
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][:,2]),NACA_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')
plt.ylabel(r'$\zeta_{i,added}$')
plt.plot([0,60],[0,0], c='black')
plt.plot([0,60],[NACA_exp[0,3],NACA_exp[0,3]], c='black')
plt.text(45, 1.25*NACA_exp[0,3], r'$\zeta_s$')
import matplotlib.patches as patches
p1 = patches.FancyArrowPatch((45, 0), (45, NACA_exp[0,3]), arrowstyle='<->', mutation_scale=20)

#plt.figure('F0', figsize = (6,8))



plt.subplot(2, 3, 5)
plt.title("b)",y=-0.32)
plt.xlim(0,max(F0_nas[0][:,1]))
plt.scatter(F0_exp[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_exp[:,3], marker="o", c='black', label="Experimental")
plt.scatter(F0_shifted[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_shifted[:,3], marker="^", c='black', label="Experimental - shifted")
plt.arrow(F0_shifted[1,1]+0.1, F0_shifted[1,3]-0.001,0.3,-0.003, width=0.0000005)
plt.text(F0_shifted[1,1]+0.4, F0_shifted[1,3]-0.005,"Experimental - shifted")
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1]*1/(F0_nas[0][0,2]/F0_nas[i][:,2]),F0_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')
plt.plot([0,3.7],[0,0], c='black')
plt.plot([0,3.7],[F0_exp[1,3],F0_exp[1,3]], c='black')
plt.text(3, 1.25*F0_exp[1,3], r'$\zeta_s$')

#plt.figure('F1', figsize = (6,8))


plt.subplot(2, 3, 6)
plt.title("c)",y=-0.32)
plt.xlim(0,max(F1_nas[0][:,1]))
plt.scatter(F1_exp[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_exp[:,3], marker="s", c='black', label="Experimental")
plt.scatter(F1_shifted[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_shifted[:,3], marker="P", c='black', label="Experimental - shifted")
plt.arrow(F1_shifted[1,1]+0.1, F1_shifted[1,3]-0.001,0.3,-0.003, width=0.0000005)
plt.text(F1_shifted[1,1]+0.4, F1_shifted[1,3]-0.005,"Experimental - shifted")
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1]*1/(F1_nas[0][0,2]/F1_nas[i][:,2]),F1_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')
plt.plot([0,5.2],[0,0], c='black')
plt.plot([0,5.2],[F1_exp[1,3],F1_exp[1,3]], c='black')
plt.text(4.5, 1.25*F1_exp[1,3], r'$\zeta_s$')

plt.show()
exit()

plt.figure("Added damping")
plt.ylabel(r'$\zeta_{added}$')
plt.xlabel(r'$V$ [m/s]')
plt.title("Adimensional added damping according to velocity")
plt.grid(True)
plt.scatter(F1_exp[:,0],F1_exp[:,1], marker="o", c = 'C0', label = 'F1 Experimental')
plt.scatter(CFX_num[:,0],CFX_num[:,1],marker="s", c = 'C1', label = 'NACA0003 CFX')
plt.plot(F1_nas[:,0],F1_nas[:,1],linestyle = 'solid', c = 'C0', label = 'F1 NASTRAN')
plt.plot(CFX_nas[:,0],CFX_nas[:,1],linestyle = '--', c = 'C1', label = 'NACA0003 NASTRAN')
plt.legend(loc='best')
plt.xlim([0,30])
plt.show()