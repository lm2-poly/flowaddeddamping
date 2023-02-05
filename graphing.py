#Plotting the important data for the article
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-09-05

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

matplotlib.rcParams.update({'font.size': 18})

nmodes = 3
lines = ['-','--','-.']
widths = [1.0, 1.0, 1.0]

#Structure typique d'une ligne: [U, UR, omega, zeta, Re(lambda), Im(lambda)]


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

#Plots

#Part 1: Multiple modes

plt.figure("Multiple modes")

ax = plt.subplot(2, 3, 1)
plt.xlim(0,max(NACA_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][:,2]),NACA_nas[i][:,2]/NACA_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(NACA_nas[i][0,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][0,2]),NACA_nas[i][0,2]/NACA_nas[0][0,2]+0.05,"Hydroelastic mode "+str(i+1))
plt.ylabel(r'$\Omega_i$')
plt.ylim([0.6,4])

ax = plt.subplot(2, 3, 2)
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1]*1/(F0_nas[0][0,2]/F0_nas[i][:,2]),F0_nas[i][:,2]/F0_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(F0_nas[i][0,1]*1/(F0_nas[0][0,2]/F0_nas[i][0,2]),F0_nas[i][0,2]/F0_nas[0][0,2]+0.1,"Hydroelastic mode "+str(i+1))
plt.xlim(0,max(F0_nas[0][:,1]))
plt.ylim([0.4,4])

ax = plt.subplot(2, 3, 3)
plt.xlim(0,max(F1_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1]*1/(F1_nas[0][0,2]/F1_nas[i][:,2]),F1_nas[i][:,2]/F1_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(F1_nas[i][0,1]*1/(F1_nas[0][0,2]/F1_nas[i][0,2]),F1_nas[i][0,2]/F1_nas[0][0,2]+0.1,"Hydroelastic mode "+str(i+1))
plt.ylim([0.4,4])




plt.subplot(2, 3, 4)
plt.xlim(0,max(NACA_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][:,2]),NACA_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')
plt.ylabel(r'$\zeta_{i,added}$')



plt.subplot(2, 3, 5)
plt.xlim(0,max(F0_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1]*1/(F0_nas[0][0,2]/F0_nas[i][:,2]),F0_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')


plt.subplot(2, 3, 6)
plt.xlim(0,max(F1_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1]*1/(F1_nas[0][0,2]/F1_nas[i][:,2]),F1_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')


#With and without added mass vs experimental/numerical

nmodes = 1

plt.figure("Comparison")

ax = plt.subplot(2, 3, 1)
plt.xlim(0,max(NACA_nas[0][:,1]))
plt.scatter(NACA_exp[:,1], NACA_exp[:,2]/NACA_exp[0,2], marker="x", c='black', label="Numerical")
plt.arrow(1, 0.97,3,-0.1)
plt.text(5, 0.8,"Numerical")
plt.scatter([None],[None], marker = "v", c="black", label = "Numerical - shifted")
plt.plot(NACA_nas[0][:,1]*1/(NACA_nas[0][0,2]/NACA_nas[0][:,2]),NACA_nas[0][:,2]/NACA_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(NACA_nas[0][0,1], 1.05,"Hydroelastic mode 1")
plt.ylabel(r'$\Omega_i$')
plt.ylim([0.6,1.5])

ax = plt.subplot(2, 3, 2)
plt.scatter(F0_exp[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_exp[:,2]/F0_exp[0,2], marker="o", c='black', label="Experimental")
plt.arrow(0.1, 1.05,0.2, 0.1)
plt.text(0.3, 1.2,"Experimental")
plt.scatter([None],[None], marker = "^", c="black", label = "Experimental - shifted")
fviv = 11/(0.012*F0_exp[0,2])
plt.plot([0,fviv,2*fviv,3*fviv],[0,1,2,3], c='black', linestyle = ':', linewidth = 1.0)
plt.text(0.85, 0.5, "Vortex shedding frequency")
plt.plot(F0_nas[0][:,1]*F0_nas[0][:,2]/F0_nas[0][0,2],F0_nas[0][:,2]/F0_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(1.5,0.83,'Hydroelastic mode 1')
plt.xlim(0,max(F0_nas[0][:,1]))
plt.ylim([0.4,1.5])

ax = plt.subplot(2, 3, 3)
plt.xlim(0,max(F1_nas[0][:,1]))
plt.scatter(F1_exp[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_exp[:,2]/F1_exp[0,2], marker="s", c='black', label="Experimental")
plt.arrow(0.1, 1.03,0.175,0.1)
plt.text(F1_exp[0,1]+0.15, 1.21,"Experimental")
plt.scatter([None],[None], marker = "P", c="black", label = "Experimental - shifted")
fviv = 8/(0.012*F1_exp[0,2])
plt.plot([0,fviv,2*fviv,3*fviv,4*fviv],[0,1,2,3,4], c='black', linestyle = ':', linewidth = 1.0)
plt.text(0.8, 0.5, "Vortex shedding frequency")
plt.plot(F1_nas[0][:,1]*F1_nas[0][:,2]/F1_nas[0][0,2],F1_nas[0][:,2]/F1_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(2, 0.83, "Hydroelastic mode 1")
plt.ylim([0.4,1.5])




plt.subplot(2, 3, 4)
plt.xlim(0,max(NACA_nas[0][:,1]))
plt.scatter(NACA_exp[:,1], NACA_exp[:,3], marker="x", c='black', label="Numerical")
plt.scatter(NACA_shifted[:,1], NACA_shifted[:,3], marker="v", c='black', label="Numerical - shifted")
plt.arrow(NACA_shifted[0,1]+1, NACA_shifted[0,3]-0.01,3,-0.03)
plt.text(NACA_shifted[0,1]+5, NACA_shifted[0,3]-0.05,"Numerical - shifted")
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][:,2]),NACA_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')
plt.ylabel(r'$\zeta_{i,added}$')



plt.subplot(2, 3, 5)
plt.xlim(0,max(F0_nas[0][:,1]))
plt.scatter(F0_exp[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_exp[:,3], marker="o", c='black', label="Experimental")
plt.scatter(F0_shifted[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_shifted[:,3], marker="^", c='black', label="Experimental - shifted")
plt.arrow(F0_shifted[1,1]+0.1, F0_shifted[1,3]-0.001,0.3,-0.003, width=0.0000005)
plt.text(F0_shifted[1,1]+0.4, F0_shifted[1,3]-0.005,"Experimental - shifted")
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1]*1/(F0_nas[0][0,2]/F0_nas[i][:,2]),F0_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')


plt.subplot(2, 3, 6)
plt.xlim(0,max(F1_nas[0][:,1]))
plt.scatter(F1_exp[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_exp[:,3], marker="s", c='black', label="Experimental")
plt.scatter(F1_shifted[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_shifted[:,3], marker="P", c='black', label="Experimental - shifted")
plt.arrow(F1_shifted[1,1]+0.1, F1_shifted[1,3]-0.001,0.3,-0.003, width=0.0000005)
plt.text(F1_shifted[1,1]+0.4, F1_shifted[1,3]-0.005,"Experimental - shifted")
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1]*1/(F1_nas[0][0,2]/F1_nas[i][:,2]),F1_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$')


plt.show()