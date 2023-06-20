#Plotting the important data for the article
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-09-05

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

matplotlib.rcParams.update({'font.size': 10})

nmodes = 3
lines = ['-','--','-.','-.',':']
widths = [1.0, 1.0, 1.0,1.0,1.0]

#Structure typique d'une ligne: [U, UR, omega, zeta, Re(lambda), Im(lambda)]


NACA_exp = r'TestFiles\Results\NACA0003_exp.csv'
NACA_exp = pd.read_csv(NACA_exp, sep=';', header = None).to_numpy()
NACA_exp[:,1] = NACA_exp[:,1]*95/3.4
NACA_shifted = NACA_exp.copy()
NACA_shifted[:,3] -= NACA_shifted[0,3]*np.ones(len(NACA_shifted[:,3]))
NACA_nas = []
for i in range(0,nmodes):
    NACA_nas.append(r'TestFiles\Results\naca_finer_hydroelastic_aeroelastic'+str(i)+'.csv')
    NACA_nas[i] = pd.read_csv(NACA_nas[i], sep=',', header = None).to_numpy()
    NACA_nas[i][:,1] = NACA_nas[i][:,1]*95/3.4

F0_exp = r'TestFiles\Results\F0_exp.csv'
F0_exp = pd.read_csv(F0_exp, sep=';', header = None).to_numpy()
F0_exp[:,1] = F0_exp[:,1]*250/12
F0_shifted = F0_exp.copy()
F0_shifted[:,3] -= F0_shifted[1,3]*np.ones(len(F0_shifted[:,3]))
F0_nas = []
for i in range(0,nmodes):
    F0_nas.append(r'TestFiles\Results\francis99_finer_hydroelastic_aeroelastic'+str(i)+'.csv')
    F0_nas[i] = pd.read_csv(F0_nas[i], sep=',', header = None).to_numpy()    
    F0_nas[i][:,1] = F0_nas[i][:,1]*250/12
  

F1_exp = r'TestFiles\Results\F1_exp.csv'
F1_exp = pd.read_csv(F1_exp, sep=';', header = None).to_numpy()
F1_exp[:,1] = F1_exp[:,1]*250/12
F1_shifted = F1_exp.copy()
F1_shifted[:,3] -= F1_shifted[1,3]*np.ones(len(F1_shifted[:,3]))
F1_nas = []
for i in range(0,nmodes):
    F1_nas.append(r'TestFiles\Results\bergan_precis_hydroelastic_aeroelastic'+str(i)+'.csv')
    F1_nas[i] = pd.read_csv(F1_nas[i], sep=',', header = None).to_numpy()
    F1_nas[i][:,1] = F1_nas[i][:,1]*250/12




blunt_exp = r'TestFiles\Results\Blunt_exp.csv'
blunt_exp = pd.read_csv(blunt_exp, sep=';', header = None).to_numpy()
blunt_exp[:,1] = blunt_exp[:,1]*10
blunt_shifted = blunt_exp.copy()
blunt_shifted[:,3] -= blunt_shifted[1,3]*np.ones(len(blunt_shifted[:,3]))
blunt_nas = []
for i in range(0,nmodes):
    blunt_nas.append(r'TestFiles\Results\blunt_precis_hydroelastic_aeroelastic'+str(i)+'.csv')
    blunt_nas[i] = pd.read_csv(blunt_nas[i], sep=',', header = None).to_numpy()
    blunt_nas[i][:,1] = blunt_nas[i][:,1]*100/10

blunt_simu = r'TestFiles\Results\Blunt_simu.csv'
blunt_simu = pd.read_csv(blunt_simu, sep=';', header = None).to_numpy()
blunt_simu[:,1] = blunt_simu[:,1]*100/10
blunt_simu_shifted = blunt_simu.copy()
blunt_simu_shifted[:,3] -= blunt_simu_shifted[1,3]*np.ones(len(blunt_simu_shifted[:,3]))

dona_nas = []
for i in range(0,nmodes):
    dona_nas.append(r'TestFiles\Results\dona_finer_hydroelastic_aeroelastic'+str(i)+'.csv')
    dona_nas[i] = pd.read_csv(dona_nas[i], sep=',', header = None).to_numpy()
    dona_nas[i][:,1] = dona_nas[i][:,1]*100/10

dona_simu = r'TestFiles\Results\Donaldson_simu.csv'
dona_simu = pd.read_csv(dona_simu, sep=';', header = None).to_numpy()
dona_simu[:,1] = dona_simu[:,1]*100/10
dona_shifted = dona_simu.copy()
dona_shifted[:,3] -= dona_shifted[1,3]*np.ones(len(dona_shifted[:,3]))
#Plots

#Part 1: Multiple modes

plt.figure("nacamodes",figsize=(6,9.5))


ax = plt.subplot(2, 1, 1)
plt.xlim(0,max(NACA_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1],NACA_nas[i][:,2]/NACA_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(NACA_nas[i][0,1]*1/(NACA_nas[0][0,2]/NACA_nas[i][0,2]),NACA_nas[i][0,2]/NACA_nas[0][0,2]+0.05,"Hydroelastic mode "+str(i+1),fontsize = 18)
plt.ylabel(r'$\Omega_i$',fontsize = 18)
plt.ylim([0.6,4])
plt.xticks([])

plt.subplot(2, 1, 2)
plt.xlim(0,max(NACA_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1],NACA_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.ylabel(r'$\zeta_{i,added}$',fontsize = 18)
plt.savefig('NACA_modes.pdf')

plt.figure("f0modes",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1],F0_nas[i][:,2]/F0_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(F0_nas[i][0,1]*1/(F0_nas[0][0,2]/F0_nas[i][0,2]),F0_nas[i][0,2]/F0_nas[0][0,2]+0.1,"Hydroelastic mode "+str(i+1),fontsize = 18)
plt.xlim(0,max(F0_nas[0][:,1]))
plt.ylim([0.4,4])
plt.xticks([])

plt.subplot(2, 1, 2)
plt.xlim(0,max(F0_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1],F0_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)

plt.savefig('F0_modes.pdf')


plt.figure("f1modes",figsize=(6,9.5))

ax = plt.subplot(2, 1, 1)
plt.xlim(0,max(F1_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1],F1_nas[i][:,2]/F1_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(F1_nas[i][0,1]*1/(F1_nas[0][0,2]/F1_nas[i][0,2]),F1_nas[i][0,2]/F1_nas[0][0,2]+0.1,"Hydroelastic mode "+str(i+1),fontsize = 18)
plt.ylim([0.4,4])
plt.xticks([])

ax = plt.subplot(2, 1, 2)
plt.xlim(0,max(F1_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1],F1_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)

plt.savefig('F1_modes.pdf')



plt.figure("bluntmodes",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
plt.xlim(0,max(blunt_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(blunt_nas[i][:,1],blunt_nas[i][:,2]/blunt_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(blunt_nas[i][0,1]*1/(blunt_nas[0][0,2]/blunt_nas[i][0,2]),blunt_nas[i][0,2]/blunt_nas[0][0,2]+0.25,"Hydroelastic mode "+str(i+1),fontsize = 18)

plt.ylim([0.6,8])
plt.xticks([])

plt.subplot(2, 1, 2)
plt.xlim(0,max(blunt_nas[0][:,1]))
for i in range(0,nmodes):
    plt.plot(blunt_nas[i][:,1],blunt_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.savefig('blunt_modes.pdf')

plt.figure("donamodes",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
for i in range(0,nmodes):
    plt.plot(dona_nas[i][:,1],dona_nas[i][:,2]/dona_nas[0][0,2], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
    plt.text(dona_nas[i][0,1]*1/(dona_nas[0][0,2]/dona_nas[i][0,2]),dona_nas[i][0,2]/dona_nas[0][0,2]+0.25,"Hydroelastic mode "+str(i+1),fontsize = 18)
plt.ylabel(r'$\Omega_i$',fontsize = 18)
plt.xlim(0,max(dona_nas[0][:,1]))
plt.ylim([0.4,8])
plt.xticks([])





plt.subplot(2,1 , 2)
plt.xlim(0,max(dona_nas[0][:,1]))
plt.ylim([-0.005,0.06])
for i in range(0,nmodes):
    plt.plot(dona_nas[i][:,1],dona_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.ylabel(r'$\zeta_{i,added}$',fontsize = 18)
plt.savefig('dona_modes.pdf')


nmodes = 1


plt.figure("nacacompa",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)

plt.xlim(0,max(NACA_nas[0][:,1]))
plt.scatter(NACA_exp[:,1], NACA_exp[:,2]/NACA_exp[0,2], marker="x", c='black', label="Numerical")
plt.arrow(1, 0.97,4.5,-0.25)
plt.text(6, 0.70,"Numerical",fontsize = 18)
plt.scatter([None],[None], marker = "v", c="black", label = "Numerical - shifted")
plt.plot(NACA_nas[0][:,1],NACA_nas[0][:,2]/NACA_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(NACA_nas[0][8,1], 0.85,"Hydroelastic mode 1",fontsize = 18)
plt.ylabel(r'$\Omega_i$',fontsize = 18)
plt.ylim([0.6,1.2])
plt.xticks([])

plt.subplot(2, 1, 2)
plt.xlim(0,max(NACA_nas[0][:,1]))

plt.scatter(NACA_shifted[:,1], NACA_shifted[:,3], marker="v", c='black', label="Numerical - shifted")
plt.arrow(NACA_shifted[1,1]+1, NACA_shifted[1,3]-0.01,3,-0.03)
plt.text(NACA_shifted[1,1]+4, NACA_shifted[1,3]-0.05,"Numerical - shifted",fontsize = 18)
for i in range(0,nmodes):
    plt.plot(NACA_nas[i][:,1],NACA_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.ylabel(r'$\zeta_{i,added}$',fontsize = 18)
plt.savefig('NACA_comparison.pdf')

plt.figure("f0compa",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
plt.scatter(F0_exp[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_exp[:,2]/F0_exp[0,2], marker="o", c='black', label="Experimental")
plt.arrow(0.1, 1.05,0.2, 0.1)
plt.text(0.3, 1.2,"Experimental",fontsize = 18)
plt.scatter([None],[None], marker = "^", c="black", label = "Experimental - shifted")
fviv = 11/(0.012*F0_exp[0,2])
plt.plot([0,fviv,2*fviv,3*fviv],[0,1,2,3], c='black', linestyle = ':', linewidth = 1.0)
plt.text(0.85, 0.5, "Vortex shedding frequency",fontsize = 18)
plt.plot(F0_nas[0][:,1],F0_nas[0][:,2]/F0_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(1.5,0.83,'Hydroelastic mode 1',fontsize = 18)
plt.xlim(0,max(F0_nas[0][:,1]))
plt.ylim([0.4,1.5])
plt.xticks([])

plt.subplot(2, 1, 2)
plt.xlim(0,max(F0_nas[0][:,1]))
plt.scatter(F0_shifted[:,1]*1/(F0_exp[0,2]/F0_exp[:,2]), F0_shifted[:,3], marker="^", c='black', label="Experimental - shifted")
plt.arrow(F0_shifted[1,1]+0.1, F0_shifted[1,3]-0.001,0.3,-0.003, width=0.0000005)
plt.text(F0_shifted[1,1]+0.4, F0_shifted[1,3]-0.005,"Experimental - shifted",fontsize = 18)
for i in range(0,nmodes):
    plt.plot(F0_nas[i][:,1],F0_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.savefig('F0_comparison.pdf')

plt.figure("f1compa",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
plt.xlim(0,max(F1_nas[0][:,1]))
plt.scatter(F1_exp[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_exp[:,2]/F1_exp[0,2], marker="s", c='black', label="Experimental")
plt.arrow(2.5, 0.97,0.19,0.25)
plt.text(F1_exp[6,1]+0.15, 1.21,"Experimental",fontsize = 18)
plt.scatter([None],[None], marker = "P", c="black", label = "Experimental - shifted")
fviv = 8/(0.012*F1_exp[0,2])
plt.plot([0,fviv,2*fviv,3*fviv,4*fviv],[0,1,2,3,4], c='black', linestyle = ':', linewidth = 1.0)
plt.text(0.8, 0.5, "Vortex shedding frequency",fontsize = 18)
plt.plot(F1_nas[0][:,1],F1_nas[0][:,2]/F1_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(2, 0.83, "Hydroelastic mode 1",fontsize = 18)
plt.ylim([0.4,1.5])
plt.xticks([])


plt.subplot(2, 1, 2)
plt.xlim(0,max(F1_nas[0][:,1]))
plt.scatter(F1_shifted[:,1]*1/(F1_exp[0,2]/F1_exp[:,2]), F1_shifted[:,3], marker="P", c='black', label="Experimental - shifted")
plt.arrow(F1_shifted[1,1]+0.1, F1_shifted[1,3]-0.001,0.3,-0.003, width=0.0000005)
plt.text(F1_shifted[1,1]+0.4, F1_shifted[1,3]-0.005,"Experimental - shifted",fontsize = 18)
for i in range(0,nmodes):
    plt.plot(F1_nas[i][:,1],F1_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.savefig('F1_comparison.pdf')


plt.figure("donacompa",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
plt.xlim(0,max(dona_nas[0][:,1]))
plt.scatter(dona_simu[:,1], dona_simu[:,2]/dona_simu[0,2], marker="o", c='black', label="Numerical")
plt.arrow(5.5, 1.05,2.7,0.18)
plt.text(dona_simu[4,1]+0.15, 1.21,"Numerical",fontsize = 18)
plt.scatter([None],[None], marker = "P", c="black", label = "Experimental - shifted")
plt.plot(dona_nas[0][:,1],dona_nas[0][:,2]/dona_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(2, 0.8, "Hydroelastic mode 1",fontsize = 18)
plt.ylim([0.4,1.5])
plt.ylabel(r'$\Omega_i$',fontsize = 18)
plt.xticks([])

ax = plt.subplot(2, 1, 2)
plt.xlim(0,max(dona_nas[0][:,1]))
plt.scatter(dona_shifted[:,1], dona_shifted[:,3], marker="^", c='black', label="Numerical - shifted")
plt.arrow(dona_shifted[2,1]+0.1, dona_shifted[2,3]-0.001,2.5,-0.01, width=0.0000005)
plt.text(dona_shifted[2,1]+2.8, dona_shifted[2,3]-0.015,"Numerical - shifted",fontsize = 18)
for i in range(0,nmodes):
    plt.plot(dona_nas[i][:,1],dona_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)
plt.ylabel(r'$\zeta_{i,added}$',fontsize = 18)
plt.savefig('dona_comparison.pdf')

plt.figure("bluntcompa",figsize=(6,9.5))
ax = plt.subplot(2, 1, 1)
plt.xlim(0,max(blunt_nas[0][:,1]))

plt.scatter(blunt_exp[:,1], blunt_exp[:,2]/blunt_exp[0,2], marker="s", c='black', label="Experimental")
plt.scatter(blunt_simu[:,1], blunt_simu[:,2]/blunt_simu[0,2], marker="x", c='black', label="Numerical")
plt.arrow(0.1, 1.03,0.8,0.15)
plt.text(blunt_exp[0,1]+1, 1.21,"Experimental",fontsize = 18)
plt.arrow(5.5, 1,2.7,0.25)
plt.text(blunt_simu[4,1]+0.15, 1.21,"Numerical",fontsize = 18)
plt.scatter([None],[None], marker = "P", c="black", label = "Experimental - shifted")
plt.plot(blunt_nas[0][:,1],blunt_nas[0][:,2]/blunt_nas[0][0,2], linestyle = lines[0], linewidth = widths[0], c='black', label = "Hydroelastic mode 1")
plt.text(2, 0.8, "Hydroelastic mode 1",fontsize = 18)
plt.ylim([0.4,1.5])
plt.xticks([])



plt.subplot(2, 1,2)
plt.xlim(0,max(blunt_nas[0][:,1]))
plt.ylim(-0.015,0.07)
plt.scatter(blunt_shifted[:,1], blunt_shifted[:,3], marker="P", c='black', label="Experimental - shifted")
plt.scatter(blunt_simu_shifted[:,1], blunt_simu_shifted[:,3], marker="v", c='black', label="Experimental - shifted")
plt.arrow(blunt_shifted[2,1]+0.1, blunt_shifted[2,3]-0.001,1.5,-0.01, width=0.0000005)
plt.text(blunt_shifted[2,1]+1.7, blunt_shifted[2,3]-0.015,"Experimental - shifted",fontsize = 18)
plt.arrow(blunt_simu_shifted[3,1]+0.1, blunt_simu_shifted[3,3]-0.001,0.7,-0.02, width=0.0000005)
plt.text(blunt_simu_shifted[3,1]+0.9, blunt_simu_shifted[3,3]-0.025,"Numerical - shifted",fontsize = 18)
for i in range(0,nmodes):
    plt.plot(blunt_nas[i][:,1],blunt_nas[i][:,3], linestyle = lines[i], linewidth = widths[i], c='black', label = "Hydroelastic mode "+str(i+1))
plt.xlabel(r'$U_R$',fontsize = 18)

plt.savefig('blunt_comparaison.pdf')


plt.show()