#Plotting the important data for the article
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-09-05

import matplotlib.pyplot as plt
import pandas as pd

F1_exp = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\results\F1_Exp.csv'
F1_nas = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\results\F1_NAS.csv'

CFX_num = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\results\CFX_num.csv'
CFX_nas = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\results\CFX_NAS.csv'


F1_exp = pd.read_csv(F1_exp, sep=';', header = None).to_numpy()
F1_nas = pd.read_csv(F1_nas, sep=';', header = None).to_numpy()
f_F1 = 286.35
L_F1 = 250.0E-3
F1_exp[:,0] = F1_exp[:,0]*f_F1*L_F1
F1_exp[:,1] = F1_exp[:,1] - F1_exp[0,1]
F1_nas[:,0] = F1_nas[:,0]*f_F1*L_F1

CFX_num = pd.read_csv(CFX_num, sep=';', header = None).to_numpy()
CFX_nas = pd.read_csv(CFX_nas, sep=';', header = None).to_numpy()
f_CFX = 173.71
L_CFX = 95E-3
#CFX_num[:,0] = CFX_num[:,0]/(f_CFX*L_CFX)
#CFX_num[:,1] = CFX_num[:,1] - CFX_num[0,1]
#CFX_nas[:,0] = CFX_nas[:,0]/(f_CFX*L_CFX)

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