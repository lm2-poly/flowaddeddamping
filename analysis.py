#Functions to analyze the different outputs
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-12

#Obsolete

from running import *
import numpy as np
import matplotlib.pyplot as plt

def analyse(filename, analysis_type):
    file_in, file_out, file_out_txt, file_verif, file_run, filename = naming(filename, analysis_type)
    
    if analysis_type == "modes" or analysis_type == "acoustic" or analysis_type == "static":
        print("TODO: modes analysis")
        exit()

    elif analysis_type == "aeroelastic":
        f = open(file_out_txt,'r')
        rows = f.readlines()
        f.close()
        
        row2 = []
        for i in range(len(rows)):
            row = rows[i].split(' ')
            row = [item for item in row if item != '']
            row2.append(row)
            
        flutter = False
        for i in range(len(row2)):
            if len(row2[i])>1 and flutter == False:
                if row2[i][1] == 'FLUTTER':
                    flutter = True
                    i_start = i
                
            if flutter == True:
                if row2[i][0] == 'END-OF-DATA':
                    i_end = i
        rows = row2[i_start:i_end]
        
        index = []
        for i in range(len(rows)):
            if rows[i][0] == 'KFREQ' or rows[i] == " " or rows[i][0] == 'CONFIGURATION' or rows[i][0] == '0' or\
                    rows[i][0] == '1' or rows[i][0] == '\n':
                index.append(i)
        for i in range(len(index)):
            ind = index[i]
            del rows[ind]
            for j in range(i,len(index)):
                index[j] = index[j] - 1
        modes = []
        machs = []
        densities = []
        velocities = []
        for i in range(len(rows)):
            if rows[i][0] == "POINT":
                if len(modes) == 0:
                    modes.append(int(rows[i][2]))
                    machs.append(float(rows[i][6]))
                    densities.append(float(rows[i][10]))
                else:
                    if int(rows[i][2]) != modes[-1]:
                        modes.append(int(rows[i][2]))
                    if float(rows[i][6]) != machs[-1]:
                        machs.append(float(rows[i][6]))
                    if float(rows[i][10]) != densities[-1]:
                        densities.append(float(rows[i][10]))

        modes = np.asarray(modes)
        machs = np.asarray(machs)
        densities = np.asarray(densities)


        met = False
        for i in range(len(rows)):
            if rows[i][0] == "POINT":
                if int(rows[i][2]) == modes[-1] and float(rows[i][6]) == machs[-1] and\
                        float(rows[i][10]) == densities[-1] and met == False:
                    met = True
            else:
                if met == True and is_float(rows[i][0]) == False:
                    last_index = i
                    break
        rows = rows[0:last_index]

        for i in range(len(rows)):
            if rows[i][0] != "POINT":
                if len(velocities) == 0:
                    velocities.append(float(rows[i][2]))
                else:
                    if float(rows[i][2]) == velocities[0]:
                        break
                    if float(rows[i][2]) != velocities[-1]:
                        velocities.append(float(rows[i][2]))
        velocities = np.asarray(velocities)

        states = np.zeros([len(modes), len(machs), len(densities), len(velocities),7])
        for i in range(len(rows)):
            if rows[i][0] == 'POINT':
                mode = int(rows[i][2])
                mach = float(rows[i][6])
                density = float(rows[i][10])
                mode_i = np.where(modes == mode)
                mach_i = np.where(machs == mach)
                density_i = np.where(densities == density)
            if is_float(rows[i][0]) == True:
                velocity = float(rows[i][2])
                velocity_i = np.where(velocities == velocity)
                for j in range(7):
                    #See F06 file to understand what each column represent
                    states[mode_i, mach_i, density_i, velocity_i, j] = rows[i][j]

        wn = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
        wd = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
        zeta = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
        for i in range(len(modes)):
            for j in range(len(machs)):
                for k in range(len(densities)):
                    for m in range(len(velocities)):
                        wn[i,j,k,m] = states[i,j,k,m,4]*2*np.pi
                        wd[i,j,k,m] = states[i,j,k,m,6]
                        zetawn = states[i,j,k,m,5]
                        zeta[i,j,k,m] = zetawn/wn[i,j,k,m]

        plt.figure('Damping')
        plt.xlabel('Velocity [m/s]')
        plt.ylabel(r'$\zeta$')
        plt.title('Adimensional damping according to velocity')
        plt.grid(True)
        for i in range(len(modes)):
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.plot(velocities/1000, zeta[i,j,k,:], label='mode = '+str(modes[i])+', mach = '+str(machs[j])+\
                                                                 ', density = '+str(densities[k]))

        plt.figure('Frequency')
        plt.xlabel('Velocity [m/s]')
        plt.ylabel(r'$f_n$')
        plt.title('Frequency according to velocity')
        plt.grid(True)
        for i in range(len(modes)):
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.plot(velocities, wn[i, j, k, :]/(2*np.pi),
                             label='mode = ' + str(modes[i]) + ', mach = ' + str(machs[j]) + \
                                   ', density = ' + str(densities[k]))

        plt.show()

        return states

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

if __name__ == "__main__":
    filename = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\NACA4412\naca4412'
    analysis_type = 'aeroelastic'
    states = analyse(filename, analysis_type)