#Functions to analyze the different outputs
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-12


from running import *
import numpy as np
import matplotlib.pyplot as plt

def aero_analysis(filename, file_out_txt, ref_length):
    f = open(file_out_txt, 'r')
    rows = f.readlines()
    f.close()

    row2 = []
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    flutter = False
    fn = [0, 0, 0, 0, 0]
    for i in range(len(row2)):
        if len(row2[i]) > 1 and flutter == False:
            if row2[i][1] == 'FLUTTER' and row2[i][0] != "PK" and row2[i][2] != "NOLIST":
                flutter = True
                i_start = i
            if row2[i][0] == 'MODE' and row2[i][2] == 'EIGENVALUE':
                fn[0] = np.sqrt(float(row2[i + 2][2])) / (2 * np.pi)
                fn[1] = np.sqrt(float(row2[i + 3][2])) / (2 * np.pi)
                fn[2] = np.sqrt(float(row2[i + 4][2])) / (2 * np.pi)
                fn[3] = np.sqrt(float(row2[i + 5][2])) / (2 * np.pi)
                fn[4] = np.sqrt(float(row2[i + 6][2])) / (2 * np.pi)

        if flutter == True:
            if row2[i][0] == 'END-OF-DATA':
                i_end = i
    rows = row2[i_start:i_end]

    index = []
    for i in range(len(rows)):
        if rows[i][0] == 'KFREQ' or rows[i] == " " or rows[i][0] == 'CONFIGURATION' or rows[i][0] == '0' or \
                rows[i][0] == '1' or rows[i][0] == '\n':
            index.append(i)
    for i in range(len(index)):
        ind = index[i]
        del rows[ind]
        for j in range(i, len(index)):
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
            if int(rows[i][2]) == modes[-1] and float(rows[i][6]) == machs[-1] and \
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

    states = np.zeros([len(modes), len(machs), len(densities), len(velocities), 7])
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
                # See F06 file to understand what each column represent
                states[mode_i, mach_i, density_i, velocity_i, j] = rows[i][j]

    wn = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
    wd = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
    zeta = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
    vstar = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
    Re = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
    Im = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
    from scipy.optimize import fsolve
    import cmath
    for i in range(len(modes)):
        for j in range(len(machs)):
            for k in range(len(densities)):
                for m in range(len(velocities)):
                    # Re[i,j,k,m] = states[i,j,k,m,6]
                    # Im[i,j,k,m] = states[i,j,k,m,5]
                    wd[i, j, k, m] = states[i, j, k, m, 6]
                    # wn[i, k, k, m] = states[i, j, k, m, 4] * 2 * np.pi
                    wn[i, k, k, m] = states[i, j, k, m, 4] * 2 * np.pi
                    zeta[i, j, k, m] = -states[i, j, k, m, 3] / 2
                    # elif int(states[i,j,k,m,6])==0 and int(states[i,j,k,m,4])==0:
                    # wn[i, k, k, m] = -states[i,j,k,m,5]/(zeta[i,j,k,m]+np.sqrt((zeta[i,j,k,m])**2-1))
                    # zetawn = -states[i, j, k, m, 5]
                    # root = fsolve(lambda x: eigenvals(x,wd[i,j,k,m],zetawn), [wd[i,j,k,m], zetawn/wd[i,j,k,m]])
                    # wn[i,j,k,m] = root[0]
                    # zeta[i,j,k,m] = root[1]
                    coeff = [1, 2 * zeta[i, j, k, m] * wn[i, j, k, m], wn[i, j, k, m] ** 2]
                    eigenvalue = np.roots(coeff)
                    if len(eigenvalue) == 2:
                        eigenvalue = eigenvalue[0]
                    Re[i, j, k, m] = eigenvalue.real
                    Im[i, j, k, m] = eigenvalue.imag
                    with np.errstate(divide='ignore'):
                        warnings.simplefilter("ignore")
                        # zeta[i, j, k, m] = zetawn / wn[i, j, k, m]
                        vstar[i, j, k, m] = velocities[m] / (ref_length * wn[i, j, k, m] / (2 * np.pi))
    data1 = np.array(
        [velocities, vstar[0, 0, 0, :], wn[0, 0, 0, :] / (2 * np.pi), zeta[0, 0, 0, :], Re[0, 0, 0, :],
         Im[0, 0, 0, :]]).T
    data2 = np.array(
        [velocities, vstar[1, 0, 0, :], wn[1, 0, 0, :] / (2 * np.pi), zeta[1, 0, 0, :], Re[1, 0, 0, :],
         Im[1, 0, 0, :]]).T
    data3 = np.array(
        [velocities, vstar[2, 0, 0, :], wn[2, 0, 0, :] / (2 * np.pi), zeta[2, 0, 0, :], Re[2, 0, 0, :],
         Im[2, 0, 0, :]]).T
    data = [data1, data2, data3]
    for i in range(3):
        file_num = list(filename)
        file_num.append(str(i) + ".csv")
        file_num = "".join(file_num)
        np.savetxt(file_num, data[i], delimiter=",")

    # Velocity
    plt.figure('Frequency and Damping - Velocity - Constant', figsize=(12, 6))
    plt.suptitle("Frequency and damping according to velocity")
    # subplot 1
    plt.subplot(1, 2, 2)
    x = [0, 5, 10, 15, 20, 25, 28]
    y = [0.02, 0.06, 0.095, 0.14, 0.18, 0.225, 0.27] - 0.02 * np.ones(7)
    # plt.scatter(x,y, color = 'black', marker = 'o', label = 'CFX')
    # exp_data(self.results_filename)
    plt.ylabel(r'$\zeta$')
    plt.title("Adimensional damping according to velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            if True:
                for j in range(len(machs)):
                    for k in range(len(densities)):
                        plt.xlabel(r'$U$ [m/s]')
                        plt.plot(velocities / 1000, zeta[i, j, k, :],
                                 label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                        # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
    plt.legend(loc='best')
    plt.tight_layout()

    # subplot 2
    plt.subplot(1, 2, 1)
    plt.ylabel(r'$f_n$')
    plt.title("Frequency according to velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$U$ [m/s]')
                    plt.plot(velocities / 1000, wn[i, j, k, :] / (2 * np.pi),
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Frequency and Damping - Velocity - Constant.png')

    plt.figure('Eigenvalues - Velocity - Constant', figsize=(12, 6))
    plt.suptitle("Eigenvalues according to velocity")
    # Subplot 1
    plt.subplot(2, 1, 1)
    plt.ylabel(r'$Re(\omega)$')
    plt.title("Real part of the eigenvalue according to velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$U$ [m/s]')
                    plt.plot(velocities / 1000, Re[i, j, k, :],
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()

    # subplot 2
    plt.subplot(2, 1, 2)
    plt.ylabel(r'$Im(\omega)$')
    plt.title("Imaginary part of the eigenvalue according to velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$U$ [m/s]')
                    plt.plot(velocities / 1000, Im[i, j, k, :],
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('Eigenvalues - Velocity - Constant.png')

    # Reduced velocity
    plt.figure('Frequency and Damping - Reduced Velocity - Constant', figsize=(12, 6))
    plt.suptitle("Frequency and damping according to reduced velocity")
    # subplot 1
    plt.subplot(1, 2, 2)
    # exp_data(self.results_filename)
    plt.ylabel(r'Flow added damping $\zeta_f$')
    plt.title("Adimensional added damping according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            if True:
                for j in range(len(machs)):
                    for k in range(len(densities)):
                        plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                        plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)),
                                 zeta[i, j, k, :],
                                 label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                        # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Added Damping.pdf')

    # subplot 2
    plt.subplot(1, 2, 1)
    plt.ylabel(r'$f_n$')
    plt.title("Frequency according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                    plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)),
                             wn[i, j, k, :] / (2 * np.pi),
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Frequency and Damping - Reduced Velocity - Constant.png')

    plt.figure('Eigenvalues - Reduced Velocity - Constant', figsize=(12, 6))
    plt.suptitle("Eigenvalues according to reduced velocity")
    # Subplot 1
    plt.subplot(2, 1, 1)
    plt.ylabel(r'$Re(\omega)$')
    plt.title("Real part of the eigenvalue according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                    plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)), Re[i, j, k, :],
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()

    # subplot 2
    plt.subplot(2, 1, 2)
    plt.ylabel(r'$Im(\omega)$')
    plt.title("Imaginary part of the eigenvalue according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                    plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)), Im[i, j, k, :],
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Eigenvalues - Reduced Velocity - Constant.png')

    # Variable frequencies
    # Reduced velocity
    plt.figure('Frequency and Damping - Reduced Velocity - Variable', figsize=(12, 6))
    plt.suptitle("Frequency and damping according to reduced velocity")
    # subplot 1
    plt.subplot(1, 2, 2)
    # exp_data(self.results_filename)
    plt.ylabel(r'$\zeta$')
    plt.title("Adimensional damping according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            if True:
                for j in range(len(machs)):
                    for k in range(len(densities)):
                        plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                        plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)),
                                 zeta[i, j, k, :],
                                 label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                        # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Added Damping.pdf')

    # subplot 2
    plt.subplot(1, 2, 1)
    plt.ylabel(r'$f_n$')
    plt.title("Frequency according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                    plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)),
                             wn[i, j, k, :] / (2 * np.pi),
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Frequency and Damping - Reduced Velocity - Variable.png')

    plt.figure('Eigenvalues - Reduced Velocity - Variable', figsize=(12, 6))
    plt.suptitle("Eigenvalues according to reduced velocity")
    # Subplot 1
    plt.subplot(2, 1, 1)
    plt.ylabel(r'$Re(\omega)$')
    plt.title("Real part of the eigenvalue according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                    plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)),
                             Re[i, j, k, :],
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()

    # subplot 2
    plt.subplot(2, 1, 2)
    plt.ylabel(r'$Im(\omega)$')
    plt.title("Imaginary part of the eigenvalue according to reduced velocity")
    plt.grid(True)
    for i in range(min(len(modes), 5)):
        if True:
            for j in range(len(machs)):
                for k in range(len(densities)):
                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                    plt.plot(velocities / (wn[i, j, k, :] * ref_length / (2 * np.pi)),
                             Im[i, j, k, :],
                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Eigenvalues - Reduced Velocity - Variable.png')

    plt.figure("Argand", figsize=(8, 6))
    plt.xlabel(r'$Re(\omega)$')
    plt.ylabel(r'$Im(\omega)$')
    plt.title("Argand's Diagram")
    plt.grid(True)
    plt.scatter(Re[0, 0, 0, 0], Im[0, 0, 0, 0], marker='o', color='black', label="Start")
    plt.scatter(Re[0, 0, 0, -1], Im[0, 0, 0, -1], marker='X', color='black', label="End")
    for i in range(min(len(modes), 5)):
        if True:
            if True:
                for j in range(len(machs)):
                    for k in range(len(densities)):
                        plt.scatter(Re[i, j, k, 0], Im[i, j, k, 0], marker='o', color='black')
                        plt.scatter(Re[i, j, k, -1], Im[i, j, k, -1], marker='X', color='black')
                        plt.plot(Re[i, j, k, :], Im[i, j, k, :],
                                 label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                        # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
    plt.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('Argand.png')

    plt.show(block=True)

    return states

if __name__ == "__main__":
    filename = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\flowaddeddamping\TestFiles\Bergan\Bergan'
    file_out_txt = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\flowaddeddamping\TestFiles\Bergan\Bergan_hydroelastic\Bergan_hydroelastic_aeroelastic\Bergan_hydroelastic_aeroelastic.f06'
    states = aero_analysis(filename = filename, file_out_txt = file_out_txt, ref_length = 250)