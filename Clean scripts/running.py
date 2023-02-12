#Functions to run the complete code
#Author: Danick Lamoureux
#Project under Frédérick Gosselin and Sébastien Houde's supervision
#Date: 2022-05-11

import os
import time
import subprocess
from pyNastran.op2.op2 import read_op2
from nodes import *
from loads import *
from elements import *
from headings import *
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

def naming(filename, analysis_type = None, create_folder = False):
    """Naming the different files generated

    Args:
        filename (string): General filename and location
        analysis_type (string, optional): Analysis type to perform. Defaults to None.
        create_folder (bool, optional): Defines if a new folder is created. Defaults to False.

    Returns:
        strings: all the named files
    """
    if analysis_type != None:
        #Extracting the non-extension part of the filename and adding the analysis type to it
        filename = list(filename)

        for i in range(len(filename)):
            j = -1-i
            if filename[j] == "\\" or filename[j] == "/":
                path = filename[:j]
                folder_name = filename[j:]
                break

        analysis_type = list(analysis_type)
        folder_name.append("_")
        for letter in analysis_type:
            folder_name.append(letter)
        for letter in folder_name:
            path.append(letter)
        folder_path = "".join(path)
        Path(folder_path).mkdir(parents = True, exist_ok = True)
        for letter in folder_name:
            path.append(letter)
        filename = "".join(path)
    filename = os.path.realpath(filename)
    #Naming the input file with bdf extension
    file_in = list(filename)
    file_in.append(".bdf")
    file_in = "".join(file_in)

    #Naming the output file to find it with op2 extension
    file_out = list(filename)
    file_out.append(".op2")
    file_out = "".join(file_out)

    #Naming the output file to find it with f06 extension to verify the analysis is done
    file_out_txt = list(filename)
    file_out_txt.append(".f06")
    file_out_txt = "".join(file_out_txt)

    #Naming a file that should be deleted once done to find it
    file_verif = list(filename)
    file_verif.append(".plt")
    file_verif = "".join(file_verif)

    #Naming the bat file to run
    file_run = list(filename)
    file_run.append(".bat")
    file_run = "".join(file_run)

    #Numerical results file
    file_num = list(filename)
    file_num.append(".csv")
    file_num = "".join(file_num)

    return file_in, file_out, file_out_txt, file_verif, file_run, file_num, filename

def run(filename, nastran_location = r'"C:\Program Files\Siemens\Simcenter3D_2020.2\NXNASTRAN\bin\nastranw.exe"'):
    """Running the NASTRAN simulations automatically. Based on Olivier Duchesne's code for launching NASTRAN's simulations automatically    

    Args:
        filename (string): General filename and location
        nastran_location (regexp, optional): NASTRAN's location. Defaults to r'"C:\\Program Files\\Siemens\Simcenter3D_2020.2\\NXNASTRAN\\bin\\nastranw.exe"'.
    """
    #Based on Olivier Duchesne's code for launching NASTRAN's simulations automatically    
    
    #Naming the files
    file_in, file_out, file_out_txt, file_verif, file_run, file_num, filename = naming(filename)

    # Creating the bat file
    bat = open(file_run, 'w')
    bat.write(nastran_location+' ')
    bat.write('"' + file_in + '" scr=yes old=no delete=f04,log,xdb')
    bat.close()

    # Changing working directory
    path = os.path.dirname(os.path.realpath(file_run))
    os.chdir(path)
    
    # Running the bat file
    status = subprocess.call(file_run)

    # Waiting for the solution to be done
    Ready_to_continue = False
    while not Ready_to_continue:
        time.sleep(5)
        if os.path.exists(file_out_txt) and not os.path.exists(file_verif):
            Ready_to_continue = True

def read(filename):
    """Reading the OP2 analysis file

    Args:
        filename (string): General filename and location

    Returns:
        OP2Model: OP2 model
    """
    file_in, file_out, file_out_txt, file_verif, file_run, file_num, filename = naming(filename)
    results = read_op2(file_out, debug = False)
    return results

#Class of the different analyses:
class model:
    """Model class for analysis
    """
    def __init__(self, filename, profile, analysis_type):
        #Defining the filenames and other important parameters
        self.filename = filename
        self.profile = profile
        #Making sure the analysis type exists
        if analysis_type == "static" or analysis_type == "modes" or analysis_type == "aeroelastic" or\
                analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics" or\
                analysis_type == 'hydroelastic':
            self.analysis_type = analysis_type
        else:
            print("Analysis type unknown")
            exit()
        file_in, file_out, file_out_txt, file_verif, file_run, file_num, filename = naming(self.filename, self.analysis_type)
        self.filename = filename
        self.file_in = file_in
        self.file_out_txt = file_out_txt
        self.file_out = file_out
        self.file_num = file_num
        self.addedmass = False
        self.fn = np.zeros(1000)


    def setup(self, geometry_object, mesh_object, solid_object, fluid_object, flow_object = None, show = False):
        """Saving the various objects used for the analysis

        Args:
            geometry_object: Geometry
            mesh_object: Mesh parameters
            solid_object: Solid
            fluid_object: Acoustic fluid
            flow_object: Aeroelastic flow. Defaults to None.
        """
        if self.analysis_type == "hydroelastic":
            # For an hydroelastic analysis, define the different objects
            self.geom = geometry_object
            self.mesh = mesh_object
            self.solid = solid_object
            self.fluid = fluid_object
            self.flow = flow_object
            self.hydroelastic(show)
        else:
            # Nodes generation (GRID)
            nodelist = nodes()

            # Loading the profile
            nodelist.profile_load(self.profile, n = geometry_object.n_hydrofoils, spacing = geometry_object.hydrofoil_spacing)
            
            #To add cascades (work in progress)
            #nodelist.cascade(n = geometry_object.n_hydrofoils, spacing = geometry_object.hydrofoils_spacing)
            
            nodelist.envelope_load(envelope_chord = geometry_object.envelope_chord, thick=geometry_object.thick,
                                   n = geometry_object.n_hydrofoils, spacing = geometry_object.hydrofoil_spacing)

            # Meshing the profile using gmsh
            # To show the mesh, show = True to use gmsh software
            nodelist.profile_mesh(mesh_size_normalized=mesh_object.solid_mesh_size, n=geometry_object.n_hydrofoils, show=show)
            if self.analysis_type == "uncoupled_acoustics" or self.analysis_type == "coupled_acoustics" or self.addedmass:
                nodelist.envelope_mesh(mesh_size_normalized=mesh_object.fluid_mesh_size, n = geometry_object.n_hydrofoils, show=show)

            # Meshing the wing in 3D
            nodelist.wing_def(rootchord=geometry_object.rootchord,
                              tipchord=geometry_object.tipchord,
                              span=geometry_object.span,
                              roottwist=geometry_object.roottwist,
                              tiptwist=geometry_object.tiptwist,
                              sweep=geometry_object.sweep,
                              dihedral=geometry_object.dihedral,
                              nspan=mesh_object.nspan,
                              nx=mesh_object.nx,
                              ny=mesh_object.ny,
                              n = geometry_object.n_hydrofoils,
                              spacing = geometry_object.hydrofoil_spacing)

            if self.analysis_type == "uncoupled_acoustics" or self.analysis_type == "coupled_acoustics" or self.addedmass:
                #Generating the fluid nodes if acoustics
                nodelist.flow_def()
            self.nodelist = nodelist

            self.E = solid_object.E
            self.nu = solid_object.nu
            self.rho_solid = solid_object.rho

            self.rho_flow = fluid_object.rho
            self.bulk = fluid_object.bulk

            self.n = geometry_object.n_hydrofoils
            self.spacing = geometry_object.hydrofoil_spacing
            self.chord = geometry_object.rootchord
            self.thickness = geometry_object.thickness

            if self.analysis_type == "aeroelastic":
                self.ref_velocity = flow_object.ref_velocity
                self.ref_length = flow_object.ref_length
                self.rho_flow = flow_object.rho_flow
                self.velocities = flow_object.velocities
                self.density_ratios = flow_object.density_ratios
                self.machs = flow_object.machs
                self.mach_matrix = flow_object.mach_matrix
                self.freq_matrix = flow_object.freq_matrix

    def write(self):
        """Writes the analysis file to .bdf
        """
        analysis_type = self.analysis_type
        if analysis_type == "hydroelastic":
            print("Writing the modal analysis input file")
            self.modes.write()
            print("Writing the coupled acoustics analysis input file")
            self.coupled.write()

        else:
            nodelist = self.nodelist
            file_in = self.file_in
            f = open(file_in, 'w')
            # Defining the headers from the analysis type
            print("Writing headers")
            if self.addedmass == False:
                if analysis_type == "static":
                    headers_bending(f)
                elif analysis_type == "modes":
                    headers_modes(f)
                elif analysis_type == "uncoupled_acoustics":
                    headers_acoustics_real(f)
                elif analysis_type == "coupled_acoustics":
                    headers_acoustics_complex(f)
                elif analysis_type == "aeroelastic":
                    headers_aero(f)
            elif self.addedmass == True and analysis_type == "aeroelastic":
                headers_hydro(f, P1 = self.P1, nmodes = self.nmodes)

            # Defining the materials and element properties
            print("Writing material and element properties")
            f.write('$* Materials and properties\n$*\n')
            MAT1(f, MID=1, E=self.E, nu=self.nu, rho=self.rho_solid)
            PSOLID(f, PID=1, MID=1)
            if analysis_type == "aeroelastic":
                PAERO(f, PID=2)
            f.write('$*\n')

            if analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics":
                MAT10(f, MID = 2, BULK = self.bulk, RHO = self.rho_flow)
                PFLUID(f, PID=2, MID = 2)

            # Defining the elements
            print("Writing elements")
            f.write('$* Elements\n')
            self.used_nodes = CPENTA(f, PID=1, nodes_object=nodelist)
            if analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics":
                CPENTAFluid(f, PID = 2, nodes_object = nodelist)
            if self.addedmass == True:
                f.write('$* Materials and properties for simulated nodes\n$*\n')
                MAT1(f, MID=3, E=1.0, nu=0.3, rho=1.0)
                PSOLID(f, PID=3, MID=3)

            # Writing the points using the GRID cards
            print("Writing the nodes")
            nodelist.GRIDWrite(f)
            if analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics":
                nodelist.GRIDWriteFluid(f)
            fluid_nodes_len = 0

            if analysis_type == "aeroelastic":
                print("Writing the aeroelastic model")
                CAERO(f, PID=2, nodes_object=nodelist, n = self.n, spacing = self.spacing*self.chord, IGID=1, fluid_nodes_len = fluid_nodes_len)
                SPLINE(file=f,
                       nodes_object=nodelist,
                       fluid_nodes_len = fluid_nodes_len,
                       EID=2,
                       n=self.n)

            if analysis_type == "aeroelastic":
                # Aero model
                AERO(file=f,
                     ref_length=self.ref_length,
                     rhoref=self.rho_flow)
                MKAERO(file=f,
                       SID=200,
                       mach_matrix=self.mach_matrix,
                       freq_matrix=self.freq_matrix,
                       velocities=self.velocities,
                       densities=self.density_ratios,
                       machs=self.machs)

            # Defining the constraints and loads:
            print("Writing the constraints and loads")
            SPC(f, nodelist=nodelist.nodelist, profile_nodes_len = len(nodelist.mesh_profile))
            if analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics": ACMODL(f)
            if analysis_type == "static": GRAV(f)

            # Footer and closing
            footers(f)

    def run(self, nmodes = 10):
        if self.analysis_type == "hydroelastic":
            import time
            time1 = time.time()
            timea = time.time()
            print("Vacuum Modal analysis")
            self.modes.run()
            time2 = time.time()
            vac_an = time2 - time1
            print("Required time for the vacuum modal analysis: "+str(vac_an)+"s")

            # Extracting the calculated modes without added mass
            op2 = read(self.modes.filename)
            eigenvectors = op2.eigenvectors[1].data[:, :, :].real
            eigenvalues = op2.eigenvalues[''].eigenvalues

            LAMBDAvac = np.zeros([len(eigenvalues[:])], dtype='float32')
            for i in range(len(eigenvalues)):
                LAMBDAvac[i] = eigenvalues[i]

            phi = np.zeros([len(eigenvectors[0, :, 0]) * 6, len(eigenvectors[:, 0, 0])], dtype='float32')

            for i in range(len(eigenvectors[0, :, 0])):
                DOF1 = 6 * i
                modes1 = eigenvectors[:, :, 0].T
                DOF2 = DOF1 + 1
                modes2 = eigenvectors[:, :, 1].T
                DOF3 = DOF1 + 2
                modes3 = eigenvectors[:, :, 2].T
                DOF4 = DOF1 + 3
                modes4 = eigenvectors[:, :, 3].T
                DOF5 = DOF1 + 4
                modes5 = eigenvectors[:, :, 4].T
                DOF6 = DOF1 + 5
                modes6 = eigenvectors[:, :, 5].T
                phi[DOF1, :] = modes1[i, :]
                phi[DOF2, :] = modes2[i, :]
                phi[DOF3, :] = modes3[i, :]
                phi[DOF4, :] = modes4[i, :]
                phi[DOF5, :] = modes5[i, :]
                phi[DOF6, :] = modes6[i, :]
            phi_s = phi
            
            time1 = time.time()
            print("Coupled Modal analysis")
            self.coupled.run()
            time2 = time.time()
            fr_an = time2 - time1
            print("Required time for the coupled modal analysis: " + str(fr_an) + "s")


            # Extracting the calculated modes with added mass
            op2 = read(self.coupled.filename)
            eigenvectors = op2.eigenvectors[1].data[:, :, :].real
            eigenvalues = op2.eigenvalues[''].eigenvalues

            LAMBDAcoup = np.zeros([len(eigenvalues)], dtype='float32')
            for i in range(len(eigenvalues)):
                LAMBDAcoup[i] = np.imag(eigenvalues[i]) ** 2
                

            
            phi = np.zeros([len(eigenvectors[0, :, 0]) * 6, len(eigenvectors[:, 0, 0])], dtype='float32')

            for i in range(len(eigenvectors[0, :, 0])):
                DOF1 = 6 * i
                modes1 = eigenvectors[:, :, 0].T
                DOF2 = DOF1 + 1
                modes2 = eigenvectors[:, :, 1].T
                DOF3 = DOF1 + 2
                modes3 = eigenvectors[:, :, 2].T
                DOF4 = DOF1 + 3
                modes4 = eigenvectors[:, :, 3].T
                DOF5 = DOF1 + 4
                modes5 = eigenvectors[:, :, 4].T
                DOF6 = DOF1 + 5
                modes6 = eigenvectors[:, :, 5].T
                phi[DOF1, :] = modes1[i, :]
                phi[DOF2, :] = modes2[i, :]
                phi[DOF3, :] = modes3[i, :]
                phi[DOF4, :] = modes4[i, :]
                phi[DOF5, :] = modes5[i, :]
                phi[DOF6, :] = modes6[i, :]

            phi_f = phi

            #Building the P1 matrix
            n = min(len(LAMBDAvac), len(LAMBDAcoup))
            j=0
            k=0
            for i in range(n):
                if LAMBDAvac[i] < (10*2*np.pi)**2:
                    j += 1
                elif LAMBDAcoup[i] < (10*2*np.pi)**2:
                    k += 1
            n = min(len(LAMBDAvac[j:]), len(LAMBDAcoup[k:]))
            LAMBDAvac = LAMBDAvac[:n]
            LAMBDAcoup = LAMBDAcoup[:n]
            phi_s = phi_s[:,:n]
            phi_f = phi_f[:,:n]

            n = len(LAMBDAcoup)
            self.aero.nmodes = min(n,nmodes)

            #As long as cascades are not implemented, keep following section
            MAC = np.zeros([n,n])
            from numpy.linalg import norm
            for ii in range(n):
                for jj in range(n): #Modal assurance criterion as described in the article
                    MAC[ii,jj] = norm(phi_s[:,ii].T@phi_f[:,jj])**2/(norm(phi_s[:,ii])**2*norm(phi_f[:,jj])**2)

            print("MAC = ")
            print(str(MAC))
            
            mode_index = np.zeros(self.aero.nmodes)
            for i in range(self.aero.nmodes):
                mode_index[i] = np.argmax(MAC[i,:]) #Rearranging to get the best fit between modes
                
                print("MAC("+str(i+1)+','+str(int(mode_index[i]+1))+') = '+str(MAC[i,int(mode_index[i])]))

            #Calculating first matrix to change frequencies
            invLAMBDAcoupmat = np.zeros([self.aero.nmodes,self.aero.nmodes])
            LAMBDAvacmat = np.zeros([self.aero.nmodes,self.aero.nmodes])
            LAMBDA = np.zeros([self.aero.nmodes,self.aero.nmodes])
            self.aero.fn = np.zeros(self.aero.nmodes)
            for i in range(self.aero.nmodes):
                invLAMBDAcoupmat[i,i] = np.sqrt(1/LAMBDAcoup[int(mode_index[i])])
                LAMBDAvacmat[i,i] = np.sqrt(LAMBDAvac[i])
                LAMBDA[i,i] = LAMBDAvac[i]
                self.aero.fn[i] = np.sqrt(LAMBDAcoup[int(mode_index[i])])/(2*np.pi)
            #P matrix as described in the article
            P1 = invLAMBDAcoupmat@LAMBDAvacmat
            n = len(P1)
            
            plt.figure('Compared Frequencies')
            plt.xlabel(r'$modes$')
            plt.ylabel(r'$f_n$')
            plt.title("Frequency of the different modes in vacuum and in resting water")
            plt.grid(True)
            mode_number = np.arange(1,n+1,1)
            plt.plot(mode_number, np.sqrt(LAMBDAvac[:n]), label='Vacuum')
            plt.plot(mode_number, np.sqrt(LAMBDAcoup[:n]), label='Resting fluid')
            plt.legend(loc='best')
            
            self.aero.P1 = P1

            #Flag to indicate next aeroelastic simulation is going to be hydroelastic
            self.aero.addedmass = True

            timeb = time.time()
            print("Required time for the added mass calculations: " + str(timeb - timea - fr_an - vac_an) + "s")

            print("Writing the aeroelastic input file with added mass")
            self.aero.write()
            time1 = time.time()
            print("Running the aeroelastic input file with added mass")
            self.aero.run()
            time2 = time.time()
            print("Temps requis par l'analyse aéroélastique: " + str(time2 - time1) + "s")

        else:
            print("Running")
            run(self.filename)

    def analyse(self, show=True):
        """Depending on the analysis file, performs the analysis.

        Args:
            show (bool, optional): Showing the results in pyplots graphs. Defaults to True.

        Returns:
            variable: Results object
        """
        analysis_type = self.analysis_type
        file_out_txt = self.file_out_txt
        file_out = self.file_out

        if analysis_type == "modes" or analysis_type == "coupled_acoustics":
            coupling = False
            op2 = read(self.filename)
            mode_number = []

            wn = op2.eigenvalues[''].eigenvalues
            if analysis_type == "coupled_acoustics":
                coupling = True
                for i in range(len(wn)):
                    wn[i] = abs(wn[i].imag)**2
            fn = []
            for i in range(len(wn)):
                fn_i = np.sqrt(wn[i])/(2*np.pi)
                fn.append(fn_i)
                mode_number.append(op2.eigenvectors[1].modes[i])

            if show == True:
                plt.figure('Frequency - '+str(analysis_type))
                plt.xlabel(r'$modes$')
                plt.ylabel(r'$f_n$')
                if analysis_type == "modes":
                    plt.title("Frequency of the different modes")
                else:
                    plt.title("Frequency of the different modes under resting fluid-structure coupling")
                plt.grid(True)
                plt.plot(mode_number, fn)
                if coupling == False:
                    plt.savefig('Frequency in vacuum.pdf')
                else:
                    plt.savefig('Frequency in resting fluid.pdf')

                plt.show(block = True)
            return fn

        elif analysis_type == "uncoupled_acoustics":
            print("Please open Simcenter 3D to analyze the uncoupled acoustic modes")
            
        elif analysis_type == 'hydroelastic':
            self.aero.analyse(show)
        
        elif analysis_type == "aeroelastic":
            #pyNastran cannot read aeroelastic results from OP2 file -> custom f06 reader
            # The following lines aim to read the text file until the appropriate section is found
            # It is recommended to have an aeroelastic f06 file besides to follow along
            f = open(file_out_txt, 'r')
            rows = f.readlines()
            f.close()

            row2 = []
            for i in range(len(rows)):
                row = rows[i].split(' ')
                row = [item for item in row if item != '']
                row2.append(row)

            flutter = False
            for i in range(len(row2)):
                if len(row2[i]) > 1 and flutter == False:
                    if row2[i][1] == 'FLUTTER' and row2[i][0] != "PK" and row2[i][2] != "NOLIST":
                        flutter = True
                        i_start = i

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

            Omega = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            wn = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            fn = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            wd = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            zeta = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            UR = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            Re = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            Im = np.zeros([len(modes), len(machs), len(densities), len(velocities)])
            
            #Extracting the results from f06 file
            for i in range(len(modes)):
                for j in range(len(machs)):
                    for k in range(len(densities)):
                        for m in range(len(velocities)):
                            wd[i, j, k, m] = states[i, j, k, m, 6]
                            if self.fn[0] != 0:
                                Omega[i,j,k,m] = states[i,j,k,m,4]/self.fn[0]
                            else:
                                Omega[i,j,k,m] = states[i,j,k,m,4]/states[0,j,k,0,4]
                            fn[i,j,k,m] = states[i,j,k,m,4]
                            wn[i, k, k, m] = states[i, j, k, m, 4] * 2 * np.pi
                            zeta[i,j,k,m] = -states[i,j,k,m,3]/2
                            coeff = [1, 2*zeta[i,j,k,m]*wn[i,j,k,m], wn[i,j,k,m]**2]
                            eigenvalue = np.roots(coeff)
                            if len(eigenvalue) == 2:
                                eigenvalue = eigenvalue[0]
                            Re[i,j,k,m] = eigenvalue.real
                            Im[i,j,k,m] = eigenvalue.imag
                            with np.errstate(divide='ignore'):
                                warnings.simplefilter("ignore")
                                if self.fn[0] != 0:
                                    UR[i,j,k,m] = velocities[m]/(self.thickness*self.fn[0])
                                else:
                                    UR[i,j,k,m] = velocities[m]/(self.thickness*wn[0,j,k,0]/(2*np.pi))
            
            #Saving the data to a csv file
            for i in range(5):
                data = np.array(
                [velocities, UR[i, 0, 0, :], fn[i, 0, 0, :], zeta[i, 0, 0, :], Re[i, 0, 0, :],
                 Im[i, 0, 0, :]]).T
                file_num = list(self.filename)
                file_num.append(str(i)+".csv")
                file_num = "".join(file_num)
                np.savetxt(file_num, data, delimiter = ",")

            #Displaying results
            #To obtain similar graphs to what is presented in the article, see graphing.py
            if show == True:
                plt.figure('Frequency and Damping', figsize = (4,6))
                # subplot 1
                plt.subplot(2, 1, 2)
                plt.ylabel(r'$\zeta_{i,added}$')
                plt.grid(True)
                for i in range(min(len(modes), 3)):
                    for j in range(len(machs)):
                        for k in range(len(densities)):
                            plt.xlabel(r'$U_R$')
                            plt.plot(UR[i,j,k,:], zeta[i, j, k, :],
                                        label='Mode ' + str(modes[i]))
                plt.legend(loc='best')
                plt.tight_layout()
                

                # subplot 2
                plt.subplot(2, 1, 1)
                plt.ylabel(r'$\Omega_i$')
                plt.grid(True)
                for i in range(min(len(modes), 3)):
                    for j in range(len(machs)):
                        for k in range(len(densities)):
                            plt.plot(UR[i,j,k,:], Omega[i,j,k,:],
                                        label='Mode ' + str(modes[i]))
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('FrequencyDamping.pdf')

            return states

    def hydroelastic(self, show = False):
        """Defining the analyses required for the hydroelastic analysis
        """

        if self.analysis_type == "hydroelastic":
            self.coupled = model(self.filename, self.profile, 'coupled_acoustics')
            self.aero = model(self.filename, self.profile, 'aeroelastic')
            self.modes = model(self.filename, self.profile, 'modes')

        print("Setup the coupled acoustics analysis")
        self.coupled.setup(geometry_object=self.geom,
                           mesh_object=self.mesh,
                           solid_object=self.solid,
                           fluid_object=self.fluid,
                           flow_object=self.flow,
                           show = show)

        print("Setup the modal analysis")
        self.modes.setup(geometry_object=self.geom,
                        mesh_object=self.mesh,
                        solid_object=self.solid,
                        fluid_object=self.fluid,
                        flow_object=self.flow,
                        show = show)

        print("Setup the aeroelastic analysis with added mass")
        self.aero.setup(geometry_object=self.geom,
                         mesh_object=self.mesh,
                         solid_object=self.solid,
                         fluid_object=self.fluid,
                         flow_object=self.flow,
                         show = show)

class geometry:
    """Defining the hydrofoil geometry
    """
    def __init__(self, rootchord, tipchord, thickness, span, envelope_chord, thick, roottwist = 0.0, tiptwist = 0.0,
                 sweep = 0.0, dihedral = 0.0, n_hydrofoils = 1, spacing = 10):
        self.rootchord = rootchord
        self.tipchord = tipchord
        self.thickness = thickness
        self.span = span
        self.roottwist = roottwist
        self.tiptwist = tiptwist
        self.sweep = sweep
        self.dihedral = dihedral
        self.envelope_chord = envelope_chord
        self.thick = thick
        self.n_hydrofoils = n_hydrofoils
        self.hydrofoil_spacing = spacing

class mesh:
    """Defining the mesh parameters for both solid and fluid
    """
    
    def __init__(self, solid_mesh_size, fluid_mesh_size, nspan, nx, ny):
        self.solid_mesh_size = solid_mesh_size
        self.fluid_mesh_size = fluid_mesh_size
        self.nspan = nspan
        self.ny = ny
        self.nx = nx

class solid:
    """Defining the solid for dynamic analysis
    """
    def __init__(self, E, nu, rho):
        self.E = E
        self.nu = nu
        self.rho = rho

class fluid:
    """Defining the fluid for acoustics analysis
    """
    def __init__(self, rho, bulk):
        self.rho = rho
        self.bulk = bulk

class flow:
    """Defining the flow for flutter analysis
    """
    def __init__(self, ref_velocity, ref_length, rho_flow, velocities, density_ratios, machs, mach_matrix, freq_matrix):
        self.ref_velocity = ref_velocity
        self.ref_length = ref_length
        self.rho_flow = rho_flow
        self.velocities = velocities
        self.density_ratios = density_ratios
        self.machs = machs
        self.mach_matrix = mach_matrix
        self.freq_matrix = freq_matrix

def is_float(string):
    """Defining is a string contains a float

    Args:
        string (string): String potentially containing a float

    Returns:
        bool: Is the string a float
    """
    try:
        float(string)
        return True
    except ValueError:
        return False

import pandas as pd
import math as mt

def exp_data(filename):
    """Plots the experimental data on a pyplot graph

    Args:
        filename (string): Experimental data filename
    """
    if filename != None:
        #From an excel file, plot the damping according to reduced velocity
        data = pd.read_csv(filename, sep=';').to_numpy()
        test_data = np.zeros([int(len(data[0,:])/2), 2, len(data[:,0])])
        markers = ["o"]
        labels = ["F1"]
        for i in range(len(data[0,:])):
            mode = (i)//2
            value = i%2
            for j in range(len(data[:,i])):
                if not mt.isnan(data[j,i]):
                    test_data[mode, value, j] = data[j,i]
        #plt.figure('Damping')
        plt.subplot(1,2,2)
        plt.xlim([min(test_data[:,0,:].flatten()),max(test_data[:,0,:].flatten())])
        plt.ylim([min(test_data[:, 1, :].flatten()), max(test_data[:, 1, :].flatten())])
        for i in range(len(test_data[:,0,0])):
                plt.scatter(test_data[i,0,:], test_data[i,1,:], marker=markers[i], c = "black", label = labels[i])