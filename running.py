#Functions to run the code
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-11

import os
import time
import subprocess
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import read_op2
from nodes import *
from loads import *
from elements import *
from headings import *
from analysis import *
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
from coupling import *
import pickle
import scipy.linalg as scpl
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def naming(filename, analysis_type = None, create_folder = False):
    if analysis_type != None:
        #Extracting the non-extension part of the filename and adding the analysis type to it
        filename = list(filename)

        for i in range(len(filename)):
            j = -1-i
            if filename[j] == "\\":
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
    file_in, file_out, file_out_txt, file_verif, file_run, file_num, filename = naming(filename)
    results = read_op2(file_out, debug = False)
    return results

#Class of the different analyses:
class model:
    def __init__(self, filename, profile, analysis_type):
        #Defining the filenames and other important parameters
        self.filename = filename
        self.profile = profile
        #Making sure the analysis type exit
        if analysis_type == "static" or analysis_type == "modes" or analysis_type == "aeroelastic" or\
                analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics" or\
                analysis_type == 'hydroelastic' or analysis_type == 'simulated_modes':
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
        self.Kf = None
        self.AGG = None
        self.addedmass = False
        self.fn = np.zeros(1000)


    def setup(self, geometry_object, mesh_object, solid_object, fluid_object, flow_object = None):
        if self.analysis_type == "hydroelastic" or self.analysis_type == "simulated_modes":
            # For an hydroelastic analysis, define the different objects
            self.geom = geometry_object
            self.mesh = mesh_object
            self.solid = solid_object
            self.fluid = fluid_object
            self.flow = flow_object
            self.hydroelastic()
        else:
            # Nodes generation (GRID)
            nodelist = nodes()

            # Loading the profile
            nodelist.profile_load(self.profile, n = geometry_object.n_hydrofoils, spacing = geometry_object.hydrofoil_spacing)
            #Adding cascades
            #nodelist.cascade(n = geometry_object.n_hydrofoils, spacing = geometry_object.hydrofoils_spacing)
            nodelist.envelope_load(envelope_chord = geometry_object.envelope_chord, thick=geometry_object.thick,
                                   n = geometry_object.n_hydrofoils, spacing = geometry_object.hydrofoil_spacing)

            # Meshing the profile using gmsh
            # To show the mesh, show = True to use gmsh software
            nodelist.profile_mesh(mesh_size_normalized=mesh_object.solid_mesh_size, n=geometry_object.n_hydrofoils, show=False)
            if self.analysis_type == "uncoupled_acoustics" or self.analysis_type == "coupled_acoustics" or self.addedmass:
                nodelist.envelope_mesh(mesh_size_normalized=mesh_object.fluid_mesh_size, n = geometry_object.n_hydrofoils, show=False)

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


            if self.analysis_type == "aeroelastic":
                self.ref_velocity = flow_object.ref_velocity
                self.ref_length = flow_object.ref_length
                self.rho_flow = flow_object.rho_flow
                self.velocities = flow_object.velocities
                self.density_ratios = flow_object.density_ratios
                self.machs = flow_object.machs
                self.mach_matrix = flow_object.mach_matrix
                self.freq_matrix = flow_object.freq_matrix
                self.results_filename = flow_object.results_filename
                self.velocity = flow_object.velocity

    def show(self):
        #Selecting the modes accordingly
        if self.analysis_type != "hydroelastic":
            nodelist = self.nodelist
        else:
            nodelist = self.modes.nodelist
        # Showing the imported profile points
        nodelist.show_profile()

        # Showing the wing points
        nodelist.show_wing()

        # Showing the aero panels points for the wing
        nodelist.show_aero()

        plt.show(block = True)  # Show figures if necessary
        exit()

    def write(self, EXTOUT = False):
        analysis_type = self.analysis_type
        if analysis_type == "hydroelastic":
            print("Writing the modal analysis input file")
            self.modes.write(EXTOUT = False)
            print("Writing the coupled acoustics analysis input file")
            self.coupled.write(EXTOUT = False)

        else:
            nodelist = self.nodelist
            file_in = self.file_in
            f = open(file_in, 'w')
            # Defining the headers from the analysis type
            print("Writing headers")
            if self.addedmass == False:
                if analysis_type == "static":
                    headers_bending(f)
                elif analysis_type == "modes" and EXTOUT == False:
                    headers_modes(f)
                elif analysis_type == "modes" and EXTOUT == True:
                    headers_modes_EXTOUT(f)
                elif analysis_type == "uncoupled_acoustics":
                    headers_acoustics_real(f)
                elif analysis_type == "coupled_acoustics" and EXTOUT == False:
                    headers_acoustics_complex(f)
                elif analysis_type == "coupled_acoustics" and EXTOUT == True:
                    self.used_nodes = nodes_to_use(nodelist)
                    corresp_solid_nodes, corresp_fluid_nodes = nodelist.corresponding(self.used_nodes)
                    nsol = len(corresp_solid_nodes) * 6
                    nflu = len(corresp_fluid_nodes)
                    headers_acoustics_complex_EXTOUT(f, nsol, nflu)
                elif analysis_type == "aeroelastic":
                    headers_aero(f)
            elif self.addedmass == True and analysis_type == "aeroelastic" and EXTOUT == False:
                headers_hydro(f, P1 = self.P1, nmodes = self.nmodes)
            elif self.addedmass == True and analysis_type == 'modes' and EXTOUT == False:
                headers_simulated_modes(f)
            elif self.addedmass == True and analysis_type == "modes" and EXTOUT == True:
                headers_addedmass_modes_EXTOUT(f)


            # Writing the points using the GRID cards

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

            print("Writing the nodes")
            self.corresp_nodes = nodelist.GRIDWrite(f, self.used_nodes)
            if analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics":
                self.corresp_fluid_nodes = nodelist.GRIDWriteFluid(f)
            if self.addedmass == True and EXTOUT == True:
                print("Writing the simulating nodes")
                nodelist.simulateWrite(f, self.simNodes)
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
                     nodes_object=nodelist,
                     velocity=self.ref_velocity,
                     ref_length=self.ref_length,
                     rhoref=self.rho_flow)
                MKAERO(file=f,
                       SID=200,
                       mach_matrix=self.mach_matrix,
                       freq_matrix=self.freq_matrix,
                       velocities=self.velocities,
                       densities=self.density_ratios,
                       machs=self.machs,
                       velocity = self.velocity)

            # Defining the constraints and loads:
            print("Writing the constraints and loads")
            SPC(f, nodelist=nodelist.nodelist, profile_nodes_len = len(nodelist.mesh_profile))
            if analysis_type == "uncoupled_acoustics" or analysis_type == "coupled_acoustics":
                ACMODL(f)
                #SPCF(f, nodelist_object=nodelist)
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
            print("Temps requis par l'analyse modale dans le vide: "+str(time2-time1)+"s")

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
            print("Temps requis par l'analyse modale vibro-acoustique: " + str(time2 - time1) + "s")


            # Extracting the calculated modes with added mass
            op2 = read(self.coupled.filename)
            eigenvectors = op2.eigenvectors[1].data[:, :, :].real
            eigenvalues = op2.eigenvalues[''].eigenvalues
            modes_s = eigenvectors[:,:len(np.unique(self.modes.used_nodes)),:]

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

            n = min(len(LAMBDAvac), len(LAMBDAcoup))
            j=0
            k=0
            for i in range(n):
                if LAMBDAvac[i] < (10*2*np.pi)**2:
                    j += 1
                elif LAMBDAcoup[i] < (10*2*np.pi)**2:
                    k += 1
            n = min(len(LAMBDAvac[j:]), len(LAMBDAcoup[k:]),nmodes)
            LAMBDAvac = LAMBDAvac[j:n+j]
            LAMBDAcoup = LAMBDAcoup[k:n+k]
            phi_s = phi_s[:,j:n+j]
            phi_f = phi_f[:,k:n+k]

            n = len(LAMBDAcoup)
            self.aero.nmodes = n

            plt.figure('Compared Frequencies')
            plt.xlabel(r'$modes$')
            plt.ylabel(r'$f_n$')
            plt.title("Frequency of the different modes in vacuum and in resting water")
            plt.grid(True)
            mode_number = np.arange(1,n+1,1)
            plt.plot(mode_number, np.sqrt(LAMBDAvac), label='Vacuum')
            plt.plot(mode_number, np.sqrt(LAMBDAcoup), label='Resting fluid')
            plt.legend(loc='best')


            #Tant que le calcul des cascades n'est pas implémenté, garder cette section
            MAC = np.zeros([n,n])
            from numpy.linalg import norm
            for ii in range(n):
                for jj in range(n):
                    MAC[ii,jj] = norm(phi_s[:,ii].T@phi_f[:,jj])**2/(norm(phi_s[:,ii])**2*norm(phi_f[:,jj])**2)

            mode_index = np.zeros(n)
            for i in range(n):
                mode_index[i] = np.argmax(MAC[i,:])

            #Calcul de la première matrice de passage en fréquences
            invLAMBDAcoupmat = np.zeros([len(LAMBDAcoup),len(LAMBDAcoup)])
            LAMBDAvacmat = np.zeros([len(LAMBDAvac),len(LAMBDAvac)])
            LAMBDA = np.zeros([len(LAMBDAcoup),len(LAMBDAcoup)])
            self.aero.fn = np.zeros(n)
            for i in range(len(LAMBDAcoup)):
                invLAMBDAcoupmat[i,i] = np.sqrt(1/LAMBDAcoup[int(mode_index[i])])
                LAMBDAvacmat[i,i] = np.sqrt(LAMBDAvac[i])
                LAMBDA[i,i] = LAMBDAvac[i]
                self.aero.fn[i] = np.sqrt(LAMBDAcoup[i])/(2*np.pi)
            P1 = invLAMBDAcoupmat@LAMBDAvacmat

            #Calcul de la matrice de passage en modes (implémentation du calcul en cascades)
            #P2 = np.linalg.pinv(phi_s)@phi_f
            for i in range(n):
                phi_s[:,i] = phi_s[:,i]/norm(phi_s[:,i])
                phi_f[:, i] = phi_f[:, i] / norm(phi_f[:, i])
            phi_hh = np.linalg.pinv(phi_f)@phi_s
            for i in range(n):
                phi_hh[:,i] = phi_hh[:,i]/norm(phi_hh[:,i])

            # from numpy.linalg import inv
            # P1 = phi_hh.T
            # P2 = phi_hh
            #
            # Mhh = np.zeros([n,n])
            # Khh = np.zeros([n,n])
            # for i in range(n):
            #     mat = np.array([P1[:,i]]).T@np.array([P2[:,i]])
            #     Mhh += mat
            #     Khh += LAMBDA[i,i]*mat
            #
            # for i in range(n):
            #     diagM = Mhh[i,i]
            #     diagK = Khh[i,i]
            #     for j in range(n):
            #         if Mhh[i,j]/diagM < 10**-2:
            #             Mhh[i,j] = 0.0
            #         if Khh[i,j]/diagK < 10**-2:
            #             Khh[i,j] = 0.0
            #
            # Mhh2 = P1@P2
            # Khh2 = P1@LAMBDA@P2
            #
            # for i in range(n):
            #     diagM = Mhh2[i,i]
            #     diagK = Khh2[i,i]
            #     for j in range(n):
            #         if Mhh2[i,j]/diagM < 10**-2:
            #             Mhh2[i,j] = 0.0
            #         if Khh2[i,j]/diagK < 10**-2:
            #             Khh2[i,j] = 0.0
            #
            # w1, vl1 = scpl.eig(Khh, Mhh)
            # w2, vl2 = scpl.eig(Khh2, Mhh2)
            #
            # print('Mhh: '+str(np.allclose(Mhh,Mhh2)))
            # print('Khh: ' + str(np.allclose(Khh, Khh2)))

            #Fin de l'implémentation en cascades

            #Attribution des matrices de passage
            self.aero.P1 = P1

            self.aero.addedmass = True


            timeb = time.time()
            print("Temps requis par le calcul de la masse ajoutée: " + str(timeb - timea - fr_an - vac_an) + "s")

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

    def analyse(self, scale = 1, show=True, fn_vac = None, fn_coup = None):
        analysis_type = self.analysis_type
        file_out_txt = self.file_out_txt
        file_out = self.file_out

        if analysis_type == "modes" or analysis_type == "coupled_acoustics":
            coupling = False
            op2 = read(self.filename)
            modes = op2.eigenvectors[1].data
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
            
            used_nodes = np.unique(self.used_nodes)
            undeformed = self.nodelist.nodelist[used_nodes]
            d = len(undeformed)
            i=0
            p = 0
            if show == True:
                for mode_i in range(min(5,len(mode_number))):
                    mode = mode_number[mode_i]
                    deformation = modes[mode - 1,:,:]
                    Tx = deformation[:,0].real
                    Ty = deformation[:,1].real
                    Tz = deformation[:,2].real
                    deformation = np.array([Tx,Ty,Tz]).T
                    deformed = undeformed + deformation*scale
                    show_deformation(freq = fn[mode-1],
                                     nodes = deformed,
                                     mode = mode,
                                     scale = scale,
                                     coupling = coupling)

                plt.show(block=True)
            return fn, modes

        elif analysis_type == "uncoupled_acoustics":
            print("Please open Simcenter 3D to analyze the uncoupled acoustic modes")


        elif analysis_type == "hydroelastic":
            self.aero.analyse()#fn_vac = self.fn_vac, fn_coup = self.fn_coup)
        
        elif analysis_type == "aeroelastic":
            f = open(file_out_txt, 'r')
            rows = f.readlines()
            f.close()

            row2 = []
            for i in range(len(rows)):
                row = rows[i].split(' ')
                row = [item for item in row if item != '']
                row2.append(row)

            flutter = False
            fn=[0,0,0,0,0]
            for i in range(len(row2)):
                if len(row2[i]) > 1 and flutter == False:
                    if row2[i][1] == 'FLUTTER' and row2[i][0] != "PK" and row2[i][2] != "NOLIST":
                        flutter = True
                        i_start = i
                    if row2[i][0] == 'MODE' and row2[i][2] == 'EIGENVALUE':
                        fn[0] = np.sqrt(float(row2[i+2][2]))/(2*np.pi)
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
                            #Re[i,j,k,m] = states[i,j,k,m,6]
                            #Im[i,j,k,m] = states[i,j,k,m,5]
                            wd[i, j, k, m] = states[i, j, k, m, 6]
                            #wn[i, k, k, m] = states[i, j, k, m, 4] * 2 * np.pi
                            if fn_vac != None and fn_coup != None:
                                wn[i, k, k, m] = states[i, j, k, m, 4] * 2 * np.pi#*fn_coup[i+1]/fn_vac[i]
                            else:
                                wn[i, k, k, m] = states[i, j, k, m, 4] * 2 * np.pi
                            zeta[i,j,k,m] = -states[i,j,k,m,3]/2
                            #elif int(states[i,j,k,m,6])==0 and int(states[i,j,k,m,4])==0:
                                #wn[i, k, k, m] = -states[i,j,k,m,5]/(zeta[i,j,k,m]+np.sqrt((zeta[i,j,k,m])**2-1))
                            #zetawn = -states[i, j, k, m, 5]
                            #root = fsolve(lambda x: eigenvals(x,wd[i,j,k,m],zetawn), [wd[i,j,k,m], zetawn/wd[i,j,k,m]])
                            #wn[i,j,k,m] = root[0]
                            #zeta[i,j,k,m] = root[1]
                            coeff = [1, 2*zeta[i,j,k,m]*wn[i,j,k,m], wn[i,j,k,m]**2]
                            eigenvalue = np.roots(coeff)
                            if len(eigenvalue) == 2:
                                eigenvalue = eigenvalue[0]
                            Re[i,j,k,m] = eigenvalue.real
                            Im[i,j,k,m] = eigenvalue.imag
                            with np.errstate(divide='ignore'):
                                warnings.simplefilter("ignore")
                                #zeta[i, j, k, m] = zetawn / wn[i, j, k, m]
                                vstar[i,j,k,m] = velocities[m]/(self.ref_length*wn[i,j,k,m]/(2*np.pi))
            for i in range(5):
                data = np.array(
                [velocities, vstar[i, 0, 0, :], wn[i, 0, 0, :] / (2 * np.pi), zeta[i, 0, 0, :], Re[i, 0, 0, :],
                 Im[i, 0, 0, :]]).T
                file_num = list(self.filename)
                file_num.append(str(i)+".csv")
                file_num = "".join(file_num)
                np.savetxt(file_num, data, delimiter = ",")

            if show == True:
                #Velocity
                plt.figure('Frequency and Damping - Velocity - Constant', figsize = (12,6))
                plt.suptitle("Frequency and damping according to velocity")
                # subplot 1
                plt.subplot(1, 2, 2)
                x = [0,5,10,15,20,25,28]
                y = [0.02, 0.06, 0.095, 0.14, 0.18, 0.225, 0.27] - 0.02*np.ones(7)
                #plt.scatter(x,y, color = 'black', marker = 'o', label = 'CFX')
                #exp_data(self.results_filename)
                plt.ylabel(r'$\zeta$')
                plt.title("Adimensional damping according to velocity")
                plt.grid(True)
                for i in range(min(len(modes), 5)):
                    if True:
                        if True:
                            for j in range(len(machs)):
                                for k in range(len(densities)):
                                    plt.xlabel(r'$U$ [m/s]')
                                    plt.plot(velocities/1000, zeta[i, j, k, :],
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
                                plt.plot(velocities/1000, wn[i, j, k, :] / (2 * np.pi),
                                         label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Frequency and Damping - Velocity - Constant.png')

                plt.figure('Eigenvalues - Velocity - Constant', figsize = (12,6))
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
                                plt.plot(velocities/1000, Re[i, j, k, :],
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
                                plt.plot(velocities/1000, Im[i, j, k, :],
                                         label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Eigenvalues - Velocity - Constant.png')


                #Reduced velocity
                plt.figure('Frequency and Damping - Reduced Velocity - Constant', figsize = (12,6))
                plt.suptitle("Frequency and damping according to reduced velocity")
                # subplot 1
                plt.subplot(1, 2, 2)
                #exp_data(self.results_filename)
                plt.ylabel(r'Flow added damping $\zeta_f$')
                plt.title("Adimensional added damping according to reduced velocity")
                plt.grid(True)
                for i in range(min(len(modes), 5)):
                    if True:
                        if True:
                            for j in range(len(machs)):
                                for k in range(len(densities)):
                                    if self.fn[i] != 0:
                                        plt.xlabel(r'$v^{*} = \frac{U}{f_nL}$')
                                        plt.plot(velocities / (self.fn[i]*self.ref_length), zeta[i, j, k, :],
                                                 label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                        # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                                    else:
                                        plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                                        plt.plot(velocities / (wn[i, j, k, :]*self.ref_length / (2 * np.pi)), zeta[i, j, k, :],
                                                 label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                        # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Added Damping.pdf')

                # subplot 2
                plt.subplot(1, 2, 1)
                plt.ylabel(r'$f_n$')
                plt.title("Frequency according to reduced velocity")
                plt.grid(True)
                for i in range(min(len(modes), 5)):
                    if True:
                        for j in range(len(machs)):
                            for k in range(len(densities)):
                                if self.fn[i] != 0:
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_nL}$')
                                    plt.plot(velocities / (self.fn[i]*self.ref_length), wn[i, j, k, :] / (2 * np.pi),
                                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                                else:
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                                    plt.plot(velocities / (wn[i, j, k, :]*self.ref_length / (2 * np.pi)),
                                             wn[i, j, k, :] / (2 * np.pi),
                                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Frequency and Damping - Reduced Velocity - Constant.png')

                plt.figure('Eigenvalues - Reduced Velocity - Constant', figsize = (12,6))
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
                                if self.fn[i] != 0:
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_nL}$ [m]')
                                    plt.plot(velocities / (self.fn[i]*self.ref_length), Re[i, j, k, :],
                                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                                else:
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                                    plt.plot(velocities / (wn[i, j, k, :]*self.ref_length / (2 * np.pi)), Re[i, j, k, :],
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
                                if self.fn[i] != 0:
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_nL}$ [m]')
                                    plt.plot(velocities / (self.fn[i]*self.ref_length), Im[i, j, k, :],
                                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                                else:
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                                    plt.plot(velocities / (wn[i, j, k, :]*self.ref_length / (2 * np.pi)), Im[i, j, k, :],
                                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Eigenvalues - Reduced Velocity - Constant.png')

                # Variable frequencies
                # Reduced velocity
                plt.figure('Frequency and Damping - Reduced Velocity - Variable', figsize = (12,6))
                plt.suptitle("Frequency and damping according to reduced velocity")
                # subplot 1
                plt.subplot(1, 2, 2)
                #exp_data(self.results_filename)
                plt.ylabel(r'$\zeta$')
                plt.title("Adimensional damping according to reduced velocity")
                plt.grid(True)
                for i in range(min(len(modes), 5)):
                    if True:
                        if True:
                            for j in range(len(machs)):
                                for k in range(len(densities)):
                                    plt.xlabel(r'$v^{*} = \frac{U}{f_n(U)L}$ [m]')
                                    plt.plot(velocities / (wn[i, j, k, :] * self.ref_length / (2 * np.pi)),
                                             zeta[i, j, k, :],
                                             label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                    # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Added Damping.pdf')

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
                                plt.plot(velocities / (wn[i, j, k, :] * self.ref_length / (2 * np.pi)),
                                         wn[i, j, k, :] / (2 * np.pi),
                                         label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Frequency and Damping - Reduced Velocity - Variable.png')

                plt.figure('Eigenvalues - Reduced Velocity - Variable', figsize = (12,6))
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
                                plt.plot(velocities / (wn[i, j, k, :] * self.ref_length / (2 * np.pi)),
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
                                plt.plot(velocities / (wn[i, j, k, :] * self.ref_length / (2 * np.pi)),
                                         Im[i, j, k, :],
                                         label='Mode ' + str(modes[i]))  # + ', mach = ' + str(machs[j]) + \
                                # ', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')

                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Eigenvalues - Reduced Velocity - Variable.png')

                plt.figure("Argand", figsize = (8,6))
                plt.xlabel(r'$Re(\omega)$')
                plt.ylabel(r'$Im(\omega)$')
                plt.title("Argand's Diagram")
                plt.grid(True)
                plt.scatter(Re[0, 0, 0, 0], Im[0, 0, 0, 0], marker='o', color='black', label="Start")
                plt.scatter(Re[0, 0, 0, -1], Im[0, 0, 0, -1], marker='X', color='black', label="End")
                for i in range(min(len(modes),5)):
                    if True:
                        if True:
                            for j in range(len(machs)):
                                for k in range(len(densities)):
                                    plt.scatter(Re[i,j,k,0],Im[i,j,k,0],marker = 'o',color='black')
                                    plt.scatter(Re[i, j, k, -1], Im[i, j, k, -1], marker = 'X',color='black')
                                    plt.plot(Re[i,j,k,:],Im[i,j,k,:],
                                             label='Mode ' + str(modes[i]))# + ', mach = ' + str(machs[j]) + \
                                                   #', density = ' + str(densities[k]*self.rho_flow)+r'$kg/mm^3$')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig('Argand.png')

                plt.show(block = True)

            return states

    def hydroelastic(self):
        #Defining the analyses required for the hydroelastic analysis

        if self.analysis_type == "hydroelastic":
            # self.modes = model(self.filename,self.profile,'modes')
            self.coupled = model(self.filename, self.profile, 'coupled_acoustics')
            self.aero = model(self.filename, self.profile, 'aeroelastic')
            self.modes = model(self.filename, self.profile, 'modes')
        elif self.analysis_type == "simulated_modes":
            #self.modes = model(self.filename,self.profile,'modes')
            self.coupled = model(self.filename,self.profile,'coupled_acoustics')
            self.aero = model(self.filename, self.profile, 'modes')
            self.modes = model(self.filename, self.profile, 'modes')
            self.analysis_type = "hydroelastic"

        print("Setup the coupled acoustics analysis")
        self.coupled.setup(geometry_object=self.geom,
                           mesh_object=self.mesh,
                           solid_object=self.solid,
                           fluid_object=self.fluid,
                           flow_object=self.flow)

        print("Setup the modal analysis")
        self.modes.setup(geometry_object=self.geom,
                        mesh_object=self.mesh,
                        solid_object=self.solid,
                        fluid_object=self.fluid,
                        flow_object=self.flow)

        print("Setup the aeroelastic analysis with added mass")
        self.aero.setup(geometry_object=self.geom,
                         mesh_object=self.mesh,
                         solid_object=self.solid,
                         fluid_object=self.fluid,
                         flow_object=self.flow)

class geometry:
    #Defining the hydrofoil geometry
    def __init__(self, rootchord, tipchord, span, envelope_chord, thick, roottwist = 0.0, tiptwist = 0.0,
                 sweep = 0.0, dihedral = 0.0, n_hydrofoils = 1, spacing = 10):
        self.rootchord = rootchord
        self.tipchord = tipchord
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
    #Defining the mesh parameters for both solid and fluid
    def __init__(self, solid_mesh_size, fluid_mesh_size, nspan, nx, ny):
        self.solid_mesh_size = solid_mesh_size
        self.fluid_mesh_size = fluid_mesh_size
        self.nspan = nspan
        self.ny = ny
        self.nx = nx

class solid:
    #Defining the solid for dynamic analysis
    def __init__(self, E, nu, rho):
        self.E = E
        self.nu = nu
        self.rho = rho

class fluid:
    #Defining the fluid for acoustics analysis
    def __init__(self, rho, bulk):
        self.rho = rho
        self.bulk = bulk

class flow:
    #Defining the flow for flutter analysis
    def __init__(self, ref_velocity, ref_length, rho_flow, velocities, density_ratios, machs, mach_matrix, freq_matrix, results_file, velocity = None):
        self.ref_velocity = ref_velocity
        self.ref_length = ref_length
        self.rho_flow = rho_flow
        self.velocities = velocities
        self.density_ratios = density_ratios
        self.machs = machs
        self.mach_matrix = mach_matrix
        self.freq_matrix = freq_matrix
        self.results_filename = results_file
        self.velocity = velocity

def is_float(string):
    #Defining is a string contains a float
    try:
        float(string)
        return True
    except ValueError:
        return False
    
def show_deformation(freq, nodes, mode, scale, coupling):
    #Showing the wing structural nodes after applied scale deformation
    #If coupling, specify it
    if coupling == False:
        plt.figure('Mode '+str(mode)+' at scale = '+str(scale*100)+'% for fn = '+str(round(freq,2))+'Hz')
    else:
        plt.figure('Mode '+str(mode)+' at scale = '+str(scale*100)+'% for fn = '+str(round(freq,2))+'Hz under resting fluid-structure coupling')
    ax = plt.axes(projection='3d')
    ax.scatter3D(nodes[:,0],nodes[:,1],nodes[:,2])
    equal_axes(ax)
    ax.set_xlabel('chord')
    ax.set_ylabel('span')
    ax.set_zlabel('z')
    if coupling == False:
        ax.set_title(r'Mode ' + str(mode) + ' at scale = ' + str(scale * 100) + '%, $f_n = $' + str(
            round(freq, 2)) + 'Hz')
    else:
        ax.set_title(r'Mode '+str(mode)+' at scale = '+str(scale*100)+'%, $f_n = $'+str(round(freq,2))+'Hz under resting fluid-structure coupling')
    if coupling == False:
        plt.savefig('Mode '+str(mode)+' in vacuum.pdf')
    else:
        plt.savefig('Mode '+str(mode)+' in resting fluid.pdf')

def eigenvals(params, wd, zetawn):
    wn = params[0]
    zeta = params[1]
    alpha = zeta*wn
    beta = wn*np.sqrt(1-zeta**2)
    return alpha-zetawn, beta-wd

import pandas as pd
import math as mt
# def exp_data(filename):
#     if filename != None:
#         #From an excel file, plot the damping according to reduced velocity
#         data = pd.read_csv(filename, sep=';').to_numpy()
#         test_data = np.zeros([int(len(data[0,:])/2), 2, len(data[:,0])])
#         markers = ["1", "x","s", "D", "^", "v", "*", "+", "o"]
#         labels = ["M1", "M2", "M3", "M4", "H0", "H1", "H3", "F0", "F1"]
#         for i in range(len(data[0,:])):
#             mode = (i)//2
#             value = i%2
#             for j in range(len(data[:,i])):
#                 if not mt.isnan(data[j,i]):
#                     test_data[mode, value, j] = data[j,i]
#         #plt.figure('Damping')
#         plt.subplot(1,2,2)
#         plt.xlim([min(test_data[:,0,:].flatten()),max(test_data[:,0,:].flatten())])
#         plt.ylim([min(test_data[:, 1, :].flatten()), max(test_data[:, 1, :].flatten())])
#         for i in range(len(test_data[:,0,0])):
#             if i!=6:
#                 plt.scatter(test_data[i,0,:], test_data[i,1,:]-test_data[i,1,0], marker=markers[i], c = "black", label = labels[i])

def exp_data(filename):
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

if __name__ == '__main__':
    exp_data(r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\ExpData.csv')
    plt.legend(loc='best')
    plt.show()