#Evaluating the added damping of a fluid using Nastran AeroElastic module
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-04

import numpy as np
import gmsh
from nodes import *
from elements import *
from loads import *
from headings import *
from running import *
from analysis import *
import pyNastran
import os
import subprocess
import time
from coupling import *

if __name__ == '__main__':
    #Input data (mm-kg-mN-s)

    #File names
    #Put the file names used for the analysis files (no extensions) and the profile (with extension)
    filename = r'TestFiles\Bergan\F1'
    profile = r'TestFiles\Bergan\F1.dat'
    #results_filename = r'C:\Users\danic\Documents\3 - STG-CRSNG_E2022\HydroElasticNastran\TestFiles\NACA0003.csv'
    results_filename = None

    # Available analysis types
    # static: Gravity load (mainly debugging oriented)
    # modes: Modal analysis in vacuum
    # uncoupled_acoustics: Modal analysis in resting fluid (real eigenvalues, no damping, no fluid-structure coupling))
    # coupled_acoustics: Modal analysis in resting fluid (complex eigenvalues, no damping, added mass)
    # aeroelastic: Modal analysis in fluid (with damping, no added mass or rigidity)
    # simulated_modes: Modal analysis with added mass without acoustic fluid
    # hydroelastic: Modal analysis in fluid (with added mass, damping and rigidity) (In progress)
    analysis_type = "hydroelastic"

    #Dimensions: mm
    rootchord = 250.0
    tipchord = 250.0
    span = 150.0
    roottwist = 0.0
    tiptwist = 0.0
    sweep = 0.0
    dihedral = 0.0

    #Cascade parameters
    n = 1
    spacing = 39/250

    #Normalized by the chord
    #Vérifier la grandeur du volume fluide, ça peut affecter les résultats de manière assez importante
    envelope_chord = 2.0
    envelope_thickness = 150/rootchord

    #Mesh parameters
    mesh_size_solid = 5/rootchord #Mesh size for the profile, normalized by the chord
    mesh_size_fluid = 7/rootchord
    SpanDensity = 41 #Number of solid elements spanwise, must be odd
    nchord = 16 #Number of panels chordwise
    nspan = 16 #Number of panels spanwise

    #Material properties: N
    #Solid
    E = 68.890E6 #Young's modulus in kPa (mN/mm^2)
    nu = 0.33 #Poisson's ratio
    rho_solid = 2711.0E-9 #Solid density in kg/mm^3

    #Other material
    #E = 115.0E6
    #nu = 0.33
    #rho_solid = 7800.0E-9

    #Fluid properties and characteristics

    #Water
    rho_flow = 997.00E-9  # Density for response analysis kg/mm^3
    rho_flow_aero = rho_flow # Density for the aeroelastic analysis
    bulk = 2.21E6 #Bulk modulus in kPa
    soundspeed = np.sqrt(bulk/rho_flow)

    # # Air
    # rho_flow = 1.225E-9
    # rho_flow_aero = rho_flow  # Density for the aeroelastic analysis
    # bulk = 142.0
    # soundspeed = 340.0E3

    #Minimum water velocity to test, mm/s
    Urealmin = 5.0E3

    # Maximum water velocity to test, mm/s
    Urealmax = 30.0E3

    #Detailed velocity
    velocity = None

    #Minimum natural frequency to be expected
    fmin = 100 #Hz

    #Maximum natural frequency to be expected
    fmax = 1000 #Hz

########################################################################################################################
####################################################   Running   #######################################################
########################################################################################################################
    #Further calculations
    densities = [1.0]  # Density

    kfreqmax = fmax*rootchord/Urealmin
    kfreqmin = fmin*rootchord/Urealmax

    print("Velocities used go from "+str(Urealmin/1000)+" to "+str(Urealmax/1000)+" m/s")
    print("Mach numbers used go from " + str(Urealmin/soundspeed) + " to " + str(Urealmax/soundspeed))
    print("Reduced frequencies used go from " + str(kfreqmin) + " to " + str(kfreqmax))

    densities = [1.0]  # Density
    velocities = np.array([2.326783868, 5.015511892, 7.445708376, 9.927611169, 12.40951396, 14.94312306, 17.42502585, 20.06412691,
                  20.99276112, 22.44053775, 24.0434333, 27.92140641])*1000
    #velocities = np.array([5, 10, 15, 20, 25, 28])*1000
    #velocities = np.array([5, 8, 10, 15, 20, 25])*1000
    #velocities = np.linspace(Urealmin,Urealmax,10)
    machs = [np.mean(velocities) / soundspeed]

    #Flow characteristics (reference values for reduced frequency)
    ref_velocity = (Urealmax+Urealmin)/2 #Reference velocity for response analysis in mm/s
    ref_length = rootchord #Length for response analysis

    # Aerodynamic matrix input
    mach_matrix = np.linspace(10**-3, 0.99, 30) #velocities / soundspeed
    freq_matrix = np.logspace(-2, 2, 30, base = 10)

    start_time = time.time()

    #You should not have to modify anything below unless you want to do modifications to the script
    
    wing = model(filename, profile, analysis_type)

    geom = geometry(rootchord = rootchord,
                    tipchord = tipchord,
                    span = span,
                    envelope_chord = envelope_chord,
                    thick = envelope_thickness,
                    roottwist = roottwist,
                    tiptwist = tiptwist,
                    sweep = sweep,
                    dihedral = dihedral,
                    n_hydrofoils=n,
                    spacing=spacing)

    mesh_params = mesh(solid_mesh_size = mesh_size_solid,
                       fluid_mesh_size = mesh_size_fluid,
                       nspan = SpanDensity,
                       nx = nchord,
                       ny = nspan)

    alu = solid(E = E,
                nu = nu,
                rho = rho_solid)

    water = fluid(rho = rho_flow,
                  bulk = bulk)

    #Flow model: Optional depending on the analysis type
    fluid = flow(ref_velocity = ref_velocity,
                 ref_length = ref_length,
                 rho_flow = rho_flow_aero,
                 velocities = velocities,
                 density_ratios = densities,
                 machs = machs,
                 mach_matrix = mach_matrix,
                 freq_matrix = freq_matrix,
                 results_file = results_filename,
                 velocity = velocity)
    
    wing.setup(geometry_object = geom,
               mesh_object = mesh_params,
               solid_object = alu,
               fluid_object = water,
               flow_object = fluid)

    #wing.show()
    #exit()
    wing.write()

    wing.run(nmodes = 10)

    end_time = time.time()

    execution_time = end_time - start_time

    results = wing.analyse(scale = 10, show = True)

    print("Total execution time: "+str(execution_time)+' seconds')