#Evaluating the added damping of a fluid using Nastran AeroElastic module
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-04

import numpy as np
from nodes import *
from elements import *
from loads import *
from headings import *
from running import *
import time

if __name__ == '__main__':
    #Input data (mm-kg-mN-s)

    #File names
    #Put the file names used for the analysis files (no extensions) and the profile (with extension)
    filename = r'TestFiles\Francis99\F0'
    profile = r'TestFiles\Francis99\F0.dat'

    # Available analysis types
    # static: Gravity load (mainly debugging oriented)
    # modes: Modal analysis in vacuum
    # uncoupled_acoustics: Modal analysis in resting fluid (real eigenvalues, no damping, no fluid-structure coupling))
    # coupled_acoustics: Modal analysis in resting fluid (complex eigenvalues, no damping, added mass)
    # aeroelastic: Modal analysis in fluid (with damping, no added mass or rigidity)
    # hydroelastic: Modal analysis in fluid (with added mass, damping and rigidity) (In progress)
    analysis_type = "aeroelastic"

    #Dimensions: mm
    rootchord = 250.0
    tipchord = 250.0
    span = 150.0
    roottwist = 0.0
    tiptwist = 0.0
    sweep = 0.0
    dihedral = 0.0

    #Cascade parameters (work in progress)
    n = 1
    spacing = 39/250

    #Normalized by the chord
    envelope_chord = 6
    envelope_thickness = 100/rootchord

    #Mesh parameters
    mesh_size_solid = 10/rootchord #Mesh size for the profile, normalized by the chord
    mesh_size_fluid = 10/rootchord
    SpanDensity = 31 #Number of solid elements spanwise, must be odd
    nchord = 16 #Number of panels chordwise
    nspan = 16 #Number of panels spanwise

    #Material properties: N
    #Aluminum
    E = 68.890E6 #Young's modulus in kPa (mN/mm^2)
    nu = 0.33 #Poisson's ratio
    rho_solid = 2711.0E-9 #Solid density in kg/mm^3

    #Bronze
    #E = 115.0E6
    #nu = 0.33
    #rho_solid = 7800.0E-9

    #Fluid properties and characteristics
    #Water
    rho_flow = 997.00E-9  # Density for response analysis kg/mm^3
    rho_flow_aero = rho_flow # Density for the aeroelastic analysis
    bulk = 2.21E6 #Bulk modulus in kPa
    soundspeed = np.sqrt(bulk/rho_flow)

    # Air
    #rho_flow = 1.225E-9
    #rho_flow_aero = rho_flow  # Density for the aeroelastic analysis
    #bulk = 142.0
    #soundspeed = 340.0E3

    #Minimum water velocity to test, mm/s
    Urealmin = 5.0E3

    # Maximum water velocity to test, mm/s
    Urealmax = 30.0E3

    #Detailed velocity
    velocity = None


########################################################################################################################
####################################################   Running   #######################################################
########################################################################################################################

    #You should not have to modify anything below unless you want to do modifications to the script

    #Further calculations
    densities = [1.0]  # Density
    velocities = np.linspace(Urealmin,Urealmax,30)
    machs = [np.mean(velocities) / soundspeed]

    #Flow characteristics (reference values for reduced frequency)
    ref_velocity = (Urealmax+Urealmin)/2 #Reference velocity for aeroelastic analysis in mm/s
    ref_length = rootchord #Length for aeroelastic analysis

    # Aerodynamic matrix input
    mach_matrix = np.linspace(10**-3, 0.99, 30) #velocities / soundspeed
    freq_matrix = np.logspace(-2, 2, 30, base = 10)

    start_time = time.time()
    
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

    solid1 = solid(E = E,
                nu = nu,
                rho = rho_solid)

    fluid1 = fluid(rho = rho_flow,
                  bulk = bulk)

    #Flow model: Optional depending on the analysis type
    flow1 = flow(ref_velocity = ref_velocity,
                 ref_length = ref_length,
                 rho_flow = rho_flow_aero,
                 velocities = velocities,
                 density_ratios = densities,
                 machs = machs,
                 mach_matrix = mach_matrix,
                 freq_matrix = freq_matrix,
                 velocity = velocity)
    
    wing.setup(geometry_object = geom,
               mesh_object = mesh_params,
               solid_object = solid1,
               fluid_object = fluid1,
               flow_object = flow1)
    
    wing.write()

    wing.run(nmodes = 10)

    end_time = time.time()

    execution_time = end_time - start_time

    results = wing.analyse(show = True)
    print("Other results can be analyzed in Simcenter 3D using the OP2 file produced")
    print("Total execution time: "+str(execution_time)+' seconds')