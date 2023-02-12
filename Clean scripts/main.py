#Evaluating the added damping of a fluid using Nastran aeroelastic module
#Author: Danick Lamoureux
#Project under Frédérick Gosselin and Sébastien Houde's supervision
#Date: 2022-05-04

import numpy as np
from nodes import *
from elements import *
from loads import *
from headings import *
from running import *
import time

if __name__ == '__main__':
    # Input data (mm-kg-mN-s)
    # Please use decimals to differenciate floats and integers

    # File names
    # Put the file names used for the analysis files (no extensions) and the profile (with extension)
    filename = r'TestFiles\NACA0003\NACA0003'
    profile = r'TestFiles\NACA0003\NACA0003.dat'

    # Available analysis types
    # static: Gravity load (mainly debugging oriented)
    # modes: Modal analysis in vacuum
    # uncoupled_acoustics: Modal analysis in resting fluid (real eigenvalues, no damping, no fluid-structure coupling))
    # coupled_acoustics: Modal analysis in resting fluid (complex eigenvalues, no damping, added mass)
    # aeroelastic: Modal analysis in fluid (with damping, no added mass or rigidity)
    # hydroelastic: Modal analysis in fluid (with added mass, damping and rigidity)
    analysis_type = "hydroelastic"

    # Dimensions: mm
    rootchord = 95.0
    tipchord = rootchord # Is currently required, future contributions will add support of tapered hydrofoils
    thickness = 3.4 #Maximum thickness of the hydrofoil
    span = 150.0
    roottwist = 0.0 # Future contributions will add support for twist
    tiptwist = 0.0 # Future contributions will add support for twist
    sweep = 0.0 # Future contributions will add support for sweep
    dihedral = 0.0 # Future contributions will add support for dihedral

    # Cascade parameters (Work in progress - future contributions will add support for cascades)
    n = 1 # Is currently required to be 1
    spacing = 39/250

    # The following dimensions are normalized by the chord
    envelope_chord = 7.5
    envelope_thickness = 100/rootchord

    # Mesh parameters
    mesh_size_solid = 6/rootchord # Mesh size for the profile, normalized by the chord
    mesh_size_fluid = 10/rootchord # Mesh size for the fluid envelope, normalized by the chord
    SpanDensity = 41 # Number of solid elements spanwise, must be odd
    nchord = 16 # Number of aerodynamic panels chordwise
    nspan = 16 # Number of aerodynamic panels spanwise

    # Material properties:
    # Aluminum
    E = 68.890E6 #Young's modulus in kPa (mN/mm^2)
    nu = 0.33 #Poisson's ratio
    rho_solid = 2711.0E-9 #Solid density in kg/mm^3

    # Uncomment the following lines to use Bronze instead    
    # Bronze
    # E = 115.0E6
    # nu = 0.33
    # rho_solid = 7800.0E-9

    # Fluid properties and characteristics
    # Water
    rho_flow = 997.00E-9  # Density for response analysis kg/mm^3
    rho_flow_aero = rho_flow # Density for the aeroelastic analysis - should be the same as rho_flow
    bulk = 2.21E6 # Bulk modulus in kPa
    soundspeed = np.sqrt(bulk/rho_flow)

    # Uncomment the following lines if an aeroelastic analysis is performed to use air instead of water
    # Air
    # rho_flow = 1.225E-9 # Density for response analysis kg/mm^3
    # rho_flow_aero = rho_flow  # Density for the aeroelastic analysis - should be the same as rho_flow
    # bulk = 142.0 # Bulk modulus in kPa
    # soundspeed = np.sqrt(bulk/rho_flow)

    # Minimum flow velocity to test in mm/s
    Urealmin = 5.0E3

    # Maximum flow velocity to test in mm/s
    Urealmax = 30.0E3
    
    # Number of modes to transfer from vibro-acoustic analysis to aeroelastic analysis
    nmodes = 10


########################################################################################################################
####################################################   Running   #######################################################
########################################################################################################################

    # You should not have to modify anything below unless you want to do modifications to the underlying script

    # Further calculations
    densities = [1.0]  # Density
    velocities = np.linspace(Urealmin,Urealmax,30)
    machs = [np.mean(velocities) / soundspeed]

    # Flow characteristics (reference values for reduced frequency)
    ref_velocity = (Urealmax+Urealmin)/2 #Reference velocity for aeroelastic analysis in mm/s
    ref_length = rootchord #Length for aeroelastic analysis

    # Aerodynamic matrix input
    mach_matrix = np.linspace(10**-3, 0.99, 30)
    freq_matrix = np.logspace(-2, 2, 30, base = 10)

    start_time = time.time()
    
    # Overall model
    wing = model(filename, profile, analysis_type)

    # Geometry of the wing
    geom = geometry(rootchord = rootchord,
                    tipchord = tipchord,
                    thickness = thickness,
                    span = span,
                    envelope_chord = envelope_chord,
                    thick = envelope_thickness,
                    roottwist = roottwist,
                    tiptwist = tiptwist,
                    sweep = sweep,
                    dihedral = dihedral,
                    n_hydrofoils=n,
                    spacing=spacing)

    # Mesh parameters
    mesh_params = mesh(solid_mesh_size = mesh_size_solid,
                       fluid_mesh_size = mesh_size_fluid,
                       nspan = SpanDensity,
                       nx = nchord,
                       ny = nspan)

    # Definition of the hydrofoil
    solid1 = solid(E = E,
                nu = nu,
                rho = rho_solid)

    # Definition of the acoustic fluid
    fluid1 = fluid(rho = rho_flow,
                  bulk = bulk)
    
    # Definition of the flow in the aeroelastic analysis
    flow1 = flow(ref_velocity = ref_velocity,
                 ref_length = ref_length,
                 rho_flow = rho_flow_aero,
                 velocities = velocities,
                 density_ratios = densities,
                 machs = machs,
                 mach_matrix = mach_matrix,
                 freq_matrix = freq_matrix)
    
    # Setting up the analysis
    wing.setup(geometry_object = geom,
               mesh_object = mesh_params,
               solid_object = solid1,
               fluid_object = fluid1,
               flow_object = flow1,
               show = False) # show = True to show GMSH software
    
    # Writing the analysis to .bdf file
    wing.write()

    # Run the .bdf file using NASTRAN
    wing.run(nmodes = nmodes)

    end_time = time.time()

    execution_time = end_time - start_time

    # Perform the analysis of the performed simulation. For further information, use Simcenter 3D by opening the .op2 files
    results = wing.analyse(show = True)
    print("Other results can be analyzed in Simcenter 3D using the OP2 file produced")
    print("Total execution time: "+str(execution_time)+' seconds')
    plt.show()