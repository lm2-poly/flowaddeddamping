#Class definition of the nodes and methods to create them
#Author: Danick Lamoureux
#Project under Frédérick Gosselin and Sébastien Houde's supervision
#Date: 2022-05-04

import numpy as np
from gmeshing import *


class nodes:
    """Class containing the nodes of the simulation, both solid and fluid
    """
    def __init__(self):
        # List of the structural and aerodynamic nodes
        self.nodelist = np.array([[0,0,0]])
        self.aeronodes = np.array([[0, 0, 0]])
        self.fluidmesh = np.array([[0,0,0]])

    def profile_load(self, profile_filename, n, spacing):
        """Function to read a txt file separated by space of a 2D profile and get the coordinates

        Args:
            profile_filename (string): Filename for the profile with 2 columns, first being the x coordinates, the second being the y.
                                       The points should have decimals noted by dots.
            n (int): Number of hydrofoils in cascade
            spacing (float): Spacing between the hydrofoils
        """

        # Loading the profile points and putting them in a matrix
        profile = open(profile_filename, 'r')
        raw_coords = profile.readlines()
        npoints = len(raw_coords)
        if npoints < 3:
            print("Not enough points, please try again")
        profile_coords = np.zeros([npoints, 2])
        for i in range(npoints):
            line_coords = raw_coords[i].split(" ")
            line_coords = [item for item in line_coords if item != '']
            profile_coords[i, 0] = line_coords[0]
            profile_coords[i, 1] = line_coords[1]
        
        # Normalizing the profile
        profile_coords = np.array([profile_coords[:,0]/max(profile_coords[:,0]), profile_coords[:,1]]).T
        profile_len = len(profile_coords[:,0])
        for i in range(1, n):
            new_profile_coords = profile_coords[0:profile_len,:].copy()
            new_profile_coords[:,1] = new_profile_coords[:,1] + i*spacing
            profile_coords = np.append(profile_coords, new_profile_coords, axis = 0)

        self.profile_coords = profile_coords

        # Defining the camber line
        # The real camber line is not the mean of the points, but this should give an approximation
        meanx, meany = running_mean(profile_coords[:profile_len,0],profile_coords[:profile_len,1], deg = 5)
        self.camberline_coeffs = np.polyfit(meanx, meany, deg = 5)

    def envelope_load(self, envelope_chord, thick, n, spacing):
        """Defining the envelope around the hydrofoil

        Args:
            envelope_chord (float): Chord of the required envelope
            thick (float): Thickness of the required envelope per hydrofoil
            n (int): Number of hydrofoils in the global envelope
            spacing (float): Spacing between the hydrofoils in cascade
        """
        
        # The envelope is considered perfectly centered on the hydrofoil
        
        points = [0,0,0,0]
        points[0] = [0.5 + envelope_chord / 2, -thick / 2]
        points[1] = [0.5 + envelope_chord / 2, thick / 2 + (n-1)*spacing]
        points[2] = [0.5 - envelope_chord / 2, thick / 2 + (n-1)*spacing]
        points[3] = [0.5 - envelope_chord / 2, -thick / 2]
        self.envelope = np.asarray(points)

    def profile_mesh(self, mesh_size_normalized = 0.01, n = 1, show = False):
        """Meshing the profile using GMSH

        Args:
            mesh_size_normalized (float, optional): Mesh size. Defaults to 0.01.
            n (int, optional): Number of hydrofoils in cascade. Defaults to 1.
            show (bool, optional): Bool to detemine if showing the mesh using GMSH software is required. Defaults to False.
        """
        
        profile_coords = np.array([self.profile_coords[:,0],np.zeros(len(self.profile_coords[:,0])),self.profile_coords[:,1]]).T

        self.mesh_profile, self.profile_simplices = airfoil_gmeshing(profile_coords, mesh_size = mesh_size_normalized, n = n, show = show)

    def envelope_mesh(self, mesh_size_normalized = 0.01, n = 1, show = False):
        """Using GMSH, meshing the envelope similarly to the profile

        Args:
            mesh_size_normalized (float, optional): Mesh size. Defaults to 0.01.
            n (int, optional): Number of hydrofoils in cascade. Defaults to 1.
            show (bool, optional): Bool to detemine if showing the mesh using GMSH software is required. Defaults to False.
        """
        profile_coords = np.array([self.profile_coords[:,0],np.zeros(len(self.profile_coords[:,0])),self.profile_coords[:,1]]).T
        envelope_coords = np.array(
            [self.envelope[:, 0], np.zeros(len(self.envelope[:, 0])), self.envelope[:, 1]]).T

        self.envelope_mesh, self.envelope_simplices = envelope_gmeshing(envelope_coords, profile_coords,
                                                                      mesh_size = mesh_size_normalized, n = n, show = show)

    def wing_def(self, rootchord, tipchord, span, roottwist = 0, tiptwist = 0, sweep = 0, dihedral = 0,
                 nspan = 101, nx = 8, ny = 24, x0 = 0, n = 1, spacing = 0.5):
        """Function to define the wing (extruded profile) and define its nodes

        Args:
            rootchord (float): Rootchord of the wing
            tipchord (float): Tip chord of the wing
            span (float): Span of the wing
            roottwist (float, optional): Twist at the root [°] **Quarter chord stacked. Defaults to 0.
            tiptwist (float, optional): Twist at the tip [°] **Quarter chord stacked. Defaults to 0.
            sweep (float, optional): Sweep angle [°] **Angle relative to the leading edge. Defaults to 0.
            dihedral (float, optional): Dihedral angle [°] of the wing. Defaults to 0.
            nspan (int, optional): Number of solid elements on the span, must be odd. Defaults to 100.
            nx (int, optional): Number of aero-panels along the chord. Defaults to 8.
            ny (int, optional): Number of aero-panels nodes along the span. Defaults to 24.
            x0 (float, optional): Origin of the leading edge axially. Defaults to 0.
            n (int, optional): Number of hydrofoils in cascade. Defaults to 1.
            spacing (float, optional): Spacing of the hydrofoils in cascade. Defaults to 0.5.
        """

        self.nspan = nspan
        self.span = span
        self.chord = rootchord

        print("Meshing the 3D wing structure")

        # Converting the angles to rad
        roottwist = roottwist * np.pi / 180
        tiptwist = tiptwist * np.pi / 180
        sweep = sweep * np.pi / 180
        dihedral = dihedral * np.pi / 180

        # Defining the structural nodes
        npoints = len(self.mesh_profile)
        wing_mesh = np.array([[0,0,0]])
        wing_simplices = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])
        mesh_profile=np.array([self.mesh_profile[:,0],np.zeros(len(self.mesh_profile[:,0])),self.mesh_profile[:,1]]).T
        for i in range(nspan):
            # Defining the profile's dimensions at this point on the span
            chord_i = (tipchord - rootchord) * (i / (nspan - 1)) + rootchord
            twist_i = (tiptwist - roottwist) * (i / (nspan - 1)) + roottwist
            span_i = i / (nspan - 1) * span

            # Defining the wing's leading edge position at this point
            LE_x_i = span_i * np.tan(sweep) + (1 - np.cos(twist_i)) * 0.25 * chord_i
            quarterchord_z_i = span_i * np.tan(dihedral)
            LE_z_i = quarterchord_z_i + 0.25 * chord_i * np.sin(twist_i)
            LE_i = np.array([LE_x_i, span_i, LE_z_i])

            # Rotating the profil to account for the twist and pushing it by the leading edge position
            R = np.array([[np.cos(-twist_i), 0, np.sin(-twist_i)],
                          [0, 1, 0],
                          [-np.sin(-twist_i), 0, np.cos(-twist_i)]])
            # Pushing all the airfoil points by the LE distance
            delta_pos = LE_i * np.ones([npoints,3])
            profile_mesh_i = mesh_profile * chord_i @ R + delta_pos
            wing_mesh = np.append(wing_mesh, profile_mesh_i, axis=0)

        # Adding the simplices for each pentahedron
        non_write = np.array([0])
        for i in np.arange(0, nspan-1, 2):
            wing_simplices_i = np.append(self.profile_simplices+i*(npoints),self.profile_simplices[:,0:3]+(i+1)*(npoints), axis=1)
            non_write = np.append(non_write, self.profile_simplices[:,3]+(i+1)*(npoints),axis = 0)
            non_write = np.append(non_write, self.profile_simplices[:,4]+(i+1)*(npoints),axis = 0)
            non_write = np.append(non_write, self.profile_simplices[:, 5] + (i + 1) * (npoints), axis=0)
            wing_simplices_i = np.append(wing_simplices_i, self.profile_simplices+(i+2)*(npoints), axis = 1)
            wing_simplices = np.append(wing_simplices, wing_simplices_i, axis = 0)
        self.non_write = np.unique(non_write[1:])
        # Deleting the first row of the different matrices
        wing_coords = np.delete(wing_mesh,0,0)
        self.simplices = np.delete(wing_simplices,0,0)
        deletion = False
        if len(self.nodelist)==1:
            deletion = True
        self.nodelist = np.append(self.nodelist, wing_coords, axis=0)
        if deletion == True:
            self.nodelist = np.delete(self.nodelist,0,0)
        print("Eliminating useless nodes")
        self.nodelist, self.simplices = is_on_3Delement(nodelisting = self.nodelist, simplices = self.simplices, non_write = self.non_write)

        # Defining the aero nodes
        # For each ny span position, find the camber line and place nx points along the camber line to create each box
        # **For the moment, linear placement is the only choice
        # Placing the xyz coordinates

        print("Meshing the aerodynamic model")
        
        x = np.linspace(0,1,nx)

        # Defining the camberline
        camberline = np.polyval(self.camberline_coeffs, x)
        camberline = np.array([x, np.zeros(nx), camberline]).T
        camberline_len = len(camberline[:, 0])
        for i in range(1,n):
            new_camberline = camberline[0:camberline_len, :].copy()
            new_camberline[:, 2] = new_camberline[:, 2] + i * spacing
            camberline = np.append(camberline, new_camberline, axis=0)


        aero_mesh = np.array([[0, 0, 0]])
        aero_simplices = np.array([[0, 0, 0, 0]])
        for i in range(ny): # For all planes on the span
            # Defining the profile's dimensions at this point on the span
            chord_i = (tipchord - rootchord) * (i / (ny - 1)) + rootchord
            twist_i = (tiptwist - roottwist) * (i / (ny - 1)) + roottwist
            span_i = i / (ny - 1) * span

            # Defining the wing's leading edge position at this point
            LE_x_i = span_i * np.tan(sweep) + (1 - np.cos(twist_i)) * 0.25 * chord_i
            quarterchord_z_i = span_i * np.tan(dihedral)
            LE_z_i = quarterchord_z_i + 0.25 * chord_i * np.sin(twist_i)
            LE_i = np.array([LE_x_i, span_i, LE_z_i])
            if i == 0:
                self.LE_0 = LE_i
                self.rootchord = rootchord
            elif i == ny-1 :
                self.LE_span = LE_i
                self.tipchord = tipchord
            self.nx = nx
            self.ny = ny


            # Rotating the camberline to account for the twist and pushing it by the leading edge position
            R = np.array([[np.cos(-twist_i), 0, np.sin(-twist_i)],
                          [0, 1, 0],
                          [-np.sin(-twist_i), 0, np.cos(-twist_i)]])
            # Pushing all the airfoil points by the LE distance
            delta_pos = LE_i * np.ones([nx*n,3])
            camberline_i = camberline * chord_i @ R + delta_pos
            aero_mesh = np.append(aero_mesh, camberline_i, axis=0)
        camber_simplices = np.array([np.linspace(0, nx-2, nx-1), np.linspace(1,nx-1,nx-1)]).T
        for i in range(ny-1):
            aero_simplices_i = np.append(camber_simplices+i*(nx), camber_simplices+(i+1)*(nx),axis=1)
            aero_simplices = np.append(aero_simplices, aero_simplices_i, axis = 0)

        self.aero_simplices = np.delete(aero_simplices, 0, 0)
        aero_mesh = np.delete(aero_mesh,0,0)
        deletion = False
        if len(self.aeronodes)==1:
            deletion = True
        self.aeronodes = np.append(self.aeronodes, aero_mesh, axis=0)
        if deletion == True:
            self.aeronodes = np.delete(self.aeronodes,0,0)

    def flow_def(self):
        """Meshing the 3D fluid domain. Currently only works for straight wings. Similar to how the solid is built.
        """
        nspan = self.nspan
        span = self.span

        print("Meshing the 3D fluid volume")

        envelope_mesh = self.envelope_mesh
        npoints = len(envelope_mesh)

        # Calculating where to put the fluid nodes
        delta_span = span/(nspan-1)
        fluidmesh = np.array([[0,0,0]])
        for i in range(nspan):
            fluidmesh = np.append(fluidmesh, np.array([envelope_mesh[:,0]*self.chord, i*delta_span*np.ones(len(envelope_mesh)), envelope_mesh[:,1]*self.chord]).T, axis = 0)
        # Adding the simplices for each pentahedron
        fluid_simplices = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])
        fluid_non_write = np.array([0])
        for i in np.arange(0, nspan-1, 2):
            fluid_simplices_i = np.append(self.envelope_simplices + i * (npoints),
                                         self.envelope_simplices[:, 0:3] + (i + 1) * (npoints), axis=1)
            fluid_non_write = np.append(fluid_non_write, self.envelope_simplices[:,3]+(i+1)*(npoints),axis = 0)
            fluid_non_write = np.append(fluid_non_write, self.envelope_simplices[:,4]+(i+1)*(npoints),axis = 0)
            fluid_non_write = np.append(fluid_non_write, self.envelope_simplices[:, 5] + (i + 1) * (npoints), axis=0)
            fluid_simplices_i = np.append(fluid_simplices_i, self.envelope_simplices + (i + 2) * (npoints), axis=1)
            fluid_simplices = np.append(fluid_simplices, fluid_simplices_i, axis=0)
        self.fluid_non_write = np.unique(fluid_non_write[1:])
        # Deleting the first row of the different matrices
        fluidmesh = np.delete(fluidmesh, 0, 0)
        self.fluid_simplices = np.delete(fluid_simplices, 0, 0)
        deletion = False


        if len(self.fluidmesh) == 1:
            deletion = True
        self.fluidmesh = np.append(self.fluidmesh, fluidmesh, axis=0)
        if deletion == True:
            self.fluidmesh = np.delete(self.fluidmesh, 0, 0)
        self.fluidmesh, self.fluid_simplices = is_on_3Delement(nodelisting=self.fluidmesh, simplices=self.fluid_simplices, non_write = self.fluid_non_write)

    def GRIDWrite(self, file):
        """Writing the GRID cards. Consult NASTRAN's documentation for details.

        Args:
            file: File object to write to. Obtained using file = open(string, "w")
        """
        file.write('$* GRID CARDS\n')
        file.write('$*\n')
        for i in range(len(self.nodelist)):
            x=self.nodelist[i,0]
            y=self.nodelist[i,1]
            z=self.nodelist[i,2]
            #GRID, NID, , X, Y, Z
            file.write('GRID, '+str(i+1)+', , '+str(round(x,4))+', '+str(round(y,4))+', '+str(round(z,4))+', , 456 \n')
        file.write('$*\n')

    def GRIDWriteFluid(self, file):
        """Writing the GRID cards for fluid nodes. Consult NASTRAN's documentation for details.

        Args:
            file: File object to write to. Obtained using file = open(string, "w")
        """
        file.write('$* GRID CARDS\n')
        file.write('$*\n')
        ngrid = len(self.nodelist)
        for i in range(len(self.fluidmesh)):
            x = self.fluidmesh[i, 0]
            y = self.fluidmesh[i, 1]
            z = self.fluidmesh[i, 2]
            file.write('GRID, ' + str(i + 1 + ngrid) + ', , ' + str(round(x, 4)) + ', ' + str(round(y, 4)) + ', ' + str(
                round(z, 4)) + ',-1\n')
        file.write('$*\n')

def running_mean(x, y, deg = 3):
    """Used for camber lines, currently can only be flat points. Performs a moving average on physical points

    Args:
        x (array): x points
        y (array): y points
        deg (int, optional): Degree of the moving average. Defaults to 3.

    Returns:
        array, array: averaged data points
    """
    # Used for camber lines, currently can only be flat points
    # Running mean of xy points
    length = len(x)
    # Sorting the points according to the x coordinates
    indices = np.argsort(x)
    x = np.array(x[indices])
    y = np.array(y[indices])
    # Number of elements not present at the beginning and end of the array
    absent_length = (deg-1)/2
    # Mean of the points
    xnew = np.zeros(length - deg + 1)
    ynew = np.zeros(length - deg + 1)
    for i in range(len(xnew)):
        true_i = int(absent_length + i)
        xnew[i] = np.mean(x[int(true_i-absent_length): int(true_i+absent_length)])
        ynew[i] = np.mean(y[int(true_i-absent_length): int(true_i+absent_length)])
    return xnew, ynew

def is_on_3Delement(nodelisting, simplices, non_write):
    """Determine if a node is used by an element, if not, delete it

    Args:
        nodelisting (nodes): List of nodes
        simplices (list): List of nodes on elements
        non_write (_type_): Determines if nodes do not have to be written

    Returns:
        list, list: List of existing nodes on 3D elements and simplices
    """
    i = 0
    non_write = np.unique(non_write)
    p = 0
    nodelisting2 = np.zeros((len(nodelisting)-len(non_write),3), dtype = 'float32')
    r=0
    done = False
    for i in range(len(nodelisting)):
        used = True
        if i == non_write[p] and done == False:
            p += 1
            if p == len(non_write):
                p -= 1
                done = True
            used = False
            indices = np.argwhere(simplices > (i - p))
            row = indices[:,0]
            col = indices[:,1]
            simplices[row,col] -= 1
        if used == True:
            nodelisting2[r,:] = nodelisting[i, :]
            r = r + 1
    return nodelisting2, simplices