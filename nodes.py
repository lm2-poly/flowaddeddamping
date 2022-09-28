#Class definition of the nodes and methods to create them
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-04

import numpy as np
import matplotlib.pyplot as plt
from gmeshing import *
from numba import njit


class nodes:
    def __init__(self):
        #List of the structural and aerodynamic nodes
        self.nodelist = np.array([[0,0,0]])
        self.aeronodes = np.array([[0, 0, 0]])
        self.fluidmesh = np.array([[0,0,0]])

    def profile_load(self, profile_filename, n, spacing):
        # Function to read a txt file separated by tabulations of a 2D profile and get the coordinates
        # Inputs:
        # profile_filename: Filename for the profile with 2 columns, first being the x coordinates, the second being the y.
        #                   The points should have decimals noted by dots.
        # Outputs: Adds the profile nodes to the nodelist

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

        profile_coords = np.array([profile_coords[:,0]/max(profile_coords[:,0]), profile_coords[:,1]]).T
        profile_len = len(profile_coords[:,0])
        for i in range(1, n):
            new_profile_coords = profile_coords[0:profile_len,:].copy()
            new_profile_coords[:,1] = new_profile_coords[:,1] + i*spacing
            profile_coords = np.append(profile_coords, new_profile_coords, axis = 0)

        self.profile_coords = profile_coords

        #Defining the camber line
        # The real camber line is not the mean of the points, but this should give an approximation
        meanx, meany = running_mean(profile_coords[:profile_len,0],profile_coords[:profile_len,1], deg = 5)
        self.camberline_coeffs = np.polyfit(meanx, meany, deg = 5)

    def envelope_load(self, envelope_chord, thick, n, spacing):
        #Defining the envelope around the hydrofoil
        points = [0,0,0,0]
        points[0] = [0.5 + envelope_chord / 2, -thick / 2]
        points[1] = [0.5 + envelope_chord / 2, thick / 2 + (n-1)*spacing]
        points[2] = [0.5 - envelope_chord / 2, thick / 2 + (n-1)*spacing]
        points[3] = [0.5 - envelope_chord / 2, -thick / 2]
        self.envelope = np.asarray(points)

    def profile_mesh(self, mesh_size_normalized = 0.01, n = 1, show = False):
        #Meshing the profile using GMSH
        profile_coords = np.array([self.profile_coords[:,0],np.zeros(len(self.profile_coords[:,0])),self.profile_coords[:,1]]).T

        self.mesh_profile, self.profile_simplices = airfoil_gmeshing(profile_coords, mesh_size = mesh_size_normalized, n = n, show = show)

    def envelope_mesh(self, mesh_size_normalized = 0.01, n = 1, show = False):
        #Using GMSH, meshing the envelope similarly to the profile
        profile_coords = np.array([self.profile_coords[:,0],np.zeros(len(self.profile_coords[:,0])),self.profile_coords[:,1]]).T
        envelope_coords = np.array(
            [self.envelope[:, 0], np.zeros(len(self.envelope[:, 0])), self.envelope[:, 1]]).T

        self.envelope_mesh, self.envelope_simplices = envelope_gmeshing(envelope_coords, profile_coords,
                                                                      mesh_size = mesh_size_normalized, n = n, show = show)

    def show_profile(self): #Obsolete
        #Showing the profile points and the camber line, used for debugging
        plt.figure('Profile')
        plt.scatter(self.profile_coords[:, 0], self.profile_coords[:, 1])
        x = np.linspace(0,1,50)
        y = np.polyval(self.camberline_coeffs, x)
        plt.plot(x,y)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.title('Profile of the structure')
        plt.grid(True)
        plt.xlim(min(self.profile_coords[:,0])-0.1,max(self.profile_coords[:,0])+0.1)
        plt.ylim(min(self.profile_coords[:, 1])-0.1, max(self.profile_coords[:, 1])+0.1)
        plt.gca().set_aspect('equal', adjustable='box')

    def show_profile_mesh(self): #Obsolete
        #Showing the profile mesh with the triangular elements, used for debugging
        plt.figure('Profile Mesh')
        plt.scatter(self.mesh_profile[:, 0], self.mesh_profile[:, 1])
        plt.scatter(self.profile_coords[:,0], self.profile_coords[:,1],color='red')
        plt.triplot(self.mesh_profile[:, 0], self.mesh_profile[:, 1], self.profile_simplices, linewidth = 2)
        plt.xlabel('x')
        plt.ylabel('z')
        plt.title('Mesh of the profile of the structure')
        plt.grid(True)
        plt.xlim(min(self.mesh_profile[:,0])-0.1,max(self.mesh_profile[:,0])+0.1)
        plt.ylim(min(self.mesh_profile[:, 1])-0.1, max(self.mesh_profile[:, 1])+0.1)
        plt.gca().set_aspect('equal', adjustable='box')

    def wing_def(self, rootchord, tipchord, span, roottwist = 0, tiptwist = 0, sweep = 0, dihedral = 0,
                 nspan = 100, nx = 8, ny = 24, x0 = 0, n = 1, spacing = 0.5):
        # Function to define the wing (extruded profile) and define its nodes
        # Inputs:
        #       rootchord: Chord at the root
        #       tipchord: Chord at the tip
        #       span: Length of the wing
        #       roottwist: Twist at the root [°] **Quarter chord stacked
        #       tiptwist: Twist at the tip [°] **Quarter chord stacked
        #       sweep: Sweep angle [°] **Angle relative to the leading edge
        #       dihedral: Dihedral angle [°]

        #       nspan: Number of elements along the span for the structural nodes
        #       nx: Number of aero-panels along the chord
        #       ny: Number of aero-panels nodes along the span

        self.nspan = nspan
        self.span = span
        self.chord = rootchord

        print("Meshing the 3D wing structure")

        #Converting the angles to rad
        roottwist = roottwist * np.pi / 180
        tiptwist = tiptwist * np.pi / 180
        sweep = sweep * np.pi / 180
        dihedral = dihedral * np.pi / 180

        #Defining the structural nodes
        npoints = len(self.mesh_profile)
        wing_mesh = np.array([[0,0,0]])
        wing_simplices = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])
        mesh_profile=np.array([self.mesh_profile[:,0],np.zeros(len(self.mesh_profile[:,0])),self.mesh_profile[:,1]]).T
        for i in range(nspan):
            #Defining the profile's dimensions at this point on the span
            chord_i = (tipchord - rootchord) * (i / (nspan - 1)) + rootchord
            twist_i = (tiptwist - roottwist) * (i / (nspan - 1)) + roottwist
            span_i = i / (nspan - 1) * span

            #Defining the wing's leading edge position at this point
            LE_x_i = span_i * np.tan(sweep) + (1 - np.cos(twist_i)) * 0.25 * chord_i
            quarterchord_z_i = span_i * np.tan(dihedral)
            LE_z_i = quarterchord_z_i + 0.25 * chord_i * np.sin(twist_i)
            LE_i = np.array([LE_x_i, span_i, LE_z_i])

            #Rotating the profil to account for the twist and pushing it by the leading edge position
            R = np.array([[np.cos(-twist_i), 0, np.sin(-twist_i)],
                          [0, 1, 0],
                          [-np.sin(-twist_i), 0, np.cos(-twist_i)]])
            # Pushing all the airfoil points by the LE distance
            delta_pos = LE_i * np.ones([npoints,3])
            profile_mesh_i = mesh_profile * chord_i @ R + delta_pos
            wing_mesh = np.append(wing_mesh, profile_mesh_i, axis=0)

        #Adding the simplices for each pentahedron
        non_write = np.array([0])
        for i in np.arange(0, nspan-1, 2):
            wing_simplices_i = np.append(self.profile_simplices+i*(npoints),self.profile_simplices[:,0:3]+(i+1)*(npoints), axis=1)
            non_write = np.append(non_write, self.profile_simplices[:,3]+(i+1)*(npoints),axis = 0)
            non_write = np.append(non_write, self.profile_simplices[:,4]+(i+1)*(npoints),axis = 0)
            non_write = np.append(non_write, self.profile_simplices[:, 5] + (i + 1) * (npoints), axis=0)
            wing_simplices_i = np.append(wing_simplices_i, self.profile_simplices+(i+2)*(npoints), axis = 1)
            wing_simplices = np.append(wing_simplices, wing_simplices_i, axis = 0)
        self.non_write = np.unique(non_write[1:])
        #Deleting the first row of the different matrices
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

        #Defining the aero nodes
        #For each ny span position, find the camber line and place nx points along the camber line to create each box
        #**For the moment, linear placement is the only choice
        #Placing the xyz coordinates

        print("Meshing the aerodynamic model")

        aero_nodes = np.zeros([3, nx, ny])
        x = np.linspace(0,1,nx)

        #Defining the camberline
        camberline = np.polyval(self.camberline_coeffs, x)
        camberline = np.array([x, np.zeros(nx), camberline]).T
        camberline_len = len(camberline[:, 0])
        for i in range(1,n):
            new_camberline = camberline[0:camberline_len, :].copy()
            new_camberline[:, 2] = new_camberline[:, 2] + i * spacing
            camberline = np.append(camberline, new_camberline, axis=0)


        aero_mesh = np.array([[0, 0, 0]])
        aero_simplices = np.array([[0, 0, 0, 0]])
        for i in range(ny): #For all planes on the span
            #Defining the profile's dimensions at this point on the span
            chord_i = (tipchord - rootchord) * (i / (ny - 1)) + rootchord
            twist_i = (tiptwist - roottwist) * (i / (ny - 1)) + roottwist
            span_i = i / (ny - 1) * span

            #Defining the wing's leading edge position at this point
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


            #Rotating the camberline to account for the twist and pushing it by the leading edge position
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
        #Currently only works for straight wings
        nspan = self.nspan
        span = self.span

        print("Meshing the 3D fluid volume")

        envelope_mesh = self.envelope_mesh
        npoints = len(envelope_mesh)

        #Calculating where to put the fluid nodes
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

    def show_aero(self): #Used for debugging, mostly obsolete
        #Showing the aero panels nodes
        plt.figure('Aero mesh')
        ax = plt.axes(projection='3d')
        ax.scatter3D(self.aeronodes[:,0],self.aeronodes[:,1],self.aeronodes[:,2])
        equal_axes(ax)
        ax.set_xlabel('chord')
        ax.set_ylabel('span')
        ax.set_zlabel('z')
        ax.set_title('Aero mesh')
        plt.show()

    def show_wing(self): #Used for debugging, now obsolete
        #Showing the wing structural nodes (gives an idea of the shape of the wing
        plt.figure('Wing structure')
        ax = plt.axes(projection='3d')
        ax.scatter3D(self.nodelist[:,0],self.nodelist[:,1],self.nodelist[:,2])
        equal_axes(ax)
        ax.set_xlabel('chord')
        ax.set_ylabel('span')
        ax.set_zlabel('z')
        ax.set_title('Wing structure')

    def show_mesh(self): #Obsolete, does not work anymore
        plt.figure("Mesh")
        ax=plt.axes(projection='3d')
        ax.plot3D(self.nodelist[:,0], self.nodelist[:,1], self.nodelist[:,2], 'o')
        for i in range(len(self.simplices)):
            j=np.linspace(0,len(self.simplices[i])-1,len(self.simplices[i]))
            j=np.append(j,0).astype(int)
            for k in range(len(j)-1):
                x2=self.nodelist[self.simplices[i][j[k+1]],0]
                x1=self.nodelist[self.simplices[i][j[k]],0]
                y2=self.nodelist[self.simplices[i][j[k+1]],1]
                y1=self.nodelist[self.simplices[i][j[k]],1]
                z2=self.nodelist[self.simplices[i][j[k+1]],2]
                z1=self.nodelist[self.simplices[i][j[k]],2]
                x=np.array([x1,x2])
                y=np.array([y1,y2])
                z=np.array([z1,z2])
                ax.plot3D(x, y, z, 'gray')
        equal_axes(ax)
        ax.set_xlabel('chord')
        ax.set_ylabel('span')
        ax.set_zlabel('z')
        ax.set_title('Mesh')
        plt.show()

    def GRIDWrite(self, file, used_nodes):
        #Writing the GRID cards
        file.write('$* GRID CARDS\n')
        file.write('$*\n')
        k=0
        corresp_nodes = []
        for i in range(len(self.nodelist)):
            x=self.nodelist[i,0]
            y=self.nodelist[i,1]
            z=self.nodelist[i,2]
            #GRID, NID, , X, Y, Z
            corresp_nodes.append(i)
            file.write('GRID, '+str(i+1)+', , '+str(round(x,4))+', '+str(round(y,4))+', '+str(round(z,4))+', , 456 \n')
        file.write('$*\n')
        return corresp_nodes

    def GRIDWriteFluid(self, file, addedmass = False):
        # Writing the GRID cards for fluid nodes
        file.write('$* GRID CARDS\n')
        file.write('$*\n')
        ngrid = len(self.nodelist)
        corresp_nodes = []
        k=0
        for i in range(len(self.fluidmesh)):
            x = self.fluidmesh[i, 0]
            y = self.fluidmesh[i, 1]
            z = self.fluidmesh[i, 2]
            corresp_nodes.append(i)
            file.write('GRID, ' + str(i + 1 + ngrid) + ', , ' + str(round(x, 4)) + ', ' + str(round(y, 4)) + ', ' + str(
                round(z, 4)) + ',-1\n')
        file.write('$*\n')
        return corresp_nodes

    def simulateWrite(self, file, fluidmesh):
        # Writing the GRID cards for fluid nodes
        file.write('$* FLUID SIMULATING GRID CARDS\n')
        file.write('$*\n')
        ngrid = len(self.nodelist)
        corresp_nodes = []
        for i in range(len(fluidmesh)):
            x = fluidmesh[i, 0]
            y = fluidmesh[i, 1]
            z = fluidmesh[i, 2]
            corresp_nodes.append(i)
            # GRID, NID, , X, Y, Z
            file.write('GRID, ' + str(i + 1 + ngrid) + ', , ' + str(round(x, 4)) + ', ' + str(round(y, 4)) + ', ' + str(
                round(z, 4)) + ', ,23456\n')
        file.write('$*\n')

    def simulateWrite2(self, file, fluidmesh, j = 1):
        # Writing the GRID cards for fluid nodes
        file.write('$* FLUID SIMULATING GRID CARDS\n')
        file.write('$*\n')
        ngrid = len(self.nodelist)
        corresp_nodes = []
        for i in range(len(fluidmesh)):
            nsimul = len(fluidmesh)
            x = fluidmesh[i, 0]
            y = fluidmesh[i, 1]
            z = fluidmesh[i, 2]
            corresp_nodes.append(i)
            file.write('GRID, ' + str(i + 1 + ngrid+nsimul) + ', , ' + str(round(x, 4)) + ', ' + str(round(y, 4)) + ', ' + str(
                round(z, 4)) + ', ,23456\n')
        file.write('$*\n')

    def simulateWriteSolidFluid(self, file, nodelist, fluidmesh):
        # Writing the GRID cards for fluid nodes
        file.write('$* FLUID SIMULATING GRID CARDS\n')
        file.write('$*\n')
        ngrid = len(self.nodelist)
        nsimul = len(fluidmesh)
        corresp_nodes = []
        for i in range(len(fluidmesh)):
            x = fluidmesh[i, 0]
            y = fluidmesh[i, 1]
            z = fluidmesh[i, 2]
            corresp_nodes.append(i)
            file.write('GRID, ' + str(i + 1 + ngrid) + ', , ' + str(round(x, 4)) + ', ' + str(
                round(y, 4)) + ', ' + str(
                round(z, 4)) + ', ,23456\n')
        for i in range(len(self.nodelist)):
            x = self.nodelist[i, 0]
            y = self.nodelist[i, 1]
            z = self.nodelist[i, 2]
            # GRID, NID, , X, Y, Z
            corresp_nodes.append(i)
            file.write('GRID, ' + str(i + ngrid + 1 + nsimul) + ', , ' + str(round(x, 4)) + ', ' + str(round(y, 4)) + ', ' + str(
                round(z, 4)) + ', , 456 \n')
        for i in range(len(fluidmesh)):
            x = fluidmesh[i, 0]
            y = fluidmesh[i, 1]
            z = fluidmesh[i, 2]
            corresp_nodes.append(i)
            file.write('GRID, ' + str(i + 1 + ngrid*2 + nsimul) + ', , ' + str(round(x, 4)) + ', ' + str(
                round(y, 4)) + ', ' + str(
                round(z, 4)) + ', ,23456\n')
        file.write('$*\n')

    def corresponding(self, used_nodes):
        #Similarly to what is done in the GRIDWrite functions, selecting how the written nodes correspond to existing nodes
        sol_corresp_nodes = []
        flu_corresp_nodes = []
        for i in range(len(self.nodelist)):
            if i in used_nodes:
                sol_corresp_nodes.append(i)
        for i in range(len(self.fluidmesh)):
            if all(self.fluid_non_write != i):
                flu_corresp_nodes.append(i)
        return sol_corresp_nodes, flu_corresp_nodes

def equal_axes(plot):
    #Doing what matlab does when calling "axis(equal)"

    #Getting the plot's limits
    x_limits = plot.get_xlim3d()
    y_limits = plot.get_ylim3d()
    z_limits = plot.get_zlim3d()

    #Calculating the ranges and mean values for each axis
    x_range = abs(x_limits[1] - x_limits[0])
    x_mean = (x_limits[1]+x_limits[0])/2
    y_range = abs(y_limits[1] - y_limits[0])
    y_mean = (y_limits[1]+y_limits[0])/2
    z_range = abs(z_limits[1] - z_limits[0])
    z_mean = (z_limits[1]+z_limits[0])/2

    #Calculating half of the maximum range
    max_half = 0.5 * max(x_range, y_range, z_range)

    #The ranges, to be equal, need to be all the same length, centered to their respective means
    plot.set_xlim3d([x_mean - max_half, x_mean + max_half])
    plot.set_ylim3d([y_mean - max_half, y_mean + max_half])
    plot.set_zlim3d([z_mean - max_half, z_mean + max_half])

def is_on_segment(p,q,r): #Obsolete
    # Does a line going right from the point p go through the segment joining points q and r?
    #Points are defined as [x coordinate, y coordinate]
    if q[0] != r[0]:
        m = (q[1] - r[1])/(q[0] - r[0])
        b = q[1] - m*q[0]
        yfit = m*p[0] + b
        if m>=0 and p[0]<=max(q[0],r[0]):
            if p[0]>=min(q[0],r[0]) and p[1]>=yfit and p[1]<=max(q[1],r[1]):
                return True
            elif p[0]<min(q[0],r[0]) and p[1]>=min(q[1],r[1]) and p[1]<=max(q[1],r[1]):
                return True
            else:
                return False

        if m<0 and p[0]<=max(q[0],r[0]):
            if p[0]>=min(q[0],r[0]) and p[1]<=yfit and p[1]>=min(q[1],r[1]):
                return True
            elif p[0]<min(q[0],r[0]) and p[1]>=min(q[1],r[1]) and p[1]<=max(q[1],r[1]):
                return True
            else:
                return False
    else:
        if p[0] <= q[0] and p[1]<=max(q[1],r[1]) and p[1]>=min(q[1],r[1]):
            return True
        else:
            return False

def running_mean(x, y, deg = 3): #Used for camber lines, currently obsolete
    #Running mean of xy points
    length = len(x)
    #Sorting the points according to the x coordinates
    x, y = zip(*sorted(zip(x,y)))
    x = np.array(x)
    y = np.array(y)
    #Number of elements not present at the beginning and end of the array
    absent_length = (deg-1)/2
    #Mean of the points
    xnew = np.zeros(length - deg + 1)
    ynew = np.zeros(length - deg + 1)
    for i in range(len(xnew)):
        true_i = int(absent_length + i)
        for j in range(deg):
            xnew[i] += x[int(true_i-absent_length+j)]
            ynew[i] += y[int(true_i-absent_length+j)]
        xnew[i] = xnew[i]/deg
        ynew[i] = ynew[i]/deg
    return xnew, ynew

def is_on_3Delement(nodelisting, simplices, non_write):
    # Determine if a node is used by an element, if not, delete it
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