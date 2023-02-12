#Using gmsh to define the points on the profile
#Author: Danick Lamoureux and Clément Audefroy, based on Olivier Duchesne and David Lessard's codes
#Project under Frédérick Gosselin and Sébastien Houde's supervision
#Date: 2022-05-04

import gmsh
import numpy as np
import os

def airfoil_gmeshing(polygon_vertices, mesh_size, n=1, show = False):
    """Meshing the airfoil using GMSH

    Args:
        polygon_vertices (list): Vertices of the airfoil
        mesh_size (float): Mesh size.
        n (int, optional): Number of hydrofoils in cascade. Defaults to 1.
        show (bool, optional): Bool to define if GMSH software is used to show the mesh. Defaults to False.

    Returns:
        list, list: List of nodes and elements
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("modele_test")

    npoints = int(len(polygon_vertices)/n)

    for vertex_index in range(len(polygon_vertices)): # Adding all the points of the profile
        gmsh.model.geo.addPoint(polygon_vertices[vertex_index,0],polygon_vertices[vertex_index,1],
                                polygon_vertices[vertex_index,2], mesh_size, vertex_index)

    # Adding the overall spline and creating the surface
    for i in range(n):
        gmsh.model.geo.addSpline(np.linspace(i*npoints,(i+1)*npoints-1,npoints),2*i)
        gmsh.model.geo.addLine((i+1)*npoints-1,i*npoints,2*i+1)
        gmsh.model.geo.addCurveLoop([2*i,2*i+1],tag = 10*i + 2*n + 1)
        gmsh.model.geo.addPlaneSurface(wireTags= [10*i + 2*n + 1],tag = i)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.Smoothing", 10)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)

    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size)
    gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber("Mesh.Format", 31)

    # Gmsh plotting
    if show:
        gmsh.fltk.run()

    # Getting all the nodes
    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes(includeBoundary = True)
    nNodes = int(max(nodeTags))
    nodelisting = np.zeros([nNodes,2])
    for i in range(len(nodeTags)):
        index = int(nodeTags[i]-1)

        nodelisting[index, 0] = coord[int(3*i)]
        nodelisting[index, 1] = coord[int(3*i+2)]

    filename = 'Polygon.bdf'
    gmsh.write(filename)

    # Finding the simplices from the written file
    f = open(filename,'r')
    rows = f.readlines()
    simplices = np.array([[0,0,0,0,0,0]])
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == 'CTRIA6':
            new_simplice = np.array([int(row[3])-1,int(row[4])-1,int(row[5])-1,int(row[6])-1,int(row[7])-1,int(row[8])-1])
            simplices = np.append(simplices,[new_simplice],axis=0)
    simplices = np.delete(simplices,0,0)
    f.close()

    nodelisting, simplices = is_on_element(nodelisting, simplices)

    gmsh.finalize()

    os.remove(filename)

    return nodelisting, simplices

def envelope_gmeshing (envelope_vertices, airfoil_vertices, mesh_size, n = 1, show = False):
    """Meshing the airfoil using GMSH

    Args:
        envelope_vertices (list): Vertices of the envelope
        airfoil_vertices (list): Vertices of the airfoil
        mesh_size (float): Mesh size.
        n (int, optional): Number of hydrofoils in cascade. Defaults to 1.
        show (bool, optional): Bool to define if GMSH software is used to show the mesh. Defaults to False.

    Returns:
        list, list: List of nodes and elements
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("modele_test")

    npoints_profile = int(len(airfoil_vertices)/n)

    small_envelope_vertices = np.array(([[1.5,0,-0.25],[1.5,0,0.25],[-0.5,0,0.25],[-0.5,0,-0.25]]))

    for vertex_index in range(len(airfoil_vertices)): # Adding all the points of the profile
        gmsh.model.geo.addPoint(airfoil_vertices[vertex_index,0],airfoil_vertices[vertex_index,1],
                                airfoil_vertices[vertex_index,2], mesh_size, vertex_index)
        last_index = vertex_index

    for vertex_index in range(len(small_envelope_vertices)): # Adding all the points of the small envelope
        gmsh.model.geo.addPoint(small_envelope_vertices[vertex_index,0],small_envelope_vertices[vertex_index,1],
                                small_envelope_vertices[vertex_index,2], mesh_size*4, vertex_index+last_index+1)
        new_last_index = vertex_index+last_index+1

    # Adding the overall spline  of the small enveloppe
    gmsh.model.geo.addLine(new_last_index, new_last_index-1, 0)
    gmsh.model.geo.addLine(new_last_index-1, new_last_index-2, 1)
    gmsh.model.geo.addLine(new_last_index-2, new_last_index-3, 2)
    gmsh.model.geo.addLine(new_last_index-3, new_last_index, 3)
    gmsh.model.geo.addCurveLoop([0, 1, 2, 3], tag=10+2*n)
    wiretags = [10+2*n]

    # Adding the overall spline of the profile and creating the surface
    for i in range(n):
        gmsh.model.geo.addSpline(np.linspace(i * npoints_profile, (i + 1) * npoints_profile - 1, npoints_profile), 2 * i+4)
        gmsh.model.geo.addLine((i + 1) * npoints_profile - 1, i * npoints_profile, 2 * i + 5)
        gmsh.model.geo.addCurveLoop([2 * i+4, 2 * i + 5], tag=10 * i + 2 * n + 20)
        wiretags.append(10 * i + 2 * n + 20)

    gmsh.model.geo.addPlaneSurface(wireTags= wiretags,tag = 0)


    big_envelope_vertices = envelope_vertices

    for vertex_index in range(len(small_envelope_vertices)): # Adding all the points of the small envelope
        gmsh.model.geo.addPoint(small_envelope_vertices[vertex_index,0],small_envelope_vertices[vertex_index,1],
                                small_envelope_vertices[vertex_index,2], mesh_size*4,900+ vertex_index)
        

    for vertex_index in range(len(big_envelope_vertices)): # Adding all the points of the big envelope
        gmsh.model.geo.addPoint(big_envelope_vertices[vertex_index,0],big_envelope_vertices[vertex_index,1],
                                big_envelope_vertices[vertex_index,2], mesh_size*12, 910+vertex_index)
        

    #Adding the line of the big envelope
    gmsh.model.geo.addLine(913, 912, 6)
    gmsh.model.geo.addLine(912, 911, 7)
    gmsh.model.geo.addLine(911, 910, 8)
    gmsh.model.geo.addLine(910, 913, 9)
 
    gmsh.model.geo.addCurveLoop([6,7,8,9], tag=30+2*n)
    wiretags2 = [30+2*n]

    # Adding the line of the small envelope and creating the surface
    
    gmsh.model.geo.addLine(903, 902, 10)
    gmsh.model.geo.addLine(902, 901, 11)
    gmsh.model.geo.addLine(901, 900, 12)
    gmsh.model.geo.addLine(900, 903, 13)

    
    gmsh.model.geo.addCurveLoop([10,11,12,13], tag=  2 * n + 40)
    wiretags2.append( 2 * n + 40)

    gmsh.model.geo.addPlaneSurface(wireTags= wiretags2,tag = 1)

    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.Smoothing", 10)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.model.mesh.generate(2)

    gmsh.option.setNumber("Mesh.Format", 31)

    # GMSH VISUALISATION (plot of the meshing)
    if show:
        gmsh.fltk.run()

    #Getting all the nodes
    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes(includeBoundary = True)
    nNodes = int(max(nodeTags))
    nodelisting = np.zeros([nNodes,2])
    for i in range(len(nodeTags)):
        index = int(nodeTags[i]-1)

        nodelisting[index, 0] = coord[int(3*i)]
        nodelisting[index, 1] = coord[int(3*i+2)]

    filename = 'Polygon.bdf'
    gmsh.write(filename)

    #Finding the simplices from the written file
    f = open(filename,'r')
    rows = f.readlines()
    simplices = np.array([[0,0,0,0,0,0]])
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == 'CTRIA6':
            new_simplice = np.array([int(row[3])-1,int(row[4])-1,int(row[5])-1,int(row[6])-1,int(row[7])-1,int(row[8])-1])
            simplices = np.append(simplices,[new_simplice],axis=0)
    simplices = np.delete(simplices,0,0)
    f.close()
    nodelisting, simplices = is_on_element(nodelisting, simplices)



    gmsh.finalize()

    os.remove(filename)

    return nodelisting, simplices

def is_on_element(nodelisting, simplices):
    """Determine if a node is used by an element, if not, delete it

    Args:
        nodelisting (nodes): List of nodes
        simplices (list): List of nodes on elements

    Returns:
        list, list: List of existing nodes on 2D elements and simplices
    """
    i = 0
    d = len(nodelisting)
    while i<d:
        used = False
        for j in range(len(simplices)):
            if i == simplices[j,0] or i == simplices[j,1] or i == simplices[j,2] or i == simplices[j,3] \
                    or i == simplices[j,4] or i == simplices[j,5]:
                used = True
                break
        if used == False:
            nodelisting = np.delete(nodelisting, i, axis = 0)
            for j in range(len(simplices)):
                for k in range(6):
                    if simplices[j,k]>i:
                        simplices[j,k] = simplices[j,k] - 1
            d = len(nodelisting)
            i = i - 1
        i = i + 1
    return nodelisting, simplices
