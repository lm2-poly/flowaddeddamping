#Using gmsh to define the points on the profile
#Author: Danick Lamoureux, based on Olivier Duchesne and David Lessard's code
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-04

import gmsh
import numpy as np
import os

def airfoil_gmeshing(polygon_vertices, mesh_size, n=1, show = False):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("modele_test")

    npoints = int(len(polygon_vertices)/n)

    geom = gmsh.model.geo

    for vertex_index in range(len(polygon_vertices)): #Adding all the points of the profile
        gmsh.model.geo.addPoint(polygon_vertices[vertex_index,0],polygon_vertices[vertex_index,1],
                                polygon_vertices[vertex_index,2], mesh_size, vertex_index)

    #Adding the overall spline and creating the surface
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

def envelope_gmeshing(envelope_vertices, airfoil_vertices, mesh_size, n = 1, spacing = 0, show = False):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("modele_test")

    npoints_profile = int(len(airfoil_vertices)/n)
    npoints_envelope = len(envelope_vertices) - 1
    print(npoints_envelope)
    print(npoints_profile)
    npoints = len(envelope_vertices)+len(airfoil_vertices)-1

    geom = gmsh.model.geo

    for vertex_index in range(len(airfoil_vertices)): #Adding all the points of the profile
        gmsh.model.geo.addPoint(airfoil_vertices[vertex_index,0],airfoil_vertices[vertex_index,1],
                                airfoil_vertices[vertex_index,2], mesh_size, vertex_index)
        last_index = vertex_index

    for vertex_index in range(len(envelope_vertices)): #Adding all the points of the envelope
        gmsh.model.geo.addPoint(envelope_vertices[vertex_index,0],envelope_vertices[vertex_index,1],
                                envelope_vertices[vertex_index,2], mesh_size*8, vertex_index+last_index+1)
        new_last_index = vertex_index+last_index+1
    print(new_last_index)
    #Adding the overall spline and creating the surface
    gmsh.model.geo.addLine(new_last_index, new_last_index-1, 0)
    gmsh.model.geo.addLine(new_last_index-1, new_last_index-2, 1)
    gmsh.model.geo.addLine(new_last_index-2, new_last_index-3, 2)
    gmsh.model.geo.addLine(new_last_index-3, new_last_index, 3)
    gmsh.model.geo.addCurveLoop([0, 1, 2, 3], tag=10+2*n)
    
    # Adding the overall spline and creating the surface
    small_envelope_vertices = np.array(([[1.5,0,-0.25 ],[1.5,0,0.25+ (n-1) * spacing],[-0.5,0,0.25+ (n-1) * spacing],[-0.5,0,-0.25 ]]))
    for x in range(1,len(small_envelope_vertices)+1): #Adding all the points of the small envelope
            gmsh.model.geo.addPoint(small_envelope_vertices[x-1,0],small_envelope_vertices[x-1,1], small_envelope_vertices[x-1,2], mesh_size*4, x + new_last_index )
    
    gmsh.model.geo.addLine(new_last_index +  1, new_last_index +  2,  6 )
    gmsh.model.geo.addLine(new_last_index +  2, new_last_index +  3,  7 )
    gmsh.model.geo.addLine(new_last_index +  3, new_last_index +  4,  8 )
    gmsh.model.geo.addLine(new_last_index +  4, new_last_index +  1,  9 )
    gmsh.model.geo.addCurveLoop([6,  7, 8,  9], tag= 2 * n + 25)
    wiretags = [2 * n + 10]
    for i in range(n):


        gmsh.model.geo.addSpline(np.linspace(i * npoints_profile, (i + 1) * npoints_profile - 1, npoints_profile), 6 * i+4)
        gmsh.model.geo.addLine((i + 1) * npoints_profile - 1, i * npoints_profile, 6 * i + 5)
        gmsh.model.geo.addCurveLoop([6 * i + 4, 6 * i + 5], tag=10 * i + 2 * n + 20)
        wiretags.append(10 * i + 2 * n + 20)

        
    gmsh.model.geo.addPlaneSurface(wireTags= wiretags,tag = 0)
    #gmsh.model.geo.addPlaneSurface([ 2 * n + 10,  2 * n + 25 ],tag = 1)

    

    

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
    # Determine if a node is used by an element, if not, delete it
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
