#From the nodes and simplices, generating the elements for Nastran
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-05

import numpy as np
from random import uniform

def MAT1(file, MID, E, nu, rho):
    #This function assumes an isotropic material
    #See documentation for the MAT1 command in Nastran
    file.write('$* Material used\n')
    file.write('MAT1, '+str(MID)+', '+str(E)+', , '+str(nu)+', '+str(rho)+'\n')

def MAT10(file, MID, BULK, RHO):
    #See documentation for the MAT10 command in Nastran
    file.write('$* Fluid used\n')
    file.write('MAT10, '+str(MID)+', '+str(BULK)+', '+str(RHO)+'\n')

def PSOLID(file, PID, MID):
    #See documentation for PSOLID
    file.write('$* Solid property\n')
    file.write('PSOLID, '+str(PID)+', '+str(MID)+'\n')

def PFLUID(file, PID, MID):
    #See documentation for PSOLID
    file.write('$* Fluid property\n')
    file.write('PSOLID, '+str(PID)+', '+str(MID)+', , , , , PFLUID, \n')

def PAERO(file, PID):
    #See documentation for PAERO1
    file.write('$* Aerodynamic panels properties\n')
    file.write('PAERO1, '+str(PID)+'\n')

def PAABSF(file, PID, RHOC): #OBSOLETE
    #See documentation for PAABSF
    #Currently unused
    file.write('$* Acoustic absorbers properties\n')
    file.write('PAABSF, '+str(PID)+', , , , , '+str(RHOC)+', , '+str(RHOC)+'\n')

def CPENTA(file, PID, nodes_object):
    #See documentation for CPENTA, using 15 nodes
    file.write('$* Element CARDS: Pentahedrons\n')
    file.write('$*\n')
    used_nodes = []
    for i in range(len(nodes_object.simplices)):
        EID = i
        G1 = nodes_object.simplices[i,0]
        G2 = nodes_object.simplices[i, 1]
        G3 = nodes_object.simplices[i, 2]
        G7 = nodes_object.simplices[i, 3]
        G8 = nodes_object.simplices[i,4]
        G9 = nodes_object.simplices[i,5]
        G10 = nodes_object.simplices[i,6]
        G11 = nodes_object.simplices[i,7]
        G12 = nodes_object.simplices[i,8]
        G4 = nodes_object.simplices[i,9]
        G5 = nodes_object.simplices[i,10]
        G6 = nodes_object.simplices[i,11]
        G13 = nodes_object.simplices[i,12]
        G14 = nodes_object.simplices[i,13]
        G15 = nodes_object.simplices[i,14]
        file.write('CPENTA, '+str(EID+1)+', '+str(PID)+', '+str(G1+1)+', '+str(G2+1)+', '+str(G3+1)+', '+str(G4+1)+', '+\
                    str(G5+1)+', '+str(G6+1)+'\n, '+str(G7+1)+', '+str(G8+1)+', '+str(G9+1)+', '+str(G10+1)+', '\
                   +str(G11+1)+', '+str(G12+1)+', '+str(G13+1)+', '+str(G14+1)+'\n, '+str(G15+1)+'\n')
        #Not all nodes are used, we save the nodes that are used for writing the used nodes
        used_nodes.append(G1)
        used_nodes.append(G2)
        used_nodes.append(G3)
        used_nodes.append(G4)
        used_nodes.append(G5)
        used_nodes.append(G6)
        used_nodes.append(G7)
        used_nodes.append(G8)
        used_nodes.append(G9)
        used_nodes.append(G10)
        used_nodes.append(G11)
        used_nodes.append(G12)
        used_nodes.append(G13)
        used_nodes.append(G14)
        used_nodes.append(G15)
    
    file.write('$*\n')
    return used_nodes

def CAABSF(file, PID, nodes_object): #OBSOLETE
    front_nodes = []
    back_nodes = []
    ngrid = len(nodes_object.nodelist)
    for i in range(len(nodes_object.fluidmesh)):
        x = nodes_object.fluidmesh[i, 0]
        if x == min(nodes_object.fluidmesh[:,0]):
            front_nodes.append(i)
        elif x == max(nodes_object.fluidmesh[:,0]):
            back_nodes.append(i)
    front_coords = nodes_object.fluidmesh[front_nodes,1:]
    back_coords = nodes_object.fluidmesh[back_nodes, 1:]
    from scipy.spatial import Delaunay
    tri_front = Delaunay(front_coords)
    tri_back = Delaunay(back_coords)
    for i in range(len(tri_front.simplices)):
        element = tri_front.simplices[i]
        G1 = ngrid+1+front_nodes[element[0]]
        G2 = ngrid + 1 + front_nodes[element[1]]
        G3 = ngrid + 1 + front_nodes[element[2]]
        EID = i + len(nodes_object.simplices) + len(nodes_object.fluid_simplices) + 1
        file.write('CAABSF, '+str(EID)+', '+str(PID)+', '+str(G1)+', '+str(G2)+', '+str(G3)+'\n')

    for i in range(len(tri_back.simplices)):
        element = tri_back.simplices[i]
        G1 = ngrid + 1 + back_nodes[element[0]]
        G2 = ngrid + 1 + back_nodes[element[1]]
        G3 = ngrid + 1 + back_nodes[element[2]]
        EID = i + len(nodes_object.simplices) + len(nodes_object.fluid_simplices) + len(tri_front.simplices) + 1
        file.write('CAABSF, ' + str(EID) + ', ' + str(PID) + ', ' + str(G1) + ', ' + str(G2) + ', ' + str(G3) + '\n')

def nodes_to_use(nodes_object):
    #We save the used nodes in order to only write those and not useless nodes
    used_nodes = []
    for i in range(len(nodes_object.simplices)):
        EID = i
        G1 = nodes_object.simplices[i, 0]
        G2 = nodes_object.simplices[i, 1]
        G3 = nodes_object.simplices[i, 2]
        G7 = nodes_object.simplices[i, 3]
        G8 = nodes_object.simplices[i, 4]
        G9 = nodes_object.simplices[i, 5]
        G10 = nodes_object.simplices[i, 6]
        G11 = nodes_object.simplices[i, 7]
        G12 = nodes_object.simplices[i, 8]
        G4 = nodes_object.simplices[i, 9]
        G5 = nodes_object.simplices[i, 10]
        G6 = nodes_object.simplices[i, 11]
        G13 = nodes_object.simplices[i, 12]
        G14 = nodes_object.simplices[i, 13]
        G15 = nodes_object.simplices[i, 14]
        used_nodes.append(G1)
        used_nodes.append(G2)
        used_nodes.append(G3)
        used_nodes.append(G4)
        used_nodes.append(G5)
        used_nodes.append(G6)
        used_nodes.append(G7)
        used_nodes.append(G8)
        used_nodes.append(G9)
        used_nodes.append(G10)
        used_nodes.append(G11)
        used_nodes.append(G12)
        used_nodes.append(G13)
        used_nodes.append(G14)
        used_nodes.append(G15)
    return used_nodes

def CPENTAFluid(file, PID, nodes_object):
    # See documentation for CPENTA for fluids, using 15 nodes
    # Same process as the CPENTA function
    file.write('$* Fluid Element CARDS: Pentahedrons\n')
    file.write('$*\n')
    ngrid = len(nodes_object.nodelist)
    used_nodes = []
    for i in range(len(nodes_object.fluid_simplices)):
        EID = i + len(nodes_object.simplices)

        G1 = nodes_object.fluid_simplices[i,0] + ngrid
        G2 = nodes_object.fluid_simplices[i, 1] + ngrid
        G3 = nodes_object.fluid_simplices[i, 2] + ngrid
        G7 = nodes_object.fluid_simplices[i, 3] + ngrid
        G8 = nodes_object.fluid_simplices[i,4] + ngrid
        G9 = nodes_object.fluid_simplices[i,5] + ngrid
        G10 = nodes_object.fluid_simplices[i,6] + ngrid
        G11 = nodes_object.fluid_simplices[i,7] + ngrid
        G12 = nodes_object.fluid_simplices[i,8] + ngrid
        G4 = nodes_object.fluid_simplices[i,9] + ngrid
        G5 = nodes_object.fluid_simplices[i,10] + ngrid
        G6 = nodes_object.fluid_simplices[i,11] + ngrid
        G13 = nodes_object.fluid_simplices[i,12] + ngrid
        G14 = nodes_object.fluid_simplices[i,13] + ngrid
        G15 = nodes_object.fluid_simplices[i,14] + ngrid

        file.write('CPENTA, '+str(EID+1)+', '+str(PID)+', '+str(G1+1)+', '+str(G2+1)+', '+str(G3+1)+', '+str(G4+1)+', '+\
                    str(G5+1)+', '+str(G6+1)+'\n, '+str(G7+1)+', '+str(G8+1)+', '+str(G9+1)+', '+str(G10+1)+', '\
                   +str(G11+1)+', '+str(G12+1)+', '+str(G13+1)+', '+str(G14+1)+'\n, '+str(G15+1)+'\n')
    file.write('$*\n')
    used_nodes.append(G1)
    used_nodes.append(G2)
    used_nodes.append(G3)
    used_nodes.append(G4)
    used_nodes.append(G5)
    used_nodes.append(G6)
    used_nodes.append(G7)
    used_nodes.append(G8)
    used_nodes.append(G9)
    used_nodes.append(G10)
    used_nodes.append(G11)
    used_nodes.append(G12)
    used_nodes.append(G13)
    used_nodes.append(G14)
    used_nodes.append(G15)

def CAERO(file, PID, nodes_object, n= 1, spacing = 50, IGID = 1, fluid_nodes_len = 0):
    #Only uncambered hydrofoils are considered for the moment
    file.write('$* Element CARDS: Aero panels\n')
    nelm = len(nodes_object.simplices) + fluid_nodes_len*2 + 1
    nboxes = nodes_object.ny*nodes_object.nx
    #Writing the CAERO1 entries
    for i in range(n):
        file.write('CAERO1, '+str(nelm+i*nboxes)+', '+str(PID)+', , '+str(nodes_object.ny)+', , , 1, '+\
                   str(IGID)+',\n')
        file.write(', '+str(round(nodes_object.LE_0[0],3))+', '+str(round(nodes_object.LE_0[1],3))+', '+\
                   str(round(nodes_object.LE_0[2]+i*spacing,3))+', '+str(round(nodes_object.rootchord,3))+', '+\
                   str(round(nodes_object.LE_span[0],3))+', '+str(round(nodes_object.LE_span[1],3))+', '+\
                   str(round(nodes_object.LE_span[2]+i*spacing,3))+', '+str(round(nodes_object.tipchord,3))+'\n')

    #Writing the chord distribution (cosine)
    x = np.linspace(np.pi, 0, (nodes_object.nx+1))
    AEFACT2 = np.cos(x)
    AEFACT = [(i - min(AEFACT2)) / (max(AEFACT2) - min(AEFACT2)) for i in AEFACT2]
    file.write("""AEFACT, 1, """)
    for i in range(len(AEFACT)):
        if i<7:
            file.write(str(round(AEFACT[i],4))+', ')
            if i==6:
                file.write('\n, ')
        if i>=7:
            file.write(str(round(AEFACT[i],4)) + ', ')
            if (i-6)%8 == 0 and i>10 and i != len(AEFACT):
                file.write('\n, ')
    file.write('\n$*\n')

def SPLINE(file, nodes_object, fluid_nodes_len, EID, n=1, DZ = .1, METH = "IPS", USAGE = "BOTH"):
    #Only good for straight wings for the moment
    file.write("$ SPLINE OBJECTS\n")
    #Calcuting the end box
    BOX2 = nodes_object.nx*nodes_object.ny
    nelm = len(nodes_object.simplices) + 1
    nnodes = len(nodes_object.nodelist)
    nnodes_profiles = int(1.5*len(nodes_object.mesh_profile))
    nnodes_profile = int(len(nodes_object.mesh_profile)/n)
    nspan = int((nodes_object.nspan-1)/2)
    # Finding the splined nodes
    # Sorting the nodes in xyz
    x = nodes_object.nodelist[:, 0]
    y = nodes_object.nodelist[:, 1]
    z = nodes_object.nodelist[:, 2]
    ind = np.argsort(x)
    for i in range(n):
        height = max(z) - min(z)
        minheight = min(z) + i * height / n
        maxheight = min(z) + (i + 1) * height / n
        file.write("SPLINE1, "+str(EID+i)+', '+str(nelm+BOX2*i)+', '+str(nelm+BOX2*i)+', '+str(nelm+BOX2*(i+1)-1)+', '+str(i+1)+', '+str(DZ)+', '+\
                   METH+', '+USAGE+'\n')
        SET = []
        for k in range(nnodes):
            if z[k] <= maxheight and z[k] >= minheight and k%2 == 0:
                SET.append(k)
        SET1 = [(str(x + 1) + ', ') for x in SET]
        SET = ["SET1, ", str(i+1), ", "]
        for ID in SET1:
            SET.append(ID)
        index = np.arange(10, len(SET), 8)
        file.write("".join(SET[0:10]) + "\n")
        for k in range(len(index)):
            if k != len(index) - 1:
                file.write(", " + "".join(SET[index[k]:index[k + 1]]) + "\n")
            else:
                file.write(", " + "".join(SET[index[k]:]))
        file.write('\n$*\n')


