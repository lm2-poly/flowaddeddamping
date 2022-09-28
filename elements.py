#From the nodes and simplices, generating the elements for Nastran
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-05

import numpy as np
from random import uniform

def MAT1(file, MID, E, nu, rho, GE = 0):
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

def simulatedCPENTA(file, PID, simplices, nsol, ngrid):
    #See documentation for CPENTA, using 15 nodes
    file.write('$* Element CARDS: Simulated Pentahedrons\n')
    file.write('$*\n')
    used_nodes = []
    for i in range(len(simplices)):
        EID = i + nsol
        G1 = simplices[i,0] + ngrid
        G2 = simplices[i, 1] + ngrid
        G3 = simplices[i, 2] + ngrid
        G7 = simplices[i, 3] + ngrid
        G8 = simplices[i,4] + ngrid
        G9 = simplices[i,5] + ngrid
        G10 = simplices[i,6] + ngrid
        G11 = simplices[i,7] + ngrid
        G12 = simplices[i,8] + ngrid
        G4 = simplices[i,9] + ngrid
        G5 = simplices[i,10] + ngrid
        G6 = simplices[i,11] + ngrid
        G13 = simplices[i,12] + ngrid
        G14 = simplices[i,13] + ngrid
        G15 = simplices[i,14] + ngrid
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

def simulatedCPENTA2(file, PID, simplices, nsol, nflu, ngrid):
    #See documentation for CPENTA, using 15 nodes
    file.write('$* Element CARDS: Simulated Pentahedrons\n')
    file.write('$*\n')
    used_nodes = []
    for i in range(len(simplices)):
        EID = i + nsol + nflu
        G1 = simplices[i,0] + ngrid
        G2 = simplices[i, 1] + ngrid
        G3 = simplices[i, 2] + ngrid
        G7 = simplices[i, 3] + ngrid
        G8 = simplices[i,4] + ngrid
        G9 = simplices[i,5] + ngrid
        G10 = simplices[i,6] + ngrid
        G11 = simplices[i,7] + ngrid
        G12 = simplices[i,8] + ngrid
        G4 = simplices[i,9] + ngrid
        G5 = simplices[i,10] + ngrid
        G6 = simplices[i,11] + ngrid
        G13 = simplices[i,12] + ngrid
        G14 = simplices[i,13] + ngrid
        G15 = simplices[i,14] + ngrid
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

def nodes_to_use(nodes_object):
    #Similarly to what is done in the previous function
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
    #Considering one flat plate for the moment
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

def CAEROS(file, PID, nodes_object, IGID = 1): #Obsolete
    #Elements of a CAERO function are all coupled with the same IGID (interaction between them)
    #MULTIPLE CAERO ENTRIES FOR CURVED CAMBER PROFILES, CURRENTLY USELESS
    file.write('$* Element CARDS: Aero panels\n')
    length = len(nodes_object.aero_simplices)
    aeronodes = nodes_object.aeronodes
    simplices = nodes_object.aero_simplices
    nodes = nodes_object.nodelist
    for i in range(length):
        #Identifying the 4 nodes of each panel
        node1 = aeronodes[int(simplices[i,0])]
        node2 = aeronodes[int(simplices[i,1])]
        node4 = aeronodes[int(simplices[i,2])]
        node3 = aeronodes[int(simplices[i,3])]
        nelm = len(nodes_object.simplices)+1
        chord12 = np.linalg.norm(node2-node1)
        chord43 = np.linalg.norm(node4-node3)
        #PROBLEM: All CAERO1 elements are aligned with the aerodynamic CSYS... No camber can be represented by this
        # implementation...
        #Writing the CAERO1 entries
        file.write('CAERO1, '+str(nelm+i)+', '+str(PID)+', , 1, 1, , , '+str(IGID)+',\n')
        file.write(', '+str(round(node1[0],3))+', '+str(round(node1[1],3))+', '+\
                   str(round(node1[2],3))+', '+str(round(chord12,3))+', '+str(round(node4[0],3))+\
                   ', '+str(round(node4[1],3))+', '+str(round(node4[2],3))+', '+\
                   str(round(chord43,3))+'\n')
    file.write('$*\n')


    #AEFACT is not used as it only specifies the coordinates of the divisions, for the moment only
    #uniform division is allowed

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


