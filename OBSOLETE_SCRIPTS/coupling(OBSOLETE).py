#Coupling the coupled acoustics coupling matrix to the structural mass matrix
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-17

#Completely obsolete

import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib.pyplot as plt
import pickle
import time
import multiprocessing as mp
import concurrent.futures
import os
from numba import njit
import cProfile

def AGG_extract(AGGpch_file, structural_nodes, fluid_nodes, corresp_fluid_nodes, corresp_solid_nodes):
    print("Extracting the coupling matrix")
    f = open(AGGpch_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the punch file and only keep the mass ones
    row2 = []
    first_encounter = False
    line_start = 0
    line_break = 0
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG" and first_encounter == False:
            if row[1] == "AGGMAT":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG" and first_encounter == True:
            if row[1] != "AGGMAT":
                line_break = i
                break
    rows = rows[line_start:line_break]

    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    #Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    #Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)


    #Create the A matrix
    #nrows = nddl solid
    nrows = len(structural_nodes)*6
    #ncols = nddl fluid
    ncols = len(fluid_nodes)*6

    # Struct node, Struct DOF, Fluid node, Fluid DOF
    #AGG = np.zeros([nrows*6, ncols*1], dtype = 'float16')
    #AGG = sp.lil_array((nrows, ncols), dtype = "float32")
    AGG = sp.lil_array((len(corresp_solid_nodes) * 6, len(corresp_fluid_nodes)*6), dtype="float32")
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        fluid_node = int(rows[entry][2])-len(structural_nodes)-1
        used_fluid_node = fluid_node
        fluid_DOF = 0

        if i != len(DMIG_ind) - 1:
            for j in range(entry+1, DMIG_ind[i+1]):
                solid_DOF = int(rows[j][2]) - 1
                solid_node = int(rows[j][1]) - 1
                used_solid_node = solid_node
                row = used_solid_node*6 + solid_DOF
                col = used_fluid_node*6 + fluid_DOF
                AGG[row,col] = float(rows[j][3])
        else:
            for j in range(entry+1, len(rows)):
                solid_DOF = int(rows[j][2]) - 1
                solid_node = int(rows[j][1]) - 1
                used_solid_node = solid_node
                row = used_solid_node * 6 + solid_DOF
                col = used_fluid_node * 6 + fluid_DOF
                AGG[row,col] = float(rows[j][3])
    outfile = open("AGG", 'wb')
    pickle.dump(AGG, outfile)
    outfile.close()
    return AGG.tocsc()

def Mf_extract(Kf_file, structural_nodes, fluid_nodes, corresp_fluid_nodes):
    print("Extracting the fluid mass matrix")
    f = open(Kf_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the punch file and only keep the mass ones
    row2 = []
    first_encounter = False
    line_start = 0
    line_break = 0
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG" and first_encounter == False:
            if row[1] == "MAAX":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG" and first_encounter == True:
            if row[1] != "MAAX":
                line_break = i
                break
    rows = rows[line_start:line_break]

    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    #The fluid elements IDs are at least len(struct)
    nelm = len(structural_nodes)

    # Create the fluid mass matrix
    #nrows = nddl fluid = number of nodes
    #ncols = same thing
    #Mf matrix is symmetric!
    nelm_fluid = len(fluid_nodes)

    Mf = sp.lil_array((6*len(corresp_fluid_nodes), 6*len(corresp_fluid_nodes)), dtype="float32")

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        fluid_node = int(rows[entry][2]) - nelm - 1
        if fluid_node >= 0:
            col = 6*fluid_node
            if i != len(DMIG_ind) - 1:
                for j in range(entry + 1, DMIG_ind[i+1]):
                    fluid_node = int(rows[j][1]) - nelm - 1
                    row = 6*fluid_node
                    Mf[row, col] = float(rows[j][3])
                    Mf[col, row] = float(rows[j][3])
            else:
                for j in range(entry + 1, len(rows)):
                    fluid_node = int(rows[j][1]) - nelm - 1
                    row = 6*fluid_node
                    Mf[row, col] = float(rows[j][3])
                    Mf[col, row] = float(rows[j][3])
    outfile = open("Mf", 'wb')
    pickle.dump(Mf, outfile)
    outfile.close()
    return Mf.tocsc()

def Kf_extract(Kf_file, structural_nodes, fluid_nodes, corresp_fluid_nodes):
    print("Extracting the fluid rigidity matrix")
    f = open(Kf_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the punch file and only keep the mass ones
    row2 = []
    first_encounter = False
    line_start = 0
    line_break = 0
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG" and first_encounter == False:
            if row[1] == "KAAX":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG" and first_encounter == True:
            if row[1] != "KAAX":
                line_break = i
                break
    rows = rows[line_start:line_break]

    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    #The fluid elements IDs are at least len(struct)
    nelm = len(structural_nodes)

    # Create the fluid mass matrix
    #nrows = nddl fluid = number of nodes
    #ncols = same thing
    #Mf matrix is symmetric!
    nelm_fluid = len(fluid_nodes)

    Kf = sp.lil_array((6*len(corresp_fluid_nodes), 6*len(corresp_fluid_nodes)), dtype="float32")

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        fluid_node = int(rows[entry][2]) - nelm - 1
        if fluid_node >= 0:
            col = 6*fluid_node
            if i != len(DMIG_ind) - 1:
                for j in range(entry + 1, DMIG_ind[i+1]):
                    fluid_node = int(rows[j][1]) - nelm - 1
                    row = 6*fluid_node
                    Kf[row, col] = float(rows[j][3])
                    Kf[col, row] = float(rows[j][3])
            else:
                for j in range(entry + 1, len(rows)):
                    fluid_node = int(rows[j][1]) - nelm - 1
                    row = 6*fluid_node
                    Kf[row, col] = float(rows[j][3])
                    Kf[col, row] = float(rows[j][3])
    outfile = open("Kf", 'wb')
    pickle.dump(Kf, outfile)
    outfile.close()
    return Kf.tocsc()

def calc_addedmass(Kf, AGG, corresp_fluid_nodes, corresp_solid_nodes):
    print("Assembling the added mass matrix")
    Kf = Kf.tocsc()
    A_indices = AGG.nonzero()
    AGG = AGG.tocsc()
    n = AGG.shape[0]
    Ma = sp.lil_array((n, n), dtype='float32')

    print("Solving for addedmass matrix")
    invKf = spl.spilu(Kf)
    start = time.time()
    Ma = AGG@invKf.solve(AGG.transpose().todense())
    print("Time taken = "+str(time.time()-start))
    Ma = sp.csc_matrix(Ma)
    outfile = open("Ma", 'wb')
    pickle.dump(Ma, outfile)
    outfile.close()
    return Ma

def addedmass_calcwrite(file, Kf, AGG): #Obsolete
    print("Assembling the added mass matrix")
    Kf = Kf.tocsc()
    A_indices = AGG.nonzero()
    AGG = AGG.tocsc()
    n = AGG.shape[0]
    Ma = sp.lil_array((n, n), dtype='float32')

    print("Solving for addedmass matrix")
    invKf = spl.spilu(Kf)
    print("Matrix inversed")
    Ma = AGG@invKf.solve(AGG.transpose().todense())
    """
    ms = range(len(A_indices[0]))
    addedmass = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        #for m in range(len(A_indices[0])):
            #print("Pourcentage calculé: " + str(m * 100 / len(A_indices[0])) + '%')
        results = [executor.submit(AKfAT, (AGG, A_indices, Kf, m)) for m in ms]
            #addedmass.append(executor.submit(AKfAT, (AGG, A_indices, Kf, m)))
        
        running = True
        while running == True:
            running = False
            done = 0
            print("Calculations done = "+str(done*100/len(addedmass))+'%')
            for action in addedmass:
                if action.done() == False:
                    running = True
                else:
                    done += 1
        print("Running done")
        
        print("Running")
        concurrent.futures.wait(results)
        print("Parallel computations done")
        for f in concurrent.futures.as_completed(results):
            #print("Pourcentage assigné: " + str(m * 100 / len(A_indices[0])) + '%')
            #print("Pourcentage assigné: " + str(i * 100 / len(results)) + '%')
            print(results[f].result())
            #Ma[:, f.result[1]] = f.result[0]
    """
    print("Done")
    outfile = open("Ma", 'wb')
    pickle.dump(Ma, outfile)
    outfile.close()
    indices = np.nonzero(Ma)

    print("Writing the added mass matrix")
    file.write("$*\n$* Coupled mass matrix input\n$*\n")
    file.write("DMIG, MAAX, 0, 6, 1, 0, , , " + str(n//6) + '\n')

    for i in range(len(indices[0])):
        row = indices[0][i]
        col = indices[1][i]
        node_row = (row//6) + 1
        node_col = (col//6) + 1
        DOF_row = (row%6) + 1
        DOF_col = (col%6) + 1
        if row >= col:
            file.write("DMIG, MAAX, " + str(node_col) + ', ' + str(DOF_col) + ', , '+str(node_row)+', '+str(DOF_row)+', '+\
                       str(format(Ma[row,col],'.8E'))+'\n')

def write_addedmass(file, Ma, corresp_nodes):
    print("Indexing added mass matrix")
    indices1 = np.nonzero(Ma)
    used_indices = np.nonzero(indices1[0] <= indices1[1])
    indices = np.array([indices1[0][used_indices[0]], indices1[1][used_indices[0]]])
    n = Ma.shape[0]

    #row = [str(index) for index in indices[0]]
    #col = [str(index) for index in indices[1]]
    #addedmass = [str(mass) for mass in Ma[indices]]

    addedmass = Ma[indices[0,:],indices[1,:]]

    del Ma
    del used_indices
    del indices1

    print("Writing the added mass matrix")
    file.write("$*\n$* Coupled mass matrix input\n$*\n")
    file.write("DMIG, MAAX, 0, 6, 1, 0, , , " + str(n // 6) + '\n')
    #write_DMIG(row, col, addedmass, corresp_nodes)
    """
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(write_DMIG, (file, Ma, corresp_nodes, indices, i)) for i in range(len(indices[0]))]

        for f in concurrent.futures.as_completed(results):
            print(f.result()*100/len(indices[0]))"""
    #Idée: Utiliser numba en batch de 1000?
    DMIG = ['$\n']
    """
    j=0
    for i in range(1000, len(indices[0]), 10000):
        node_row, node_col, DOF_row, DOF_col, masses = matrix_indexing(indices, addedmass, i, j, corresp_nodes)
        node_row = [str(index) for index in node_row]
        node_col = [str(index) for index in node_col]
        DOF_row = [str(index) for index in DOF_row]
        DOF_col = [str(index) for index in DOF_col]
        # node_row = corresp_nodes[node_row]
        # node_col = corresp_nodes[node_col]
        masses = [str(format(mass,'.8E')) for mass in np.asarray(addedmass[0, j:i]).flatten()]
        DMIG = write_DMIG(node_row, node_col, DOF_row, DOF_col, masses, corresp_nodes, DMIG)
        j = i
        print(j*100/len(indices[0]))
    i = len(indices[0])
    node_row, node_col, DOF_row, DOF_col, masses = matrix_indexing(indices, addedmass, i, j, corresp_nodes)
    node_row = [str(index) for index in node_row]
    node_col = [str(index) for index in node_col]
    DOF_row = [str(index) for index in DOF_row]
    DOF_col = [str(index) for index in DOF_col]
    # node_row = corresp_nodes[node_row]
    # node_col = corresp_nodes[node_col]
    masses = [str(format(mass, '.8E')) for mass in np.asarray(addedmass[0, j:i]).flatten()]
    DMIG = write_DMIG(node_row, node_col, DOF_row, DOF_col, masses, corresp_nodes, DMIG)
    j = i
    print(j * 100 / len(indices[0]))
    print(DMIG)
    exit()
    file.write(DMIG+'$*\n')

    """
    for i in range(min(len(indices[0]), 50)):
        #node_row = (row // 6) + 1
        #node_row = corresp_nodes[node_row]
        #node_col = (col // 6) + 1
        #node_col = corresp_nodes[node_col]
        #DOF_row = (row % 6) + 1
        #DOF_col = (col % 6) + 1
        #mass = addedmass[0,i]
        start = time.time()
        node_row, node_col, DOF_row, DOF_col, mass = matrix_indexing(indices, addedmass, i, corresp_nodes)
        DMIG.append("DMIG, MAAX, " + str(node_col) + ', ' + str(DOF_col) + ', , ' + str(node_row) + ', ' + str(DOF_row) + ', ' + str(format(mass, '.8E')) + '\n')
        print(str(i*100/len(indices[0]))+"%")

    file.write("".join(DMIG)+'$*\n')

    print("Done")

def matrix_indexing(indices, addedmass, index, corresp_nodes):
    node_row = (index // 6 + 1)
    node_col = (index // 6 + 1)
    DOF_row = (index % 6 + 1)
    DOF_col = (index % 6 + 1)
    # node_row = corresp_nodes[node_row]
    # node_col = corresp_nodes[node_col]
    masses = addedmass[0, index]
    return node_row, node_col, DOF_row, DOF_col, masses

def write_DMIG(node_rows, node_cols, DOF_rows, DOF_cols, addedmass, corresp_nodes, DMIG):
    for i in range(len(node_rows)):
        node_row = node_rows[i]
        node_col = node_cols[i]
        DOF_row = DOF_rows[i]
        DOF_col = DOF_cols[i]
        Ma = addedmass[i]
        DMIG.append("DMIG, MAAX, "+node_col+', '+DOF_col+', , '+node_row+', '+DOF_row+', '+Ma+'\n')
        """file.write("DMIG, MAAX, " + str(node_col) + ', ' + str(DOF_col) + ', , ' + str(node_row) + ', ' + str(
            DOF_row) + ', ' + \
                    str(format(Ma[row, col], '.8E')) + '\n')"""
    return DMIG

def AKfAT(AGG, A_indices, Kf, m): #Obsolete
    #print("Pourcentage fait: " + str(m * 100 / len(A_indices[0])) + '%')
    col = A_indices[0][m]
    #col = A_indices[0][m]
    #p = invKf.solve(AGG.transpose()[:, [col]].todense())
    b = AGG.T[:,[col]].todense()
    p = spl.minres(Kf, b)
    #p = np.linalg.solve(Kf, AGG.T[:, col])
    addedmass = AGG@p[0]
    return addedmass, col

def mass_extract(M2GG_file, npoints):
    f = open(M2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the punch file and only keep the mass ones
    row2 = []
    first_encounter = False
    line_start = 0
    line_break = 0
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG" and first_encounter == False:
            if row[1] == "MAAX":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG" and first_encounter == True:
            if row[1] != "MAAX":
                line_break = i
                break
    rows = rows[line_start:line_break]

    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    # Create the mass matrix
    lenM2GG = npoints

    Ms = sp.lil_array((lenM2GG*6, lenM2GG*6), dtype="float32")

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        col_node = int(rows[entry][2])
        col_DOF = int(rows[entry][3])

        if i != len(DMIG_ind) - 1:
            for j in range(entry + 1, DMIG_ind[i + 1]):
                row_DOF = int(rows[j][2])
                if row_DOF == 0:
                    row_DOF = 1
                row_node = int(rows[j][1])
                row = (row_node - 1) * 6 + row_DOF - 1
                col = (col_node - 1) * 6 + col_DOF - 1
                Ms[row, col] = float(rows[j][3])
                Ms[col, row] = float(rows[j][3])
        else:
            for j in range(entry + 1, len(rows)):
                row_DOF = int(rows[j][2])
                if row_DOF == 0:
                    row_DOF = 1
                row_node = int(rows[j][1])
                row = (row_node - 1) * 6 + row_DOF - 1
                col = (col_node - 1) * 6 + col_DOF - 1
                Ms[row, col] = float(rows[j][3])
                Ms[col, row] = float(rows[j][3])
    outfile = open("Ms", 'wb')
    pickle.dump(Ms, outfile)
    outfile.close()
    return Ms.tocsc()

def rigidity_extract(K2GG_file, npoints):
    f = open(K2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the punch file and only keep the mass ones
    row2 = []
    first_encounter = False
    line_start = 0
    line_break = 0
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG" and first_encounter == False:
            if row[1] == "KAAX":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG" and first_encounter == True:
            if row[1] != "KAAX":
                line_break = i
                break
    rows = rows[line_start:line_break]

    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    # Create the mass matrix
    lenM2GG = npoints

    Ks = sp.lil_array((lenM2GG*6, lenM2GG*6), dtype="float32")

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        col_node = int(rows[entry][2])
        col_DOF = int(rows[entry][3])

        if i != len(DMIG_ind) - 1:
            for j in range(entry + 1, DMIG_ind[i + 1]):
                row_DOF = int(rows[j][2])
                if row_DOF == 0:
                    row_DOF = 1
                row_node = int(rows[j][1])
                row = (row_node - 1) * 6 + row_DOF - 1
                col = (col_node - 1) * 6 + col_DOF - 1
                Ks[row, col] = float(rows[j][3])
                Ks[col, row] = float(rows[j][3])
        else:
            for j in range(entry + 1, len(rows)):
                row_DOF = int(rows[j][2])
                if row_DOF == 0:
                    row_DOF = 1
                row_node = int(rows[j][1])
                row = (row_node - 1) * 6 + row_DOF - 1
                col = (col_node - 1) * 6 + col_DOF - 1
                Ks[row, col] = float(rows[j][3])
                Ks[col, row] = float(rows[j][3])
    outfile = open("Ks", 'wb')
    pickle.dump(Ks, outfile)
    outfile.close()
    return Ks.tocsc()

def damping_extract(B2GG_file, npoints):
    print("Extracting the fluid rigidity matrix")
    f = open(B2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the punch file and only keep the mass ones
    row2 = []
    first_encounter = False
    line_start =0
    line_break = 0
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG" and first_encounter == False:
            if row[1] == "BAAX":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG" and first_encounter == True:
            if row[1] != "BAAX":
                line_break = i
                break
    rows = rows[line_start:line_break]

    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    #The fluid elements IDs are at least len(struct)
    nelm = npoints

    Bs = sp.lil_array((6*nelm, 6*nelm), dtype="float32")

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        col_node = int(rows[entry][2])
        col_DOF = int(rows[entry][3])

        if i != len(DMIG_ind) - 1:
            for j in range(entry + 1, DMIG_ind[i + 1]):
                row_DOF = int(rows[j][2])
                if row_DOF == 0:
                    row_DOF = 1
                row_node = int(rows[j][1])
                row = (row_node - 1) * 6 + row_DOF - 1
                col = (col_node - 1) * 6 + col_DOF - 1
                Bs[row, col] = float(rows[j][3])
                Bs[col, row] = float(rows[j][3])
        else:
            for j in range(entry + 1, len(rows)):
                row_DOF = int(rows[j][2])
                if row_DOF == 0:
                    row_DOF = 1
                row_node = int(rows[j][1])
                row = (row_node - 1) * 6 + row_DOF - 1
                col = (col_node - 1) * 6 + col_DOF - 1
                Bs[row, col] = float(rows[j][3])
                Bs[col, row] = float(rows[j][3])
    outfile = open("Kf", 'wb')
    pickle.dump(Bs, outfile)
    outfile.close()
    return Bs.tocsc()

def damping_extract2(B2GG_file, npoints):
    f = open(B2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the rigidity file and only keep the rigidity ones
    row2 = []
    for i in range(1, len(rows)):
        row = rows[i].split(' ')
        rowa = rows[i-1].split(' ')
        row = [item for item in row if item != '']
        rowa = [item for item in rowa if item != '']
        if row[0] == "DMIG*" or rowa[0] == "DMIG*":
            if row[1] == "BAAX" or rowa[1] == "BAAX":
                row2.append(row)
    rows = row2

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    # Create the mass matrix
    lenM2GG = npoints

    Bs = sp.lil_array((lenM2GG*6, lenM2GG*6), dtype="float32")

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        col_node = int(rows[entry][2]) - 1
        col_DOF = int(rows[entry][3]) - 1

        if i != len(DMIG_ind) - 1:
            for j in range(entry + 1, DMIG_ind[i + 1]):
                row_DOF = int(rows[j][2]) - 1
                row_node = int(rows[j][1]) - 1
                row = row_node * 6 + row_DOF
                col = col_node * 6 + col_DOF
                Bs[row, col] = float(rows[j][3])
                Bs[col, row] = float(rows[j][3])
        else:
            for j in range(entry + 1, len(rows)):
                row_DOF = int(rows[j][2]) - 1
                row_node = int(rows[j][1]) - 1
                row = row_node * 6 + row_DOF
                col = col_node * 6 + col_DOF
                Bs[row, col] = float(rows[j][3])
                Bs[col, row] = float(rows[j][3])
    outfile = open("Bs", 'wb')
    pickle.dump(Bs, outfile)
    outfile.close()
    return Bs.tocsc()

def M2GG_extract(M2GG_file, structural_nodes):
    f = open(M2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the mass file and only keep the mass ones
    row2 = []
    for i in range(1, len(rows)):
        row = rows[i].split(' ')
        rowa = rows[i-1].split(' ')
        row = [item for item in row if item != '']
        rowa = [item for item in rowa if item != '']
        if row[0] == "DMIG*" or rowa[0] == "DMIG*":
            if row[1] == "MAAX" or rowa[1] == "MAAX":
                row2.append(row)
    rows = row2

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    # Create the mass matrix
    lenM2GG = len(structural_nodes)

    M2GG = np.zeros([3, lenM2GG, lenM2GG], dtype='float16')

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        node = int(rows[entry][2]) - 1
        if i != len(DMIG_ind) - 1:
            for j in range(entry + 1, DMIG_ind[i + 1]):
                DOF = int(rows[j][2]) - 1
                M2GG[DOF, node, node] = float(rows[j][3])
        else:
            for j in range(entry + 1, len(rows)):
                DOF = int(rows[j][2]) - 1
                M2GG[DOF, node, node] = float(rows[j][3])
    return M2GG

def K2GG_extract(K2GG_file, structural_nodes):
    f = open(K2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the mass file and only keep the rigidity ones
    row2 = []
    for i in range(1, len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG*":
            if row[1] == "MAAX":
                line_break = i
                break
    rows = rows[:line_break-1]
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        row2.append(row)

    #Divide negative rows in two sections
    rows = row2
    for i in range(len(rows)):
        if rows[i][0] != "DMIG*":
            if len(rows[i]) == 3:
                last_data = [char for char in rows[i][2]]
                del rows[i][2]
                rows[i].append(last_data[0])
                rows[i].append("".join(last_data[1:]))

    # Delete the \n at the end of the rows
    for i in range(len(rows)):
        last_data = [char for char in rows[i][-1]]
        del rows[i][-1]
        rows[i].append("".join(last_data[:-1]))

    # Find all the rows with DMIG cards
    DMIG_ind = []
    for i in range(len(rows)):
        if rows[i][0] == "DMIG*":
            DMIG_ind.append(i)

    # Create the mass matrix
    lenK2GG = len(structural_nodes)

    K2GG = np.zeros([3, 3, lenK2GG, lenK2GG], dtype='float16')

    # Get the components in the matrix
    for i in range(len(DMIG_ind)):
        entry = DMIG_ind[i]
        col = int(rows[entry][2]) - 1
        DOF1 = int(rows[entry][3]) - 1
        if i != len(DMIG_ind) - 1:
            for j in range(entry + 1, DMIG_ind[i + 1]):
                row = int(rows[j][1]) - 1
                DOF2 = int(rows[j][2]) - 1
                K2GG[DOF1, DOF2, row, col] = float(rows[j][3])
        else:
            for j in range(entry + 1, len(rows)):
                row = int(rows[j][1]) - 1
                DOF2 = int(rows[j][2]) - 1
                K2GG[DOF1, DOF2, row, col] = float(rows[j][3])
    return K2GG

def write_P1(file, P1):
    indices = np.nonzero(P1)
    used_indices = np.nonzero(indices[0] <= indices[1])
    indices = np.array([indices[0][used_indices[0]], indices[1][used_indices[0]]])
    file.write("$*\n$* Coupled mass matrix input\n$*\n")
    file.write("DMI, P1P, 0, 3, 1, 0, , " + str(int(P1.shape[0])) + '\n')
    for i in range(len(indices[0])):
        col = indices[0][i]
        row = indices[1][i]
        file.write(
            "DMI, P1P, " + str(col) + ', ' + str(row) + ', ' + str(P1[row, col]) + '\n')
    file.write('$*\n')

def write_mass(file, M2GG):
    indices = np.nonzero(M2GG)
    used_indices = np.nonzero(indices[0] <= indices[1])
    indices = np.array([indices[0][used_indices[0]], indices[1][used_indices[0]]])
    file.write("$*\n$* Coupled mass matrix input\n$*\n")
    file.write("DMIG, MAJ, 0, 6, 1, 0, , , "+str(int(M2GG.shape[0]/6))+'\n')
    for i in range(len(indices[0])):
        col = indices[0][i]
        row = indices[1][i]
        node_row = (row // 6) + 1
        node_col = (col // 6) + 1
        DOF_row = (row % 6) + 1
        DOF_col = (col % 6) + 1
        file.write("DMIG, MAJ, " + str(node_col) + ', ' + str(DOF_col) + ', , '+str(node_row)+', '+str(DOF_row)+', '+\
                   str(M2GG[indices[0][i],indices[1][i]])+'\n')
    """
    for j in range(len(M2GG[0,:,:])):
        for k in range(3):
            if np.any(M2GG[k,:,j] != 0):
                file.write("DMIG, M2GG, "+str(j+1)+', '+str(k+1)+',\n')
                for i in range(j, len(M2GG[k,:,:])):
                    if M2GG[k,i,j] != 0.0:
                        file.write(", "+str(i+1)+', '+str(k+1)+', '+str(round(M2GG[k,i,j],4))+',\n')"""
    file.write('$*\n')

def write_rigidity(file, K2GG):
    indices = np.nonzero(K2GG)
    used_indices = np.nonzero(indices[0] <= indices[1])
    indices = np.array([indices[0][used_indices[0]], indices[1][used_indices[0]]])
    file.write("$*\n$* Rigidity matrix input\n$*\n")
    file.write("DMIG, KAJ, 0, 6, 1, 0, , , "+str(int(K2GG.shape[0]/6))+'\n')

    for i in range(len(indices[0])):
        col = indices[0][i]
        row = indices[1][i]
        node_row = (row // 6) + 1
        node_col = (col // 6) + 1
        DOF_row = (row % 6) + 1
        DOF_col = (col % 6) + 1
        file.write("DMIG, KAJ, " + str(node_col) + ', ' + str(DOF_col) + ', , '+str(node_row)+', '+str(DOF_row)+', '+\
                   str(K2GG[indices[0][i],indices[1][i]])+'\n')
    file.write('$*\n')

def write_CDAMP(file, B2GG, nelm):
    file.write('$* CDAMP2 entries to represent fluid structure coupling\n')
    indices = np.nonzero(B2GG)
    used_indices = np.nonzero(indices[0] <= indices[1])
    indices = np.array([indices[0][used_indices[0]], indices[1][used_indices[0]]])
    for i in range(len(indices[0])):
        col = indices[0][i]
        row = indices[1][i]
        node_row = (row // 6) + 1
        node_col = (col // 6) + 1
        DOF_row = (row % 6) + 1
        DOF_col = (col % 6) + 1
        file.write("CDAMP2, " +str(nelm+i)+', ' +str(format(B2GG[indices[0][i],indices[1][i]],'.6E'))+', '+\
                    str(node_col) + ', ' + str(DOF_col) + ', '+str(node_row)+', '+str(DOF_row)+'\n')
    file.write('$*')

def write_damping(file, B2GG):
    indices = np.nonzero(B2GG)
    used_indices = np.nonzero(indices[0] <= indices[1])
    indices = np.array([indices[0][used_indices[0]], indices[1][used_indices[0]]])
    file.write("$*\n$* Rigidity matrix input\n$*\n")
    file.write("DMIG, BAJ, 0, 6, 1, 0, , , "+str(int(B2GG.shape[0]/6))+'\n')
    for i in range(len(indices[0])):
        col = indices[0][i]
        row = indices[1][i]
        node_row = row // 6 + 1
        node_col = (col // 6) + 1
        DOF_row = (row % 6) + 1
        DOF_col = (col % 6) + 1
        file.write("DMIG, BAJ, " + str(node_col) + ', ' + str(DOF_col) + ', , '+str(node_row)+', '+str(DOF_row)+', '+\
                   str(B2GG[indices[0][i],indices[1][i]])+'\n')
    file.write('$*\n')

def M2GG_write(file, M2GG_file):
    f = open(M2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the mass file and only keep the mass ones
    row2 = []
    first_encounter = False
    for i in range(1, len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG*" and first_encounter == False:
            if row[1] == "MAAX":
                line_start = i
                first_encounter = True
        if row[0] == "DMIG*":
            if row[1] == "VAX":
                line_break = i
                break
    rows = rows[line_start-1:line_break - 1]
    rows = "".join(rows)
    file.write("$ Structural mass matrix\n")
    file.write(rows)
    file.write("$\n")

def K2GG_write(file, K2GG_file):
    f = open(K2GG_file, 'r')
    rows = f.readlines()
    f.close()

    # Get all the rows from the mass file and only keep the rigidity ones
    row2 = []
    for i in range(1, len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == "DMIG*":
            if row[1] == "MAAX":
                line_break = i
                break
    rows = rows[:line_break - 1]
    rows = "".join(rows)
    file.write("$ Structural rigidity matrix\n")
    file.write(rows)
    file.write("$\n")

def CMASS(file, AGG, nodelist):
    file.write('$* CMASS2 and CONM1 entries to represent Fluid-Structure Coupling\n')
    nelm = len(nodelist.simplices)+len(nodelist.aero_simplices)
    indices = np.nonzero(AGG)
    M = np.zeros([len(AGG[0,:,:]),3])
    Gs = []
    for i in range(len(indices[0])):
        DOF = indices[0][i]
        G1 = indices[1][i]
        C1 = 1
        G2 = indices[2][i]
        M = -AGG[DOF,G1,G2]
        #M[G2,DOF] += -AGG[DOF,G1,G2]
        #Gs.append(G2+1)
        if G1 == G2:
            file.write('CMASS2, ' + str(i + nelm) + ', ' + str(M) + ', ' + str(G1 + 1) + ', '+ str(DOF+1) + '\n')
        elif G1>G2:
            file.write('CMASS2, ' + str(i + nelm) + ', ' + str(M) + ', ' + str(G1 + 1) + ', ' + str(DOF + 1) + ', ' + str(
                G2 + 1) + ', ' + str(DOF + 1) + '\n')
        """
        for i in range(len(M)):
        file.write('CONM1, ' + str(i + nelm) + ', ' + str(Gs[i]) + ', , ' + str(M[i, 0]) + ', , ' + str(
            M[i, 1]) + ', , \n, ' + str(M[i, 2]) + '\n')
        if G1 == G2:
            if DOF == 0:
                file.write('CONM1, ' + str(i + nelm) + ', ' + str(G1 + 1) + ', , ' +str(M) + '\n')
            elif DOF == 1:
                file.write('CONM1, ' + str(i + nelm) + ', ' + str(G1 + 1) + ', , , , ' + str(M) + '\n')
            elif DOF == 2:
                file.write('CONM1, ' + str(i + nelm) + ', ' + str(G1 + 1) + ', , , , , , \n, ' + str(M) + '\n')"""
    file.write('$*\n')

if __name__ == '__main__':
    damping = pickle.load(open("DAMP", "rb"))
    f = open("essai.dat", 'w')
    write_damping(file = f, B2GG = damping)
    newDamp = damping_extract(B2GG_file = 'essai.dat', npoints = 19253)
    plt.figure("DesiredDamping")
    plt.spy(damping)
    plt.title("Desired damping matrix")
    plt.figure("ObtainedDamping")
    plt.spy(newDamp)
    plt.title("Obtained damping matrix")

    plt.show()