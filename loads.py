#Defining the constraints and loads in the BDF file
#Author: Danick Lamoureux
#LM2 project under Frédérick Gosselin's supervision
#Date: 2022-05-05

#It is assumed that the nodes at the root and tip are fixed, and that's it

import numpy as np

def SPC(file, nodelist, profile_nodes_len):
    #See documentation for SPC1
    file.write('$* Constraints and loads\n')
    file.write('$* SPCs\n$*\n')
    fixed_nodes = []
    for i in range(profile_nodes_len):
        file.write('SPC1, 100, 123456, '+str(i+1)+'\n')
        file.write('SPC1, 100, 123456, ' + str(len(nodelist) - i) + '\n')

def ACMODL(file):
    #See documentation for ACMODL
    file.write("$*\n")
    file.write('ACMODL, , , , , 1.0, , , REL,\n')
    file.write(", 0.5, 0, , , STRONG\n")

def AML(file, nodelist): #Currently useless but can be useful for direct response analysis
    #See documentation for the automatically matched layer in the acoustics documentation
    #Applies as a farfield condition
    file.write("$*\n")
    file.write("AMLREG, 1, 1, AML REGION 1,\n")
    file.write(", 5, AML, 0, 0, 0\n")


    file.write('BSURFS, 1, , , , ')

def GRAV(file):
    #See documentation for GRAV
    #Used for debugging the structural part of this program
    file.write('GRAV, 100, , 9810.00, 0.0000, 0.0000, -1.0000\n')

def DOWNWASH(file, flow_angle, nodes_object): #To complete
    #TO CONTINUE ONCE NEEDED when the straight wing supposition is reversed
    #DMI CARDs
    flow_angle = flow_angle*np.pi/180
    #Kutta's condition states that the flow leaving a panel is parallel to its trailing edge
    alpha = flow_angle
    it = np.zeros(len(nodes_object.aero_simplices))
    downwash = np.zeros(len(nodes_object.aero_simplices))
    w = np.zeros(len(nodes_object.aero_simplices))
    for i in range(len(nodes_object.aero_simplices)):
        node1 = nodes_object.aero_simplices[i,0]
        node2 = nodes_object.aero_simplces[i,1]

    file.write('DMI, W2GJ, 0, 2, 1, 0, '+str(len(nodes_object.aero_simplices))+', 1\n')
    for i in range(len(nodes_object.aero_simplices)):
        file.write('DMI, W2GJ, 1, '+str(i)+', '+str(w[i])+'\n')

def AERO(file, nodes_object, velocity = 100, ref_length = 1, rhoref = 1.225):
    #See documentation for AERO
    file.write('AERO, , , '+str(ref_length/2)+', '+str(rhoref)+',\n$*\n')

def MKAERO(file, SID, mach_matrix, freq_matrix, kfreq, densities, machs, velocity = None):
    #See documentation for MKAERO1, FLFACT, FLUTTER
    #This is not the most optimal way of writing it in terms of space but it works
    #For each mach number, write every reduced frequency
    file.write('$*')
    for i in range(len(mach_matrix)):
        file.write('\nMKAERO1, '+str(format(mach_matrix[i],".8E"))+'\n,')
        for j in range(len(freq_matrix)):
            if j%7 ==0 and j!=0:
                file.write('\nMKAERO1, ' + str(format(mach_matrix[i],".8E")) + '\n,')
            file.write(str(round(freq_matrix[j],4))+", ")
    file.write('\n$*\n')

    #Writing the FLFACT cards for values to test
    file.write('$ FLFACT and FLUTTER Cards\n')
    file.write('$ Density\n')
    file.write('FLFACT, 10')
    for i in range(len(densities)):
        file.write(', '+str(densities[i]))
    file.write('\n')

    file.write('$ Mach numbers\n')
    file.write('FLFACT, 20')
    for i in range(len(machs)):
        file.write(', '+str(format(machs[i],".8E")))
    file.write('\n')

    #If a velocity is negative, it is studied in depth
    #Therefore, if there is a negative velocity to check, add it to the test list
    if velocity != None:
        kfreq = np.append(kfreq,abs(velocity))
        kfreq = np.unique(kfreq)
    file.write('$ Velocities\n')
    file.write('FLFACT, 30')
    skip = 0
    for i in range(len(kfreq)):
        if i == 7:
            file.write("\n")
            skip = i
        elif i - skip == 8 and i > 8:
            file.write(", \n")
            skip = i
        if kfreq[i] == velocity:
            file.write(', ' + str(format(-kfreq[i], '.2E')))
        else:
            file.write(', ' + str(format(kfreq[i], '.2E')))

    file.write('\n')

    file.write('$*\n')

    #Flutter card
    file.write('FLUTTER, '+str(SID)+', PK, 10, 20, 30\n$*\n')