#Defining the constraints and loads in the BDF file
#Author: Danick Lamoureux
#Project under Frédérick Gosselin and Sébastien Houde's supervision
#Date: 2022-05-05

import numpy as np

def SPC(file, nodelist, profile_nodes_len):
    """See documentation for SPC1. To be modified to have one free end

    Args:
        file: File object obtained with file = open(string, 'w')
        nodelist (list): List of nodes
        profile_nodes_len (int): Number of nodes for the profile of the hydrofoils
    """
    file.write('$* Constraints and loads\n')
    file.write('$* SPCs\n$*\n')
    for i in range(profile_nodes_len):
        file.write('SPC1, 100, 123456, '+str(i+1)+'\n')
        file.write('SPC1, 100, 123456, ' + str(len(nodelist) - i) + '\n') # comment this line if you want the hydrofoil to be free at one side

def ACMODL(file):
    """See documentation for ACMODL

    Args:
        file: File object obtained with file = open(string, 'w')
    """
    file.write("$*\n")
    file.write('ACMODL, , , , , 1.0, , , REL,\n')
    file.write(", 0.5, 0, , , STRONG\n")

def GRAV(file):
    """See documentation for GRAV. Used for debugging the structural part of this program

    Args:
        file: File object obtained with file = open(string, 'w')
    """
    file.write('GRAV, 100, , 9810.00, 0.0000, 0.0000, -1.0000\n')

def AERO(file, ref_length = 1, rhoref = 1.225):
    """See documentation for AERO

    Args:
        file: File object obtained with file = open(string, 'w')
        ref_length (float, optional): Reference length of the hydrofoil. Defaults to 1.
        rhoref (float, optional): Reference density. Defaults to 1.225.
    """
    file.write('AERO, , , '+str(ref_length)+', '+str(rhoref)+',\n$*\n')

def MKAERO(file, SID, mach_matrix, freq_matrix, velocities, densities, machs):
    """See documentation for MKAERO1, FLFACT, FLUTTER

    Args:
        file: File object obtained with file = open(string, 'w')
        SID: Identifier
        mach_matrix: Mach matrix for interpolation
        freq_matrix: Reduced frequency matrix for interpolation
        velocities: Velocities to test
        densities: Densities to test
        machs: Machs to test
    """
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
    file.write('$ Velocities\n')
    file.write('FLFACT, 30')
    skip = 0
    for i in range(len(velocities)):
        if i == 7:
            file.write("\n")
            skip = i
        elif i - skip == 8 and i > 8:
            file.write(", \n")
            skip = i
        file.write(', ' + str(format(velocities[i], '.2E')))

    file.write('\n')

    file.write('$*\n')

    #Flutter card
    file.write('FLUTTER, '+str(SID)+', PK, 10, 20, 30\n$*\n')