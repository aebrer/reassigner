#!/usr/bin/env python3

import numpy as np
from urllib.request import urlopen
from Bio.PDB.PDBParser import PDBParser
import argparse
import os
import sys
import shutil
import glob
import re

#silence warnings from numpy when doing gap_checks (and all others but it's all dealt with internally [hopefully])
old_settings = np.seterr(all='ignore') 


# Oregon State 2014

# Andrew E. Brereton

# In collaboration with:

# Dmitriy Bobrovnikov
# Christopher King
# Robby Blizzard
# Dr. Victor Hsu


#checks if there is a valid file at a specified location
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return (open(arg, 'r'))  # return an open file handle
# importing arguments from the user
parser=argparse.ArgumentParser(
    description='''Assigns secondary structure without using hydrogen bonds. One method uses virtual dihedrals and bond angles, the other using phi/psi2 motifs to identify secondary structure. ''',
    epilog="""IMPORTANT NOTES:\nExpect there to be strange ? residues at the end of the output if there are any ligands or other non-water HETATOMS present. This can be safely ignored."""
    )
parser.add_argument('-i',dest="input",metavar="FILE",type=lambda x: is_valid_file(parser, x), help='This should be a pdb file. Do not use in combination with the "-c" option.')
parser.add_argument('-c',dest="code",type=str, help='This should be a four letter pdb code. Only use this option if you want to download directly from the PDB.')
parser.add_argument('-o',dest="output",metavar="FILE", help='Output file, a table using tabs as seperators. If not specified, default is to output to STDOUT as human readable output.')
parser.add_argument('--legend',dest='legend', action='store_true', help='Option to print a legend for the Secondary Structure codes.')
parser.add_argument('--verbose',dest='verbose', action='store_true', help='Option to print a all the output, including the behind the scenes methods for structure assignment.')
parser.set_defaults(legend=False, verbose=False)
args=parser.parse_args()
if args.legend == True:
    print ("\nLEGEND:\n_:\tBreak in the chain (as detected by checking bond distance)\n-:\tUnassigned trans-residue\n=:\tUnassigned cis-residue\nP:\tPII-Helix\nt:\tTurn defined by CA to CA Distance\nN:\tNon-Hydrogen Bonded Turn (Similar to 's')\nT:\tHydrogen Bonded Turn T\nE:\tExtended or Beta-strand conformation\nH:\tAlpha Helical Conformation\nG:\t3,10 Helix\nBb:\tBeta-bulge\nU:\tPi-helical bulge\nX:\tThe Stig\n")
# Dictionary and function that converts triple letter codes to single letter
# Any non-conventional three letter code becomes "?"
# Victor made the first version of this function
def to_single(single,triple):

    return (single[triple])
one_letter = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', \
          'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'BLA':'-'}
#function for getting pdb files online
def fetch_pdb(id):
    pdb = "%s.pdb" % str(id.lower())
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id.lower()
    tmp = urlopen(url)
    fh = open(pdb,'wb')
    fh.write(tmp.read())
    fh.close()
if args.code == None:
    try:
        pdb = args.input
        print ("\nWorking...\n")
    except: 
        pass

elif args.input == None:
    try:
        fetch_pdb(args.code)
        pdb = open('%s.pdb' %args.code.lower(), 'r')
        print ("\nWorking...\n")
    except:
        print ("\n\n\nPlease enter a valid code or path.\n\n\n")

else:
    print ("\n\n\nPlease only choose option -i or option -c.\n\n\n")





# function to get indicies, matches the index from resnum list with the index for the atomic coords,  which is why the resnum has to be in each
def index_getter(resnum):
    indices = []    
    i = 0
    for atom in atom_list:
        try:
            index = atom.index(resnum)
            if index == 0:
                indices.append(i)
        except:
            pass
        i += 1
        
    return (indices)

# checks for certain atom types and grabs just those coordinates

#this is for correctly sorting atom types
def atom_getter(index,atom):
    if atom_list[index][1] == atom:
        return (atom_xyz[index])
    else:
        return ('no')


#### GAP CHECK FUNCTION ######

def gap_check(resnum):

    indices = index_getter(resnum)
    prev_indices = []
    next_indices = []
    for i in indices:
        prev_indices.append(i-4)
        next_indices.append(i+4)


    atom_types = ['C','N']

    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                iC = atom_getter(i,atom)
            elif atom == 'N':
                iN = atom_getter(i,atom)
    for i in prev_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                prevC = atom_getter(i,atom)
    for i in next_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                nextN = atom_getter(i,atom)
    try:
        ahead = np.subtract(nextN, iC)
        behind = np.subtract(iN, prevC)
        
        ahead_mag = np.sqrt(ahead.dot(ahead))
        behind_mag = np.sqrt(behind.dot(behind))

        if ((ahead_mag > 1.5) and (behind_mag > 1.5)):
            return('isolated')
        elif(ahead_mag > 1.5):
            return('ahead')
        elif(behind_mag > 1.5):
            return('behind')
        else:
            return('no')
    except:
        return("Fatal Error in Gap Check")


#### GAP CHECK FUNCTION for dison3, discn3, discaca3, and dison4 ######

def long_gap_check(resnum):

    indices = index_getter(resnum)
    next_indices = []
    plus_two_indices = []
    plus_three_indices =[]
    plus_four_indices =[]

    for i in indices:
        plus_two_indices.append(i+8)
        next_indices.append(i+4)
        plus_three_indices.append(i+12)
        plus_four_indices.append(i+16)



    atom_types = ['C','N']

    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                iC = atom_getter(i,atom)
    for i in plus_two_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                plus_twoC = atom_getter(i,atom)
            elif atom == 'N':
                plus_twoN = atom_getter(i,atom)
    for i in next_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                nextC = atom_getter(i,atom)
            elif atom == 'N':
                nextN = atom_getter(i,atom)
    for i in plus_three_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                plus_threeC = atom_getter(i,atom)
            elif atom == 'N':
                plus_threeN = atom_getter(i,atom)
    for i in plus_four_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                plus_fourN = atom_getter(i,atom)
    try:
        ahead = np.subtract(nextN, iC)
        two_ahead = np.subtract(plus_twoN, nextC)
        three_ahead = np.subtract(plus_threeN, plus_twoC)
        four_ahead = np.subtract(plus_fourN, plus_threeC)
        
        ahead_mag = np.sqrt(ahead.dot(ahead))
        two_ahead_mag = np.sqrt(two_ahead.dot(ahead))
        three_ahead_mag = np.sqrt(three_ahead.dot(ahead))
        four_ahead_mag = np.sqrt(four_ahead.dot(ahead))


        if ((ahead_mag > 1.5) or (two_ahead_mag > 1.5) or (three_ahead_mag > 1.5)):
            return('threegap')
        elif ((ahead_mag > 1.5) or (two_ahead_mag > 1.5) or (three_ahead_mag > 1.5) or (four_ahead_mag > 1.5)):
            return('fourgap')
        else:
            return('no')
    except:
        return("Fatal Error in Gap Check")

# ZETA FUNCTION
# returns a single value, zeta for the resnum entered
# There are functions within functions here: I'm sorry. I was new to python
def zeta_calc(resnum):
    # The order of atoms in atom_list is N, CA, C, O
    
    #this returns the index values of each of the atoms at this residue number
    

    indices = index_getter(resnum)
    prev_indices = []
    for i in indices:
        prev_indices.append(i-4)
    if ((gap_check(resnum) == 'behind') or (gap_check(resnum) == 'isolated')):
        return('ERROR')

    atom_types = ['C','O']
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                iC = atom_getter(i,atom)
            elif atom == 'O':
                iO = atom_getter(i,atom)
                
    for i in prev_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                prevC = atom_getter(i,atom)
            elif atom == 'O':
                prevO = atom_getter(i,atom)
    
    
    
    v1 = np.subtract(iC, iO)
    v2 = np.subtract(prevC, iC)
    v3 = np.subtract(prevO, prevC)
    n1 = np.cross(v1,v2)
    n2 = np.cross(v2,v3)

    dot1 = np.dot(n1,n2)
    n1mag = np.sqrt(n1.dot(n1))
    n2mag = np.sqrt(n2.dot(n2))
    cos_zeta = dot1/(n1mag*n2mag)

    zeta = np.degrees(np.arccos(cos_zeta))


    #testing for direction
    cross = np.cross(n1,n2)
    direction = np.dot(cross, v2)
    if direction < 0:
        zeta = -1 * zeta
    
    return (zeta)

# Omega FUNCTION
# returns a single value, ome for the resnum entered
# There are functions within functions here: I'm sorry. I was new to python
def ome_calc(resnum):
    # The order of atoms in atom_list is N, CA, C, O
    


    indices = index_getter(resnum)
    prev_indices = []
    for i in indices:
        prev_indices.append(i-4)
    if ((gap_check(resnum) == 'behind') or (gap_check(resnum) == 'isolated')):
        return('ERROR')

    atom_types = ['C','N','CA']
    

    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                iN = atom_getter(i,atom)
            elif atom == 'CA':
                iCA = atom_getter(i,atom)
                
    for i in prev_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                prevC = atom_getter(i,atom)
            elif atom == 'CA':
                prevCA = atom_getter(i,atom)
    
    
    
    v1 = np.subtract(iN, iCA)
    v2 = np.subtract(prevC, iN)
    v3 = np.subtract(prevCA, prevC)
    n1 = np.cross(v1,v2)
    n2 = np.cross(v2,v3)

    dot1 = np.dot(n1,n2)
    n1mag = np.sqrt(n1.dot(n1))
    n2mag = np.sqrt(n2.dot(n2))
    cos_ome = dot1/(n1mag*n2mag)

    ome = np.degrees(np.arccos(cos_ome))


    #testing for direction
    cross = np.cross(n1,n2)
    direction = np.dot(cross, v2)
    if direction < 0:
        ome = -1 * ome
    
    return (ome)
# PHI FUNCTION
# returns a single value, phi for the resnum entered
# There are functions within functions here: I'm sorry. I was new to python
def phi_calc(resnum):
    

    indices = index_getter(resnum)
    prev_indices = []
    for i in indices:
        prev_indices.append(i-4)
    if ((gap_check(resnum) == 'behind') or (gap_check(resnum) == 'isolated')):
        return('ERROR')

    atom_types = ['C','N','CA']
    

    
    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                iC = atom_getter(i,atom)
            elif atom == 'CA':
                iCA = atom_getter(i,atom)
            elif atom == 'N':
                iN = atom_getter(i,atom)
                
    for i in prev_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                prevC = atom_getter(i,atom)
    
    
    
    v1 = np.subtract(iCA, iC)
    v2 = np.subtract(iN, iCA)
    v3 = np.subtract(prevC, iN)
    n1 = np.cross(v1,v2)
    n2 = np.cross(v2,v3)

    dot1 = np.dot(n1,n2)
    n1mag = np.sqrt(n1.dot(n1))
    n2mag = np.sqrt(n2.dot(n2))
    cos_phi = dot1/(n1mag*n2mag)

    phi = np.degrees(np.arccos(cos_phi))


    #testing for direction
    cross = np.cross(n1,n2)
    direction = np.dot(cross, v2)
    if direction < 0:
        phi = -1 * phi
    
    return (phi)
# PSI FUNCTION
# returns a single value, phi for the resnum entered
# There are functions within functions here: I'm sorry. I was new to python
def psi_calc(resnum):
    
    

    indices = index_getter(resnum)
    next_indices = []
    for i in indices:
        next_indices.append(i+4)
    if ((gap_check(resnum) == 'ahead') or (gap_check(resnum) == 'isolated')):
        return('ERROR')

    atom_types = ['C','N','CA']
    

    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                iC = atom_getter(i,atom)
            elif atom == 'CA':
                iCA = atom_getter(i,atom)
            elif atom == 'N':
                iN = atom_getter(i,atom)
                
    for i in next_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                nextN = atom_getter(i,atom)
    
    
    
    v1 = np.subtract(iC, nextN)
    v2 = np.subtract(iCA, iC)
    v3 = np.subtract(iN, iCA)
    n1 = np.cross(v1,v2)
    n2 = np.cross(v2,v3)

    dot1 = np.dot(n1,n2)
    n1mag = np.sqrt(n1.dot(n1))
    n2mag = np.sqrt(n2.dot(n2))
    cos_psi = dot1/(n1mag*n2mag)

    psi = np.degrees(np.arccos(cos_psi))


    #testing for direction
    cross = np.cross(n1,n2)
    direction = np.dot(cross, v2)
    if direction < 0:
        psi = -1 * psi
    
    return (psi)
def tau_calc(resnum):
    # The order of atoms in atom_list is N, CA, C, O
    

    indices = index_getter(resnum)
    prev_indices = []
    next_indices = []
    for i in indices:
        prev_indices.append(i-4)
        next_indices.append(i+4)
    if ((gap_check(resnum) == 'behind') or (gap_check(resnum) == 'isolated') or (gap_check(resnum) == 'ahead')):
        return('ERROR')
      
    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in prev_indices:
        if atom_getter(i,'CA') == 'no':
            pass
        else:
            prev_CA = atom_getter(i,'CA')
    
    for i in indices:
        if atom_getter(i,'CA') == 'no':
            pass
        else:
            iCA = atom_getter(i,'CA')
    
    for i in next_indices:
        if atom_getter(i,'CA') == 'no':
            pass
        else:
            next_CA = atom_getter(i,'CA')

    
    
    v1 = np.subtract(prev_CA, iCA)
    v2 = np.subtract(next_CA, iCA)
    d1 = np.dot(v1,v2)

    v1mag = np.sqrt(v1.dot(v1))
    v2mag = np.sqrt(v2.dot(v2))
    cos_tau = d1/(v1mag*v2mag)

    tau = np.degrees(np.arccos(cos_tau))

    return (tau)
def dison3_calc(resnum):

        
    indices = index_getter(resnum)
    atom_types = ['N','O']
    
    plus_three_indices = []
    for i in indices:
        plus_three_indices.append(i+12)
    if ((gap_check(resnum) == 'ahead') or (gap_check(resnum) == 'isolated')):
        return('ERROR')
    if ((long_gap_check(resnum) == 'threegap')):
        return('ERROR')
    
    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'O':
                iO = atom_getter(i,atom)
                
    for i in plus_three_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                plus_three_N = atom_getter(i,atom)

    
    dison3 = np.sqrt( np.sum( ( np.array(plus_three_N) - np.array(iO) ) **2 ) )
    
    return (dison3)
def dison4_calc(resnum):
    # The order of atoms in atom_list is N, CA, C, O
    #this returns the index values of each of the atoms at this residue number
    

        
    indices = index_getter(resnum)
    plus_four_indices = []
    for i in indices:
        plus_four_indices.append(i+16)
    if ((gap_check(resnum) == 'ahead') or (gap_check(resnum) == 'isolated')):
        return('ERROR')
    if ((long_gap_check(resnum) == 'fourgap') or (long_gap_check(resnum) == 'threegap')):
        return('ERROR')

    atom_types = ['N','O']
    
    
    
    
    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'O':
                iO = atom_getter(i,atom)
                
    for i in plus_four_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                plus_four_N = atom_getter(i,atom)

    
    dison4 = np.sqrt( np.sum( ( np.array(plus_four_N) - np.array(iO) ) **2 ) )
    
    return (dison4)
def discn3_calc(resnum):
    # The order of atoms in atom_list is N, CA, C, O
    #this returns the index values of each of the atoms at this residue number
    

        
    indices = index_getter(resnum)

    plus_three_indices = []
    for i in indices:
        plus_three_indices.append(i+12)
    if ((gap_check(resnum) == 'ahead') or (gap_check(resnum) == 'isolated')):
        return('ERROR')
    if ((long_gap_check(resnum) == 'threegap')):
        return('ERROR')


    atom_types = ['N','C']
    

    
    
    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    for i in indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'C':
                iC = atom_getter(i,atom)
                
    for i in plus_three_indices:
        for atom in atom_types:
            if atom_getter(i,atom) == 'no':
                pass
            elif atom == 'N':
                plus_three_N = atom_getter(i,atom)

    
    discn3 = np.sqrt( np.sum( ( np.array(plus_three_N) - np.array(iC) ) **2 ) )
    
    return (discn3)
def disnn1_calc(resnum):


    indices = index_getter(resnum)

    next_indices = []
    for i in indices:
        next_indices.append(i+4)
    if ((gap_check(resnum) == 'ahead') or (gap_check(resnum) == 'isolated')):
        return('ERROR')
    
    
    
    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    
    for i in indices:
        if atom_getter(i,'N') == 'no':
            pass
        else:
            iN = atom_getter(i,'N')
    
    for i in next_indices:
        if atom_getter(i,'N') == 'no':
            pass
        else:
            next_N = atom_getter(i,'N')

    
    
    disnn1 = np.sqrt( np.sum( ( np.array(next_N) - np.array(iN) ) **2 ) )

    return (disnn1)
def discaca3_calc(resnum):

    indices = index_getter(resnum)

    plus_three_indices = []
    
    for i in indices:
        plus_three_indices.append(i+12)
    
    if ((gap_check(resnum) == 'ahead') or (gap_check(resnum) == 'isolated')):
        return('ERROR')
    if ((long_gap_check(resnum) == 'threegap')):
        return('ERROR')
    

    
    # gets coords for each atom and creates a variable for each, this makes atom assignment pdbfile order independant
    
    for i in indices:
        if atom_getter(i,'CA') == 'no':
            pass
        else:
            iCA = atom_getter(i,'CA')
    
    for i in plus_three_indices:
        if atom_getter(i,'CA') == 'no':
            pass
        else:
            next_CA = atom_getter(i,'CA')

    
    
    discaca3 = np.sqrt( np.sum( ( np.array(next_CA) - np.array(iCA) ) **2 ) )

    return (discaca3)
###############################################################################################
######### Calculators above this line   #######################################################
###############################################################################################

# helixa, helix10, and helix1 used to contain work by Chris, but needed to be entirely redone
def helixa_test(count):
    if zeta_dict[resnum_count_dict[count]] > -10 \
    and zeta_dict[resnum_count_dict[count]] < 45 \
    and dison4_dict[resnum_count_dict[count-1]] < 3.5 \
    and (\
    (dison4_dict[resnum_count_dict[count-1]] < (dison3_dict[resnum_count_dict[count-1]] + 0.1))\
    or (\
    (dison4_dict[resnum_count_dict[count]] < (dison3_dict[resnum_count_dict[count]] + 0.1))\
    ) \
    and dison4_dict[resnum_count_dict[count+1]] < (dison3_dict[resnum_count_dict[count+1]] + 0.1)\
    ):
        return (1)
def helixa_calc(count):

    if helixa_test(count) == 1 \
    and (\
    (\
    (helixa_test(count-2) == 1) and (helixa_test(count-1) == 1) and (helixa_test(count-3) == 1)\
    ) or (\
    (helixa_test(count+2) == 1) and (helixa_test(count+1) == 1) and (helixa_test(count+3) == 1)\
    ) or (\
    (helixa_test(count-1) == 1) and (helixa_test(count+1) == 1) and (helixa_test(count+2) == 1)\
    ) or (\
    (helixa_test(count-1) == 1) and (helixa_test(count-2) == 1) and (helixa_test(count+1) == 1)\
    )\
    ) :
        return (1)
    
    else:
        return (0)
def helix10_test(count):
    if zeta_dict[resnum_count_dict[count]] > -10 \
    and zeta_dict[resnum_count_dict[count]] < 45 \
    and dison3_dict[resnum_count_dict[count-1]] < 3.5 \
    and (\
    (dison3_dict[resnum_count_dict[count-1]] < (dison4_dict[resnum_count_dict[count-1]] + 0.1)) \
    or (dison3_dict[resnum_count_dict[count]] < (dison4_dict[resnum_count_dict[count]] + 0.1)) \
    ):
        return (1)                    
def helix10_calc(count):

    if helix10_test(count) == 1 \
    and (\
    (\
    (helix10_test(count-2) == 1) and (helix10_test(count-1) == 1) and (helix10_test(count-3) == 1)\
    ) or (\
    (helix10_test(count+2) == 1) and (helix10_test(count+1) == 1) and (helix10_test(count+3) == 1)\
    ) or (\
    (helix10_test(count-1) == 1) and (helix10_test(count+1) == 1) and (helix10_test(count+2) == 1)\
    ) or (\
    (helix10_test(count-1) == 1) and (helix10_test(count-2) == 1) and (helix10_test(count+1) == 1)\
    )\
    ) :
        return (1)
    else:
        return (0)
def helix1_test_1(count):
    if (zeta_dict[resnum_count_dict[count]] > -140 and zeta_dict[resnum_count_dict[count]] < -90) \
    and (tau_dict[resnum_count_dict[count]] > 105 and tau_dict[resnum_count_dict[count]] < 135):
        return (1)
    else:
        return(0)
def helix1_test_2(count):
    if (helix1_test_1(count-1) == 1) \
    and (\
    (zeta_dict[resnum_count_dict[count]] > -150 and zeta_dict[resnum_count_dict[count]] < -80) \
    and (tau_dict[resnum_count_dict[count]] > 100 and tau_dict[resnum_count_dict[count]] < 140)\
    ):
        return (1)
    else:
        return(0)
def helix1_test_3(count):
    if (\
    (helix1_test_1(count-1) == 1) or (helix1_test_2(count-1) == 1)\
    ) and (\
    (helix1_test_1(count+1) == 1 or (helix1_test_2(count+1) == 1)\
    ) and (\
    (zeta_dict[resnum_count_dict[count]] < 0) and (tau_dict[resnum_count_dict[count]] > 100 and tau_dict[resnum_count_dict[count]] < 140)\
    )\
    ):
        return (1)
    else:
        return(0)
def helix1_short_calc(count):
    # for short PP helicies (ie without checking length) return 1
    if ((helix1_test_1(count) == 1) or (helix1_test_2(count) == 1) or (helix1_test_3(count) == 1)):
        return (1)
    else:
        return (0)
def helix1_long_calc(count):
    if helix1_short_calc(count) == 1 \
    and (\
    (\
    (helix1_short_calc(count-2) == 1) and (helix1_short_calc(count-1) == 1) and (helix1_short_calc(count-3) == 1)\
    ) or (\
    (helix1_short_calc(count+2) == 1) and (helix1_short_calc(count+1) == 1) and (helix1_short_calc(count+3) == 1)\
    ) or (\
    (helix1_short_calc(count-1) == 1) and (helix1_short_calc(count+1) == 1) and (helix1_short_calc(count+2) == 1)\
    ) or (\
    (helix1_short_calc(count-1) == 1) and (helix1_short_calc(count-2) == 1) and (helix1_short_calc(count+1) == 1)\
    )\
    ) :
        return (1)
    
    else:
        return (0)
def betan_test(count):
    if discn3_dict[resnum_count_dict[count]] < 4.7:
        return (1)
def betat_test(count):
    if dison3_dict[resnum_count_dict[count]] < 3.5:
        return (1)
def sheet_test(count):
    if ((zeta_dict[resnum_count_dict[count]] < -140.0 ) or (zeta_dict[resnum_count_dict[count]] > 162.0)):
        return (1)
    else:
        return (0)
def sheet_calc(count):

    if sheet_test(count) == 1 \
    and (\
    (\
    (sheet_test(count-2) == 1) and (sheet_test(count-1) == 1)\
    or \
    (sheet_test(count+2) == 1) and (sheet_test(count+1) == 1)\
    or \
    (sheet_test(count-1) == 1) and (sheet_test(count+1) == 1)\
    )\
    ) :
        return (1)
    
    else:
        return (0)
def bulge_calc(count):
    
    if (\
    (disnn1_dict[resnum_count_dict[count]] < 3.01)\
    and (disnn1_dict[resnum_count_dict[count+1]] > 3.01)\
    and (disnn1_dict[resnum_count_dict[count-1]] > 3.01)\
    ) and (\
    ( (zeta_dict[resnum_count_dict[count]] > -55.0)\
    and (zeta_dict[resnum_count_dict[count]] < 19.0)\
    ) and ( (tau_dict[resnum_count_dict[count]] > 87.0)\
    and (tau_dict[resnum_count_dict[count]] < 120.0)\
    ) ) and (\
    ( (sheet_test(count-1) == 1) and (sheet_test(count-2) == 1) )\
    or ( (sheet_test(count-2) == 1) and (sheet_test(count-3) == 1) )\
    or ( (sheet_test(count+1) == 1) and (sheet_test(count+2) == 1) )\
    or ( (sheet_test(count+2) == 1) and (sheet_test(count+3) == 1) )\
    ) and (\
    (dison3_dict[resnum_count_dict[count-2]] > 3.1)\
    ):
    
        return (1)
    
    else:
        
        return (0)
###############################################################################################
######### Above this line goes all the functions for assigning the secondary structure ########
########################### based on zeta and such, below is based on phipsi2 #################
###############################################################################################



### Just a function I made to check if a value is within a range, in order to save having to type it out each time
def range_check(bound,target,check):
    target_up = target + bound
    target_down = target - bound
    if (check < target_up and check > target_down):
        return(1)
    else:
        return(0)

def pp_helixa_calc(count):
    #range = 10 degrees
    if (\
    #alpha helix condition one - alpha.alpha.1        
    (\
    (phi_dict[resnum_count_dict[count]] < -53 and phi_dict[resnum_count_dict[count]] > -73)\
    and (psi_dict[resnum_count_dict[count]] < -32 and psi_dict[resnum_count_dict[count]] > -52)\
    and (phi_dict[resnum_count_dict[count+1]] < -53 and phi_dict[resnum_count_dict[count+1]] >  -73)\
    and (psi_dict[resnum_count_dict[count+1]] < -32 and psi_dict[resnum_count_dict[count+1]] > -52)\
    )\
    ):
        return(4)
    elif (\
    #alpha helix condition two - P.alpha.1        
    (\
    (phi_dict[resnum_count_dict[count]] < -58 and phi_dict[resnum_count_dict[count]] > -78)\
    and (psi_dict[resnum_count_dict[count]] < 169 and psi_dict[resnum_count_dict[count]] > 149)\
    and (phi_dict[resnum_count_dict[count+1]] < -47 and phi_dict[resnum_count_dict[count+1]] >  -67)\
    and (psi_dict[resnum_count_dict[count+1]] < -30 and psi_dict[resnum_count_dict[count+1]] > -50)\
    )or \
    #alpha helix condition three - beta.alpha.1        
    (\
    (phi_dict[resnum_count_dict[count]] < -126 and phi_dict[resnum_count_dict[count]] > -146)\
    and (psi_dict[resnum_count_dict[count]] < 171 and psi_dict[resnum_count_dict[count]] > 151)\
    and (phi_dict[resnum_count_dict[count+1]] < -47 and phi_dict[resnum_count_dict[count+1]] >  -67)\
    and (psi_dict[resnum_count_dict[count+1]] < -31 and psi_dict[resnum_count_dict[count+1]] > -51)\
    )or \
    #alpha helix condition four - beta.alpha.3        
    (\
    (phi_dict[resnum_count_dict[count]] < -150 and phi_dict[resnum_count_dict[count]] > -170)\
    and (psi_dict[resnum_count_dict[count]] < 176 and psi_dict[resnum_count_dict[count]] > 156)\
    and (phi_dict[resnum_count_dict[count+1]] < -52 and phi_dict[resnum_count_dict[count+1]] >  -72)\
    and (psi_dict[resnum_count_dict[count+1]] < -27 and psi_dict[resnum_count_dict[count+1]] > -47)\
    ) or \
    #alpha helix condition five - beta.alpha.4        
    (\
    (phi_dict[resnum_count_dict[count]] < -105 and phi_dict[resnum_count_dict[count]] > -125)\
    and (psi_dict[resnum_count_dict[count]] < 111 and psi_dict[resnum_count_dict[count]] > 91)\
    and (phi_dict[resnum_count_dict[count+1]] < -49 and phi_dict[resnum_count_dict[count+1]] >  -69)\
    and (psi_dict[resnum_count_dict[count+1]] < -26 and psi_dict[resnum_count_dict[count+1]] > -46)\
    )\
    ):
        return(1)
    else:
        return(0)

def pp_sheet_calc(count):
    # returns a different number depending on how many neibghoring residues are E
        #this is the first case in which I implemented the number system so they 
        #don't really make sense individually. In pp_helix_calc() return(4) is for 
        #the case with four consecutive helix residues
    # range = 10 degrees unless otherwise specified
    if (\
    #beta.beta.1 range = 15
    (\
    (phi_dict[resnum_count_dict[count]] < -100 and phi_dict[resnum_count_dict[count]] > -130)\
    and (psi_dict[resnum_count_dict[count]] < 143 and psi_dict[resnum_count_dict[count]] > 113)\
    and (phi_dict[resnum_count_dict[count+1]] < -103 and phi_dict[resnum_count_dict[count+1]] >  -133)\
    and (psi_dict[resnum_count_dict[count+1]] < 141 and psi_dict[resnum_count_dict[count+1]] > 111)\
    ) or \
    #beta.beta.2 range = 15
    (\
    (phi_dict[resnum_count_dict[count]] < -132 and phi_dict[resnum_count_dict[count]] > -162)\
    and (psi_dict[resnum_count_dict[count]] < 174 and psi_dict[resnum_count_dict[count]] > 144)\
    and (phi_dict[resnum_count_dict[count+1]] < -132 and phi_dict[resnum_count_dict[count+1]] >  -162)\
    and (psi_dict[resnum_count_dict[count+1]] < 171 and psi_dict[resnum_count_dict[count+1]] > 141)\
    )or \
    #P.beta.2 range = 15      
    (\
    (phi_dict[resnum_count_dict[count]] < -64 and phi_dict[resnum_count_dict[count]] > -94)\
    and (psi_dict[resnum_count_dict[count]] < 150 and psi_dict[resnum_count_dict[count]] > 120)\
    and (phi_dict[resnum_count_dict[count+1]] < -118 and phi_dict[resnum_count_dict[count+1]] >  -148)\
    and (psi_dict[resnum_count_dict[count+1]] < 175 and psi_dict[resnum_count_dict[count+1]] > 145)\
    )\
    ):
        return(1)
    elif (\
    #beta.P.2
    (\
    (phi_dict[resnum_count_dict[count]] < -108 and phi_dict[resnum_count_dict[count]] > -128)\
    and (psi_dict[resnum_count_dict[count]] < 135 and psi_dict[resnum_count_dict[count]] > 115)\
    and (phi_dict[resnum_count_dict[count+1]] < -56 and phi_dict[resnum_count_dict[count+1]] >  -76)\
    and (psi_dict[resnum_count_dict[count+1]] < 143 and psi_dict[resnum_count_dict[count+1]] > 123)\
    )\
    ):
        return(3)
    elif (\
    #beta.P.3        
    (\
    (phi_dict[resnum_count_dict[count]] < -80 and phi_dict[resnum_count_dict[count]] > -100)\
    and (psi_dict[resnum_count_dict[count]] < 133 and psi_dict[resnum_count_dict[count]] > 113)\
    and (phi_dict[resnum_count_dict[count+1]] < -57 and phi_dict[resnum_count_dict[count+1]] >  -77)\
    and (psi_dict[resnum_count_dict[count+1]] < 148 and psi_dict[resnum_count_dict[count+1]] > 128)\
    ) or \
    #P.beta.1        
    (\
    (phi_dict[resnum_count_dict[count]] < -53 and phi_dict[resnum_count_dict[count]] > -73)\
    and (psi_dict[resnum_count_dict[count]] < 155 and psi_dict[resnum_count_dict[count]] > 135)\
    and (phi_dict[resnum_count_dict[count+1]] < -126 and phi_dict[resnum_count_dict[count+1]] >  -146)\
    and (psi_dict[resnum_count_dict[count+1]] < 165 and psi_dict[resnum_count_dict[count+1]] > 145)\
    )\
    ):
        return(4)
    elif (\
    #P.beta.3      
    (\
    (phi_dict[resnum_count_dict[count]] < -56 and phi_dict[resnum_count_dict[count]] > -76)\
    and (psi_dict[resnum_count_dict[count]] < 154 and psi_dict[resnum_count_dict[count]] > 134)\
    and (phi_dict[resnum_count_dict[count+1]] < -108 and phi_dict[resnum_count_dict[count+1]] >  -128)\
    and (psi_dict[resnum_count_dict[count+1]] < 138 and psi_dict[resnum_count_dict[count+1]] > 118)\
    )):
        return(5)
    elif (\
    #epsilon.beta.1      
    (\
    (phi_dict[resnum_count_dict[count]] < 178 and phi_dict[resnum_count_dict[count]] > 158)\
    and (psi_dict[resnum_count_dict[count]] < -142 and psi_dict[resnum_count_dict[count]] > -162)\
    and (phi_dict[resnum_count_dict[count+1]] < -125 and phi_dict[resnum_count_dict[count+1]] >  -145)\
    and (psi_dict[resnum_count_dict[count+1]] < 161 and psi_dict[resnum_count_dict[count+1]] > 141)\
    )\
    ):
        return(2)
    else:
        return(0)

def pp_bulge_calc(count):
    # returns a different number depending on how many neibghoring residues are E
    if (\
    (\
    (range_check(10, -153,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 173, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 55, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 48, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -103,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 136, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 60, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 36, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(1)
    elif (\
    (\
    (range_check(10, -89,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -29, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -147, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 151, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -68,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 133, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -117, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -14, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -75,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 124, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -94, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -35, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(2)
    elif (\
    (\
    (range_check(10, -115,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -16, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -155, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 147, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -77,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -43, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -161, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 158, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(3)
    elif (\
    (\
    (range_check(10, 80,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 4, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -114, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 140, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, 77,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 8, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -20, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 154, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(4)
    else:
        return(0)

def pp_turn_calc(count):
    # returns a different number depending on how many neibghoring residues are E
    if (\
    (\
    (range_check(10, -64,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -23, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -87, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -2, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(1)
    elif (\
    (\
    (range_check(10, -55,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 133, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 87, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -7, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -67,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 159, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 61, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 35, psi_dict[resnum_count_dict[count+1]]) == 1)
    )
    ):
        return(2)
    elif (\
    (\
    (range_check(10, -64,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -40, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -78, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 125, psi_dict[resnum_count_dict[count+1]]) == 1)\
    ) or (\
    (range_check(10, -65,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -41, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -78, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 153, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -80,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -22, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -159, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 165, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -62,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -42, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -115, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 131, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(3)
    elif (\
    (\
    (range_check(10, -86,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -2, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 75, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 22, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(4)
    
    elif (\
    (\
    (range_check(10, 53,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 37, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 81, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 2, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, 59,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -132, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -93, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 8, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(5)   
    elif (\
    (\
    (range_check(10, 91,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -9, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -70, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 145, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, 80,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 4, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -114, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 140, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, 77,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 8, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -120, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 154, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, 72,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 22, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -105, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 151, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(6)
    elif (\
    (\
    (range_check(10, -129,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 124, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 51, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 41, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -134,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 113, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 59, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -127, psi_dict[resnum_count_dict[count+1]]) == 1)
    ) or (\
    (range_check(10, -127,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, 99, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, 59, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -124, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(7)        
    else:
        return(0)
def pp_pi_calc(count):
    # returns a different number depending on how many neibghoring residues are E
    if (\
    (\
    (range_check(10, -102,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -24, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -116, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, -59, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(1)
    else:
        return(0)
def pp_stig_calc(count):
    # returns a different number depending on how many neibghoring residues are E
    if (\
    (\
    (range_check(10, 75,phi_dict[resnum_count_dict[count]]) == 1) \
    and (range_check(10, -173, psi_dict[resnum_count_dict[count]]) == 1)\
    and (range_check(10, -63, phi_dict[resnum_count_dict[count+1]]) == 1)\
    and (range_check(10, 143, psi_dict[resnum_count_dict[count+1]]) == 1)
    )\
    ):
        return(1)
    else:
        return(0)
        
def cis_calc(count):
    # checks for cis peptides
    if (range_check(40, 0, ome_dict[resnum_count_dict[count]]) == 1):
        return(1)
    else:
        return(0)






# These three lists contain everything needed from PDB file for all calculations
# don't touch them
# maybe this doesn't make sense anymore?
# they are WAY required though, so don't touch them

atom_list = []
atom_xyz = []
resnum_list_unsort = []
pdb_reader = PDBParser(PERMISSIVE = 1, QUIET = True)
struc = pdb_reader.get_structure("temp", pdb)
model = struc[0]


########## NOTE TO THE READER! ############
########## CAUTION CAUTION CATUION ########
# There is a loop below, "for chain in model:"
# all the code needs to be in that loop. This
# is the crazy solution we came up with that 
# allows the program to work on each chain.
# Judge not lest ye be judged.


# Reads in the atom coordinates to make these lists
for chain in model:
    #### EVERYTHING MUST BE WITHIN THIS CHAIN LOOP
    for residue in chain:
        for atom in residue:
            if atom.get_name() == 'CA' or atom.get_name() == 'C' or atom.get_name() == 'N' or atom.get_name() == 'O' and residue.get_resname() != 'HOH':
                    atom_list.append([str(residue.get_id()[1])+str(residue.get_id()[2]), atom.get_id(), residue.get_resname(), chain.get_id() ])
                    atom_xyz.append(atom.get_coord())
                    resnum_list_unsort.append(float(str(residue.get_id()[1])+"."+str(ord(str(residue.get_id()[2])))))    

    pdb.close()



    #this makes the resnum_list unique, and deals correctly with insertion codes
    #ie. it gets unique values, sorts them numerically, which works because codes (A, B, C, etc) get converted to ascii values ie. 101A becomes 101.52
    # then after sorting it undoes this ascii conversion and regenerates the string
    res_set = set(resnum_list_unsort)
    resnum_list_float = list(res_set)
    resnum_list = []
    for resnum in resnum_list_float:
        num = re.sub(r"([0-9]+)(\.)([0-9]+)", r"\1", str(resnum))
        ascii = re.sub(r"([0-9]+)(\.)([0-9]+)", r"\3", str(resnum))
        if int(ascii) < 0:
            ascii = int(ascii) * -1
        insert = chr(int(ascii))
        resnum_list.append(str(num)+str(insert))
                
    # this is necessicary so that we can call the next or previous residue by using count +1 or count -1, but we need to pass the actual resnum in resnum_list to the calculators
    resnum_count = range(1,len(resnum_list) + 1)
    count = 1
    resnum_count_dict = {}
    for resnum in resnum_list:
        resnum_count_dict[count] = resnum
        count = count + 1

    ###############################################################################################
    ######### Below this line is all the looping and such that calls the calulators then ##########
    ######### assigns secondary structure. ########################################################
    ###############################################################################################

    ### loop that makes a dictionary of residue numbers and their zeta values

    zeta_dict = {}
    tau_dict = {}
    dison3_dict = {}
    dison4_dict = {}
    discn3_dict = {}
    disnn1_dict = {}
    discaca3_dict = {}
    phi_dict = {}
    psi_dict = {}
    ome_dict = {}

    for resnum in resnum_list:
        
        #using try means that we will not be thrown off by gaps
        try:
            phi_dict[resnum] = phi_calc(resnum)
        except:
            pass            #print "Phi not calculated for residue number %i\n" % resnum
        
        try:
            psi_dict[resnum] = psi_calc(resnum)
        except:
            pass            #print "psi not calculated for residue number %i\n" % resnum

        try:
            zeta_dict[resnum] = zeta_calc(resnum)
        except:
            pass            #print "Zeta not calculated for residue number %i\n" % resnum
        try:
            tau_dict[resnum] = tau_calc(resnum)
        except:
            pass            #print "Tau not calculated for residue number %i\n" % resnum
        
        try:
            dison3_dict[resnum] = dison3_calc(resnum)
        except:
            pass            #print "Dison3 not calculated for residue number %i\n" % resnum
        
        try:
            dison4_dict[resnum] = dison4_calc(resnum)
        except:
            pass            #print "Dison4 not calculated for residue number %i\n" % resnum
        
        try:
            discn3_dict[resnum] = discn3_calc(resnum)
        except:
            pass            #print "Discn3 not calculated for residue number %i\n" % resnum
        
        try:
            disnn1_dict[resnum] = disnn1_calc(resnum)
        except:
            pass

                
        try:
            discaca3_dict[resnum] = discaca3_calc(resnum)
        except:
            pass
        
        try:
            ome_dict[resnum] = ome_calc(resnum)
        except:
            pass
           

    ###############################################################################################
    ######### Above this line is all the looping and such that calls the calulators ###############
    ################## Below is the acutal assignment #############################################
    ###############################################################################################

    ## dictionary that contains all secondary structure assignments
    ss_dict = {}
    pp_dict = {}
    final_dict = {}

    ##dictionary for residue types based on residue number
    res_type_dict = {}
    chain_id_dict = {}
    aa_dict = {}

    def atom_get(i,atom):
        if atom_list[i][1] == atom:
            return (atom_list[i][2])
        else:
            return ('no')
    def chain_get(i,atom):
        if atom_list[i][1] == atom:
            return (atom_list[i][3])
        else:
            return ('no') 

    for resnum in resnum_list:
        indices = index_getter(resnum)    
        atom_types = 'CA'
        
        for i in indices:
            if atom_getter(i,'CA') == 'no':
                pass
            else:
                res_type_dict[resnum] = atom_get(i,'CA')
                chain_id_dict[resnum] = chain_get(i,'CA')
        try:
            aa_dict[resnum] = to_single(one_letter,res_type_dict[resnum])
        except:
            aa_dict[resnum] = '?'

    ### setting all the SS to blank
    for resnum in resnum_list:
        ss_dict[resnum] = '-'
        pp_dict[resnum] = '-'
        

    # assigns SS to all residues
    # ensures that there is a minimum of 4 residues in a row for a helix
    #PRIORITY: B, P(4+), H, G, T, E, N, P(2-3)
    for count in resnum_count:
        try:
            if helix1_short_calc(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'P'
        except:
            pass

        try:
            if betan_test(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'N'
        except:
            pass
                
        try:
            if sheet_calc(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'E'
        except:
            pass

        try:
            if betat_test(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'T'
        except:
            pass
            
            
        try:
            if helix10_calc(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'G'
        except:
            pass

        try:
            if helixa_calc(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'H'
        except:
            pass

        #assigns beta-bulges
        try:
            if bulge_calc(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'B'
            if ss_dict[resnum_count_dict[count-1]] == 'B':
                ss_dict[resnum_count_dict[count]] = 'b'
            if ss_dict[resnum_count_dict[count]] == 'B' and sheet_test(count-1) == 1 and ss_dict[resnum_count_dict[count - 1]] == '-':
                ss_dict[resnum_count_dict[count-1]] = 'e'
            if ss_dict[resnum_count_dict[count]] == 'b' and sheet_test(count+1) == 1 and ss_dict[resnum_count_dict[count + 1]] == '-':
                ss_dict[resnum_count_dict[count+1]] = 'e'

        except:
            pass

        try:
            if helix1_long_calc(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'P'
        except:
            pass
        

        ### PhiPsi2 motif assignment    
        ### t turns ### Way way loose, even calls helixes turns
        try:
            if (discaca3_dict[resnum_count_dict[count]] < 6.5 and pp_dict[resnum_count_dict[count]] == '-'):
                pp_dict[resnum_count_dict[count]] = 't'
        except:
            pass
        try:
            if (discaca3_dict[resnum_count_dict[count-3]] < 6.5 and pp_dict[resnum_count_dict[count]] == '-'):
                pp_dict[resnum_count_dict[count]] = 't'
        except:
            pass

        ### PhiPsi2 motif assignment    
        ### This one needs to go first because it's somewhat relaxed compared to the others, and if not first, some turns and bulges will be masked
        try:
            if pp_sheet_calc(count) == 1:
                pp_dict[resnum_count_dict[count]] = 'E'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'E'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'E'
                #except:
                #    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
            elif pp_sheet_calc(count) == 2:
                pp_dict[resnum_count_dict[count]] = 'E'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'E'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
            elif pp_sheet_calc(count) == 3:
                pp_dict[resnum_count_dict[count]] = 'E'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'E'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'E'
                #except:
                #    pass
            elif pp_sheet_calc(count) == 4:
                pp_dict[resnum_count_dict[count]] = 'E'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'E'
                except:
                    pass
            elif pp_sheet_calc(count) == 5:
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'E'
                except:
                    pass
        except:
            pass

       
        ### PhiPsi2 motif assignment    
        try:
            if pp_helixa_calc(count) == 4:
                pp_dict[resnum_count_dict[count]] = 'H'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'H'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'H'
                #except:
                #    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'H'
                #except:
                #    pass
            elif pp_helixa_calc(count) == 1:
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'H'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'H'
                #except:
                #    pass
        except:
            pass



        ### PhiPsi2 motif assignment    
            ### for turns
        try:
            if pp_turn_calc(count) == 1:
                pp_dict[resnum_count_dict[count+1]] = 'T'
            elif pp_turn_calc(count) == 2:
                pp_dict[resnum_count_dict[count]] = 'T'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'T'
                except:
                    pass
            elif pp_turn_calc(count) == 3:
                pp_dict[resnum_count_dict[count]] = 'N'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'N'
                except:
                    pass
            elif pp_turn_calc(count) == 4:
                pp_dict[resnum_count_dict[count]] = 'T'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'T'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'H'
                #except:
                #    pass
            elif pp_turn_calc(count) == 5:
                pp_dict[resnum_count_dict[count]] = 'T'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'T'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'E'
                #except:
                #    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
            elif pp_turn_calc(count) == 6:
                pp_dict[resnum_count_dict[count]] = 'T'
                try:
                    pp_dict[resnum_count_dict[count-1]] = 'T'
                except:
                    pass
            elif pp_turn_calc(count) == 7:
                pp_dict[resnum_count_dict[count]] = 'E'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'T'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'E'
                #except:
                #    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'T'
                #except:
                #    pass
        except:
            pass
        ### PhiPsi2 motif assignment    
            ### for bulges
        try:
            if pp_bulge_calc(count) == 1:
                pp_dict[resnum_count_dict[count]] = 'B'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'b'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'E'
                #except:
                #    pass
            elif pp_bulge_calc(count) == 2:
                pp_dict[resnum_count_dict[count]] = 'B'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'b'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
            elif pp_bulge_calc(count) == 3:
                pp_dict[resnum_count_dict[count]] = 'B'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'b'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'E'
                #except:
                #    pass
            elif pp_bulge_calc(count) == 4:
                pp_dict[resnum_count_dict[count]] = 'B'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'b'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'T'
                #except:
                #    pass
                #try:
                #    pp_dict[resnum_count_dict[count+2]] = 'E'
                #except:
                #    pass
        except:
            pass
        ### PhiPsi2 motif assignment    
            ### for pi-helicies
        try:
            if pp_pi_calc(count) == 1:
                pp_dict[resnum_count_dict[count]] = 'U'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'U'
                except:
                    pass
                #try:
                #    pp_dict[resnum_count_dict[count-1]] = 'H'
                #except:
                #    pass
        except:
            pass      
        ### PhiPsi2 motif assignment    
            ### for the Stig
        try:
            if pp_stig_calc(count) == 1:
                pp_dict[resnum_count_dict[count]] = 'X'
                try:
                    pp_dict[resnum_count_dict[count+1]] = 'X'
                except:
                    pass
        except:
            pass   
        #extends helicies
        try:
            if ss_dict[resnum_count_dict[count-1]] == 'G' and ss_dict[resnum_count_dict[count]] != 'G':
                ss_dict[resnum_count_dict[count]] = 'g'
            elif ss_dict[resnum_count_dict[count-1]] == 'H' and ss_dict[resnum_count_dict[count]] != 'H':
                ss_dict[resnum_count_dict[count]] = 'h'
        except:
            pass

        #extends beta-sheets based on nearby beta-sheets and passing sheet_test, and not being assigned another type of SS
        try:
            if ss_dict[resnum_count_dict[count]] == 'E' and ss_dict[resnum_count_dict[count-1]] == '-' and ss_dict[resnum_count_dict[count+1]] == '-':
                ss_dict[resnum_count_dict[count-1]] = 'e'
                ss_dict[resnum_count_dict[count+1]] = 'e'
            elif ss_dict[resnum_count_dict[count+1]] == 'E' and ss_dict[resnum_count_dict[count+2]] == 'E' and ss_dict[resnum_count_dict[count+3]] == 'E' and sheet_test(count) == 1:
                ss_dict[resnum_count_dict[count]] = 'e'
            elif ss_dict[resnum_count_dict[count-1]] == 'E' and ss_dict[resnum_count_dict[count-2]] == 'E' and ss_dict[resnum_count_dict[count-3]] == 'E' and sheet_test(count) == 1 and ss_dict[resnum_count_dict[count]] == '-':
                ss_dict[resnum_count_dict[count]] = 'e'
            elif ss_dict[resnum_count_dict[count+2]] == 'E' and ss_dict[resnum_count_dict[count+3]] == 'E' and ss_dict[resnum_count_dict[count+4]] == 'E' and sheet_test(count) == 1 and ss_dict[resnum_count_dict[count]] == '-':
                ss_dict[resnum_count_dict[count]] = 'e'
            elif ss_dict[resnum_count_dict[count-2]] == 'E' and ss_dict[resnum_count_dict[count-3]] == 'E' and ss_dict[resnum_count_dict[count-4]] == 'E' and sheet_test(count) == 1 and ss_dict[resnum_count_dict[count]] == '-':
                ss_dict[resnum_count_dict[count]] = 'e'
        except:
            pass


    # logic for assigning the final secondary structures
    for resnum in resnum_list:
        try:
            final_dict[resnum] = pp_dict[resnum]
        except:
            pass
        try:
            if(final_dict[resnum] == '-'):
                final_dict[resnum] = ss_dict[resnum]
        except:
            pass
        try:
            if((final_dict[resnum] == 't') and (ss_dict[resnum] != 'h')and (ss_dict[resnum] != 'E')and (ss_dict[resnum] != 'e')and (ss_dict[resnum] != '-')):
                final_dict[resnum] = ss_dict[resnum]
        except:
            pass
        try:
            if((pp_dict[resnum] == 'E') and ((ss_dict[resnum] == 'T') or (ss_dict[resnum] == 'N'))):
                final_dict[resnum] = ss_dict[resnum]
        except:
            pass
        try:
            if(((pp_dict[resnum] == 'E') or (pp_dict[resnum] == 'b'))and ((ss_dict[resnum] == 'B') or (ss_dict[resnum] == 'b'))):
                final_dict[resnum] = ss_dict[resnum]
        except:
            pass
        try:
            if(((pp_dict[resnum] == 'H') and ((ss_dict[resnum] == 'G')))):
                final_dict[resnum] = ss_dict[resnum]
        except:
            pass
        try:
            if( (ss_dict[resnum] == 'P') and (pp_dict[resnum] != 'T') and (pp_dict[resnum] != 'N')):
                final_dict[resnum] = ss_dict[resnum]
        except:
            pass


                    ### UNASSIGNED CIS PEPTIDES ASSIGNED HERE
    for count in resnum_count:

        try:
            if((cis_calc(count) == 1) and final_dict[resnum_count_dict[count]] == '-'):
                final_dict[resnum_count_dict[count]] = '='
        except:
            pass


    ############# DONE ASSIGNING #############

    # only prints if there is no output file specified
    if args.output == None:
        
        # Dmitriy helped out a lot here, though the code has been heavily heavily modified by this point
        def chunks(l,n):
            for i in range(0,len(l),n):
                yield l[i:i+n]

        
        # We only use the first value in this dictionary. We need it though still
        # It's frustrating but way way too late to change
        
        # all this if and elifs and such are just to get human readable output
        print("Chain:", chain_id_dict[resnum_list[0]])
        
        for chunk in chunks(resnum_count,40):
        
            j = 0
            for count in chunk:
                j += 1
                
                if j == 1:
                    sys.stdout.write("Final Struc:\t")
                try:
                    if gap_check(resnum_count_dict[count]) == 'behind' or gap_check(resnum_count_dict[count]) == 'isolated':
                        sys.stdout.write("_")
                except:
                    pass
                if j == 1:
                    sys.stdout.write(final_dict[resnum_count_dict[count]])
                    
                elif j < len(chunk) and (j % 10) != 0:
                    sys.stdout.write(final_dict[resnum_count_dict[count]])
                    
                elif j < len(chunk) and (j % 10) == 0:
                    sys.stdout.write(final_dict[resnum_count_dict[count]] + '\t')
                elif j == len(chunk):
                    print(final_dict[resnum_count_dict[count]])
                


            if args.verbose == True:
                j = 0
                for count in chunk:
                    j += 1
                    
                    if j == 1:
                        sys.stdout.write("PhiPsi2 Struc:\t")
                    try:
                        if gap_check(resnum_count_dict[count]) == 'behind' or gap_check(resnum_count_dict[count]) == 'isolated':
                            sys.stdout.write("_")
                    except:
                        pass
                    if j == 1:
                        sys.stdout.write(pp_dict[resnum_count_dict[count]])
                        
                    elif j < len(chunk) and (j % 10) != 0:
                        sys.stdout.write(pp_dict[resnum_count_dict[count]])
                        
                    elif j < len(chunk) and (j % 10) == 0:
                        sys.stdout.write(pp_dict[resnum_count_dict[count]] + '\t')
                    elif j == len(chunk):
                        print(pp_dict[resnum_count_dict[count]])

                
                j = 0
                for count in chunk:
                    j += 1
                    
                    if j == 1:
                        sys.stdout.write("Sec. Structure\t")
                    try:
                        if gap_check(resnum_count_dict[count]) == 'behind' or gap_check(resnum_count_dict[count]) == 'isolated':
                            sys.stdout.write("_")
                    except:
                        pass
                        
                    if j == 1:
                        sys.stdout.write(ss_dict[resnum_count_dict[count]])
                        
                    elif j < len(chunk) and (j % 10) != 0:
                        sys.stdout.write(ss_dict[resnum_count_dict[count]])
                        
                    elif j < len(chunk) and (j % 10) == 0:
                        sys.stdout.write(ss_dict[resnum_count_dict[count]] + '\t')
                    elif j == len(chunk):
                        print(ss_dict[resnum_count_dict[count]])



            j = 0
            for count in chunk:
                j += 1
                
                if j == 1:
                    sys.stdout.write("Residue Type\t")
                try:
                    if gap_check(resnum_count_dict[count]) == 'behind' or gap_check(resnum_count_dict[count]) == 'isolated':
                        sys.stdout.write("_")
                except:
                    pass
                if j == 1:
                    sys.stdout.write(aa_dict[resnum_count_dict[count]])
                elif j < len(chunk) and (j % 10) != 0:
                    sys.stdout.write(aa_dict[resnum_count_dict[count]])
                elif j < len(chunk) and (j % 10) == 0:
                    sys.stdout.write(aa_dict[resnum_count_dict[count]] + '\t')
                elif j == len(chunk):
                    print(aa_dict[resnum_count_dict[count]])

            j = 0
            for count in chunk:
                j += 1
                
                if j == 1:
                    sys.stdout.write("Residue Number\t")
                if j == 1:
                    sys.stdout.write(resnum_count_dict[count] + '\t\t')
                elif j < len(chunk) and ((j - 1) % 10) == 0:
                    sys.stdout.write(resnum_count_dict[count] + '\t\t')
                elif j == len(chunk):
                    print('\n')
                
    # presumably there is an output file specified, so this makes a temporary output for each chain
    else:
        
        temp = open(args.output+".temp.safe_to_delete"+".chain"+chain_id_dict[resnum_list[0]], "w")

        for resnum in resnum_list:
            if args.verbose == True:
                try:
                    temp.write(str(resnum)+"\t"+str(aa_dict[resnum])+"\t"+str(chain_id_dict[resnum])+"\t"+str(phi_dict[resnum])+"\t"+str(psi_dict[resnum])+"\t"+str(ss_dict[resnum])+"\t"+str(pp_dict[resnum])+"\t"+str(final_dict[resnum])+"\n")
                except:
                    pass    
            else:
                try:
                    temp.write(str(resnum)+"\t"+str(aa_dict[resnum])+"\t"+str(chain_id_dict[resnum])+"\t"+str(final_dict[resnum])+"\n")
                except:
                    pass     
        temp.close()

# should be obvious by now, but checking to see if the user asked for output or not
if args.output == None:
    pass

# concatenates each output chain into one output file, and removes temp files as it goes.
else:
    output = open(args.output,'wb')
    if args.verbose == True:
        output.write(bytes("Residue_ID"+"\t"+"Residue_Type"+"\t"+"Chain_ID"+"\t"+"Phi"+"\t"+"Psi"+"\t"+"Virtual_Secondary_Structure"+"\t"+"PhiPsi2_Motif_Structure"+"\t"+"Best_Secondary_Structure"+"\n", 'UTF-8'))
    else:
        output.write(bytes("Residue_ID"+"\t"+"Residue_Type"+"\t"+"Chain_ID"+"\t"+"Secondary_Structure"+"\n", 'UTF-8'))        
    for filename in glob.glob('*temp.safe*'):
        shutil.copyfileobj(open(filename, 'rb'), output)
        os.remove(filename)
    output.close()

# delete downloaded files
if args.input == None:
    code = str(args.code)
    pdb = "%s.pdb" % str(code.lower())
    os.remove(pdb)
