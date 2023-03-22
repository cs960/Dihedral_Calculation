import os

import sys

import numpy as np
import math
import glob
import re
import csv
  

#Calculate dihedral angle
def calc_dihedral(A, B, C, D):
    Vector_AB = [(B[0]-A[0]),(B[1]-A[1]),(B[2]-A[2])]
    Vector_BC = [(C[0]-B[0]),(C[1]-B[1]),(C[2]-B[2])]
    Vector_CD = [(D[0]-C[0]),(D[1]-C[1]),(D[2]-C[2])]
    v1 = np.cross(Vector_AB, Vector_BC)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(Vector_BC, Vector_CD)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * Vector_CD).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm
    deg = (rad*180)/np.pi
    return deg
    

#Add amide conformation to list
def amide_conf(angle):
    if -90 <angle<90:
        amide_conform = 'Cis'
    else:
        amide_conform = 'Trans'
    return amide_conform

 
    
    
#Add internal re or si face assignment to list
def re_si_conf(angle):
    if 0 < angle < 180:
        internal_face = 'Si'
    else:
        internal_face = 'Re'
    return internal_face


#Parseing files
def get_molecule_from_input(lines, filename, energy, indices):
    molecule = {}
    molecule['filename'] = filename
    molecule['energy'] = float(lines[1])
    xyz = []
    
#    for h in indices:
    atomid = 0
    for l in lines:
            #l = re.split(' +', l)
        l = l.split()
        if len(l) > 4 and l[3] in "CNHO":
            atomid += 1

#            if atomid == indices(h):
            for j in l[:3]:
                xyz.append(float(j))
                
    A_xyz = xyz[0:3]
    B_xyz = xyz[3:6]
    C_xyz = xyz[6:9]
    D_xyz = xyz[9:12]
    molecule['A_xyz'] = A_xyz
    molecule['B_xyz'] = B_xyz
    molecule['C_xyz'] = C_xyz
    molecule['D_xyz'] = D_xyz
    molecule['dihedral'] = calc_dihedral(A_xyz, B_xyz, C_xyz, D_xyz)
    dihedral = molecule['dihedral']
    molecule['Amide_cis-or-trans'] = amide_conf(dihedral)
    #molecule['Re-or-Si'] = re_si_conf(dihedral)
    #This way of assigning re or si doesn't really work. Need a more robust way.
    return molecule
    

def get_extension(filename):
    ff = filename.split('.')
    if not len(ff) == 2:
        return None
    return ff[1]




    
    
    


molecules = []
energy = []
dihedrals = []
folder = sys.argv[1]
indices = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])]
#error message to tell you what arguements to put in . 
print('sys.argv =', sys.argv) #the indices input comes in as string, 
if len(sys.argv) != 6:
 print ('ERROR - type "python [program.py [= sys.argv[0]]] [folder location [= sys.argv[1]]] [atom numbers of dihedral [= sys.argv[2-5]]]')


for file in os.listdir(folder):
    filename = os.fsdecode(file)
    ext = get_extension(filename)
    if ext != 'sdf':
        continue
    with open(f"{folder}/{file}","r") as infile:
        lines = infile.readlines()
#    print(lines)
    molecules.append(get_molecule_from_input(lines, filename, energy, indices))

#gets rid of "ref" molecule
molecules = molecules[0:-1]


# the array molecules will contain the x,y,z coordinates of each
# "molecule" (equivalent to the input file, as I am assuming this is one molecule)
# as a sub-array of all the atoms in it
# so this is a list, with each itemm in the list (the molecule)
# being a collection of 3 properties (filename, number, molecule)
# the molecule property is a list of the atom coordinates
# your subsequent program should be able to process this
#  [
#    {
#      "filename": <name>,
#      "number_of_atoms": <num>,
#      "atoms": 
#          [
#            [x1, y1, z1],
#            [x2, y2, z3],
#            ...
#            [xn, yn, zn]    # where n = number_of_atoms
#          ]
#    }
#  ]

# if you just want a long file full of atoms coordinates, one file 
# per molecule, then use this:

#Each molecule file reduced to relevant 4 atoms and the energy
#for molecule in molecules:
#    xyz = molecule['A_xyz']
#    f = molecule['filename'].split('.')
#    new_filename = f[0]+'.out'
#    with open(f"{folder}/{new_filename}", 'w') as outfile:
#        outfile.write(f'energy = {str(molecule["energy"])}\nA_xyz = {str(molecule["A_xyz"])}\nB_xyz = {str(molecule["B_xyz"])}\nC_xyz = {str(molecule["C_xyz"])}\nD_xyz = {str(molecule["D_xyz"])}\ndihedral = {str(molecule["dihedral"])}')

#Write file with dihedral outputs
#new_filename = 'results.out'
#with open(f"{folder}/{new_filename}", 'w') as outfile:
#    for molecule in molecules:
#       outfile.write(f'filename = {str(molecule["filename"])}\nenergy = {str(molecule["energy"])}\nA_xyz = {str(molecule["A_xyz"])}\nB_xyz = {str(molecule["B_xyz"])}\nC_xyz = {str(molecule["C_xyz"])}\nD_xyz = {str(molecule["D_xyz"])}\ndihedral = {str(molecule["dihedral"])}\n\n')
       
#Write .csv file with dihedral outputs
#new_filename = 'results.csv'
#with open(f"{folder}/{new_filename}", 'w') as outfile:
#        outfile.write(molecules)

field_names = ['filename', 'energy', 'A_xyz', 'B_xyz', 'C_xyz', 'D_xyz', 'dihedral', 'Amide_cis-or-trans']        
with open(f'{folder}/output.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames = field_names)
    writer.writeheader()
    writer.writerows(molecules)
       