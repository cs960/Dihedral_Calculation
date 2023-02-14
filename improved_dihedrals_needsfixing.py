""" Dihedral angles from coordinates. Get dihedral angles from 
    protein backbone coordinates.
    This script takes 4 vectors representing 4 points as input
    and returns the dihedral angle between them.
    
    The script was originally developed to calculate dihedral angles (phi,psi)
    from protein backbone atoms' coordinates in 3D (N-term, C-alhpa, Cterm).
"""

import numpy as np
import math
import os
import glob
import numpy
import re

energies = []
dihedrals = []
amide_conformations =[]
internal_faces = []

# Parse files
sdf_folder = r"C:\Users\s1951009\FW-Avogadro-search\Confbuster results\Unzipped\AdA-no-methyl-charged_Si-start\AdA-no-methyl-charged_Si-start"
print('Working directory = ',sdf_folder,'\n')

for file in os.listdir(sdf_folder):
     filename = os.fsdecode(file)
     if filename.startswith("conf-"):
      print ('filename = ',filename)
      continue
     else:
      break

for file in os.listdir(sdf_folder):
     filename = os.fsdecode(file)
     if filename.startswith("conf-"):
      print ('filename = ',filename)
      with open(file,"r") as outfile:
       sdf = outfile.readlines()
      
      
      
      
      
      
      continue
     else:
      break



# Enter atom index for amide bond of a set of molecules
      print ('Enter ABCD atom indexes (order of C, C(=O), N, C) (separated by spaces): ')
      atomlist = input()
      atomlist_delimited = atomlist.split()
      atomindex = []
      atomindex = list(map(int, atomlist_delimited))
      print ('Atom indexes = ',atomindex)

#Extract coordinates and atom types from index
      atomtypes_ABCD = []
      ABCD_coordinates = []
      for i in atomindex:
       atomtypes_ABCD.append(sdf[3+i][31])
       xyz_string = sdf[3+i][0:31]
       xyz_string_delimited = xyz_string.split()
      for j in xyz_string_delimited:
       ABCD_coordinates.append(float(j))
#Error message
      if atomtypes_ABCD != ['C','C','N','C']:
       print ('ERROR: Atoms types are not expected for an amide - check atom index is accurate')
#print ('ABCD coordinates = ',ABCD_coordinates)

#Add energy data to lsit
     energy = float(sdf[1])
     energies.append(int(energy))

# Processing coordinates
     atom_A = ABCD_coordinates[0:3]
     atom_B = ABCD_coordinates[3:6]
     atom_C = ABCD_coordinates[6:9]
     atom_D = ABCD_coordinates[9:12]
#print('\natom A = ',atom_A,'\natom B = ',atom_B,'\natom C = ',atom_C,'\natom D = ',atom_D)


#Calculate dihedral angle
def calc_dihedral(u1, u2, u3, u4):
    """ Calculate dihedral angle method. From bioPython.PDB
    (adapted to np.array)
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    [-pi, pi].
    """

    Vector_AB = [(atom_B[0]-atom_A[0]),(atom_B[1]-atom_A[1]),(atom_B[2]-atom_A[2])]
    Vector_BC = [(atom_C[0]-atom_B[0]),(atom_C[1]-atom_B[1]),(atom_C[2]-atom_B[2])]
    Vector_CD = [(atom_D[0]-atom_C[0]),(atom_D[1]-atom_C[1]),(atom_D[2]-atom_C[2])]

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
    
#Add dihedral angle to list
alpha_degrees = calc_dihedral(atom_A,atom_B,atom_C,atom_D)
dihedrals.append(int(alpha_degrees))

#Add amide conformation to list
if -90 <alpha_degrees <90:
 amide_conform = 'Cis'
else:
 amide_conform = 'Trans'
amide_conformations.append(amide_conform)
 
    
    
#Add internal re or si face assignment to list
if 0 < alpha_degrees < 180:
 internal_face = 'Si'
else:
 internal_face = 'Trans'
internal_faces.append(internal_face)
    
sdf_files = []
confbust_output = []
print (os.listdir(sdf_folder),'\n')
confbust_output = sorted(os.listdir(sdf_folder))
print(confbust_output)
for i in confbust_output:
 if i.startswith('conf-'):
  sdf_files.append(i)
print ('sdf files = ',sdf_files)

print ('List of Energies: ',energies)
print ('List of Dihedrals: ',dihedrals)
print ('List of Amide conformations: ',amide_conformations)
print ('List of Internal Re or Si faces: ', internal_faces)