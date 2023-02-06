#Program to detect cis or trans amides from C-N-C(=O)-C angle.

print ('\nDIHEDRAL CALCULATION & CIS/TRANS ASSIGNMENT\n\n')
print ('           D')
print ('          /')
print ('   B ――― C')
print ('  /')
print (' A\n\n\n')

import math

# Processing coordinates

xyz_input = [1,2,1,4.5,4,1,4.8,6,1,7,8,2]

dihedral_atom_coordinates = [['xa','ya','za'],['xb','yb','zb'],['xc','yb','zc'],['xd','yd','zd']]

i = 0
atom = 0
while atom < 4:
 xyz = 0
 while xyz < 3:
  dihedral_atom_coordinates[atom][xyz] = xyz_input[i]
  xyz += 1
  i += 1
 atom += 1
 
atom_A = dihedral_atom_coordinates[0]
atom_B = dihedral_atom_coordinates[1]
atom_C = dihedral_atom_coordinates[2]
atom_D = dihedral_atom_coordinates[3]

# Calculating vectors

# These are points along a line extended from the vector B-C, completing a right angle between A and B [Bp] and C and D [Cp], and thus forming two triangles A_tri and D_tri
print ('Dcol    -     -     -      -      D')
print ('|                                /|')
print ('                                /  ')
print ('|                              /  |')
print ('                              /    ')
print ('|  Right triangle D_tri -->  /     ')
print ('                            /     |')
print ('                           /       ')
print ('|                         /       |')
print ('Bp - - - B ――――――――――――― C - - - Cp')
print ('|       /')
print ('       /')
print ('|     /')
print ('     /  <--  Right triangle A_tri')
print ('    /')
print ('|  /')
print ('  /')
print ('|/')
print ('A')
print ()

# Since Bp is on the extended line C-B, angle A-B-Bp = (180 - angle A-B-C)')
# The angles A-B-C and B-C-D are calculated using the dot product of the adjacent vectors and their magnitudes (a.k.a.bond lengths)
# To calculate the angle α between two vectors in 2D space:
# α = arccos[(xa · xb + ya · yb) / (√(xa² + ya²) · √(xb² + yb²))]
# (divide the dot product of the vectors by the product of their magnitudes)
# I think it's the same for vectors in a 3D space.

# Bond Lengths
Vector_AB = [(atom_B[0]-atom_A[0]),(atom_B[1]-atom_A[1]),(atom_B[2]-atom_A[2])]
Vector_BC = [(atom_C[0]-atom_B[0]),(atom_C[1]-atom_B[1]),(atom_C[2]-atom_B[2])]
Vector_CD = [(atom_D[0]-atom_C[0]),(atom_D[1]-atom_C[1]),(atom_D[2]-atom_C[2])]
bondlength_AB = math.sqrt(Vector_AB[0]**2 + Vector_AB[1]**2 + Vector_AB[2]**2)
bondlength_BC = math.sqrt(Vector_BC[0]**2 + Vector_BC[1]**2 + Vector_BC[2]**2)
bondlength_CD = math.sqrt(Vector_CD[0]**2 + Vector_CD[1]**2 + Vector_CD[2]**2)
print ('Vector AB = ',Vector_AB,'\nVector BC = ',Vector_BC,'\nVector CD = ',Vector_CD)
print ('bondlength AB = ',bondlength_AB,'\nbondlength BC = ',bondlength_BC,'\nbondlength CD = ',bondlength_CD)

# Bond Angles
dotproduct_AB_BC = ((Vector_AB[0]*Vector_BC[0]) + (Vector_AB[1]*Vector_BC[1]) + (Vector_AB[2]*Vector_BC[2]))
dotproduct_BC_CD = ((Vector_BC[0]*Vector_CD[0]) + (Vector_BC[1]*Vector_CD[1]) + (Vector_BC[2]*Vector_CD[2]))
angle_ABC = (math.pi) - math.acos((dotproduct_AB_BC)/((bondlength_AB)*(bondlength_BC)))
angle_BCD = (math.pi) - math.acos((dotproduct_BC_CD)/((bondlength_BC)*(bondlength_CD)))
print ('dotproduct AB BC = ', dotproduct_AB_BC,'\ndotproduct BC CD = ', dotproduct_BC_CD,'\nangle A-B-C = ', angle_ABC*(180/math.pi), 'degrees.\nangle B-C-D is ', angle_BCD*(180/math.pi), 'degrees.')

# Because the vectors are directional, it's all too easy to get it mixed up due to the wrong side of a 180' you're looking at.
# From angle A-B-Bp, the known 90 degree angle B-Bp-A and the known bond length AB we can use the cos relationship to establish the distance B-Bp, an extension of the line C-B. From this distance, and the equation of the line C-B, we can calculate the coordinates of the point Bp, becaus the magnitude of the vector BC is equal to the square root of the sum of the squared differences in x and y, while at the same time the point must sit on the line BC so we can use simultaneous equations to figure out the position of the line.
# Then we need to get the lengths of the other sides in the extended out right angle triangle, the right angle of which sits on the new point P2i, and the hypotenuse the the bond P1-P2.
# (180 minus - angle_ABC) is adjacent to the new line P2_P2i and opposite to the new line P1_P2i which is perpendicular to the line P2_P3

# cos(180-angle_P1_P2-P3) = length[P2_P2i]/length[P1_P2]. So the cos function can give us the lengths B_Bp and C_Cp, and we can very simply extract - from the direction and magnitude(bondlength) of the vector BC - values for axial "displacement", telling us how the x, y and z components change per unit length along the line BC.
# So, we can use this information to identify the x, y and z coordinates of the points Bp and Cp from the lengths B_Bp and C_Cp.

# Lengths B to Bp and C to Cp
length_B_Bp = math.cos(math.pi-angle_ABC)*bondlength_AB
length_C_Cp = math.cos(math.pi-angle_BCD)*bondlength_CD
print('length B-Bp = ',length_B_Bp,'\nLength C-Cp = ',length_C_Cp,'\n')

# Coordinates of Bp and Cp
displacement_BC_x = Vector_BC[0] / bondlength_BC
displacement_BC_y = Vector_BC[1] / bondlength_BC
displacement_BC_z = Vector_BC[2] / bondlength_BC

xbp = atom_B[0] - displacement_BC_x * length_B_Bp
ybp = atom_B[1] - displacement_BC_y * length_B_Bp
zbp = atom_B[2] - displacement_BC_z * length_B_Bp
point_Bp = [xbp,ybp,zbp]
xcp = atom_C[0] + displacement_BC_x * length_C_Cp
ycp = atom_C[1] + displacement_BC_y * length_C_Cp
zcp = atom_C[2] + displacement_BC_z * length_C_Cp
point_Cp = [xcp,ycp,zcp]
print('atom A = ',atom_A,'\natom B = ',atom_B,'\natom C = ',atom_C,'\natom D = ',atom_D)
print ('Point Bp = ',point_Bp,'\nPoint Cp = ',point_Cp)

# Angle between vector Bp-A and Bp-Dcol: we need to collapse point D along the central bond so that we create the new point Dcol. 
xdcol = atom_D[0] - (point_Cp[0]-point_Bp[0])
ydcol = atom_D[1] - (point_Cp[1]-point_Bp[1])
zdcol = atom_D[2] - (point_Cp[2]-point_Bp[2])
point_Dcol = [xdcol,ydcol,zdcol]
print ('point Dcol = ',point_Dcol,'\n')

# Now angle between them is calculated the same as we calculated it before using the dot product, magnitude of vectors (from simple pythagoras theorem) and cosine equation.

# Vectors in "Dihedral Triangle"
Vector_ABp = [point_Bp[0] - atom_A[0], point_Bp[1] - atom_A[1], point_Bp[2] - atom_A[2]]
print ()
Vector_BpDcol = [point_Dcol[0] - point_Bp[0], point_Dcol[1] - point_Bp[1], point_Dcol[2] - point_Bp[2]]
print ('Vector ABp = ',Vector_ABp,'\nVector BpDcol = ',Vector_BpDcol)

# Lengths in "Dihedral Triangle"
length_A_Bp = math.sqrt((bondlength_AB**2) / (length_B_Bp)**2)
length_Bp_Dcol = math.sqrt((bondlength_CD**2) / (length_C_Cp)**2)
print('length_A_Bp = ',length_A_Bp,'\nlength_Bp_Dcol = ',length_Bp_Dcol)

# Dihedral Calculation
dotproduct_ABp_BpDcol = ((Vector_ABp[0]*Vector_BpDcol[0]) + (Vector_ABp[1]*Vector_BpDcol[1]) + (Vector_ABp[2]*Vector_BpDcol[2]))
print('dotproduct_ABp_BpDcol = ',dotproduct_ABp_BpDcol,'\n')

print (dotproduct_ABp_BpDcol/(length_A_Bp*length_Bp_Dcol))

#There's an extra math.pi in there (360 degrees) to get the value above down below pi. That is wrong and I should fix it.
dihedral_ABCD = (math.pi) - math.acos(((math.pi)-dotproduct_ABp_BpDcol/(length_A_Bp*length_Bp_Dcol)))
print ('Dihedral ABCD =', dihedral_ABCD*(180/math.pi), 'degrees.')
