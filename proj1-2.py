# Project 1-2
#
# Code by Jerry Yan

import imp_xyz as ix
import m_tools as mt
import time

file_name = "propanediol-2"

_, _, mol = ix.import_xyz("molecule_models/{}.xyz".format(file_name))
dmol = ix.select_by_name(mol,'C') + ix.select_by_name(mol,'O')

# De-signing dmol (taking absolute value)
for i, m in enumerate(dmol):
    dmol[i]['rvec'] = list(map(abs, dmol[i]['rvec']))

mol_C = ix.select_by_name(dmol, 'C')
mol_O = ix.select_by_name(dmol, 'O')

CO_BOND_LENGTH = 1.43
CC_BOND_LENGTH = 1.54

TOLERANCE_LEVEL = 0.15

# Fix the first O (O1)
# We choose the first fixed (and the only fixed) atom to be O because it should be easy to find the C adjacent to it. Other two C-s are further than its neighboring C.
O1 = mol_O[1]
print("O1 has been fixed.")

mol_Or = mol_O
mol_Or.remove(O1)

# Find the neighboring C of O1 (C2)
C2_possibles = mt.find_possibles(mol_C)

C2_af = filter(lambda m: mt.bond_length_filter(m['rvec'], O1['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL), C2_possibles) # Change the tolerance range here to see how this filter works

if len(C2_af) == 1:
    print("C2 has been found.")
    C2 = C2_af[0]
else:
    mt.bond_length_error(C2_af, 'C2')

# Find the possible neighboring C(s) of C2 (C3)
C3_possibles = filter(lambda m: m['index'] != C2['index'], C2_possibles)

C3_af = filter(lambda m: mt.bond_length_filter(m['rvec'], C2['rvec'], CC_BOND_LENGTH, TOLERANCE_LEVEL), C3_possibles) # Change the selection of O1 to see how the code below works

if len(C3_af) == 1: # C2 is adjacent to one C. Choose that C to be C3.
    print("C2 is adjacent to one C (C3).\nC3 has been found.")
    C3 = C3_af[0]
    C2_adj = 1
elif len(C3_af) == 2 and C3_af[0]['index'] != C3_af[1]['index']: # C2 is adjacent to two C-s. Let those C-s be C3 and C4 respectively.
    print("C2 is adjacent to two C-s.\nC3 has been found. \nC4 has been found.")
    C3 = C3_af[0]
    C4 = C3_af[1]
    C2_adj = 2
else:
    mt.bond_length_error(C3_af, 'C3')

# Find possible O5
O5_possibles = mt.find_possibles(mol_Or)

if C2_adj == 1: # If C2 is adjacent to one C (C3), then C3 must be connected to an O (O5).
    O5_af = filter(lambda m: mt.bond_length_filter(m['rvec'], C3['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL), O5_possibles)

    if len(O5_af) == 1:
        print("O5 has been found.")
        O5 = O5_af[0]
    else:
        mt.bond_length_error(O5_af, 'O5')

    # Find the possible C4
    C4_possibles = filter(lambda m: m['index'] != C3['index'], C3_possibles)

    C4_af = filter(lambda m: mt.bond_length_filter(m['rvec'], C3['rvec'], CC_BOND_LENGTH, TOLERANCE_LEVEL), C4_possibles)
    if len(C4_af) == 1:
        print("C4 has been found.")
        C4 = C4_af[0]
    else: # Make sure C4 is not close to O5 (which is adjacent to C3)
        C4_af2 = []
        d_C3O5 = mt.distance(C3['rvec'], O5['rvec'])
        for m in C4_af:
            if mt.distance(m['rvec'], O5['rvec']) > d_C3O5:
                C4_af2.append(m)
        if len(C4_af2) == 1:
            print("C4 has been found.")
            C4 = C4_af2[0]
        else:
            mt.bond_length_error(C4_af2, 'C4')

elif C2_adj == 2: # If C2 is adjacent to two C-s (C3, C4), then O5 must be adjacent to one of C3 and C4.
    O5_afC3 = filter(lambda m: mt.bond_length_filter(m['rvec'], C3['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL), O5_possibles)
    O5_afC4 = filter(lambda m: mt.bond_length_filter(m['rvec'], C4['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL), O5_possibles)

    O5_af = O5_afC3 + O5_afC4
    if len(O5_af) == 1:
        print("O5 has been found.")
        O5 = O5_af[0]
    else:
        mt.bond_length_error(O5_af, 'O5')

print("All atoms have been found.")
mols = [O1, C2, C3, C4, O5]
#ix.export_xyz('molecule_results/{0}_{1}.xyz'.format(file_name,int(time.time())), mols)

print("Results saved to xyz file.")
