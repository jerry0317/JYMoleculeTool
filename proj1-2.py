# Project 1-2
#
# Code by Jerry Yan

import imp_xyz as ix
import m_tools as mt
import time
import itertools

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
O1 = mol_O[0]
print("O1 has been fixed.")

mol_Or = mol_O
mol_Or.remove(O1) # The remaining set of O excluding O1

def find_C2(O1, mol_C):
    # Find the neighboring C of O1 (C2)
    C2_possibles = mt.find_possibles(mol_C)
    C2_af = filter(lambda m: mt.bond_length_filter(m['rvec'], O1['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL), C2_possibles) # Change the tolerance range here to see how this filter works
    return C2_af

def find_C3C4(C2, O1, mol_C):
    # Find the possible neighboring C(s) of C2 (C3/C4)
    C3_possibles = filter(lambda m: m['index'] != C2['index'], mt.find_possibles(mol_C))
    C3_af = filter(lambda m: mt.bond_length_filter(m['rvec'], C2['rvec'], CC_BOND_LENGTH, TOLERANCE_LEVEL), C3_possibles) # Change the selection of O1 to see how the code below works
    C3_af_index = set()

    for c3 in C3_af:
        C3_af_index.add(c3['index'])

    if len(C3_af_index) == 1:
        C2_adj = 1 # C2 is adjacent to one C. Choose that C to be C3.
        C3C4_af = C3_af
    elif len(C3_af_index) == 2:
        C2_adj = 2 # C2 is adjacent to two C-s. Let those C-s be C3 and C4 respectively.
        C3_af_grps = []
        for i in C3_af_index:
            C3_af_grps.append(filter(lambda m: m['index'] == i, C3_af))
        C3C4_af = list(map(list, itertools.product(C3_af_grps[0], C3_af_grps[1])))
    else:
        C2_adj = 0
        C3C4_af = []

    return C2_adj, C3C4_af

def find_O5_C2adj1(C3, mol_Or): # If C2 is adjacent to one C (C3), then C3 must be connected to an O (O5).
    return find_C2(C3, mol_Or) # The same procedure as finding C2

def find_C4_C2adj1(C2, C3, O5, mol_C):
    # Find the possible C4
    C4_possibles = filter(lambda m: m['index'] not in (C2['index'], C3['index']), mt.find_possibles(mol_C))
    C4_af = filter(lambda m: mt.bond_length_filter(m['rvec'], C3['rvec'], CC_BOND_LENGTH, TOLERANCE_LEVEL), C4_possibles)

    C4_af2 = []
    d_C3O5 = mt.distance(C3['rvec'], O5['rvec']) # Make sure C4 is not close to O5 (which is adjacent to C3)
    for m in C4_af:
        if mt.distance(m['rvec'], O5['rvec']) > d_C3O5:
            C4_af2.append(m)

    return C4_af2

def find_O5_C2adj2(C3, C4, mol_Or): # If C2 is adjacent to two C-s (C3, C4), then O5 must be adjacent to one of C3 and C4.
    O5_possibles = mt.find_possibles(mol_Or)
    O5_afC3 = filter(lambda m: mt.bond_length_filter(m['rvec'], C3['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL) and mt.distance(m['rvec'], C4['rvec']) > (CO_BOND_LENGTH + TOLERANCE_LEVEL), O5_possibles)
    O5_afC4 = filter(lambda m: mt.bond_length_filter(m['rvec'], C4['rvec'], CO_BOND_LENGTH, TOLERANCE_LEVEL) and mt.distance(m['rvec'], C3['rvec']) > (CO_BOND_LENGTH + TOLERANCE_LEVEL), O5_possibles)

    O5_af = O5_afC3 + O5_afC4
    return O5_af

def save_mols(mols, icode = None):
    print("-----The structure of molecule is found as follows:-----")
    for m in mols:
        print("{0}    {1}".format(m['name'], m['rvec']))
    if icode is None:
        filename = 'molecule_results/{0}_{1}.xyz'.format(file_name,int(time.time()))
    else:
        filename = 'molecule_results/{0}_{1}_{2}.xyz'.format(file_name,int(time.time()),str(icode))
    ix.export_xyz(filename, mols)
    print("Results saved to xyz file.")


C2_af = find_C2(O1, mol_C)
for C2 in C2_af:
    C2_adj, C3C4_af = find_C3C4(C2, O1, mol_C)
    if C2_adj == 1:
        for C3C4 in C3C4_af:
            C3 = C3C4
            O5_af = find_O5_C2adj1(C3, mol_Or)
            for O5 in O5_af:
                C4_af = find_C4_C2adj1(C2, C3, O5, mol_C)
                i = 1
                for C4 in C4_af:
                    save_mols([O1, C2, C3, C4, O5], i)
                    i = i + 1
    elif C2_adj == 2:
        for C3C4 in C3C4_af:
            C3 = C3C4[0]
            C4 = C3C4[1]
            O5_af = find_O5_C2adj2(C3, C4, mol_Or)
            i = 1
            for O5 in O5_af:
                save_mols([O1, C2, C3, C4, O5], i)
                i = i + 1
