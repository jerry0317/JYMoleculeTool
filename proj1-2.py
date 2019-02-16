# Project 1-2
#
# Code by Jerry Yan

import imp_xyz as ix
import numpy as np
import itertools

INDEX_RANGE = range(0,3)

_, _, mol = ix.import_xyz("molecule_models/propanediol-1.xyz")
dmol = ix.select_by_name(mol,'C') + ix.select_by_name(mol,'O')

# De-signing dmol (taking absolute value)
for i, m in enumerate(dmol):
    dmol[i]['rvec'] = list(map(abs, dmol[i]['rvec']))

mol_C = ix.select_by_name(dmol, 'C')
mol_O = ix.select_by_name(dmol, 'O')

# Fix the first O (O1)
# We choose the first fixed (and the only fixed) atom to be O because it should be easy to find the C adjacent to it. Other two C-s are further than its neighboring C.
O1 = mol_O[0]['rvec']
print("O1 has been fixed.")

def re_sign_n(rvec):
    pl = []

    SLIST = list(map(list, itertools.product([-1,1], repeat=len(INDEX_RANGE))))

    for s in SLIST:
        r = rvec
        nr = [r[i]*s[i] for i in INDEX_RANGE]
        pl.append(nr)

    return pl

# Find the neighboring C of O1 (C2)
C2_possibles = []
for i in range(0, len(mol_C)):
    m = mol_C[i]
    pl = re_sign_n(m['rvec'])
    for p in pl:
        d = {'rvec': p, 'name': m['name'], 'index': i}
        C2_possibles.append(d)

# Filtering by the bond length
def bond_length_filter(rvec_1, rvec_2, b_length, tol_range = 0.1):
    dvec = np.array(rvec_1) - np.array(rvec_2)
    d = np.sqrt(np.dot(dvec, dvec))
    if (b_length - tol_range) < d < (b_length + tol_range):
        return True
    else:
        return False

C2_af = filter(lambda x: bond_length_filter(x['rvec'], O1, 1.43, 1.0), C2_possibles)

if len(C2_af) == 1:
    print("C2 has been found.")
    C2 = mol_C[C2_af[0]['index']]
else:
    print("Unable to narrow down to one C2 by bond-length filtering.")
    print("The most narrow results for C2 are \n {}".format(str(C2_af)))
    exit()
