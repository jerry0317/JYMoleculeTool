# Project 1-1
#
# Code by Jerry Yan

import i_tensor as it
import imp_xyz as ix

_, _, mol = ix.import_xyz("molecule_models/propanediol-1.xyz")
dmol = ix.select_by_name(mol,'C') + ix.select_by_name(mol,'O')

for m in dmol:
    mass = ix.mass_cgs_by_name(m['name'])
    rvec = ix.length_cgs_by_rvec(m['rvec'])
    i_tens = it.i_single_pm(mass, rvec)
    abc = it.abc_cgs_from_i(i_tens)
    print(abc, m['name'])

# mol_m = [ix.mass_cgs_by_name(m['name']) for m in mol]
# mol_rvec = [ix.length_cgs_by_rvec(m['rvec']) for m in mol]
#
# i_tens = it.i_mpm(mol_m, mol_rvec)
# print(i_tens)
#
# abc = it.abc_cgs_from_i(i_tens)
# print(abc)
