# Project 1-2
#
# Code by Jerry Yan

import imp_xyz as ix
import m_tools as mt
import time
import itertools

file_name = "branchedalcohola"

_, _, mol = ix.import_xyz("molecule_models/{}.xyz".format(file_name))
dmol = ix.select_by_name(mol,'C') + ix.select_by_name(mol,'O')

# De-signing dmol (taking absolute value) to obtain |x|, |y|, |z| for each oxygen or carbon atom
for i, m in enumerate(dmol):
    dmol[i]['rvec'] = list(map(abs, dmol[i]['rvec']))

mol_C = ix.select_by_name(dmol, 'C')
mol_O = ix.select_by_name(dmol, 'O')

mol_comb = mol_C + mol_O

mol_comb = mt.dict_list_to_atom_list(mol_comb)

CO_BOND_LENGTH = 1.43
CC_BOND_LENGTH = 1.54

TOLERANCE_LEVEL = 0.1 # In bond-length filters, an atom passes the filter when the distance lies between BOND_LENGTH plus/minus TOLERANCE_LEVEL

# Fix the first atom
A1 = mol_comb[0]

print("The first atom has been fixed.")

mol_combr = [m for m in mol_comb if m != A1]


possible_list = []

# TODO: DEBUGGING!!!

def rcs_constructer_atom(atom_r, st_mol, tol_range = 0.1):
    possible_atoms = mt.find_possibles_atom([atom_r])

    s_m_list = []
    for p_atom in possible_atoms:
        s_m = mt.bond_length_molecule_constructer(st_mol, atom_r, tol_range)

        if len(s_m.bond_graphs) > 0:
            s_m_list.append(s_m)

    return s_m_list

def rcs_action_atom(mol_r, m_list, tol_range = 0.1):
    if len(mol_r) > 0:
        for m in m_list:
            for a_r in mol_r:
                ml_af1 = rcs_constructer_atom(a_r, m, tol_range)
                if len(ml_af1) > 0:
                    mol_rf1 = [mol for mol in mol_r if mol != a_r]
                    rcs_action_atom(mol_rf1, ml_af1, tol_range)
    elif len(mol_r) == 0:
        for m in m_list:
            possible_list.append(m)



initial_smol = mt.Strc_molecule(atoms = {A1})
rcs_action_atom(mol_combr, [initial_smol], TOLERANCE_LEVEL)

print("Number of combinations to work with: {}".format(8**(len(mol_combr))))
print("Number of plausible results: {}".format(len(possible_list)))

def save_mols_atom(st_mol, icode = None):
    print("-----The structure of molecule is found as follows:-----")
    for a in st_mol.atoms:
        print("{0}    {1}".format(a.name, a.rvec))
    if icode is None:
        filename = 'molecule_results/{0}_{1}.xyz'.format(file_name,int(time.time()))
    else:
        filename = 'molecule_results/{0}_{1}_{2}.xyz'.format(file_name,int(time.time()),str(icode))
    # ix.export_xyz(filename, mols)
    # print("Results saved to xyz file.")

i = 0
for pl in possible_list:
    i += 1
    save_mols_atom(pl, i)
