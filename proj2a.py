# Project 1-2
#
# Code by Jerry Yan

import imp_xyz as ix
import m_tools as mt
import time
import itertools

file_name = "branchedalcohola"
SAVE_XYZ = False

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
        # print("patom rvec: {}".format(p_atom.rvec.dictvec))
        s_m = mt.bond_length_molecule_constructer(st_mol, p_atom, tol_range)

        if len(s_m.bond_graphs) > 0:
            # print("M PASSED")
            s_m_list.append(s_m)
        #else:
        #     print("M FAILED")
        #
        # print("----------")

    #print("***s_m_list:{}***".format(len(s_m_list)))

    return s_m_list

# NOTE: Problems are probably in the last step

def rcs_action_atom(mol_r, m_list, tol_range = 0.1):
    if len(mol_r) > 0:
        # print("==>0 part==")
        # print("mlist length: {}".format(len(m_list)))
        # print("mol_r length: {}".format(len(mol_r)))
        for m in m_list:
            for a_r in mol_r:
                ml_af1 = rcs_constructer_atom(a_r, m, tol_range)
                # print(ml_af1)
                if len(ml_af1) > 0:
                    mol_rf1 = [mol for mol in mol_r if mol != a_r]
                    # print("molrf1:{}".format(len(mol_rf1)))
                    rcs_action_atom(mol_rf1, ml_af1, tol_range)
                # else:
                    # print("empty ml_af1")
    elif len(mol_r) == 0:
        for m in m_list:
            duplicated = False
            for plist in possible_list:
                if plist == m:
                    duplicated = True
                    continue
            if duplicated is False:
                possible_list.append(m)
                print("The possible No. {} is found.".format(len(possible_list)))
    # print("molr length after: {}".format(len(mol_r)))
    # print("----RCS ACTION----")
    # time.sleep(0.2)

t0 = time.perf_counter()

initial_smol = mt.Strc_molecule(atoms = {A1})
rcs_action_atom(mol_combr, [initial_smol], TOLERANCE_LEVEL)

t_taken = time.perf_counter() - t0

print("-----The structure of molecule is found as follows:-----")

def save_mols_atom(st_mol, icode = None):
    # print("-----The structure of molecule is found as follows:-----")
    if icode != None:
        print("**** Molecule No.{} ****".format(icode))
    for a in st_mol.atoms:
        print("{0}    {1}".format(a.name, a.rvec.dictvec))
    if icode is None:
        filename = 'molecule_results/{0}_av_{1}.xyz'.format(file_name,int(time.time()))
    else:
        filename = 'molecule_results/{0}_av_{1}_{2}.xyz'.format(file_name,int(time.time()),str(icode))

    if SAVE_XYZ:
        ix.export_xyz(filename, mt.atom_list_to_dict_list(st_mol.atoms))
        print("Results saved to xyz file.")

i = 0
for pl in possible_list:
    i += 1
    save_mols_atom(pl, i)

print("=================")
print("Duration of computaion: {} s.".format(round(t_taken, 3)))
print("Number of combinations to work with: {}".format(8**(len(mol_combr))))
print("Number of plausible results: {}".format(len(possible_list)))
