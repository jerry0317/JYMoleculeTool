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

CO_BOND_LENGTH = 1.43
CC_BOND_LENGTH = 1.54

TOLERANCE_LEVEL = 0.1 # In bond-length filters, an atom passes the filter when the distance lies between BOND_LENGTH plus/minus TOLERANCE_LEVEL

# Fix the first atom
A1 = mol_comb[0]
print("The first atom has been fixed.")

mol_combr = [m for m in mol_comb if m != A1]

def rcs_filter(a_r, base_group, tol_range = 0.1):
    possible_a = mt.find_possibles([a_r]) # Find all possible rvec for a_r

    possible_af1 = filter(lambda m: mt.bond_length_atom_finder(base_group, m, tol_range), possible_a)

    return possible_af1

possible_list = []

num_enum = 0
num_pn = 0

def rcs_action(mol_r, p_a, b_s, tol_range = 0.1):
    if len(mol_r) > 0:
        for p in p_a: # Possible atom
            for a_r in mol_r: # Remaining atoms
                b_sf = b_s + [p]

                p_af1 = rcs_filter(a_r, b_sf, tol_range)

                if len(p_af1) > 0:
                    mol_rf1 = [m for m in mol_r if m != a_r]
                    rcs_action(mol_rf1, p_af1, b_sf, tol_range)
    elif len(mol_r) == 0:
        for p in p_a:
            mols = b_s + [p]

            duplicated = False
            rvec_set = set([tuple(m['rvec']) for m in mols])
            for pl in possible_list:
                pl_rvec_set = set([tuple(m['rvec']) for m in pl])
                if rvec_set == pl_rvec_set:
                    duplicated = True
                    continue
            if duplicated is False:
                possible_list.append(mols)
                print("The possible No. {} is found.".format(len(possible_list)))

t0 = time.clock()

rcs_action(mol_combr, [A1], [], TOLERANCE_LEVEL)

t_taken = time.clock() - t0

print("-----The structure of molecule is found as follows:-----")

def save_mols(mols, icode = None):
    if icode != None:
        print("**** Molecule No.{} ****".format(icode))
    for m in mols:
        print("{0}    {1}".format(m['name'], m['rvec']))
    if icode is None:
        filename = 'molecule_results/{0}_{1}.xyz'.format(file_name,int(time.time()))
    else:
        filename = 'molecule_results/{0}_{1}_{2}.xyz'.format(file_name,int(time.time()),str(icode))

    if SAVE_XYZ:
        ix.export_xyz(filename, mols)
        # print("Results saved to xyz file.")

    #
    # rmol_C = ix.select_by_name(mols, 'C')
    # import itertools
    # CC_lengths = []
    # possible_Cpairs = itertools.combinations(rmol_C, 2)
    # for p in possible_Cpairs:
    #     d = mt.distance(p[0]['rvec'], p[1]['rvec'])
    #     if d < 3:
    #         CC_lengths.append(d)
    # CC_lengths.sort()
    #
    # CC_lengths = [round(c, 3) for c in CC_lengths]
    # print(CC_lengths)

i = 0
for pl in possible_list:
    i += 1
    save_mols(pl, i)

print("=================")
print("Duration of computaion: {} s.".format(round(t_taken, 3)))
print("Number of combinations to work with: {}".format(8**(len(mol_combr))))
print("Number of plausible results: {}".format(len(possible_list)))
