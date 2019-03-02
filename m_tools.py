# Code by Jerry Yan

import numpy as np
import itertools

INDEX_RANGE = range(0,3)

BOND_LENGTHS = {
    "CC1": 1.54,
    "CC2": 1.34,
    "CC3": 1.20,
    "CO1": 1.43,
    "OO1": 1.48,
    "OO2": 1.21
}

# Bond (between two atoms)
class chem_bond(object):
    def __init__(self, atom_1, atom_2, bond_order = 1):
        self.atoms = [atom_1, atom_2]
        self.bond_order = bond_order

    @property
    def bdcode(self):
        code = ""
        a_st = sorted(self.atoms)
        for a in a_st:
            code = code + a
        code = code + str(self.bond_order)
        return code

    def validate(self):
        try:
            l = BOND_LENGTHS[self.bdcode]
        except KeyError:
            return False
        except Exception as e:
            print(e)
            return False
        else:
            return True

    @property
    def length(self):
        if self.validate() is True:
            return BOND_LENGTHS[self.bdcode]
        else:
            return None


# Resign the rvec based on |x|, |y|, |z|
def re_sign_n(rvec):
    pl = []

    SLIST = list(map(list, itertools.product([-1,1], repeat=len(INDEX_RANGE))))

    for s in SLIST:
        r = rvec
        nr = [r[i]*s[i] for i in INDEX_RANGE]
        pl.append(nr)

    return pl

# Resign the rvec for a list of molecules
def find_possibles(mols):
    possibles = []
    for i in range(0, len(mols)):
        m = mols[i]
        pl = re_sign_n(m['rvec'])
        for p in pl:
            d = {'rvec': p, 'name': m['name'], 'index': i}
            possibles.append(d)
    return possibles

# Calculate the distance between two points
def distance(rvec_1, rvec_2):
    dvec = np.array(rvec_1) - np.array(rvec_2)
    d = np.sqrt(np.dot(dvec, dvec))
    return d

# Filtering by the bond length
def bond_length_filter(rvec_1, rvec_2, b_length, tol_range = 0.1):
    d = distance(rvec_1, rvec_2)
    if (b_length - tol_range) < d < (b_length + tol_range):
        return True
    else:
        return False

# Error Message for bond length filters
def bond_length_error(af, name, ex = True):
    print("--Unable to narrow down to one {} by bond-length filtering.--".format(name))
    print("The most narrow results for {0} are \n {1}".format(name, str(af)))
    if ex:
        exit()

# TODO: Find possible bonds between two atoms
