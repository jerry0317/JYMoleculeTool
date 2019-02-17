# Code by Jerry Yan

import numpy as np
import itertools

INDEX_RANGE = range(0,3)

def re_sign_n(rvec):
    pl = []

    SLIST = list(map(list, itertools.product([-1,1], repeat=len(INDEX_RANGE))))

    for s in SLIST:
        r = rvec
        nr = [r[i]*s[i] for i in INDEX_RANGE]
        pl.append(nr)

    return pl

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

def bond_length_error(af, name, ex = True):
    print("--Unable to narrow down to one {} by bond-length filtering.--".format(name))
    print("The most narrow results for {0} are \n {1}".format(name, str(af)))
    if ex:
        exit()
