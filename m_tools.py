# Code by Jerry Yan

import numpy as np
import itertools
import copy

INDEX_RANGE = range(0,3)

BOND_LENGTHS = {
    "CC1": 1.54,
    "CC2": 1.34,
    "CC3": 1.20,
    "CO1": 1.43,
    "OO1": 1.48,
    "OO2": 1.21
}

ATOM_ABS_VALENCES = {
    "C": 4,
    "O": 2
}

# Bond Type (between two atoms)
class Chem_bond_type(object):
    def __init__(self, atom_1, atom_2, bond_order = 1):
        self.atoms = [atom_1, atom_2]
        self.bond_order = bond_order

    def __eq__(self, other):
        if not isinstance(other, Chem_bond_type):
            return NotImplemented

        return self.bdcode == other.bdcode

    def __hash__(self):
        return hash(self.bdcode)

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

# Atom
class Atom(object):
    def __init__(self, name, rvec = None):
        self.name = name
        self.rvec = rvec

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return NotImplemented

        return self.name == other.name and self.rvec == other.rvec

    def __hash__(self):
        return hash((self.name, self.rvec))

# Bond
class Chem_bond(object):
    def __init__(self, atom_1, atom_2, bond_type):
        self.atoms = {atom_1, atom_2}
        self.type = bond_type

    def __eq__(self, other):
        if not isinstance(other, Chem_bond):
            return NotImplemented
        return self.atoms == other.atoms and self.type == other.type

    def __hash__(self):
        return hash((frozenset(self.atoms), self.type))

    # Real distance between the two atoms (real bond length)
    @property
    def diameter(self):
        at_l = list(self.atoms)
        atom_1 = at_l[0]
        atom_2 = at_l[1]
        return distance_atom(atom_1, atom_2)

# Bond Graph
class Chem_bond_graph(object):
    def __init__(self, bonds = set()):
        self.bonds = bonds

    def __eq__(self, other):
        if not isinstance(other, Chem_bond_graph):
            return NotImplemented

        return self.bonds == other.bonds

    def __hash__(self):
        return hash(frozenset(self.bonds))

    def atom_degree(self, atom):
        deg = 0
        for b in self.bonds:
            if atom in b.atoms:
                deg += 1
        return deg

    def print_bdcodes(self):
        print([bond.type.bdcode for bond in self.bonds])

# Structural Molecule
class Strc_molecule(object):
    def __init__(self, atoms = set(), bond_graphs = set()):
        self.atoms = atoms
        self.bond_graphs = bond_graphs

    @property
    def size(self):
        return len(self.atoms)

    def __eq__(self, other):
        if not isinstance(other, Strc_molecule):
            return NotImplemented

        return self.atoms == other.atoms and self.bond_graphs == other.bond_graphs

    def add_atom(self, atom):
        self.atoms.add(atom)

# Position vector
class Vector(object):
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        if not isinstance(other, Vector):
            return NotImplemented

        return (self.x == other.x and self.y == other.y and self.z == other.z)

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    @property
    def dictvec(self):
        return [self.x, self.y, self.z]

    @dictvec.setter
    def dictvec(self, dict):
        self.x = dict[0]
        self.y = dict[1]
        self.z = dict[2]

# Resign the rvec based on |x|, |y|, |z|
def re_sign_n(rvec):
    pl = []

    SLIST = list(map(list, itertools.product([-1,1], repeat=len(INDEX_RANGE))))

    for s in SLIST:
        r = rvec
        nr = [r[i]*s[i] for i in INDEX_RANGE]
        pl.append(nr)

    return pl

# Resign the rvec based on |x|, |y|, |z| (Atom version)
def re_sign_n_atom(rvec):
    pl = []

    SLIST = list(map(list, itertools.product([-1,1], repeat=len(INDEX_RANGE))))

    for s in SLIST:
        r = rvec.dictvec
        nr = [r[i]*s[i] for i in INDEX_RANGE]
        n_vec = Vector()
        n_vec.dictvec = nr
        pl.append(n_vec)

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

# Resign the rvec for a list of atoms (Atom version)
def find_possibles_atom(atoms):
    possibles = []
    for a in atoms:
        pl = re_sign_n_atom(a.rvec)
        for p in pl:
            atm = Atom(a.name, p)
            possibles.append(atm)
    return possibles

# Calculate the distance between two points
def distance(rvec_1, rvec_2):
    dvec = np.array(rvec_1) - np.array(rvec_2)
    d = np.sqrt(np.dot(dvec, dvec))
    return d

# Calculate the distance between two atoms
def distance_atom(atom_1, atom_2):
    return distance(atom_1.rvec.dictvec, atom_2.rvec.dictvec)

# Filtering by the bond length
def bond_length_filter(rvec_1, rvec_2, b_length, tol_range = 0.1):
    d = distance(rvec_1, rvec_2)
    if (b_length - tol_range) < d < (b_length + tol_range):
        return True
    else:
        return False

# Filtering by bond length (atom version)
def bondtype_length_filter_atom(atom_1, atom_2, bondtype, tol_range = 0.1):
    return bond_length_filter(atom_1.rvec.dictvec, atom_2.rvec.dictvec, bondtype.length, tol_range)

# Error Message for bond length filters
def bond_length_error(af, name, ex = True):
    print("--Unable to narrow down to one {} by bond-length filtering.--".format(name))
    print("The most narrow results for {0} are \n {1}".format(name, str(af)))
    if ex:
        exit()

# Find possible bonds between two atom name
def possible_bonds(name_1, name_2):
    p_range = range(1, 4)
    p_bonds = []
    for ord in p_range:
        bond = Chem_bond_type(name_1, name_2, ord)
        if bond.validate() == True:
            p_bonds.append(bond)
    return p_bonds

# Find possible bonds between two atom name (Atom version)
def possible_bondtypes(atom_1, atom_2):
    return possible_bonds(atom_1.name, atom_2.name)

# Use bond length filter to verify if an atom fits in a validated group of atoms
def bond_length_atom_finder(validated_mols, m, tol_range = 0.1):
    if len(validated_mols) <= 0:
        return True
    else:
        b_pass = False
        for vm in validated_mols:
            possible_bs = possible_bonds(vm['name'], m['name'])
            for b in possible_bs:
                if bond_length_filter(vm['rvec'], m['rvec'], b.length, tol_range) is True:
                    v_mols_n = [vmol for vmol in validated_mols if vmol != vm]
                    d_pass = True
                    for vm2 in v_mols_n:
                        if distance(vm2['rvec'], m['rvec']) < (b.length - tol_range):
                            d_pass = False
                    if d_pass is True:
                        b_pass = True
        return b_pass

# (Atom version) Molecule constructor from an atom [More sophisicated]
def bond_length_molecule_constructer(st_mol, atom, tol_range = 0.1):
    m = copy.deepcopy(st_mol)
    bgs = copy.deepcopy(m.bond_graphs)
    if m.size <= 0:
        m.add_atom(atom)
    else:
        b_pass = False
        m.bond_graphs.clear()
        for vm in st_mol.atoms:
            possible_bts = possible_bondtypes(vm, atom)
            for bt in possible_bts:
                if bondtype_length_filter_atom(vm, atom, bt, tol_range) is True:
                    v_n = [a for a in m.atoms if a != vm]
                    d_pass = True
                    for vm2 in v_n:
                        if distance_atom(vm2, atom) < (bt.length - tol_range):
                            d_pass = False
                    if d_pass is True:
                        p_bond = Chem_bond(vm, atom, bt)
                        m.add_atom(atom)
                        if st_mol.size == 1:
                            m.bond_graphs.add(Chem_bond_graph({p_bond}))
                        elif st_mol.size > 1:
                            for bond_graph in bgs:
                                pbg = copy.copy(bond_graph)
                                pbg.bonds.add(p_bond)
                                m.bond_graphs.add(pbg)

    return m


# Convert the dict type to atom type
def dict_to_atom(mol):
    name = mol['name']
    try:
        rvec = Vector()
        rvec.dictvec = mol['rvec']
    except Exception as e:
        rvec = None
    at = Atom(name, rvec)
    return at

# Convert dict list to a list of atoms
def dict_list_to_atom_list(mols):
    ats = []
    for m in mols:
        ats.append(dict_to_atom(m))
    return ats
