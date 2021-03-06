# Code by Jerry Yan

# Import from XYZ file
def import_xyz(file):
    with open(file) as f:
        i = 0
        molecules = []
        for l in f:
            if i == 0:  # First line: number of molecules
                try:
                    xyz_count = int(l)
                except ValueError:
                    print("The first line is not an integer.")
            elif i == 1: # Second line: comments
                xyz_comment = str(l)
            else: # Remaining lines: molecules
                s = l.split()
                try:
                    m_name = str(s[0])
                    m_rvec = [float(s[i]) for i in range(1,4)]
                    m = {'name': m_name, 'rvec': m_rvec}
                    molecules.append(m)
                except ValueError:
                    print("Value error encountered in line {}".format(i+2))
            i += 1

        if len(molecules) != xyz_count:
            print("The xyz count does not match the actual molecules count.")

    return xyz_count, xyz_comment, molecules

# Export XYZ file
def export_xyz(file, mols, comment=""):
    with open(file, "w+") as f:
        f.write(str(len(mols)))
        f.write("\n")

        f.write(comment + "\n")

        for m in mols:
            try:
                f.write("{0}       {1}       {2}       {3}\n".format(m['name'], *m['rvec']))
            except Exception as e:
                print(e)
            else:
                continue

    print("xyz file saved.")

# Select molecule from list by name
def select_by_name(mol_list, name):
    sl = []
    try:
        for m in mol_list:
            if m['name'] == name:
                sl.append(m)
    except Exception as e:
        print(e)
    return sl

# CGS unit mass by name
def mass_cgs_by_name(name):
    AUDICT = {'C': 12.0, 'O': 15.999, 'H':1.00782}
    AU = 1.66053904*(10**(-24))
    try:
        m = AUDICT[name] * AU
    except Exception as e:
        print(e)
    else:
        return m

# CGS unit length by rvec
def length_cgs_by_rvec(rvec):
    A = 10**(-8)
    rvec = [r * A for r in rvec]
    return rvec

# Save a molecule to xyz (suffix with unix time)
def save_mols(mols, icode = None):
    print("-----The structure of molecule is found as follows:-----")
    for m in mols:
        print("{0}    {1}".format(m['name'], m['rvec']))
    if icode is None:
        filename = 'molecule_results/{0}_{1}.xyz'.format(file_name,int(time.time()))
    else:
        filename = 'molecule_results/{0}_{1}_{2}.xyz'.format(file_name,int(time.time()),str(icode))
    export_xyz(filename, mols)
    print("Results saved to xyz file.")
