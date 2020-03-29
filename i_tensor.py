# Code by Jerry Yan

import numpy as np

INDEX_RANGE = range(0,3)

PI = np.pi
C = 2.99792458 * (10**10)
H = 6.6260755 * (10**(-27))

# Initialize I tensor
def i_init():
    i_tens = []
    for i in INDEX_RANGE:
        i_tens_row = []
        for j in INDEX_RANGE:
            i_tens_row.append(0)
        i_tens.append(i_tens_row)
        i_tens_row = []
    return i_tens

# I tensor for single point mass
def i_single_pm(m, r_vec):

    i_tens = i_init()

    # Calculate I tensor
    for i in INDEX_RANGE:
        for j in INDEX_RANGE:
            if i == j:
                i_tens[i][j] = m*(np.sum([r_vec[k]**2 for k in filter(lambda l: l!= i, INDEX_RANGE)]))
            else:
                i_tens[i][j] = m*r_vec[i]*r_vec[j]

    return i_tens

def guard_list_len(list_a, list_b):
    try:
        if len(list_a) != len(list_b):
            raise Exception("List lengths don't match.")
        else:
            return len(list_a)
    except Exception as e:
        print(e)

# Center of Mass for multiple point masses
def cm_mpm(m_list, rvec_list):
    count = guard_list_len(m_list, rvec_list)
    lrange = range(0, count)

    m_tot = sum(m_list)
    t_vec = [0, 0, 0]

    for n in lrange:
        rvec = rvec_list[n]
        m = m_list[n]
        for i in INDEX_RANGE:
            t_vec[i] = t_vec[i] + rvec[i] * m

    cm_vec = [c / m_tot for c in t_vec]

    return cm_vec


# I tensor for multiple point masses
def i_mpm(m_list, rvec_list, origin = None):
    count = guard_list_len(m_list, rvec_list)
    lrange = range(0, count)
    i_tens = i_init()

    if origin is None:
        cm_vec = cm_mpm(m_list, rvec_list)
        r_list = [np.array(rvec) - np.array(cm_vec) for rvec in rvec_list]
    else:
        try:
            r_list = [np.array(rvec) - np.array(origin) for rvec in rvec_list]
        except Exception as e:
            print(e)
        else:
            pass

    #Calculate I tensor
    for i in INDEX_RANGE:
        for j in INDEX_RANGE:
            if i == j:
                i_tens[i][j] = np.sum([m_list[n]*(np.sum([r_list[n][k]**2 for k in filter(lambda l: l!= i, INDEX_RANGE)])) for n in lrange])
            else:
                i_tens[i][j] = np.sum([m_list[n]*rvec_list[n][i]*r_list[n][j] for n in lrange])

    return i_tens

# Calculate ABC from I (default in MHz)
def abc_cgs_from_i(i_tens, toMHz=True):
    abc = []

    for i in INDEX_RANGE:
        b = H/(8*(PI**2)*i_tens[i][i])
        if toMHz:
            b = b * (10**(-6))
        abc.append(b)

    abc.sort(reverse=True)
    return abc
