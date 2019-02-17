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

# I tensor for multiple point masses
def i_mpm(m_list, rvec_list):
    try:
        m_count = len(m_list)
        r_count = len(rvec_list)
        if m_count != r_count:
            raise Exception("wrong lists")
        else:
            count = m_count
    except Exception as e:
        print(e)

    i_tens = i_init()

    lrange = range(0, count)

    #Calculate I tensor
    for i in INDEX_RANGE:
        for j in INDEX_RANGE:
            if i == j:
                i_tens[i][j] = np.sum([m_list[n]*(np.sum([rvec_list[n][k]**2 for k in filter(lambda l: l!= i, INDEX_RANGE)])) for n in lrange])
            else:
                i_tens[i][j] = np.sum([m_list[n]*rvec_list[n][i]*rvec_list[n][j] for n in lrange])

    return i_tens

# Calculate ABC from I (default in MHz)
def abc_cgs_from_i(i_tens, toMHz=True):
    abc = []

    for i in INDEX_RANGE:
        b = H/(8*(PI**2)*C*i_tens[i][i])
        if toMHz:
            b = b * C * (10**(-6))
        abc.append(b)

    abc.sort(reverse=True)
    return abc
