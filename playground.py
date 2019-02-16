import numpy as np
import itertools

print(list(itertools.product([-1,1],repeat=3)))

# mC = 12.0
# rC_vec = [-2.12929, -2.15420, 2.95186]
#
# INDEX_RANGE = range(0,3)
#
# # Initialize I tensor
# i_tens = []
# for i in INDEX_RANGE:
#     i_tens_row = []
#     for j in INDEX_RANGE:
#         i_tens_row.append(0)
#     i_tens.append(i_tens_row)
#     i_tens_row = []
#
# # Calculate I tensor
# for i in INDEX_RANGE:
#     for j in INDEX_RANGE:
#         if i == j:
#             i_tens[i][j] = mC*(np.sum([rC_vec[k] for k in filter(lambda l: l!= i, INDEX_RANGE)]))
#         else:
#             i_tens[i][j] = mC*rC_vec[i]*rC_vec[j]
#
# print(i_tens)

# with open('molecule_models/propanediol-1.xyz') as f:
#     i = 0
#     molecules = []
#     for l in f:
#         if i == 0:  # First line: number of molecules
#             try:
#                 xyz_count = int(l)
#             except ValueError:
#                 print("The first line is not an integer.")
#         elif i == 1: # Second line: comments
#             xyz_comment = str(l)
#         else: # Remaining lines: molecules
#             s = l.split()
#             try:
#                 m_name = str(s[0])
#                 m_rvec = [float(s[i]) for i in range(1,4)]
#                 m = {'name': m_name, 'rvec': m_rvec}
#                 molecules.append(m)
#             except ValueError:
#                 print("Value error encountered in line {}".format(i+1))
#         i += 1
#
#     if len(molecules) != xyz_count:
#         print("The xyz count does not match the actual molecules count.")
