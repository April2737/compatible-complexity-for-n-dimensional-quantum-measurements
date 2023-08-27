#   description: Find the lower bounds of R^N_3(∞) and R^N_4(∞) using canonical parent POVMs. 
#   
#   authors: Jiaxuan Zhang, Eric Chitambar
#
#   requires: This Python code does not require additional tools
#
#   last update: August, 2023

# =====================================

import numpy as np
from numpy import linalg as LA

from itertools import combinations
from itertools import product

import sys

# Give parameters =============================================

d_POVM = 4      # dimension of quantum states and POVMs

N = 15       # number of children POVMs simulated

eta = 0.114727    # We want to find the biggest eta that makes all parent elements positive

decimal_point = 6


j = complex(0,1)


# generators for qutrit POVMs ========================================


X1 = np.array([[0,1,0],[1,0,0],[0,0,0]])
X2 = np.array([[0,0,1],[0,0,0],[1,0,0]])
X3 = np.array([[0,0,0],[0,0,1],[0,1,0]])
Y1 = np.array([[0,-j,0],[j,0,0],[0,0,0]])
Y2 = np.array([[0,0,-j],[0,0,0],[j,0,0]])
Y3 = np.array([[0,0,0],[0,0,-j],[0,j,0]])
Z1 = np.array([[1,0,0],[0,-1,0],[0,0,0]])
Z2 = np.array([[1,0,0],[0,1,0],[0,0,-2]])*(1/3)**.5

if d_POVM == 3:
  Lambda_matrices = [X1,X2,X3,Y1,Y2,Y3,Z1,Z2]
  coef_lambda = 3/2     # the coef. before lambda*\hat{lambda} in the formula of parent element


# generators for 4-dimensional POVMs ================================================================================

M_u12 = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])    # Matrix u12
M_u13 = np.array([[0,0,1,0],[0,0,0,0],[1,0,0,0],[0,0,0,0]])
M_u14 = np.array([[0,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,0]])
M_u23 = np.array([[0,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,0]])
M_u24 = np.array([[0,0,0,0],[0,0,0,1],[0,0,0,0],[0,1,0,0]])
M_u34 = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,1],[0,0,1,0]])

M_v12 = np.array([[0,-j,0,0],[j,0,0,0],[0,0,0,0],[0,0,0,0]])
M_v13 = np.array([[0,0,-j,0],[0,0,0,0],[j,0,0,0],[0,0,0,0]])
M_v14 = np.array([[0,0,0,-j],[0,0,0,0],[0,0,0,0],[j,0,0,0]])
M_v23 = np.array([[0,0,0,0],[0,0,-j,0],[0,j,0,0],[0,0,0,0]])
M_v24 = np.array([[0,0,0,0],[0,0,0,-j],[0,0,0,0],[0,j,0,0]])
M_v34 = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,-j],[0,0,j,0]])

M_w1 = np.array([[1,0,0,0],[0,-1,0,0],[0,0,0,0],[0,0,0,0]])
M_w2 = np.array([[1,0,0,0],[0,1,0,0],[0,0,-2,0],[0,0,0,0]])*(1/3**.5)
M_w3 = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-3]])*(1/6**.5)

if d_POVM == 4:
  Lambda_matrices = [M_u12, M_u13, M_u14, M_u23, M_u24, M_u34,   M_v12, M_v13, M_v14, M_v23, M_v24, M_v34,   M_w1, M_w2, M_w3]
  coef_lambda = 2


# Create the set of all combinations of N matrices chosen from the set of all generators =============================================================

chooseMatrices = list(combinations(Lambda_matrices, N))   


if chooseMatrices == []:      # Make sure N <= d_POVM
  print("chooseMatrices=",chooseMatrices)
  sys.exit()    # End the program and do not restart


# Create a list for the "+/-" sign in the formula of parent element =====================================

zero_one_list = list(product(range(2), repeat=N))     # a list of length 2^N



# Loops to find the maximal eta makes all parent element positive =======================================================================================================

all_parent_pos = 0

while all_parent_pos == 0:

  pos_combinations = 0           # Calculate the number of combinations that make corresponding parent elements positive

  M_parent_element = [0]*N                # Create a set for the N matrices used in parent element for each combination


  for j in range(len(chooseMatrices)):    # A Loop checks all combinations of N matrices

    for i_M_parent_element in range(N):
      M_parent_element[i_M_parent_element] = chooseMatrices[j][i_M_parent_element]


    pos_elements = 0    # Calculate the number of positive parent POVM elements in one parent (one combination)


    for i in range(2**N):      # Create parent element one by one
      parent_element = np.identity(d_POVM)   # We drop the overall coeffient (1/(2^N)) for parent element since it does not affect positivity

      for k in range(N):
        parent_element = parent_element + coef_lambda * eta * (zero_one_list[i][k]*2-1) * M_parent_element[k]

      w, v = LA.eig(parent_element)   # Calculate the eigenvalues and eigenvectors of each parent element

      if w[0] < -1*10**(-10) or w[1] < -1*10**(-10) or w[2] < -1*10**(-10) :     # Check if all eigenvalues of each parent element are positive
        break
      else:
        pos_elements += 1

    if pos_elements == 2**N:
      pos_combinations += 1

  print()

  if pos_combinations == len(chooseMatrices):
    print('This eta is VALID for our simulation:')
    all_parent_pos = 1
    print()
    print('eta =',eta)
  else:
    print('eta =',eta,'is NOT valid for our simulation.',pos_combinations,'parent POVMs are positive in all',len(chooseMatrices),'parent POVMs (combinations of choosing',N,'lambda matrices).')
    eta = eta - 1*10**(-1*decimal_point)
