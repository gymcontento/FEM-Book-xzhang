#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provides utilities used by FE analysis.
  1. assembly: Global stiffness matrix assembly.
  2. solvedr: Solving the stiffness equations by the reduction approach.

Created on Sat May 9 17:39:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import numpy as np
import FEData as model


def assembly(e, ke):
    """
    Assemble element stiffness matrix.
    
    Args:
        e   : (int) Element number
        ke  : (numpy(nen*ndof,nen*ndof)) element stiffness matrix
    """
    model.K[np.ix_(model.LM[:,e], model.LM[:,e])] += ke

# using penalty method to assemble stiffness matrix and force vector
def penalty_solvedr():
    """
    Assemble element stiffness matrix and force vector using penalty method.
    """
    # penalty method: add a large value to the diagonal of the stiffness matrix
    # to enforce the essential boundary condition
    # the positions of the value are positioned at the bottom of the stiffness matrix
    epsilon = 1e-10  # small value for epsilon
    penalty = 1e15  # large value for penalty
    K_P = model.K
    f_P = model.f
    nd_start = model.neq - model.nd 
    nd_end = model.neq 

    # add small value to the diagonal of K_P to avoid singularity
    # when calculating 2D problems by using 3D elements
    for i in range(K_P.shape[0]):
        if np.allclose(K_P[i, :], 0):
            K_P[i, i] += epsilon

    for i in range(model.nd):
        i = 3*model.nnp - model.nd + i
        K_P[i, i] += penalty
        f_P[i] += penalty * model.d[i]  # adjust the force vector accordingly

    # solve the system of equations
    d_P = np.linalg.solve(K_P, f_P)

    # calculate the reaction force vector
    f_R = -penalty * d_P[nd_start:nd_end]
    print('\nReaction force vector f_R =\n', f_R)

    model.d = d_P  # update the global displacement vector

    return d_P

def solvedr():
    """
    Partition and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    nd_start = model.neq - model.nd 
    nd_end = model.neq 

    epsilon = 1e-10  # small value for epsilon
    # add small value to the diagonal of K to avoid singularity
    K = model.K
    for i in range(K.shape[0]):
        if np.allclose(K[i, :], 0):
            K[i, i] += epsilon

    K_E = K[nd_start:nd_end, nd_start:nd_end]
    K_F = K[0:nd_start, 0:nd_start]
    K_EF = K[nd_start:nd_end, 0:nd_start]
    f_F = model.f[0:nd_start]
    d_E = model.d[nd_start:nd_end]

    # solve for d_F
    d_F = np.linalg.solve(K_F, f_F - K_EF.T @ d_E) 

    # reconstruct the global displacement d
    model.d = np.append(d_F,d_E)
    
    # compute the reaction r
    f_E = K_E@d_E + K_EF@d_F
    
    # write to the workspace
    print('\nsolution d with partition method =\n', model.d)
    print('\nreaction f =', f_E)
    
    return f_E