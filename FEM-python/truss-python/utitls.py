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
def penalty_assembly():
    """
    Assemble element stiffness matrix and force vector using penalty method.
    
    Args:
        e   : (int) Element number
        ke  : (numpy(nen*ndof,nen*ndof)) element stiffness matrix 
    """
    # penalty method: add a large value to the diagonal of the stiffness matrix
    # to enforce the essential boundary condition
    # the positions of the value are positioned at the bottom of the stiffness matrix 
    penalty = 1e15  # large value for penalty
    for i in range(model.nd):
        i = 2*model.nnp - model.nd + i
        model.K[i, i] += penalty
        model.f[i] += penalty * model.d[i]  # adjust the force vector accordingly
    # print("Stiffness matrix K after penalty assembly:\n", model.K)
    # print("Force vector f after penalty assembly:\n", model.f)

def solvedr_penalty():
    model.d = np.linalg.solve(model.K, model.f)
    print('\nsolution d =\n', model.d)
    return model.d

def solvedr():
    """
    Partition and solve the system of equations
        
    Returns:
        f_E : (numpy.array(nd,1)) Reaction force vector
    """
    nd = model.nd; neq=model.neq
    K_E = model.K[0:nd, 0:nd]
    K_F = model.K[nd:neq, nd:neq]
    K_EF = model. K[0:nd, nd:neq]
    f_F = model.f[nd:neq]
    d_E = model.d[0:nd]
    
    # solve for d_F
    d_F = np.linalg.solve(K_F, f_F - K_EF.T @ d_E) 

    # reconstruct the global displacement d
    model.d = np.append(d_E,d_F)
    
    # compute the reaction r
    f_E = K_E@d_E + K_EF@d_F
    
    # write to the workspace
    print('\nsolution d');  print(model.d)
    print('\nreaction f =', f_E)
    
    return f_E