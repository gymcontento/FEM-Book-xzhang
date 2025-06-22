#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide methods to setup LM matrices, create FE model for a truss from a json 
file, to plot the truss, to calculate and print stresses of every element.

Created on Sat May 9 15:43:00 2020

@author: thurcni@163.com, xzhang@tsinghua.edu.cn
"""

import FEData as model
import numpy as np
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def create_model_json(DataFile):
    """ 
    Initialize the FEM model from file DataFile (in json format)
    """

    # input data from json file
    with open(DataFile) as f_obj:
        FEData = json.load(f_obj)
    
    model.Title= FEData['Title']
    model.nsd  = FEData['nsd']
    model.ndof = FEData['ndof']
    model.nnp  = FEData['nnp']
    model.nel  = FEData['nel']
    model.nen  = FEData['nen']    
    model.neq  = model.ndof*model.nnp
    model.nd   = FEData['nd']

    # initialize K, d and f 
    model.f = np.zeros((model.neq,1))            
    model.d = np.zeros((model.neq,1))        
    model.K = np.zeros((model.neq,model.neq))

    # define the mesh
    model.x = np.array(FEData['x'])
    model.y = np.array(FEData['y'])
    if model.nsd == 3:
        model.z = np.array(FEData['z'])  
    model.IEN = np.array(FEData['IEN'], dtype=np.int32)
    # 设置单元节点自由度编号矩阵
    model.LM = np.zeros((model.nen*model.ndof, model.nel), dtype=np.int32)
    set_LM()

    # element and material data (given at the element)
    model.E     = np.array(FEData['E'])
    model.CArea = np.array(FEData['CArea'])
    # 计算每个单元的长度
    if model.nsd == 2:
        model.leng  = np.sqrt(np.power(model.x[model.IEN[:, 1]-1] - 
                                    model.x[model.IEN[:, 0]-1], 2) +
                            np.power(model.y[model.IEN[:, 1]-1] - 
                                    model.y[model.IEN[:, 0]-1], 2))
    if model.nsd == 3:
        model.leng  = np.sqrt(np.power(model.x[model.IEN[:, 1]-1] - 
                                   model.x[model.IEN[:, 0]-1], 2) +
                          np.power(model.y[model.IEN[:, 1]-1] - 
                                   model.y[model.IEN[:, 0]-1], 2) +
                          np.power(model.z[model.IEN[:, 1]-1] - 
                                   model.z[model.IEN[:, 0]-1], 2))
    model.stress= np.zeros((model.nel,))

    # prescribed forces
    fdof = FEData['fdof']
    force= FEData['force']
    # 对f进行初始化
    # 枚举，ind和value分别为fdof的索引和对应的值，ind为0开始的索引，对应fdof的第一个值
    for ind, value in enumerate(fdof):
        model.f[value-1][0] = force[ind]

    # output plots
    model.plot_truss= FEData['plot_truss']
    model.plot_node = FEData['plot_node']
    model.plot_tex  = FEData['plot_tex']
    plottruss()

    # print model parameters
    # print("Location Matrix:\n", model.LM)
    # print("Element Lengths:\n", model.leng)
    # print("Nodal Forces(Initial):\n", model.f)

def set_LM():
    '''
    set up Location Matrix
    '''
    for e in range(model.nel):
        for j in range(model.nen):
            for m in range(model.ndof):
                ind = j*model.ndof + m
                model.LM[ind, e] = model.ndof*(model.IEN[e, j] - 1) + m

def plottruss():
    '''
    plot the truss
    '''
    if model.plot_truss == "yes":
        if model.ndof == 1:
            for i in range(model.nel):
                XX = np.array([model.x[model.IEN[i, 0]-1], 
                               model.x[model.IEN[i, 1]-1]])
                YY = np.array([0.0, 0.0])
                plt.plot(XX, YY, "blue")

                if model.plot_node == "yes":
                    plt.text(XX[0], YY[0], str(model.IEN[i, 0]))
                    plt.text(XX[1], YY[1], str(model.IEN[i, 1]))
        elif model.ndof == 2:
            for i in range(model.nel):
                XX = np.array([model.x[model.IEN[i, 0]-1], 
                               model.x[model.IEN[i, 1]-1]])
                YY = np.array([model.y[model.IEN[i, 0]-1], 
                               model.y[model.IEN[i, 1]-1]])
                plt.plot(XX, YY, "blue")
                d_x_index1 = 2*model.IEN[i, 0] - 2
                d_x_index2 = 2*model.IEN[i, 1] - 2
                d_x = 10000 * np.array([model.d[d_x_index1], 
                               model.d[d_x_index2]])
                d_y_index1 = 2*model.IEN[i, 0] - 1
                d_y_index2 = 2*model.IEN[i, 1] - 1
                d_y = 10000 * np.array([model.d[d_y_index1], 
                                model.d[d_y_index2]])
                XX_d = np.array([model.x[model.IEN[i, 0]-1] + d_x[0], 
                               model.x[model.IEN[i, 1]-1] + d_x[1]])
                YY_d = np.array([model.y[model.IEN[i, 0]-1] + d_y[0], 
                               model.y[model.IEN[i, 1]-1] + d_y[1]])
                plt.plot(XX_d, YY_d, "red")

                if model.plot_node == "yes":
                    plt.text(XX[0], YY[0], str(model.IEN[i, 0]))
                    plt.text(XX[1], YY[1], str(model.IEN[i, 1]))

        elif model.ndof == 3:
            # insert your code here for 3D
            for i in range(model.nel):
                XX = np.array([model.x[model.IEN[i, 0]-1], 
                               model.x[model.IEN[i, 1]-1]])
                YY = np.array([model.y[model.IEN[i, 0]-1], 
                               model.y[model.IEN[i, 1]-1]])
                ZZ = np.array([model.z[model.IEN[i, 0]-1], 
                               model.z[model.IEN[i, 1]-1]])
                
                # Create a 3D plot if it doesn't exist
                # make sure printing in one figure
                if not plt.get_fignums():
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                else:
                    ax = plt.gca()

                ax.plot(XX, YY, ZZ, color="blue")

                # Set equal aspect ratio
                ax.set_box_aspect([1, 1, 1])
                ax.set_xlabel(r"$x$")
                ax.set_ylabel(r"$y$")
                ax.set_zlabel(r"$z$")

                # Add displacement vectors
                d_x_index1 = 3*model.IEN[i, 0] - 3
                d_x_index2 = 3*model.IEN[i, 1] - 3
                d_x = 10000 * np.array([model.d[d_x_index1], 
                               model.d[d_x_index2]])
                d_y_index1 = 3*model.IEN[i, 0] - 2
                d_y_index2 = 3*model.IEN[i, 1] - 2
                d_y = 10000 * np.array([model.d[d_y_index1], 
                                model.d[d_y_index2]])
                d_z_index1 = 3*model.IEN[i, 0] - 1
                d_z_index2 = 3*model.IEN[i, 1] - 1
                d_z = 10000 * np.array([model.d[d_z_index1], 
                                model.d[d_z_index2]])
                
                XX_d = np.array([model.x[model.IEN[i, 0]-1] + d_x[0], 
                               model.x[model.IEN[i, 1]-1] + d_x[1]])
                YY_d = np.array([model.y[model.IEN[i, 0]-1] + d_y[0], 
                               model.y[model.IEN[i, 1]-1] + d_y[1]])
                ZZ_d = np.array([model.z[model.IEN[i, 0]-1] + d_z[0], 
                               model.z[model.IEN[i, 1]-1] + d_z[1]])

                ax.plot(XX_d, YY_d, ZZ_d, "red")

                if model.plot_node == "yes":
                    ax.text(XX[0], YY[0], ZZ[0], str(model.IEN[i, 0]))
                    ax.text(XX[1], YY[1], ZZ[1], str(model.IEN[i, 1]))

        else:
            raise ValueError("The dimension (ndof = {0}) given for the \
                             plottruss is invalid".format(model.ndof))
        
        plt.title("Truss Plot")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")
        plt.savefig("truss.pdf")

        # Convert matplotlib figures into PGFPlots figures stored in a Tikz file, 
        # which can be added into your LaTex source code by "\input{fe_plot.tex}"
        if model.plot_tex == "yes":
            import tikzplotlib
            tikzplotlib.clean_figure()
            tikzplotlib.save("fe_plot.tex")
    
    print("\t2D Truss Params \n")
    print(model.Title + "\n")
    print("No. of Elements  {0}".format(model.nel))
    print("No. of Nodes     {0}".format(model.nnp))
    print("No. of Equations {0}".format(model.neq))


def print_stress():
    '''
    Calculate and print stresses of every element
    '''

    # prints the element number and corresponding stresses
    print("Element\t\t\tStress")
    # Compute stress for each element
    for e in range(model.nel):
        de = model.d[model.LM[:, e]]  # nodal displacements for each element
        const = model.E[e]/model.leng[e]

        if model.ndof == 1:
            model.stress[e] = const*(np.array([-1, 1])@de)
        elif model.ndof == 2:
            IENe = model.IEN[e] - 1
            xe = model.x[IENe]
            ye = model.y[IENe]
            s = (ye[1] - ye[0])/model.leng[e]
            c = (xe[1] - xe[0])/model.leng[e]
            model.stress[e] = const*(np.array([-c, -s, c, s])@de)
        elif model.ndof == 3:
            # insert your code here for 3D
            IENe = model.IEN[e] - 1
            xe = model.x[IENe]
            ye = model.y[IENe]
            ze = model.z[IENe]
            s = (ye[1] - ye[0])/model.leng[e]
            c = (xe[1] - xe[0])/model.leng[e]
            c_phi = (ze[1] - ze[0])/model.leng[e]
            s_phi = np.sqrt(1 - c_phi*c_phi)

            coe_1 = c*s_phi
            coe_2 = s*s_phi
            coe_3 = c_phi

            model.stress[e] = const * (np.array([-coe_1, -coe_2, -coe_3, 
                                                coe_1, coe_2, coe_3]) @ de)

            pass # delete or comment this line after your implementation for 3D
        else:
            raise ValueError("The dimension (ndof = {0}) given for the \
                             problem is invalid".format(model.ndof))

        print("{0}\t\t\t{1}".format(e+1, model.stress[e]))
        