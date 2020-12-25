#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 00:31:15 2019

@author: qinayan
"""

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.sparse import dia_matrix


def linear_diffusion_2D(nrows,ncols,dx,dy,Dx,Dy):

#------------------------------------------------------------------------
# THE KEY PARAMETERS OF MATRIX A
#------------------------------------------------------------------------
    Sx = Dx/dx**2
    Sy = Dy/dy**2
    Sxy = -2*Sx-2*Sy
    
#    Sxe = Dx/dx/2.0
#    Sye = Dy/dy/2.0
#    
#    Sxm = 1.0/dx
#    Sym = 1.0/dy
    
    Scx = Dx/1.0/dx**2
    Scy = Dy/1.0/dy**2

    Nx = nrows
    Ny = ncols

#------------------------------------------------------------------------
# BUIL MATRIX
#------------------------------------------------------------------------
    idia_u = np.arange(0,Nx*Ny,Ny)
    idia_l = np.arange(Ny-2,Nx*Ny,Ny)
  
    MainD                     = Sxy*np.ones(Nx*Ny)
    Dia_top                   = Sy*np.ones(Ny*(Nx-1))
    Dia_top[0:Ny]             = 2.0*Sy  #*np.ones(Ny)
    Dia_up                    = Sx *np.ones(Nx*Ny-1)
    Dia_up[idia_u]            = 2.0*Sx  #*np.ones(Nx)
    Dia_low                   = Sx *np.ones(Nx*Ny-1)
    Dia_low[idia_l]           = 2.0*Sx  #*np.ones(Nx)
    Dia_bot                   = Sy *np.ones(Ny*(Nx-1))
    Dia_bot[Ny*(Nx-2):Ny*(Nx-1)]  = 2.0*Sy   #*np.ones(Ny)
    
    indZero            = np.arange(Ny-1, Nx*Ny-1, Ny)
    Dia_low[indZero] = 0.0
    Dia_up[indZero] = 0.0 
    
    matrix_A = diags([Dia_bot,Dia_low,MainD,Dia_up,Dia_top], [-Ny,-1,0,1,Ny],shape=(Nx*Ny, Nx*Ny))
    matrix_A = matrix_A.tocsr()
    

#------------------------------------------------------------------------
# BUIL MATRIX for carbon and soil diffusion part Ax and Ay
#------------------------------------------------------------------------

# this is for (i, j+1)-(i,j), Axp1
    DCxp1_up            = Scx*np.ones(Nx*Ny-1)
    DCxp1_up[indZero]   = 0.0 
    indZero2            = np.arange(Ny-2, Nx*Ny-1, Ny)
    DCxp1_low           = np.zeros(Nx*Ny-1)
    DCxp1_low[indZero2] = Scx
    
    matrix_Axp1 = diags([DCxp1_low, -Scx, DCxp1_up], [-1,0,1], shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Axp1 = matrix_Axp1.tocsr()
       
# this is for (i, j)-(i,j-1), Axm1
    DCxm1_low           = -Scx*np.ones(Nx*Ny-1)
    DCxm1_low[indZero]  = 0.0 
    indZero3            = np.arange(0, Nx*Ny-1, Ny)
    DCxm1_up            = np.zeros(Nx*Ny-1)
    DCxm1_up[indZero3]  = -Scx 
    matrix_Axm1 = diags([DCxm1_low, Scx, DCxm1_up], [-1,0,1], shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Axm1 = matrix_Axm1.tocsr()
    
# this is for (i+1, j)-(i,j), Ayp1
    DCyp1_low           = np.zeros(Nx*Ny-Ny)
    indZero4            = np.arange(Ny*Nx-2*Ny, Nx*Ny-Ny)
    DCyp1_low[indZero4] = Scy
    matrix_Ayp1 = diags([DCyp1_low,-Scy,Scy], [-Ny,0,Ny],shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Ayp1 = matrix_Ayp1.tocsr()    

# this is for (i, j)-(i-1,j), Aym1
    DCym1_up            = np.zeros(Nx*Ny-Ny)
    indZero5            = np.arange(0, Ny)
    DCym1_up[indZero5]  = -Scy
    matrix_Aym1 = diags([-Scy,Scy,DCym1_up], [-Ny,0,Ny],shape=(Nx*Ny, Nx*Ny))#.toarray()
    matrix_Aym1 = matrix_Aym1.tocsr()    

    
    ############################################################
    ##### Build the boundary condition 
    ############################################################
    #------------------------------------------------------------------------
    # for landscape diffusion 
    #------------------------------------------------------------------------
    Bl     = 0.0 # -0.01/365.0
    Br     = 0.0 # 0.01/365.0
    Bu     = -2.5/365.0
    Bd     = 2.5/365.0

    matrix_BC = np.zeros(Nx*Ny)
    
    matrix_BC[0:Ny]        =  2*dy*Sy*Bu
    matrix_BC[Nx*Ny-Ny:]   = -2*dy*Sy*Bd
    matrix_BC[0:Nx*Ny:Ny]=matrix_BC[0:Nx*Ny:Ny]+2*dx*Sx*Bl
    matrix_BC[Ny-1:Nx*Ny:Ny]=matrix_BC[Ny-1:Nx*Ny:Ny]-2*dx*Sx*Br
#    
#    matrix_BC = np.zeros(Nx*Ny)
#    
#    matrix_BC[0:Ny-1]=-2*dy*Sy*Cu
#    matrix_BC[0]=matrix_BC[0]-2*dx*Sx*Cl
#    matrix_BC[Ny-1]=matrix_BC[Ny-1]-2*dx*Sx*Cr
#    
#    matrix_BC[Ny]=-2*dx*Sx*Cl
#    matrix_BC[Nx*Ny-1-Ny]=-2*dx*Sx*Cr
#    
#    matrix_BC[Nx*Ny-Ny:Nx*Ny-1]=-2*dy*Sy*Cd
#    matrix_BC[Nx*Ny-Ny]=matrix_BC[Nx*Ny-Ny]-2*dx*Sx*Cl
#    matrix_BC[Nx*Ny-1]=matrix_BC[Nx*Ny-1]-2*dx*Sx*Cr
#    
    #------------------------------------------------------------------------
    # for SOC associated with landscape diffusion
    #------------------------------------------------------------------------
    BCx_p1                = np.zeros(Nx*Ny)
    BCx_p1[Ny-1:Nx*Ny:Ny] = BCx_p1[Ny-1:Nx*Ny:Ny]-2*dx*Scx*Br
    BCx_m1                = np.zeros(Nx*Ny)
    BCx_m1[0:Nx*Ny:Ny]    = BCx_m1[0:Nx*Ny:Ny]   -2*dx*Scx*Bl
    BCy_p1                = np.zeros(Nx*Ny)
    BCy_p1[Nx*Ny-Ny:]     = BCy_p1[Nx*Ny-Ny:]    -2*dy*Scy*Bd
    BCy_m1                = np.zeros(Nx*Ny)
    BCy_m1[0:Ny]          = BCy_m1[0:Ny]         -2*dy*Scy*Bu
    
    return matrix_A, matrix_BC, matrix_Axp1, matrix_Axm1, matrix_Ayp1, matrix_Aym1,BCx_p1, BCx_m1, BCy_p1,BCy_m1
