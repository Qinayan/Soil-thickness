#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 13:47:51 2019

@author: qinayan
"""

import numpy as np
from scipy.sparse import diags


def matrix_for_rolling(nrows,ncols,dx,dy,dts):
    Nx = ncols
    My = nrows
# H_i+1,j, south, pN
    diag_pN = np.zeros(Nx*My,dtype = 'int')
    diag_pN[Nx*My-Nx :] = 1
    
    A_pN = diags([diag_pN, 1], [0,Nx],shape=(Nx*My, Nx*My))
    A_pN = A_pN.tocsr()
    
    # H_i-1,j, north value,mN
    diag_mN = np.zeros(Nx*My,dtype = 'int')
    diag_mN[0:Nx] = 1
    A_mN = diags([1,diag_mN], [-Nx,0],shape=(Nx*My, Nx*My))
    A_mN = A_mN.tocsr()
    
    # H_i,j+1, East value, p1
    diag_p1_u = np.ones(Nx*My-1,dtype = 'int')
    diag_p1_u[Nx-1:Nx*My-1:Nx] = 0
    
    diag_p1_l = np.zeros(Nx*My,dtype = 'int')
    diag_p1_l[Nx-1:Nx*My:Nx] = 1
    
    A_p1 = diags([diag_p1_l,diag_p1_u], [0,1],shape=(Nx*My, Nx*My))
    A_p1 = A_p1.tocsr()
    
    # H_i,j-1, west, m1
    diag_m1_u = np.zeros(Nx*My,dtype = 'int')
    diag_m1_u[0:Nx*My:Nx] = 1
    
    diag_m1_l = np.ones(Nx*My-1,dtype = 'int')
    diag_m1_l[Nx-1:Nx*My-1:Nx] = 0
    
    A_m1 = diags([diag_m1_u,diag_m1_l], [0,-1],shape=(Nx*My, Nx*My))
    A_m1 = A_m1.tocsr()
    
    # H_i+1,j+1, southeast, p1 and pN
    diag_p1N_u1 = np.ones(Nx*My-Nx-1,dtype = 'int')
    diag_p1N_u1[Nx-1:Nx*My-Nx-1:Nx] = 0
    
    diag_p1N_dia = np.zeros(Nx*My,dtype = 'int')
    diag_p1N_dia[Nx-1:Nx*My:Nx] = 1
    diag_p1N_dia[Nx*My-Nx :]=1
    A_p1N = diags([diag_p1N_dia, diag_p1N_u1], [0,Nx+1],shape=(Nx*My, Nx*My))
    A_p1N = A_p1N.tocsr()
    
    
    # H_i-1,j-1, northwest, pm1 and mN
    diag_m1N_dia = np.zeros(Nx*My,dtype = 'int')
    diag_m1N_dia[0:Nx] = 1
    diag_m1N_dia[0:Nx*My:Nx] = 1
    
    diag_m1N_l2 = np.ones(Nx*My-Nx-1,dtype = 'int')
    diag_m1N_l2[Nx-1:Nx*My-Nx-1:Nx] = 0
    A_m1N = diags([diag_m1N_l2,diag_m1N_dia], [-Nx-1,0],shape=(Nx*My, Nx*My))
    A_m1N = A_m1N.tocsr()
    
    # H_i-1,j+1, northeast, mNp1 
    diag_mNp1_l1 = np.ones(Nx*My-Nx+1,dtype = 'int')
    diag_mNp1_l1[0:Nx*My-Nx+1:Nx] = 0
    
    diag_mNp1_dia = np.zeros(Nx*My,dtype = 'int')
    diag_mNp1_dia[0:Nx] = 1
    diag_mNp1_dia[Nx-1:Nx*My:Nx] = 1
    A_mNp1 = diags([diag_mNp1_l1,diag_mNp1_dia], [-Nx+1,0],shape=(Nx*My, Nx*My))
    A_mNp1 = A_mNp1.tocsr()
    
    # H_i+1,j-1, southeast, m1pN 
    diag_pNm1_u2 = np.ones(Nx*My-Nx+1,dtype = 'int')
    diag_pNm1_u2[0:Nx*My-Nx+1:Nx] = 0
    
    diag_pNm1_dia = np.zeros(Nx*My,dtype = 'int')
    diag_pNm1_dia[Nx*My-Nx :] = 1
    diag_pNm1_dia[0:Nx*My:Nx] = 1
    A_pNm1 = diags([diag_pNm1_dia,diag_pNm1_u2], [0,Nx-1],shape=(Nx*My, Nx*My))
    A_pNm1 = A_pNm1.tocsr()
    
    ############################################################
    ##### Build the boundary condition for AH^n+1 = H^n+BC
    ############################################################
    # the water flux on the 4 sides (boundaries)
    Qn_we =  0.00
    Qe_we =  0.00
    Qw_we =  0.00
    Qs_we =  0.00
    
    
    BC_we                  = np.zeros(Nx*My)
    BC_we[0:Nx]            = -Qn_we*dts/dy
    BC_we[Nx*My-Nx :]      = -Qs_we*dts/dy
    BC_we[0:Nx*My:Nx]      =BC_we[0:Nx*My:Nx]- Qw_we*dts/dx
    BC_we[Nx-1:Nx*My:Nx]   =BC_we[Nx-1:Nx*My:Nx]- Qe_we*dts/dx 
    
    return A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1,BC_we