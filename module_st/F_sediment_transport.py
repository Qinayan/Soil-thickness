# -*- coding: utf-8 -*-
"""
Copyright (C) 2017, Qina Corporation
All rights reserved.

Distributed Hydrologicc and Regional Analysis (DHARA) Model
DHARA model is made available as a restricted, non-exclusive, 
non-transferable license for education and research purpose only, 
and not for commercial use. See the LICENSE.txt for more details.

Author: Qina.yan@gmail.com (Qina Yan)
"""

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from scipy.sparse import dia_matrix

"""
{ item_description }
matrix_A:
matrix_BC:
"""

def transport_limited(K, wd_vector,slope,direction,beta,nrows,ncols,dt,dx,bn_ID,validID):
    
    alpha = beta+1.0/3
    xn = [-1,0,1,-1,1,-1,0,1,0]
    yn = [1,1,1,0,0,-1,-1,-1,0]
    dn = [2.**.5,1.,2.**.5,1.,1.,2.**.5,1.,2.**.5,0]

    wd = np.reshape(wd_vector,(nrows,ncols))
    qs_in = np.zeros((nrows,ncols))
    dis  = np.ones((nrows,ncols))

    for k in range(len(validID)): #len(validID_inside)
        i = validID[k]//ncols
        j = validID[k] - i*ncols
#    for i in range (0, nrows):
#        for j in range (0,ncols):
        n=direction[i][j]
        xw = i + xn[n]
        yw = j + yn[n]      
            
        qs_in[xw][yw]= K*wd[i,j]**alpha*slope[i,j]**beta + qs_in[xw][yw] #eta = eta_vector.reshape(nrows,ncols)
        dis[i,j] = dx*dn[n]
            
    dis_vector = dis.flatten()*1.0
    qs_in_vector = qs_in.flatten()/dis_vector

    qs_out = K*wd**alpha*slope**beta 
    qs_out_vector = qs_out.flatten()/dis_vector
    
#    nabula_qs = (qs_in_vector -qs_out_vector)/dis
    qs_in_vector[bn_ID] = 0.0
    return qs_in_vector,qs_out_vector
    
#def Overland_SedimentTransport(wd,Sse,Ssw,Ssn,Sss,A_p1,A_m1,A_pN,A_mN,Ca_vector,nrows,ncols,level):
#    if level ==1:
#        Kw = Kw1
#    else:
#        Kw = Kw2
#    a = 17.0/6.0
#    b = 2.5
#    Swx = Kw*dt/dx
#    Swy = Kw*dt/dy
#    
#    wd_p1 = A_p1.dot(wd)
#    wd_m1 = A_m1.dot(wd)
#    wd_pN = A_pN.dot(wd)
#    wd_mN = A_mN.dot(wd)

##    eta_w = np.zeros((nrows,ncols))
#    eta_w_vector = Swy*(((wd_pN+wd)*0.5)**a*Sss**b-((wd_mN+wd)*0.5)**a*Ssn**b)
##Swx*(((wd_p1+wd)*0.5)**a*Sse**b-((wd_m1+wd)*0.5)**a*Ssw**b)#
#    return -eta_w_vector,np.zeros((nrows,ncols)) #eta_w_Ca_vector

def detachment_limited(wd_vector,slope,Ca_vector,nrows,ncols):

    
    wd = wd_vector.reshape(nrows,ncols)
    
    
    eta_r             = -dt*Kr/rho_b*(rho_w*g*wd*slope -tau_c)
    no_move_id        = np.where(rho_w*g*wd*slope<tau_c)
    eta_r[no_move_id] = 0.0
    eta_r_vector      = np.array(eta_r).flatten()
    eta_r_Ca_vector   = eta_r_vector*Ca_vector

    
    return eta_r_vector,eta_r_Ca_vector