#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 19:37:03 2019

@author: qinayan
"""

import numpy as np

#-------------
# Mass conservation equation
#-------------

def Mass_conservation_method(h_max,nrows,ncols,invalidID, Ls, Le, ho, d_eta_1yr,d_thre,Bp,):
    
    h_ini    = np.ones(nrows*ncols)*h_max #h_Patton*1.0 #np.ones(nrows*ncols)

    m  = Ls/ho
    
    B  = Bp/Le
    B[invalidID] = -1.0
    B[np.isnan(B)]=-1.0

    h_zero_ID    = np.where((d_eta_1yr<d_thre)&(d_eta_1yr>-100)&(-d_eta_1yr>=B))[0]
    h_non0_divID = np.where((d_eta_1yr<d_thre)&(d_eta_1yr>-100)&(-d_eta_1yr<B))[0]
    
    d_eta_non0_div = d_eta_1yr[h_non0_divID]
    
    h_div = (np.log(B[h_non0_divID])-np.log(-d_eta_non0_div))/m[h_non0_divID]
    
    h_ini[h_zero_ID]   = 0.0
    h_ini[h_non0_divID] = h_div
    
    return h_ini 