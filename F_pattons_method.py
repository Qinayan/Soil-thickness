#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:26:29 2019

@author: qinayan
"""
import numpy as np

#-------------
# Patton, Nat-comm2018
# find the soil thickness from curvature
#-------------
    
def Pattons_method(curv,bn_ID,validID,validID_hill,nrows,ncols,h_bar,slope_TMR_ssd):
    
    curv[bn_ID] = 0.0
    curv_valid = curv[validID]*1.0
#    ssd_surv = np.std(curv[validID_hill])
#
#    slope_TMR_ssd
#    slope_TMR_ssd = -446.3*ssd_surv+30.3
#    print (slope_TMR_ssd)

    h_Patton = np.zeros(nrows*ncols)
    h_Patton[validID] = slope_TMR_ssd*curv_valid+ h_bar
    h_Patton[h_Patton<0.0] = 0.0
    return h_Patton