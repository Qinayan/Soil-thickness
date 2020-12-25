#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:21:31 2019

@author: qinayan
"""
import os
import numpy as np

def find_DEM_bn(eta):
#-------------
# find the boundary index
#-------------
    nrows,ncols = eta.shape
    eta_vector  = eta.flatten()*1.0
    validID     = np.where(eta_vector>0)[0]
    
    xn = [-1,0,1,-1,1,-1, 0, 1]
    yn = [1, 1,1, 0,0,-1,-1,-1]
    bn_ID = []
    ghostbn_ID = []
    for i in validID:
        irows = i//ncols
        icols = i%ncols
        for n in range (0,8):
            if (irows+xn[n]<=nrows-1) and (irows+xn[n]>=0) and (icols+yn[n]<=ncols-1) and icols+yn[n]>=0:
                if eta[irows+xn[n],icols+yn[n]] < 0.0:
                    bn_ID.append((irows)*ncols+(icols))
                    ghostbn_ID.append((irows+xn[n])*ncols+(icols+yn[n]))
    #--------------   
    # remove duplicate numbers in bn_ID. bn_ID has duplicate number happends 
    # when the bondary grid has more than one
    # neighbours which is outsided of the dem domain
    #--------------
    
    bn_ID_unique = []
    for item in bn_ID:
        if item not in bn_ID_unique:
            bn_ID_unique.append(item)
            
    bn_ID_unique_arr = np.array(bn_ID_unique)    
    bn_ID_arr = np.array(bn_ID)    
    ghostbn_ID_arr = np.array(ghostbn_ID)
    
    current_dir = os.getcwd()
    os.chdir(current_dir+'/GeoData/')      
    np.save('DEM_Sf_bn_id',bn_ID_unique_arr)
    # save the ID without deleteing duplicate numbers
    np.save('DEM_Sf_bn_id_dup',bn_ID_arr) 
    
    np.save('DEM_Sf_ghostbn_id_dup',ghostbn_ID_arr)
    
    os.chdir(current_dir)
    
    return (bn_ID_unique_arr, bn_ID_arr,ghostbn_ID_arr)