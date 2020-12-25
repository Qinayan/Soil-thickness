#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:45:21 2019

@author: qinayan
"""

import os
import gdal
import numpy as np

def loadDEM(filename,temp_path):

#    current_dir = os.getcwd()
    os.chdir(temp_path+'/'+'/GeoData/')
    gdata          = gdal.Open(filename) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    
    geo            = gdata.GetGeoTransform() 
    
    dx = geo[1]
    dy = -geo[5]
    xmin = geo[0] #+ xres * 0.5
    ymax = geo[3] #- yres * 0.5
    
#    x_convt = (x_easting  - xmin)/xres
#    y_convt = (y_northing - ymax)/yres
#    col_samp = np.round(x_convt).astype(int)
#    row_samp = np.round(y_convt).astype(int)
    
    Topography_rec = gdata.ReadAsArray()

#    Topography_rec[115,118] =  Topography_rec[115,118]+0.99
    nrows, ncols = Topography_rec.shape
    Y = np.arange(0, nrows*dy, dy)
    X = np.arange(0, (ncols)*dx, dx) 
    X, Y = np.meshgrid(X, Y)
#    os.chdir(current_dir)
    return dx,dy, X, Y,xmin,ymax, Topography_rec


        
        
def sampling_locs(filename,sample_data_file,temp_path):
#    current_dir = os.getcwd()
    os.chdir(temp_path+'/'+'/GeoData/')
    
    x_easting  = []
    y_northing = []
    eta_sites  = []
    st_field   = []
    with open(sample_data_file,'r') as f:
        next(f)
        next(f)
        next(f)
        for line in f:
            each_line = line.split(',') 
            if float(each_line[4])>0:
                x_easting.append(float(each_line[1]))
                y_northing.append(float(each_line[2]))
                eta_sites.append(float(each_line[3]))
                st_field.append(float(each_line[4]))
        
    x_easting  = np.array(x_easting)
    y_northing = np.array(y_northing)
    eta_sites  = np.array(eta_sites)
    st_field   = np.array(st_field)
    
    
    gdata          = gdal.Open(filename) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo            = gdata.GetGeoTransform() 
    
    dx = geo[1]
    dy = -geo[5]
    xmin = geo[0] #+ xres * 0.5
    ymax = geo[3] #- yres * 0.5
    
    x_convt = (x_easting  - xmin)/dx
    y_convt = (ymax-y_northing)/dy
    col_samp = np.round(x_convt).astype(int)
    row_samp = np.round(y_convt).astype(int)
    
    Topography_rec = gdata.ReadAsArray()
    nrows, ncols = Topography_rec.shape
    sample_ID = row_samp*ncols+col_samp
    
    eta_lidar = Topography_rec[row_samp,col_samp]
#    os.chdir(current_dir)
    return col_samp,row_samp,sample_ID, eta_sites,st_field,eta_lidar
#os.chdir(current_dir)
