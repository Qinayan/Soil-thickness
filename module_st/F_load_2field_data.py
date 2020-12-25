#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:43:47 2019

@author: qinayan
"""

import os
import gdal
import numpy as np

def load_2field_data(filename,sample_data_file,temp_path):
#    current_dir = os.getcwd()
    os.chdir(temp_path+'/'+'/GeoData/')
    x_easting  = []
    y_northing = []
    eta_sites  = []
    st_auger   = []
    st_CPT     = []
    bot_type   = []
    with open(sample_data_file,'r') as f:
        next(f)
        for line in f:
            each_line = line.split(',') 
            if float(each_line[4])>0:
                x_easting.append(float(each_line[1]))
                y_northing.append(float(each_line[2]))
                eta_sites.append(float(each_line[3]))
                st_auger.append(float(each_line[5]))
                st_CPT.append(float(each_line[6]))
                bot_type.append((each_line[7]))
    
    
    x_easting  = np.array(x_easting)
    y_northing = np.array(y_northing)
    eta_sites  = np.array(eta_sites)
    st_auger   = np.array(st_auger)
    st_CPT     = np.array(st_CPT)
    bot_type   = np.array(bot_type)
    
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
    rock_ID = sample_ID[bot_type=='rock']
    sapr_ID = sample_ID[bot_type=='sapr']
    fitted_ID = sample_ID[bot_type=='fitted']
    
#    os.chdir(current_dir)
         
    return col_samp,row_samp,sample_ID,st_auger,st_CPT, rock_ID,sapr_ID,fitted_ID,bot_type
