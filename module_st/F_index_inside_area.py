#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 17:17:51 2019

@author: qinayan
"""
import os
import gdal
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def find_total_area_index_of_stream_floodplain(filename,filename_str,filename_floodp,temp_path):
#    current_dir = os.getcwd()
    os.chdir(temp_path+'/GeoData/')
    
    gdata_dem          = gdal.Open(filename) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo_dem            = gdata_dem.GetGeoTransform() 
    dx = geo_dem[1]
    dy = -geo_dem[5]
    xmin_dem = geo_dem[0]
    ymax_dem = geo_dem[3] 
    ele = gdata_dem.ReadAsArray()
    nrows, ncols = ele.shape
#    
    eta_vector  = ele.flatten()*1.0
    validID_ele = np.where(eta_vector>0)[0]
    
    
    gdata_str          = gdal.Open(filename_str) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo_str            = gdata_str.GetGeoTransform() 
    dx = geo_str[1]
    dy = -geo_str[5]
    xmin_str = geo_str[0]
    ymax_str = geo_str[3] 
    ele_str = gdata_str.ReadAsArray()
    nrows_str, ncols_str = ele_str.shape
    
    row_str_up_ID = int((ymax_dem- ymax_str)/dy)
    col_str_up_ID = int((xmin_str- xmin_dem)/dx)
    row_str_low_ID = row_str_up_ID+nrows_str
    col_str_low_ID = col_str_up_ID+ncols_str
    
    # now the id of stream its own, the extra 's' means small
    row_str_up_s_ID  = 0
    col_str_up_s_ID  = 0
    row_str_low_s_ID = nrows_str
    col_str_low_s_ID = ncols_str
    
#    print (row_str_up_ID,col_str_up_ID,row_str_low_ID,col_str_low_ID)
    
    if row_str_up_ID <0:
        row_str_up_s_ID = np.abs(row_str_up_ID)
        row_str_up_ID   = 0
        
    if col_str_up_ID <0 :
        col_str_up_s_ID = np.abs(col_str_up_ID)
        col_str_up_ID   = 0
        
    if row_str_low_ID>nrows:
        row_str_low_s_ID = nrows_str-(row_str_low_ID-nrows)
        row_str_low_ID   = nrows
        
    if col_str_low_ID>ncols:
        col_str_low_s_ID = ncols_str-(col_str_low_ID-ncols)
        col_str_low_ID   = ncols
        
    ele_str_extend = -1.0*np.ones((nrows, ncols))
    ele_str_extend[row_str_up_ID:row_str_low_ID,col_str_up_ID:col_str_low_ID]=\
    ele_str[row_str_up_s_ID:row_str_low_s_ID,col_str_up_s_ID:col_str_low_s_ID]
    
    str_vector  = ele_str_extend.flatten()*1.0
    validID_str     = np.where(str_vector>0)[0]
    
    
    gdata_fld         = gdal.Open(filename_floodp) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo_fld           = gdata_fld.GetGeoTransform() 
    dx = geo_fld[1]
    dy = -geo_fld[5]
    xmin_fld = geo_fld[0]
    ymax_fld = geo_fld[3] 
    ele_fld = gdata_fld.ReadAsArray()
    nrows_fld, ncols_fld = ele_fld.shape
    row_fld_up_ID = int((ymax_dem- ymax_fld)/dy)
    col_fld_up_ID = int((xmin_fld- xmin_dem)/dx)
    row_fld_low_ID = row_fld_up_ID+nrows_fld
    col_fld_low_ID = col_fld_up_ID+ncols_fld
    
    row_fld_up_s_ID  = 0
    col_fld_up_s_ID  = 0
    row_fld_low_s_ID = nrows_fld
    col_fld_low_s_ID = ncols_fld
    
    if row_fld_up_ID <0:
        row_fld_up_s_ID = np.abs(row_fld_up_ID)
        row_fld_up_ID   = 0
        
    if col_fld_up_ID <0 :
        col_fld_up_s_ID = np.abs(col_fld_up_ID)
        col_fld_up_ID   = 0
        
    if row_fld_low_ID>nrows:
        row_fld_low_s_ID = nrows_fld-(row_fld_low_ID-nrows)
        row_fld_low_ID   = nrows
        
    if col_fld_low_ID>ncols:
        col_fld_low_s_ID = ncols_fld-(col_fld_low_ID-ncols)
        col_fld_low_ID   = ncols
        
        
    ele_fld_extend = -1.0*np.ones((nrows, ncols))
    ele_fld_extend[row_fld_up_ID:row_fld_low_ID,col_fld_up_ID:col_fld_low_ID] =\
    ele_fld[row_fld_up_s_ID:row_fld_low_s_ID,col_fld_up_s_ID:col_fld_low_s_ID]
    
    fld_vector  = ele_fld_extend.flatten()*1.0
    validID_fld = np.where(fld_vector>0)[0]

    validID_hill = np.setdiff1d(validID_ele, validID_fld)
    
#    os.chdir(current_dir)
    
    if False:
        test_three_zones = np.nan*np.ones(nrows*ncols)
        test_three_zones[validID_hill] = 0 #1 #0
        test_three_zones[validID_fld] = 1 #2 #1
        test_three_zones[validID_str] = 2 #0 #2
        
        test2D = np.reshape(test_three_zones,(nrows,ncols))
        fig, ax = plt.subplots(figsize=(8,8), dpi=70)     
        plt.rc("font", size=18)
        figplot  = ax.matshow(test2D, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.coolwarm_r) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
        plt.tight_layout()
        myLocator = mticker.MultipleLocator(100)
        ax.xaxis.set_major_locator(myLocator)
        ax.set_xlabel('X [m]', fontsize=18)
        ax.set_ylabel('Y [m]', fontsize=18)
    return (validID_str,validID_fld,validID_hill)

def find_total_area_index_of_Nf_Sf(filename,filename_Nf,filename_Sf,temp_path):
    current_dir = os.getcwd()
    os.chdir(temp_path+'/GeoData/')
 #   os.chdir(current_dir+'/GeoData/')
    
    gdata_dem          = gdal.Open(filename) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo_dem            = gdata_dem.GetGeoTransform() 
    dx = geo_dem[1]
    dy = -geo_dem[5]
    xmin_dem = geo_dem[0]
    ymax_dem = geo_dem[3] 
    ele = gdata_dem.ReadAsArray()
    nrows, ncols = ele.shape
#    
    eta_vector  = ele.flatten()*1.0
    validID_ele     = np.where(eta_vector>0)[0]
    
    
    gdata_Nf          = gdal.Open(filename_Nf) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo_Nf            = gdata_Nf.GetGeoTransform() 
    dx = geo_Nf[1]
    dy = -geo_Nf[5]
    xmin_Nf = geo_Nf[0]
    ymax_Nf = geo_Nf[3] 
    ele_Nf = gdata_Nf.ReadAsArray()
    nrows_Nf, ncols_Nf = ele_Nf.shape
    row_Nf_up_ID = int((ymax_dem- ymax_Nf)/dy)
    col_Nf_up_ID = int((xmin_Nf- xmin_dem)/dx)
    row_Nf_low_ID = row_Nf_up_ID+nrows_Nf
    col_Nf_low_ID = col_Nf_up_ID+ncols_Nf
    #ele_Nf_vector = ele_Nf.flatten()*1.0
    #Nf_ID_orig = np.where(ele_Nf_vector>0.0)[0]
    ele_Nf_extend = -1.0*np.ones((nrows, ncols))
    ele_Nf_extend[row_Nf_up_ID:row_Nf_low_ID,col_Nf_up_ID:col_Nf_low_ID] = ele_Nf
    
    Nf_vector  = ele_Nf_extend.flatten()*1.0
    validID_Nf     = np.where(Nf_vector>0)[0]
    
    
    gdata_Sf         = gdal.Open(filename_Sf) #CCW_2mDEM_deepStr, DEM_square_nofill2 saybrook_rec3_m  rec_test rec_test3_ditch rec_test3_stream
    geo_Sf           = gdata_Sf.GetGeoTransform() 
    dx = geo_Sf[1]
    dy = -geo_Sf[5]
    xmin_Sf = geo_Sf[0]
    ymax_Sf = geo_Sf[3] 
    ele_Sf = gdata_Sf.ReadAsArray()
    nrows_Sf, ncols_Sf = ele_Sf.shape
    row_Sf_up_ID = int((ymax_dem- ymax_Sf)/dy)
    col_Sf_up_ID = int((xmin_Sf- xmin_dem)/dx)
    row_Sf_low_ID = row_Sf_up_ID+nrows_Sf
    col_Sf_low_ID = col_Sf_up_ID+ncols_Sf
    ele_Sf_extend = -1.0*np.ones((nrows, ncols))
    ele_Sf_extend[row_Sf_up_ID:row_Sf_low_ID,col_Sf_up_ID:col_Sf_low_ID] = ele_Sf
    
    Sf_vector  = ele_Sf_extend.flatten()*1.0
    validID_Sf     = np.where(Sf_vector>0)[0]

    
    os.chdir(current_dir)
    
    if False:
        test_three_zones = np.nan*np.ones(nrows*ncols)
        test_three_zones[validID_ele] = 0 #1 #0
        test_three_zones[validID_Sf] = 1 #2 #1
        test_three_zones[validID_Nf] = 2 #0 #2
        
        test2D = np.reshape(test_three_zones,(nrows,ncols))
        fig, ax = plt.subplots(figsize=(8,8), dpi=150)     
        plt.rc("font", size=18)
        figplot  = ax.matshow(test2D, extent=[0,ncols*dx,nrows*dy,0],cmap=cm.coolwarm_r) #,vmin=0,vmax=0.032 #pcolor , matshow, imshow
        plt.tight_layout()
        cbar     = fig.colorbar(figplot,orientation='vertical',shrink=0.8, pad = 0.02,aspect = 30 ) #shrink=0.8, fraction=0.1, pad = 0.02;  #if shrink, shrink=0.9,
        cbar.set_label('zone')
        myLocator = mticker.MultipleLocator(100)
        ax.xaxis.set_major_locator(myLocator)
        ax.set_xlabel('X [m]', fontsize=18)
        ax.set_ylabel('Y [m]', fontsize=18)
    return (validID_Nf,validID_Sf)
