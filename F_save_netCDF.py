#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 17:04:16 2019

@author: qinayan
"""
import os
import numpy as np
from netCDF4 import Dataset

def save2D_results(data1_in, 
                   var1_name, 
                   folder, file_name, t): #data3_in,  data10_in, var3_name,  var10_name
    datanetcdf = Dataset(folder + '/' + file_name + str(t) + '.nc', 'w', format='NETCDF4')
    datanetcdf.description = '2D surface data, nrows * ncols;  1D below groud data, soil_layer_num * (nrows * ncols)'
    
#    datanetcdf2 = Dataset(folder + '/' + file_name + str(t) + '.nc', 'w', format='NETCDF4')
#    datanetcdf2.description = '3D data but the surface 2D is flatten, nlayers*(nrows*ncols) '

    # Get dimension
    My, Nx  = data1_in.shape
    
#    Myy, Nxx  = data4_in.shape
    
#    My3, Nx3  = data9_in.shape
    
    # Set up dimension for dataset in NetCDF
    datanetcdf.createDimension('y', My)
    datanetcdf.createDimension('x', Nx)
    
#    datanetcdf.createDimension('yy', Myy)
#    datanetcdf.createDimension('xx', Nxx)
#    
#    datanetcdf.createDimension('yyy', My3)
#    datanetcdf.createDimension('xxx', Nx3)

    """
    Create variables in the netcdf file
        var = netcdf.createVariable('Var_name', 'var_type', ('dimension_type'))
    """
    data1  = datanetcdf.createVariable(var1_name, 'f8', ('y','x'))
#    data2  = datanetcdf.createVariable(var2_name, 'f8', ('y','x'))
#    data3  = datanetcdf.createVariable(var3_name, 'f8', ('yy','xx'))
#    data4  = datanetcdf.createVariable(var4_name, 'f8', ('yy','xx'))
#    data5  = datanetcdf.createVariable(var5_name, 'f8', ('yy','xx'))
#    data6  = datanetcdf.createVariable(var6_name, 'f8', ('yy','xx'))
    
#    data7  = datanetcdf.createVariable(var7_name, 'f8', ('yy','xx'))
#    data8  = datanetcdf.createVariable(var8_name, 'f8', ('yy','xx'))
#    data9  = datanetcdf.createVariable(var9_name, 'f8', ('yyy','xxx'))
#    data10 = datanetcdf.createVariable(var10_name, 'f8', ('yyy','xxx'))


    # Assign data to variables in NetCDF file
    data1[:]  =  data1_in
#    data2[:]  =  data2_in
#    data3[:]  =  data3_in
#    data4[:]  =  data4_in
#    data5[:]  =  data5_in
#    data6[:]  =  data6_in
#    data7[:]  =  data7_in
#    data8[:]  =  data8_in
#    data9[:]  =  data9_in
#    data10[:] =  data10_in
         
    # Close the file
    datanetcdf.close()
#    datanetcdf2.close()

if __name__ == "__main__":
    current_dir = os.getcwd()
    os.chdir(current_dir+'/'+'../output/Surface_evolution')
    for k in range(6,501):
        i = k*2
        data = np.load('st_linear_diffus_dt_1_K_0_005_year_'+str(i*1000)+'.npy')

        save2D_results(data,'soil_thickness_accu', 'Dirichlet_bc',"ct",(i*1000)) 
    
