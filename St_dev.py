#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 12:00:49 2019

@author: qinayan
"""

import os
import sys
import timeit
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from landlab.plot.imshow import imshow_grid
#from landlab import ModelParameterDictionary, RasterModelGrid
#from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse

current_dir = os.getcwd()
sys.path.insert(1, current_dir+"/"+'/module_st/')
from F_costheta import costheta
from F_find_DEM_bn import find_DEM_bn
#from F_save_netCDF import save2D_results
from F_loadDEM import loadDEM,sampling_locs
from F_pattons_method import Pattons_method
from F_load_2field_data import load_2field_data
from F_slope_n_direction import slope_direction
from F_linear_diffusion import linear_diffusion_2D
from F_sediment_transport import transport_limited
from F_matrix_for_rolling import matrix_for_rolling
from F_mass_conservation_st import Mass_conservation_method
from F_overland_Q import overland_setup, EstimateKOverland, ImplicitSolver_2nd_BC, FlowVelocity
from F_index_inside_area import find_total_area_index_of_stream_floodplain,find_total_area_index_of_Nf_Sf
os.chdir(current_dir)

#-------------
# beginning set up
#-------------
first_time_run       = False
Diffusion            = False
OverlandFlow_erosion = False
soil_formation       = False
Overland_flow        = False
simu_side = 'Nf' #'North-facing' or 'South-facing'

#-------------
# diffusion BC setup
#-------------
# assign value boundary condtion
Dirichlet_BC   = False

# assgin flux boundary condition
Neumann_BC     = True

# if steep slope, calcuate the cos value for soil formation rate
steep_slope    = False 

#-------------
# parameters
#-------------

start = timeit.default_timer()
#-.-.-.-.-.-.-
# North facing parameter
#-.-.-.-.-.-.-
rho_s = 1.53*1e3 # soil bulk density, unit = g/cm^3 =10^3 kg/m^3
rho_r = 2.0*1e3

# soil formation parameters
ho     = 0.18    #Nf: 0.18     #Sf:0.16,    0.34,  0.1 # 1/m
Po     = 0.0001  #Nf: 0.00006  #Sf:0.0001, 0.000168,0.0021 0.0005 #0.00019 #m/yr
Bp     =(rho_r/rho_s)*Po
h_bar  = 0.6    #Nf: 0.6;    #Sf:0.56, 0.45, 0.52
a_Pattons = 45   #Nf = 45;    #Sf:50, 20, 48

d_thre = -2e-6   #Nf:-2e-6    #Sf: -2e-6, -1e-6, -0.00002 # has to be <0; #-0.0014 #-0.0026

# soil diffusion parameters
Kl = 1e-3      #Nf: 1e-3    #Sf:1.6e-3, 1.26e-3, 3e-3 per year, or 1.0*0.036

#over land flow soil erosion
Kq   = 5e-5      #Nf: 5e-5    #Sf:8e-5, 1e-3, 1e-5
beta = 1.68 #original E-H: 2.5 or generalized E-H: 1.68

#landscape evolution parameterw
U   =  0.0
Klx =  9e-3 # per year, or 1.0*0.036
Kly =  Klx*1.0 # per year, or 1.0*0.036
Kqs = 0.9/(np.sqrt(9.81)*(0.025)**2*(1.65**1.18)*(65*1e-6)**0.18)/3600/24/365#use the manning's value below

#overland flow depth parameters
inter_loop    = 2
dt_Q          = 1*365/inter_loop #1.0/24/3600 # day #Nf:0.0001; Sf:
manning_value = 0.025/3600.0/24 # sec/m^(1/3) --> hr/m^(1/3) --> day
rain_depth    = 1.54e-3/365 # 0.3-0.75 m/yr --> m./day

#--------------
# simulation period
#--------------
total_yr = 0
dt  = 1.0
evo_steps = int(total_yr/dt)

#---------------
# load data
# start from DEM 
#---------------

filename      = "st_ct_"+simu_side+"_dem.tif" #"st_ct_dem"; "st_so_dem.tif"
filename_str     = "stream_ct.tif" #"st_nr_dem.tif"
filename_floodp  = "floodplain_ct.tif" #"st_nr_dem.tif"
filename_sample_data = 'CB-PH-SoilThickness_OCT19_'+str(simu_side)+'_hs.csv'
#sample_data_file_drill = 'CB-PH-SoilThickness_OCT19_'+simu_side+'.csv'
#sample_data_file_drill = 'CB-PH-SoilThickness_OCT19_CPT_'+simu_side+'.txt'
temp_path = current_dir+"/"
dx,dy, X, Y, xmin, ymax,eta = loadDEM(filename,temp_path)
eta[eta<0.0]     = -1.0

#dx,dy, X, Y, xmin, ymax,aspect = loadDEM("Aspect_Nf.tif")

#sample_data_file_drill = sample_data_file_drill_Nf

nrows, ncols     = eta.shape
eta_vector       = eta.flatten()*1.0
# find the valid and invalid ID address
validID          = np.where(eta_vector>0)[0]
invalidID        = np.where(eta_vector<0)[0]

validID_str,validID_fld,validID_hill = \
find_total_area_index_of_stream_floodplain(filename,filename_str,filename_floodp,temp_path)

#validID_Nf,validID_Sf = \
#find_total_area_index_of_Nf_Sf(filename,filename_Nf,filename_Sf)

#-.-.-.-.-.-.-.- 
# load sampling data
#-.-.-.-.-.-.-.- 
#col_samp_drill,row_samp_drill,sample_ID_drill, eta_site_drill,st_field_drill,\
#eta_lidar_drill = sampling_locs(filename,sample_data_file_drill,temp_path)

col_samp,row_samp,sample_ID,st_auger,st_CPT, rock_ID,sapr_ID,fitted_ID,bot_type\
= load_2field_data(filename,filename_sample_data,temp_path)

#-.-.-.-.-.-.-.- 
# if run for the first time ever, run
#-.-.-.-.-.-.-.- 
if first_time_run:
    bn_ID, bn_ID_dup,ghostbn_ID_dup= find_DEM_bn(eta)


#-.-.-.-.-.-.-.- 
#load the ID index of DEM boundary
#-.-.-.-.-.-.-.- 
else:
    os.chdir(current_dir+"/"+'/GeoData/')
    bn_ID          = np.load('DEM_'+simu_side+'_bn_id.npy') #halfm
    bn_ID_dup      = np.load('DEM_'+simu_side+'_bn_id_dup.npy') #halfm
    ghostbn_ID_dup = np.load('DEM_'+simu_side+'_ghostbn_id_dup.npy') #halfm
    os.chdir(current_dir)
#-.-.-.-.-.-.-.- 
#load the surface soil transport rate
#-.-.-.-.-.-.-.- 
os.chdir(current_dir+"/"+'/GeoData/Curv_collections/')
# this is one time step of soil diffusion, dt = 1yr
#Output_curv_5.5m_ #Output_curv_poly14.5m_  # Output_curv_1000yr_inst_
curv_filename = 'Output_curv_1750yr_inst_'+simu_side+'.csv'
#data_eta_d = np.load('st_linear_'+simu_side+'_diffus_dt_1_K_0_005_year_1500.npy')
curv    =  np.zeros(nrows*ncols)
curv_validHill = []
with open(curv_filename,'r') as f:
    next(f)
    for line in f:
        each_line = line.split(',') 
        curv_validHill.append(float(each_line[0]))

curv_validHill = np.array(curv_validHill)
curv[validID_hill] = curv_validHill*1.0
d_eta_d_1yr = curv*Kl
#d_eta_d_1yr = data_eta_d.reshape(nrows*ncols)/1500.0

os.chdir(current_dir+"/"+'/GeoData/')
# this is one time step of soil overland flow erosion, dt = 1yr
data_eta_s = np.load('delta_eta_'+simu_side+'_100Kqs_dt_1_year_3000.npy') 
#data_eta_s = np.load('st_OFE_'+simu_side+'_dt_1_year_1500.npy') 
d_eta_s_1yr = data_eta_s.reshape(nrows*ncols)/3000.0


data_wd = np.load('wd_'+simu_side+'dt_Q_1e-4sec1000.npy')
data_wd = data_wd - np.nanmin(data_wd)
data_wd[np.isnan(data_wd)] = 0.0
upper_lim  = 1.9
bottom_lim = 0.7
data_wd[data_wd>upper_lim]  = upper_lim
data_wd[data_wd<bottom_lim] = bottom_lim
data_wd = data_wd-bottom_lim

data_wd_vector = data_wd.flatten()*1.0

os.chdir(current_dir)


##================
# load matrix for the cos_theta
# calculate overland flow depth, hw, and velocity, v
##================

A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1,BC_we = \
matrix_for_rolling(nrows,ncols,dx,dy,dt_Q)

stop = timeit.default_timer()
print('Loading data time (sec): ', (stop - start)) 

if Overland_flow:
    wd_vector_save = 0.0
    wd_vector,we_vector = overland_setup(eta_vector)
    
    for i in range(inter_loop):
        we_vector[ghostbn_ID_dup]    = we_vector[bn_ID_dup]
        ke,kw,kn,ks,Sse,Ssw,Ssn,Sss  = EstimateKOverland(A_p1,A_m1,A_pN,A_mN, A_p1N,A_m1N,A_pNm1,A_mNp1, wd_vector,we_vector,dx,dy,manning_value)
        we_vector_new, wd_vector_new = ImplicitSolver_2nd_BC(we_vector,eta_vector,rain_depth*dt_Q,ke,kw,kn,ks,nrows, ncols,BC_we,dx,dy,dt_Q)
        wd_vector = wd_vector_new*1.0
        we_vector = we_vector_new*1.0
        wd_vector[invalidID] = 0.0
        wd_vector_save = wd_vector_save+wd_vector
    
    #    if (i+1)%10 == 0:
    #        os.chdir(current_dir+"/"+'/output/') #st_linear_diffus_dt_1_K_0_005_year_1000_Dirichlet
    #        np.save('wd_'+simu_side+'dt_Q_1e-4sec'+str(i+1), wd_vector.reshape((nrows,ncols)))
    #        os.chdir(current_dir)
    

'''
options of surface diffusion
the following provides three options 
'''
##================
## linear 2-D diffusion
## -K*(d^2z/dx^2)
## Am^(-mh) 
## the right hand side
##================
d2zdx2_save   = np.zeros(nrows*ncols)
nabla_qs_save = np.zeros(nrows*ncols)

matrix_Ad, matrix_BC,\
matrix_Adxp1, matrix_Adxm1, matrix_Adyp1, matrix_Adym1,\
BCx_p1, BCx_m1, BCy_p1,BCy_m1 = linear_diffusion_2D(nrows,ncols,dx,dy,Klx,Kly)

start = timeit.default_timer()
if True:
#-------------
# initial soil thickness
#-------------
    h_max     = 2.0
    d_eta_1yr = d_eta_d_1yr + Kq*d_eta_s_1yr/Kqs
    d_eta_1yr[bn_ID] = 0.0
    d_eta_1yr[invalidID] = -100.0
    
    cos_theta = costheta(eta_vector+d_eta_1yr,dx,dy,A_p1,A_m1,A_pN,A_mN) ;  np.ones(nrows*ncols)
    Le        = cos_theta*1.0 # np.ones(nrows*ncols) #cos_theta*1.0
    Ls        = cos_theta*1.0 # np.ones(nrows*ncols) #cos_theta*1.0
    
    #convergent threshod, for Patton's method
    con_ID = np.where(d_eta_1yr>=d_thre)[0]
    #divergent threshod, for Mass conservation method
    div_ID  = np.where((d_eta_1yr<d_thre)&(d_eta_1yr>-100))[0]
    
    #-------------F_pattons_method
    # Patton, Nat-comm2018
    # find the curvature, and then soil thickness
    #-------------
#    curv = d_eta_d_1yr/Klx # ( matrix_Ad.dot(eta_vector)+ matrix_BC)/Klx# d_eta_1yr/Klx
    h_Patton = Pattons_method(curv,bn_ID,validID,validID_hill,nrows,ncols,h_bar,a_Pattons)
    
    #-------------
    # Mass conservation equation
    #-------------
    h_ini = Mass_conservation_method(h_max,nrows,ncols,invalidID, Ls, Le, ho, d_eta_1yr,d_thre,Bp,)
    
    #-------------
    # Combine the two method
    #-------------
    h_ini[con_ID]      = h_Patton[con_ID]
    h_ini[h_ini>h_max] = h_max
    
#    sample_ID_drill_MC,drill_ID_MC,model_ID_MC = np.intersect1d(sample_ID_drill, div_ID,return_indices=True)
#    sample_ID_drill_PM,drill_ID_PM,model_ID_PM  = np.intersect1d(sample_ID_drill,con_ID,return_indices=True)
    
    
stop = timeit.default_timer()
print('Run the soil thickness time (sec): ', (stop - start)) 

#=============================
# smooth h_ini to compare with the sampling site
#=============================
h_sim_smooth = h_ini*1.0

#h_sim_smooth = np.reshape(h_sim_smooth,(nrows,ncols))*1.0
#resol     = 3 # in m, the smoothing range of the soil thickness
#
#h_sim_temp = np.zeros((nrows,ncols))
#grid_num  = int(resol/dx)
#hafl_grid = int((grid_num-1)/2.0)
#loc_ind   = np.arange(-hafl_grid,hafl_grid+1,1)
##loc_ind   = loc_ind[loc_ind != 0]
#neigh_num = (hafl_grid*2+1)**2
#
##----------------------
## do the smoothing
##----------------------
#for i in row_samp:
#    for j in col_samp:
#        temp_store = 0.0
#        loc_cell = 0
#        
#        for k in loc_ind:
#            for l in loc_ind:
#                if (i+k > -1) and (j+l > -1) and (i+k < nrows) and (j+l < ncols) :
#                    temp_store += h_sim_smooth[i+k,j+l]
#                    loc_cell   +=1
#        if loc_cell>0:
#            h_sim_temp[i,j] = temp_store/loc_cell
#        
#h_sim_smooth = h_sim_temp.flatten()
#

#----------------
# find the spots of very high and low h 
#--------------
#h_ini_group = np.nan*np.ones(nrows*ncols)
#h_ini_le50cm_ID  = np.where(h_ini<0.5)[0]
#h_ini_la100cm_ID = np.where(h_ini>1.0)[0]
#
#h_ini_group[h_ini_le50cm_ID] = 0.5  #h_ini[h_ini_le50cm_ID]*1.0
#h_ini_group[h_ini_la100cm_ID] = 1.0 # h_ini[h_ini_la100cm_ID]*1.0


#h_ini_le50cm_row  = h_ini_le50cm_ID//ncols
#h_ini_le50cm_col  = h_ini_le50cm_ID%ncols
#
#h_ini_la100cm_row  = h_ini_la100cm_ID//ncols
#h_ini_la100cm_col  = h_ini_la100cm_ID%ncols


#----------------
# initial value of h
#--------------
h = 0.5*np.ones(nrows*ncols) # h_ini*1.0

#==============================
# soil thickness evolution
#==============================

start = timeit.default_timer()
for i in range(evo_steps):
    
    #------------------
    # diffusion process
    #------------------
    if Diffusion:
        
    #Boundary condition
        if Dirichlet_BC:
#            d2zdx2 = matrix_Ad.dot(eta_vector)+ matrix_BC
            diff_in = -(matrix_Adxm1.dot(eta_vector)+BCx_m1)-\
                       (matrix_Adym1.dot(eta_vector)+BCy_m1)
            diff_out = -(matrix_Adxp1.dot(eta_vector)+BCx_p1)-\
                        (matrix_Adyp1.dot(eta_vector)+BCy_p1) 
            d2zdx2 = diff_in-diff_out
            
            d2zdx2[bn_ID] = 0.0
            
            
        if Neumann_BC:
            eta_vector[ghostbn_ID_dup] = eta_vector[bn_ID_dup]
#            d2zdx2 = matrix_Ad.dot(eta_vector)+ matrix_BC
            diff_in = -(matrix_Adxm1.dot(eta_vector)+BCx_m1)-\
                       (matrix_Adym1.dot(eta_vector)+BCy_m1)
            diff_out = -(matrix_Adxp1.dot(eta_vector)+BCx_p1)-\
                        (matrix_Adyp1.dot(eta_vector)+BCy_p1) 
            d2zdx2 = diff_in-diff_out
            
            
        d2zdx2[validID_str] = 0.0
    
        d2zdx2_save = d2zdx2_save+d2zdx2*dt
        
    else:
        d2zdx2 = 0.0
        
    #------------------
    # overland flow erosion
    #------------------
     
    if OverlandFlow_erosion:
#        ele_vector_new = eta_vector+d_eta_d_1yr
        eta_vector[ghostbn_ID_dup] = eta_vector[bn_ID_dup]                                         
        slope,direction            = \
        slope_direction(eta_vector,nrows,ncols,dx,dy,validID,bn_ID) 
        qs_in_vector,qs_out_vector = \
        transport_limited(Kqs,data_wd_vector,slope,direction,beta,nrows,ncols,dt,dx,bn_ID,validID) 
        nabla_qs                   = qs_out_vector-qs_in_vector
        nabla_qs[validID_str]      = 0.0
#        data_wd_vector             = data_wd_vector-nabla_qs*dt
        data_wd_vector[data_wd_vector<0.0] = 0.0
        nabla_qs_save              = nabla_qs_save-nabla_qs*dt
    #------------------
    # new elevation
    #------------------
    eta_vector = eta_vector+(d2zdx2-nabla_qs)*dt + U*dt
    
    #------------------
    # soil formation process
    #------------------
    if soil_formation:
        if steep_slope:
            eta_vector_p1=A_p1.dot(eta_vector)
            eta_vector_m1=A_m1.dot(eta_vector)
            eta_vector_pN=A_pN.dot(eta_vector)
            eta_vector_mN=A_mN.dot(eta_vector)
            
            d_eta_x = eta_vector_p1-eta_vector_m1
            d_eta_y = eta_vector_pN-eta_vector_mN
            
            cosx=2.0*dx/(4.0*dx**2+d_eta_x**2)**0.5
            cosy=2.0*dy/(4.0*dy**2+d_eta_y**2)**0.5
            L = cosx*cosy # cos theta, surface angle
    
        else: 
            L = 1.0
            
        B = Bp/L
        m = L/ho 
        
        P = B*np.exp(-m*h)
        
        hnew_temp = h+dt*(P+diff_in)-dt*diff_out
        h_zero_ID = np.where(hnew_temp<0.0)[0]
        hnew_temp[h_zero_ID] = dt*B
        h = hnew_temp*1.0
    
    if (i+1) == 100 or (i+1)%500 == 0:
        nabla_qs_save[invalidID] = np.nan
        os.chdir(current_dir+"/"+'/output/') #st_linear_diffus_dt_1_K_0_005_year_1000_Dirichlet
        np.save('delta_eta_'+simu_side+'_dt_1_year_'+str(i+1), nabla_qs_save.reshape((nrows,ncols)))
        os.chdir(current_dir)
        
#        save2D_results(d2zdx2_save.reshape((nrows,ncols)),'soil_thickness_accu',
#                       'output',"ct",(i+1)) 
        
stop = timeit.default_timer()
print('Evolution time (min): ', (stop - start)/60.0)  

#col_samp,row_samp,sample_ID,st_auger,st_CPT = load_2field_data(filename,'CB-PH-SoilThickness_OCT19_'+str(simu_side)+'_hs.csv')


#=============================
# calculate the bedrock or saproilte weathering rate of soil layer bottom
#=============================
### the soil bottom layer prodcution rate:
#cos_theta = costheta(eta_vector,dx,dy,A_p1,A_m1,A_pN,A_mN)

P = (Bp/cos_theta)*np.exp(-(cos_theta/ho)*h_ini)
P[invalidID] = -1.0    

# the weathering capacity, or the maximum weathering rate if h = 0

#----------------
# save data
#--------------
os.chdir(current_dir+"/"+'output')
#np.save('surface_transport_rate_per_year'+simu_side,d_eta_1yr)
np.save('P_rate_per_year_new'+simu_side,P)

#----------------
# plot 2D data
#--------------
##---method 1: landlab
#min_z = np.min(z[np.where(z>0)])
#max_z = np.max(z[np.where(z>0)])
#imshow_grid(mg, 'topo_elevation', limits=(min_z,max_z))

##---method 2: Qina's plot


ls = LightSource(azdeg=315 , altdeg=135) #135
cmap = plt.cm.gist_earth
ve = 2

nan_ind = np.where(eta<0.0)
eta[nan_ind] = np.nan

h_ini_all = d_eta_1yr*1.0 #P, d_eta_1yr, h_ini
h_ini[invalidID] = np.nan
h_ini[validID_str]= np.nan
h_ini[validID_fld]= np.nan

h_ini_all[invalidID] = np.nan
h_ini_all[validID_str]= np.nan

xx = np.reshape(h_ini,(nrows,ncols))
xx2D = np.reshape(h_ini_all,(nrows,ncols))

fig, ax  = plt.subplots(nrows=1, ncols=1,figsize=(6.5,8), dpi=50) #Nf: figsize=(6.5,5.3)
plt.rc("font", size=18)
figplot  = ax.matshow(xx2D, extent=[0,ncols*dx,nrows*dy,0],
                      vmin=-1e-5,vmax=1e-5,cmap=cm.seismic) #cmap=cm.YlOrBr,YlGnBu, terrain, seismic # ,vmin=0,vmax=2.0 #,vmin=-8.1,vmax=8.1, #pcolor , matshow, imshow
# Vary vertical exaggeration 
ax.imshow(ls.hillshade(eta, vert_exag=ve, dx=dx, dy=dy), extent=[0,ncols*dx,nrows*dy,0], cmap='gray',alpha = 0.5) #, dx=dx, dy=dy
# #Nf:  pad = 0.06,aspect = 25
divider = make_axes_locatable(ax)
cax     = divider.append_axes('right', size='5%', pad=0.05)
cbar    = fig.colorbar(figplot,cax=cax,orientation='vertical',shrink=0.8,aspect = 30 ) 
cbar.set_label(r'Soil thickness [m]')

#plt.text(20, 300, 'Year='+str(i)+'kyr', fontsize=18)
plt.tight_layout()
myLocator = mticker.MultipleLocator(200)
ax.xaxis.set_major_locator(myLocator)
ax.yaxis.set_major_locator(myLocator)
ax.set_xlabel('X [m]', fontsize=18)
ax.set_ylabel('Y [m]', fontsize=18)
    

##----------------
## add sampling locations
##--------------

row_samp_rock = rock_ID//ncols
col_samp_rock = rock_ID%ncols

row_samp_sapr = sapr_ID//ncols
col_samp_sapr = sapr_ID%ncols

row_samp_fitted  = fitted_ID//ncols
col_samp_fitted = fitted_ID%ncols


#ax.scatter(X[row_samp_drill,col_samp_drill],Y[row_samp_drill,col_samp_drill],s=20,marker = 'o') #,label = 'Drill'
ax.scatter(X[row_samp,col_samp],Y[row_samp,col_samp],s=20,marker = 'o') #,label = 'Drill'
ax.scatter(X[row_samp_rock,col_samp_rock],Y[row_samp_rock,col_samp_rock],s=20,marker = 'o') 
ax.scatter(X[row_samp_sapr,col_samp_sapr],Y[row_samp_sapr,col_samp_sapr],s=20,marker = 'o') 
ax.scatter(X[row_samp_fitted,col_samp_fitted],Y[row_samp_fitted,col_samp_fitted],s=20,marker = 'o') 
#
#ax.errorbar(st_field_hs[bot_type=='rock'],h_ini[rock_ID],fmt='o',label='rock')
#ax.errorbar(st_field_hs[bot_type=='sapr'],h_ini[sapr_ID], fmt='o',label='saprolite')
#ax.errorbar(st_field_hs[bot_type=='fitted'],h_ini[fitted_ID], fmt='o',label='fitted')

##ax.legend(loc=2)
#plt.tight_layout()
#
plt.tight_layout()
os.chdir(current_dir+"/"+'Figures')

#fig.savefig('Soil_thickness_2D'+str(simu_side))
#----------------
# bedrock locations
#--------------
#bedrock_row = h_zero_ID//ncols
#bedrock_col = h_zero_ID%ncols
#ax.scatter(X[bedrock_row,bedrock_col],Y[bedrock_row,bedrock_col],s=20,c='w',marker = 'o')


#---------------------
# soil thickness comparision
#---------------------
#fig, ax  = plt.subplots(figsize=(5,5), dpi=50)   
#xdig = np.array([0,1.6])
#ydig = np.array([0,1.6])
#plt.plot(xdig,ydig,'k--',alpha=0.5)
##plt.plot(st_field_drill[drill_ID_MC],h_ini[sample_ID_drill_MC],'d',label='Drill_div')
##plt.plot(st_field_drill[drill_ID_PM],h_ini[sample_ID_drill_PM],'o',label='Drill_con')
##ax.legend(loc=2,fontsize = 12)
#plt.plot(st_field_drill,h_ini[sample_ID_drill],'o')
#ax.set_xlabel('Soil thickness at field [m]', fontsize=18)
#ax.set_ylabel('Soil thickness from modeling [m]', fontsize=18)
#plt.tight_layout()
#fig.savefig('Soil_thickness_comp')


#---------------------
# Compare auger with CPT
#---------------------
#fig, ax  = plt.subplots(figsize=(5,5), dpi=50)   
#xdig = np.array([0,1.6])
#ydig = np.array([0,1.6])
#ax.plot(xdig,ydig,'k--',alpha=0.5)
##plt.plot(st_field_drill[drill_ID_MC],h_ini[sample_ID_drill_MC],'d',label='Drill_div')
##plt.plot(st_field_drill[drill_ID_PM],h_ini[sample_ID_drill_PM],'o',label='Drill_con')
##ax.legend(loc=2,fontsize = 12)
#ax.plot(st_auger,st_CPT,'o')
#ax.set_xlabel('Soil thickness from auger [m]', fontsize=18)
#ax.set_ylabel('Soil thickness from CPT [m]', fontsize=18)
#plt.tight_layout()


#---------------------
# Compare model results with auger & CPT
#---------------------

st_field_hs = (st_auger+st_CPT)/2.0
xerr = np.abs(st_auger-st_field_hs)/2+0.01

fig, ax  = plt.subplots(figsize=(5,5), dpi=50)   
xdig = np.array([0,1.6])
ydig = np.array([0,1.6])
plt.plot(xdig,ydig,'k--',alpha=0.5)
#plt.plot(st_field_drill[drill_ID_MC],h_ini[sample_ID_drill_MC],'d',label='Drill_div')
#plt.plot(st_field_drill[drill_ID_PM],h_ini[sample_ID_drill_PM],'o',label='Drill_con')
#ax.legend(loc=2,fontsize = 12)
ax.errorbar(st_field_hs,h_sim_smooth[sample_ID],xerr=xerr, fmt='o')
ax.errorbar(st_field_hs[bot_type=='rock'],h_sim_smooth[rock_ID],fmt='o',label='rock')
ax.errorbar(st_field_hs[bot_type=='sapr'],h_sim_smooth[sapr_ID], fmt='o',label='saprolite')
ax.errorbar(st_field_hs[bot_type=='fitted'],h_sim_smooth[fitted_ID], fmt='o',label='fitted')
ax.legend(loc=4,fontsize = 12)
ax.set_xlabel('Soil thickness at field [m]', fontsize=16)
ax.set_ylabel('Soil thickness from modeling [m]', fontsize=16)
plt.tight_layout()
#fig.savefig('Soil_thickness_1to1_compar_'+str(simu_side))

#---------------------
# histogram plot
#---------------------
#total = h_ini[validID]*1.0
#fig, ax  = plt.subplots(figsize=(10,4), dpi=100)     
#plt.rc("font", size=18)
#plt.hist(total, bins=100)
#ax.set_xlim([0.0,2])
##ax.set_xlabel('Surface elevation change rate [m/yr]', fontsize=18)
#ax.set_xlabel('Soil thickness [m]', fontsize=18)
#ax.set_ylabel('Counts', fontsize=18)

#---------------------
# curvature v.s. soil thickness
#---------------------

#fig, ax  = plt.subplots(figsize=(5,5), dpi=80)   
#plt.plot(curv[validID],h_ini[validID],'*',label = 'All grids')
#plt.plot(curv[sample_ID_drill],h_ini[sample_ID_drill],'*',label = 'Sample sites')
#ax.set_ylabel('Soil thickness at field [m]', fontsize=18)
#ax.set_xlabel('Curvature [1/m]', fontsize=18)
#ax.legend()


#-----------------------
# combine curvature-thickness with pdf plot
#------------------------
# Set up the axes with gridspec
#
fig = plt.figure(figsize=(6, 6),dpi = 50)
grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
main_ax = fig.add_subplot(grid[1:, :-1])
x_hist = fig.add_subplot(grid[0, :-1], yticklabels=[], sharex=main_ax)
y_hist = fig.add_subplot(grid[1:, -1], xticklabels=[], sharey=main_ax)

#d_eta_1yr; curv
## scatter points on the main axes
main_ax.plot(d_eta_1yr[validID_hill],h_ini[validID_hill],'o',alpha = 0.2, label = 'All grids')
main_ax.plot(d_eta_1yr[sample_ID],st_field_hs, '*',label = 'Sample sites') #h_ini[sample_ID_drill]
main_ax.set_ylim([0,2])
main_ax.set_ylabel('Soil thickness at field [m]', fontsize=16)
main_ax.set_xlabel('Transport rate [m/yr]', fontsize=16)
#main_ax.legend(loc=4)
## histogram on the attached axes
x_hist.hist(d_eta_1yr[validID_hill], 40, histtype='stepfilled',
            orientation='vertical', color='gray')
#x_hist.invert_yaxis()
x_hist.label_outer()

y_hist.hist(h_ini[validID], 40, histtype='stepfilled',
            orientation='horizontal', color='gray')
#y_hist.invert_xaxis()
y_hist.label_outer()
plt.gcf().subplots_adjust(left=0.15)
plt.tight_layout()
#fig.savefig('St_vs_surf_transp_'+str(simu_side))

#---------------------
# soil formation rate
# this is the rate of soil foramtion
# it is rho_r/rho_s times faster than the bedrock erosion rate
#---------------------
#h_1D = np.linspace(0,1.5,30)
#P_1D = B*np.exp(-m*h_1D)
#
#fig, ax  = plt.subplots(figsize=(5,9), dpi=50)     
#plt.rc("font", size=18)
#plt.plot(P_1D*1000, h_1D,'k')
#plt.gca().invert_yaxis()
#ax.set_ylim([1.5,0])
#ax.set_xlim([0.0,0.35])
#ax.set_xlabel('Soil formation rate [mm/yr]', fontsize=18)
#ax.set_ylabel('Soil depth [m]', fontsize=18)

#----------------------------------------
# box-plot of N-f and S-f data from drill
#----------------------------------------
#fig, ax  = plt.subplots(figsize=(10,6), dpi=50)   
#data = [st_field_drill_S,st_field_drill_N]
#ax.boxplot(data)
#ax.set_ylim([-0.1,1.5])
#ax.set_xticklabels(['South-facing', 'North-facing'])