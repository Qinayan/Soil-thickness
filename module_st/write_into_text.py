#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 11:24:28 2019

@author: qinayan
"""

#https://www.pythonforbeginners.com/files/reading-and-writing-files-in-python

#https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python


import numpy as np



ho    = 0.135#0.2 # 1/m
B_p   = 0.000595 #(rho_r/rho_s)*Po
h_bar = 0.54
curv_thre = -0.0047#-0.0014 #-0.0026

# soil diffusion parameters
Klx =  5e-3 # per year, or 1.0*0.036


file = open("input.txt","w")
file.write("paramater symbol value\n")


lines_of_text = ["formation_power ho 0.135\n", \
                 "maximum_formation_rate Bp 0.000595\n", \
                 "mean_soil_thickness h_bar 0.54\n", \
                 "curvature_threshold curv_thre -0.0047\n",\
                 "diffusion_coefficient Klx 0.005\n"] 
file.writelines(lines_of_text) 

file.close()