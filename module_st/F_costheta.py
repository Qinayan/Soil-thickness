#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 14:50:31 2019

@author: qinayan
"""
import os
def costheta(elev_vector,dx,dy,A_p1,A_m1,A_pN,A_mN):
    current_dir = os.getcwd()
    eta_vector_p1=A_p1.dot(elev_vector)
    eta_vector_m1=A_m1.dot(elev_vector)
    eta_vector_pN=A_pN.dot(elev_vector)
    eta_vector_mN=A_mN.dot(elev_vector)
    
    d_eta_x = eta_vector_p1-eta_vector_m1
    d_eta_y = eta_vector_pN-eta_vector_mN
    
    cosx=2.0*dx/(4.0*dx**2+d_eta_x**2)**0.5
    cosy=2.0*dy/(4.0*dy**2+d_eta_y**2)**0.5
    cos_theta = cosx*cosy # cos theta, surface angle
    os.chdir(current_dir)
    return (cos_theta)