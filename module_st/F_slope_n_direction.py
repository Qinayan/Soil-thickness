import numpy as np
#from LEM_initial_landscape import *

def slope_direction(eta_vector,nrows,ncols,dx,dy,validID,bn_ID):
    #neighbors
    z = eta_vector.reshape(nrows,ncols)
    xn = [-1,0,1,-1,1,-1, 0, 1]
    yn = [1, 1,1, 0,0,-1,-1,-1]
    dn = [2.**.5,1.,2.**.5,1.,1.,2.**.5,1.,2.**.5]
    
    # 530
    # 6 1
    # 742
    
    delta_z = np.zeros((nrows,ncols)) #[[0. for i in xrange(grid_cells+2)]for j in xrange(grid_cells+2)]
    slope = np.zeros((nrows,ncols))
    direction = np.zeros((nrows,ncols),dtype=np.int)
#    hole = np.zeros((nrows,ncols))
    hole=np.zeros(1)
    
#===========================
# outlet elevation drops...
#===========================     
#    o_nrow            = int(np.argmin(z)/ncols)
#    o_ncol            = int(np.argmin(z)-o_nrow*ncols)
#    z[o_nrow][o_ncol] = z[o_nrow][o_ncol]*1.0-2.0
    
#===========================
# slope inside boundary
#===========================
#    validID
    #    for i in range (1, nrows-1):
#        for j in range (1,ncols-1):
    validID_insideID = np.where((validID>ncols-1)&(validID<(nrows-1)*ncols)& (validID%ncols!=0)& ((validID+1)%(ncols)!=0))[0]
    validID_inside = validID[validID_insideID]
    
    for k in range(len(validID_inside)): #len(validID_inside)
        i = validID_inside[k]//ncols
        j = validID_inside[k] - i*ncols
            
        slope_temp = np.zeros(8)
        for n in range (0,8):
            slope_temp[n]=(z[i][j]-z[i+xn[n]][j+yn[n]])/(dx*dn[n])
        direction[i][j]=np.argmax(slope_temp)
        slope[i][j]= np.max(slope_temp) # or slope_temp[direction[x][y]]
        
        if slope[i][j]<0.0:
            n = direction[i][j]
            delta_z[i][j]= - dx*dn[n]*slope[i][j]
            z[i][j]=z[i][j]+delta_z[i][j]
            slope[i][j]=0.0	
            hole[0]=1
##                n1, n2= slope_temp[np.argpartition(slope_temp, -2)][-2:]
#                n1, n2=slope_temp.argsort()[-2:][::-1]
#                slope2 = slope_temp[n2]
#                delta_z[i][j]= - ((dx*dn[n1]*slope[i][j])+(dx*dn[n2]*slope2-dx*dn[n1]*slope[i][j])/4.0)
#                z[i][j]=z[i][j]+delta_z[i][j]
#                slope[i][j]=(z[i][j]-z[i+xn[n1]][j+yn[n1]])/(dx*dn[n1])
#                hole[0]=1

#===========================
# up boundary
#===========================    
    order =np.array([1,2,4,6,7])
    xn_u = [0,1,1, 0, 1]
    yn_u = [1,1,0,-1,-1]
    dn_u = [1.,2.**.5,1.,1.,2.**.5]
    slope_temp = np.zeros(5)
    for j in range (1,ncols-1):
        
        for n in range (0,5):
            slope_temp[n]=(z[0][j]-z[0+xn_u[n]][j+yn_u[n]])/(dx*dn_u[n])
        direction_temp=np.argmax(slope_temp)
        direction[0][j] = order[direction_temp]
        slope[0][j]= np.max(slope_temp) # or slope_temp[direction[x][y]]
            
        if slope[0][j]<0.0:
            n = direction[0][j]
            delta_z[0][j]= - dx*dn[n]*slope[0][j]
            z[0][j]=z[0][j]+delta_z[0][j]
            slope[0][j]=0.0	
            hole[0]=1
            
##===========================
## the real boundary of the dem polygon
##===========================               
    direction_vector = direction.flatten()*1
    direction_vector[bn_ID] = 4
    direction = direction_vector.reshape(nrows,ncols)
    
#=========================================
# down boundary & boundary condition
#=========================================      
    order = np.array([6,5,3,1,0])
    for j in range (1,ncols-1):
        
        for n in range (0,5):
            slope_temp[n]=(z[nrows-1][j]-z[nrows-1-xn_u[n]][j-yn_u[n]])/(dx*dn_u[n])
        direction_temp = np.argmax(slope_temp)
        direction[nrows-1][j]=order[direction_temp]
        slope[nrows-1][j]= np.max(slope_temp) # or slope_temp[direction[x][y]]
            
        if slope[nrows-1][j]<0.0:
#            direction[nrows-1][j]= 8
            n = direction[nrows-1][j]
            delta_z[nrows-1][j]= - dx*dn[n]*slope[nrows-1][j]
            z[nrows-1][j]=z[nrows-1][j]+delta_z[nrows-1][j]
            slope[nrows-1][j]=0.0	
            hole[0]=1
#            slope[nrows-1][j]= outlet_slope #0.005   
            
#            
            
#    for j in xrange (1,ncols-1):
#        
##        for n in xrange (0,5):
##            slope_temp[n]=(z[nrows-1][j]-z[nrows-1-xn_u[n]][j-yn_u[n]])/(dx*dn_u[n])
##        direction_temp = np.argmax(slope_temp)
#        direction[nrows-1][j]= 8 #order[direction_temp]
#        slope[nrows-1][j]= outlet_slope# np.max(slope_temp) # or slope_temp[direction[x][y]]
#                   
#===========================
# left boundary
#===========================  
    for i in range(1, nrows-1):
        
        for n in range (0,5):
            slope_temp[n]=(z[i][0]-z[i+xn[n]][0+yn[n]])/(dx*dn[n])
        direction[i][0]=np.argmax(slope_temp)
        slope[i][0]= np.max(slope_temp) # or slope_temp[direction[x][y]]
        
        if slope[i][0]<0.0:
            n = direction[i][0]
            delta_z[i][0]= - dx*dn[n]*slope[i][0]
            z[i][0]=z[i][0]+delta_z[i][0]
            slope[i][0]=0.0	
            hole[0]=1
#===========================
# right boundary
#===========================   
    for i in range(1, nrows-1):
        
        for n in range (0,5):
            slope_temp[n]=(z[i][ncols-1]-z[i+xn[n+3]][ncols-1+yn[n+3]])/(dx*dn[n+3])
        direction[i][ncols-1]=np.argmax(slope_temp)+3
        slope[i][ncols-1]= np.max(slope_temp) # or slope_temp[direction[x][y]]

        if slope[i][ncols-1]<0.0:
            n = direction[i][ncols-1]
            delta_z[i][ncols-1]= - dx*dn[n]*slope[i][ncols-1]
            z[i][ncols-1]=z[i][ncols-1]+delta_z[i][ncols-1]
            slope[i][ncols-1]=0.0	
            hole[0]=1    
#===========================
# Up-left boundary
#=========================== 
    slope_temp = np.zeros(3)
    order_c = np.array([1,2,4])
    xn_ul = [0,1,1]
    yn_ul = [1,1, 0]
    dn_ul = [1.,2.**.5,1.]
    for n in range (0,3):
        slope_temp[n]=(z[0][0]-z[0+xn_ul[n]][0+yn_ul[n]])/(dx*dn_ul[n])
    direction_temp = np.argmax(slope_temp)        
    direction[0][0]=order_c[direction_temp]
    slope[0][0]= np.max(slope_temp) # or slope_temp[direction[x][y]]
    
    if slope[0][0]<0.0:
        n = direction[0][0]
        delta_z[0][0]= - dx*dn[n]*slope[0][0]
        z[0][0]=z[0][0]+delta_z[0][0]
        slope[0][0]=0.0	
        hole[0]=1                            
#===========================
# Up-right boundary
#===========================  
    order_c = np.array([6,7,4])
    for n in range (0,3):
        slope_temp[n]=(z[0][ncols-1]-z[0+xn_ul[n]][ncols-1-yn_ul[n]])/(dx*dn_ul[n])
    direction_temp = np.argmax(slope_temp) 
    direction[0][ncols-1]=order_c[direction_temp]
    slope[0][ncols-1]= np.max(slope_temp) # or slope_temp[direction[x][y]]
    
    if slope[0][ncols-1]<0.0:
        n = direction[0][ncols-1]
        delta_z[0][ncols-1]= - dx*dn[n]*slope[0][ncols-1]
        z[0][ncols-1]=z[0][ncols-1]+delta_z[0][ncols-1]
        slope[0][ncols-1]=0.0	
        hole[0]=1 
#===========================
# down_left boundary
#===========================  
    order_c = np.array([1,0,3])
    for n in range (0,3):
        slope_temp[n]=(z[nrows-1][0]-z[nrows-1-xn_ul[n]][0+yn_ul[n]])/(dx*dn_ul[n])
    direction_temp = np.argmax(slope_temp) 
    direction[nrows-1][0]=order_c[direction_temp]
    slope[nrows-1][0]= np.max(slope_temp) # or slope_temp[direction[x][y]]
    
    if slope[nrows-1][0]<0.0:
        n = direction[nrows-1][0]
        delta_z[nrows-1][0]= - dx*dn[n]*slope[nrows-1][0]
        z[nrows-1][0]=z[nrows-1][0]+delta_z[nrows-1][0]
        slope[nrows-1][0]=0.0	
        hole[0]=1 
        
        
    
#===========================
# down_right boundary
#===========================   
    order_c = np.array([6,5,3])
    slope_temp = np.zeros(3)
    for n in range (0,3):
        slope_temp[n]=(z[nrows-1][ncols-1]-z[nrows-1-xn_ul[n]][ncols-1-yn_ul[n]])/(dx*dn_ul[n])
    direction_temp = np.argmax(slope_temp) 
    direction[nrows-1][ncols-1]=order_c[direction_temp]
   
    slope[nrows-1][ncols-1]= np.max(slope_temp) # or slope_temp[direction[x][y]]
    
    if slope[nrows-1][ncols-1]<0.0:
        n = direction[nrows-1][ncols-1]
        delta_z[nrows-1][ncols-1]= - dx*dn[n]*slope[nrows-1][ncols-1]
        z[nrows-1][ncols-1]=z[nrows-1][ncols-1]+delta_z[nrows-1][ncols-1]
        slope[nrows-1][ncols-1]=0.0	
        hole[0]=1       
        
#===========================
# outlet
#===========================     
#    slope[o_nrow][o_ncol] = np.max(slope)
#    direction[o_nrow][o_ncol] = 8
    
#    slope_vector = np.array(slope).flatten()
#    delta_z_vector = np.array(delta_z).flatten()

    return slope, direction


#if __name__ == "__main__":
#    X,Y,Z = ini_landscape()
#    
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
#    
#    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
#                       linewidth=0, antialiased=False)           
#    ax.zaxis.set_major_locator(LinearLocator(10))
#    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#    
#    fig.colorbar(surf, shrink=0.5, aspect=10)
#    plt.show()