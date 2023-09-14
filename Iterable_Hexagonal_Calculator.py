# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:14:46 2023

@author: repha
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 12:53:45 2023

@author: repha
"""

import numpy as np
import pandas as pd
import math
import time 

tic = time.time()
import tkinter as tk
from tkinter import filedialog
root = tk.Tk()
root.withdraw()
root.attributes("-topmost", True)
directory = filedialog.askdirectory()
root.destroy()

kansas = pd.read_csv(directory+'/Slip_Angle_Ti02_New.csv', sep = ',', header = None) #The slip band and angle of slip band from ImageJ (angle from bottom up)
vermont = pd.read_csv(directory + '/List_of_Grains.csv',header = None) #The exported grains file from OIM w/ Euler Angles
pennsylvania = pd.read_csv(directory + '/Slip Grains_Ti02_New.csv') #Slip bands and the grain they exist in

a = (1+math.sqrt(3))
b = 1+(math.sqrt(3)/2)
c = -1.59
d = math.sqrt(1+0+c**2)
e = math.sqrt(3)/2
f = 0.877678086
g = 0.760091518
h = 0.438839043
i = 0.479250642
j = 0.532389592
k = 0.266194796
l = 0.461062912
m = -0.846499452
n = -0.001020592

constants = [a,b,c,d,e,f,g,h,i,j,k,l,m,n]


basal = [[0, 0, 0, 1],[0,0,0,1],[0,0,0,1]]
prismatic = [[1,0,-1,0],[0,1,-1,0],[-1,1,0,0]]
pyramidal = [[1,0,-1,1],[0,-1,1,1],[1,-1,0,1],[1,0,-1,-1],[0,-1,1,-1],[1,-1,0,-1]]
pyramidal_ca = [[1,0,-1,1],[1,0,-1,1],[-1,0,1,1],[-1,0,1,1],[0,1,-1,1],[0,1,-1,1],[0,-1,1,1],[0,-1,1,1],[1,-1,0,1],[1,-1,0,1],[-1,1,0,1],[-1,1,0,1]]

basal_direction = [[1,1,-2],[1,-2,1],[-2,1,1]]
prismatic_direction = [[1,-2,1],[2,-1,-1],[1,1,-2]]
pyramidal_direction = [[1,-2,1],[2,-1,-1],[1,1,-2],[1,-2,1],[2,-1,-1],[1,1-2]]
pyramidal_ca_direction = [[1,-2,1],[2,-1,-1],[1,1,-2],[1,-2,1],[2,-1,-1],[1,1,-2],[1,-2,1],[2,-1,-1],[1,1,-2],[1,-2,1],[2,-1,-1],[1,1,-2]]

norme = [[0,1/3]]

basal_B = [[1,-0.5,1,0],[(1+math.sqrt(3)),(1+(math.sqrt(3)/2)),1,0],[0,0,0,1]]
prismatic_B = [[1,-0.5,1,0],[(1+math.sqrt(3)),(1+(math.sqrt(3)/2)),1,0],[0,0,0,1]]

pyramidal_B = [[1,-0.5,1,0],
               [a,b,1,0],
               [0,0,0,1],
               [1,-0.5,1,0],
               [a,b,1,0],
               [0,0,0,1]]

pyramidal_ca_B = [[1,0,c,d],[0.5,e,c,d],[-1,0,c,d],[-0.5,-e,c,d],
                  [-0.5,e,c,d],[0.5,e,c,d],[0.5,-e,c,d],
                  [-0.5,-e,c,d],[1,0,c,d],[0.5,-e,c,d],
                  [-1,0,c,d],[-0.5,e,c,d]]

basal_plan = [[0,0,1],[0,0,1],[0,0,1]]
prismatic_plan = [[0,1,0],[-e,0.5,0],[-e,-0.5,0]]
pyramidal_plan = [[0,-f,i,1],
                  [g,-h,-i,1],
                  [g,h,-i,1],
                  [0,f,i,1],
                  [-g,h,-i,1],
                  [-g,-h,-i,1]]
pyramidal_ca_plan = [[g,h,i,1],
                     [g,h,i,1],
                     [-g,-h,i,1],
                     [-g,-h,i,1],
                     [0,f,i,1],
                     [0,f,i,1],
                     [0,-f,i,1],
                     [0,-f,i,1],
                     [g,-h,i,1],
                     [g,-h,i,1],
                     [g,-h,-i,1],
                     [g,-h,-i,1]]

basal_direction2 = [[-0.5,e,0],[1,0,0],[-0.5,-e,0]]
prismatic_direction2 = [[1,0,0],[0.5,e,0],[-0.5,e,0]]
pyramidal_direction2 = [[1,0,0],
                        [0.5,e,0],
                        [-0.5,e,0],
                        [1,0,0],
                        [0.5,e,0],
                        [-0.5,e,0]]
pyramidal_ca_direction2 = [[j,0,m,n],
                           [k,l,m,n],
                           [-j,0,m,n],
                           [-k,-l,m,n],
                           [-k,l,m,n],
                           [k,l,m,n],
                           [k,-l,m,n],
                           [-k,-l,m,n],
                           [j,0,m,n],
                           [k,-l,m,n],
                           [-j,0,m,-n],
                           [-k,l,m,-n]]


penny = np.asarray(pennsylvania)
philly = np.asarray(vermont)
penny = np.asarray(pennsylvania)
philly = np.asarray(vermont)
philly = philly[:,0] #List
valentine = penny[:,1] #My grains
imdex = []
for i in range(len(penny)):
    nike = np.where(valentine[i] == philly)
    imdex.append(nike)
cherry = []
for i in range(len(imdex)):
    umdex = imdex[i][0][0]
    cherry.append(umdex)
oak = vermont.iloc[cherry]
ridge = np.array(oak)
slip_band = kansas[0][1:]
maine = []
for j in range(1,len(kansas)):
    angle_real = float(kansas[1][j])
    phasetype = float(ridge[j-1][1])
    eulerangle_deg1 = float(ridge[j-1][2])
    eulerangle_deg2 = float(ridge[j-1][3])
    eulerangle_deg3 = float(ridge[j-1][4])
    eulerangle_rad1 = eulerangle_deg1*math.pi/180
    eulerangle_rad2 = eulerangle_deg2*math.pi/180
    eulerangle_rad3 = eulerangle_deg3*math.pi/180
    erad1 = eulerangle_rad1
    erad2 = eulerangle_rad2
    erad3 = eulerangle_rad3

    matricerotation = [[math.cos(erad1)*math.cos(erad3) - math.sin(erad1)*math.cos(erad2)*math.sin(erad3),-1*math.cos(erad1)*math.sin(erad3)-math.sin(erad1)*math.cos(erad2)*math.cos(erad3),math.sin(erad1)*(math.sin(erad2))],
                       [math.sin(erad1)*math.cos(erad3)+math.cos(erad1)*math.cos(erad2)*math.sin(erad3),-1*math.sin(erad1)*math.sin(erad3)+math.cos(erad1)*math.cos(erad2)*math.cos(erad3),-1*math.cos(erad1)*math.sin(erad2)],
                       [math.sin(erad3)*math.sin(erad2),math.sin(erad2)*math.cos(erad3),math.cos(erad2)]]
    mr = matricerotation
    basal_planoriente = [[mr[0][0]*basal_plan[0][0]+mr[0][1]*basal_plan[0][1]+ mr[0][2]*basal_plan[0][2],
                          mr[1][0]*basal_plan[0][0]+mr[1][1]*basal_plan[0][1]+mr[1][2]*basal_plan[0][2],
                          mr[2][0]*basal_plan[0][0]+mr[2][1]*basal_plan[0][1]+mr[2][2]*basal_plan[0][2]],
                         [mr[0][0]*basal_plan[1][0]+mr[0][1]*basal_plan[1][1]+ mr[0][2]*basal_plan[1][2],
                          mr[1][0]*basal_plan[1][0]+mr[1][1]*basal_plan[1][1]+mr[1][2]*basal_plan[1][2],
                          mr[2][0]*basal_plan[1][0]+mr[2][1]*basal_plan[1][1]+mr[2][2]*basal_plan[1][2]],
                         [mr[0][0]*basal_plan[2][0]+mr[0][1]*basal_plan[2][1]+ mr[0][2]*basal_plan[2][2],
                          mr[1][0]*basal_plan[2][0]+mr[1][1]*basal_plan[2][1]+mr[1][2]*basal_plan[2][2],
                          mr[2][0]*basal_plan[2][0]+mr[2][1]*basal_plan[2][1]+mr[2][2]*basal_plan[2][2]]]

    prismatic_planoriente = [[mr[0][0]*prismatic_plan[0][0]+mr[0][1]*prismatic_plan[0][1]+ mr[0][2]*prismatic_plan[0][2],
                              mr[1][0]*prismatic_plan[0][0]+mr[1][1]*prismatic_plan[0][1]+mr[1][2]*prismatic_plan[0][2],
                              mr[2][0]*prismatic_plan[0][0]+mr[2][1]*prismatic_plan[0][1]+mr[2][2]*prismatic_plan[0][2]],
                             [mr[0][0]*prismatic_plan[1][0]+mr[0][1]*prismatic_plan[1][1]+ mr[0][2]*prismatic_plan[1][2],
                              mr[1][0]*prismatic_plan[1][0]+mr[1][1]*prismatic_plan[1][1]+mr[1][2]*prismatic_plan[1][2],
                              mr[2][0]*prismatic_plan[1][0]+mr[2][1]*prismatic_plan[1][1]+mr[2][2]*prismatic_plan[1][2]],
                             [mr[0][0]*prismatic_plan[2][0]+mr[0][1]*prismatic_plan[2][1]+ mr[0][2]*prismatic_plan[2][2],
                              mr[1][0]*prismatic_plan[2][0]+mr[1][1]*prismatic_plan[2][1]+mr[1][2]*prismatic_plan[2][2],
                              mr[2][0]*prismatic_plan[2][0]+mr[2][1]*prismatic_plan[2][1]+mr[2][2]*prismatic_plan[2][2]]]

    pyramidal_planoriente = [[mr[0][0]*pyramidal_plan[0][0]+mr[0][1]*pyramidal_plan[0][1]+ mr[0][2]*pyramidal_plan[0][2],
                              mr[1][0]*pyramidal_plan[0][0]+mr[1][1]*pyramidal_plan[0][1]+mr[1][2]*pyramidal_plan[0][2],
                              mr[2][0]*pyramidal_plan[0][0]+mr[2][1]*pyramidal_plan[0][1]+mr[2][2]*pyramidal_plan[0][2]],
                             [mr[0][0]*pyramidal_plan[1][0]+mr[0][1]*pyramidal_plan[1][1]+ mr[0][2]*pyramidal_plan[1][2],
                              mr[1][0]*pyramidal_plan[1][0]+mr[1][1]*pyramidal_plan[1][1]+mr[1][2]*pyramidal_plan[1][2],
                          mr[2][0]*pyramidal_plan[1][0]+mr[2][1]*pyramidal_plan[1][1]+mr[2][2]*pyramidal_plan[1][2]],
                            [mr[0][0]*pyramidal_plan[2][0]+mr[0][1]*pyramidal_plan[2][1]+ mr[0][2]*pyramidal_plan[2][2],
                             mr[1][0]*pyramidal_plan[2][0]+mr[1][1]*pyramidal_plan[2][1]+mr[1][2]*pyramidal_plan[2][2],
                             mr[2][0]*pyramidal_plan[2][0]+mr[2][1]*pyramidal_plan[2][1]+mr[2][2]*pyramidal_plan[2][2]],
                            [mr[0][0]*pyramidal_plan[3][0]+mr[0][1]*pyramidal_plan[3][1]+ mr[0][2]*pyramidal_plan[3][2],
                             mr[1][0]*pyramidal_plan[3][0]+mr[1][1]*pyramidal_plan[3][1]+mr[1][2]*pyramidal_plan[3][2],
                             mr[2][0]*pyramidal_plan[3][0]+mr[2][1]*pyramidal_plan[3][1]+mr[2][2]*pyramidal_plan[3][2]],
                            [mr[0][0]*pyramidal_plan[4][0]+mr[0][1]*pyramidal_plan[4][1]+ mr[0][2]*pyramidal_plan[4][2],
                             mr[1][0]*pyramidal_plan[4][0]+mr[1][1]*pyramidal_plan[4][1]+mr[1][2]*pyramidal_plan[4][2],
                             mr[2][0]*pyramidal_plan[4][0]+mr[2][1]*pyramidal_plan[4][1]+mr[2][2]*pyramidal_plan[4][2]],
                            [mr[0][0]*pyramidal_plan[5][0]+mr[0][1]*pyramidal_plan[5][1]+ mr[0][2]*pyramidal_plan[5][2],
                             mr[1][0]*pyramidal_plan[5][0]+mr[1][1]*pyramidal_plan[5][1]+mr[1][2]*pyramidal_plan[5][2],
                             mr[2][0]*pyramidal_plan[5][0]+mr[2][1]*pyramidal_plan[5][1]+mr[2][2]*pyramidal_plan[5][2]]]

    pyramidal_ca_planoriente = [[mr[0][0]*pyramidal_ca_plan[0][0]+mr[0][1]*pyramidal_ca_plan[0][1]+ mr[0][2]*pyramidal_ca_plan[0][2],
                                 mr[1][0]*pyramidal_ca_plan[0][0]+mr[1][1]*pyramidal_ca_plan[0][1]+mr[1][2]*pyramidal_ca_plan[0][2],
                                 mr[2][0]*pyramidal_ca_plan[0][0]+mr[2][1]*pyramidal_ca_plan[0][1]+mr[2][2]*pyramidal_ca_plan[0][2]],
                                [mr[0][0]*pyramidal_ca_plan[1][0]+mr[0][1]*pyramidal_ca_plan[1][1]+ mr[0][2]*pyramidal_ca_plan[1][2],
                                 mr[1][0]*pyramidal_ca_plan[1][0]+mr[1][1]*pyramidal_ca_plan[1][1]+mr[1][2]*pyramidal_ca_plan[1][2],
                                 mr[2][0]*pyramidal_ca_plan[1][0]+mr[2][1]*pyramidal_ca_plan[1][1]+mr[2][2]*pyramidal_ca_plan[1][2]],
                                [mr[0][0]*pyramidal_ca_plan[2][0]+mr[0][1]*pyramidal_ca_plan[2][1]+ mr[0][2]*pyramidal_ca_plan[2][2],
                                 mr[1][0]*pyramidal_ca_plan[2][0]+mr[1][1]*pyramidal_ca_plan[2][1]+mr[1][2]*pyramidal_ca_plan[2][2],
                                 mr[2][0]*pyramidal_ca_plan[2][0]+mr[2][1]*pyramidal_ca_plan[2][1]+mr[2][2]*pyramidal_ca_plan[2][2]],
                                [mr[0][0]*pyramidal_ca_plan[3][0]+mr[0][1]*pyramidal_ca_plan[3][1]+ mr[0][2]*pyramidal_ca_plan[3][2],
                                 mr[1][0]*pyramidal_ca_plan[3][0]+mr[1][1]*pyramidal_ca_plan[3][1]+mr[1][2]*pyramidal_ca_plan[3][2],
                                 mr[2][0]*pyramidal_ca_plan[3][0]+mr[2][1]*pyramidal_ca_plan[3][1]+mr[2][2]*pyramidal_ca_plan[3][2]],
                                [mr[0][0]*pyramidal_ca_plan[4][0]+mr[0][1]*pyramidal_ca_plan[4][1]+ mr[0][2]*pyramidal_ca_plan[4][2],
                                 mr[1][0]*pyramidal_ca_plan[4][0]+mr[1][1]*pyramidal_ca_plan[4][1]+mr[1][2]*pyramidal_ca_plan[4][2],
                                 mr[2][0]*pyramidal_ca_plan[4][0]+mr[2][1]*pyramidal_ca_plan[4][1]+mr[2][2]*pyramidal_ca_plan[4][2]],
                                [mr[0][0]*pyramidal_ca_plan[5][0]+mr[0][1]*pyramidal_ca_plan[5][1]+ mr[0][2]*pyramidal_ca_plan[5][2],
                                 mr[1][0]*pyramidal_ca_plan[5][0]+mr[1][1]*pyramidal_ca_plan[5][1]+mr[1][2]*pyramidal_ca_plan[5][2],
                                 mr[2][0]*pyramidal_ca_plan[5][0]+mr[2][1]*pyramidal_ca_plan[5][1]+mr[2][2]*pyramidal_ca_plan[5][2]],
                                [mr[0][0]*pyramidal_ca_plan[6][0]+mr[0][1]*pyramidal_ca_plan[6][1]+ mr[0][2]*pyramidal_ca_plan[6][2],
                                 mr[1][0]*pyramidal_ca_plan[6][0]+mr[1][1]*pyramidal_ca_plan[6][1]+mr[1][2]*pyramidal_ca_plan[6][2],
                                 mr[2][0]*pyramidal_ca_plan[6][0]+mr[2][1]*pyramidal_ca_plan[6][1]+mr[2][2]*pyramidal_ca_plan[6][2]],
                                [mr[0][0]*pyramidal_ca_plan[7][0]+mr[0][1]*pyramidal_ca_plan[7][1]+ mr[0][2]*pyramidal_ca_plan[7][2],
                                 mr[1][0]*pyramidal_ca_plan[7][0]+mr[1][1]*pyramidal_ca_plan[7][1]+mr[1][2]*pyramidal_ca_plan[7][2],
                                 mr[2][0]*pyramidal_ca_plan[7][0]+mr[2][1]*pyramidal_ca_plan[7][1]+mr[2][2]*pyramidal_ca_plan[7][2]],
                                [mr[0][0]*pyramidal_ca_plan[8][0]+mr[0][1]*pyramidal_ca_plan[8][1]+ mr[0][2]*pyramidal_ca_plan[8][2],
                                 mr[1][0]*pyramidal_ca_plan[8][0]+mr[1][1]*pyramidal_ca_plan[8][1]+mr[1][2]*pyramidal_ca_plan[8][2],
                                 mr[2][0]*pyramidal_ca_plan[8][0]+mr[2][1]*pyramidal_ca_plan[8][1]+mr[2][2]*pyramidal_ca_plan[8][2]],
                                [mr[0][0]*pyramidal_ca_plan[9][0]+mr[0][1]*pyramidal_ca_plan[9][1]+ mr[0][2]*pyramidal_ca_plan[9][2],
                                 mr[1][0]*pyramidal_ca_plan[9][0]+mr[1][1]*pyramidal_ca_plan[9][1]+mr[1][2]*pyramidal_ca_plan[9][2],
                                 mr[2][0]*pyramidal_ca_plan[9][0]+mr[2][1]*pyramidal_ca_plan[9][1]+mr[2][2]*pyramidal_ca_plan[9][2]],
                                [mr[0][0]*pyramidal_ca_plan[10][0]+mr[0][1]*pyramidal_ca_plan[10][1]+ mr[0][2]*pyramidal_ca_plan[10][2],
                                 mr[1][0]*pyramidal_ca_plan[10][0]+mr[1][1]*pyramidal_ca_plan[10][1]+mr[1][2]*pyramidal_ca_plan[10][2],
                                 mr[2][0]*pyramidal_ca_plan[10][0]+mr[2][1]*pyramidal_ca_plan[10][1]+mr[2][2]*pyramidal_ca_plan[10][2]],
                                [mr[0][0]*pyramidal_ca_plan[11][0]+mr[0][1]*pyramidal_ca_plan[11][1]+ mr[0][2]*pyramidal_ca_plan[11][2],
                                 mr[1][0]*pyramidal_ca_plan[11][0]+mr[1][1]*pyramidal_ca_plan[11][1]+mr[1][2]*pyramidal_ca_plan[11][2],
                                 mr[2][0]*pyramidal_ca_plan[11][0]+mr[2][1]*pyramidal_ca_plan[11][1]+mr[2][2]*pyramidal_ca_plan[11][2]]]
    
    basal_directionorientee = [[mr[0][0]*basal_direction2[0][0]+mr[0][1]*basal_direction2[0][1]+ mr[0][2]*basal_direction2[0][2],
                                mr[1][0]*basal_direction2[0][0]+mr[1][1]*basal_direction2[0][1]+mr[1][2]*basal_direction2[0][2],
                                mr[2][0]*basal_direction2[0][0]+mr[2][1]*basal_direction2[0][1]+mr[2][2]*basal_direction2[0][2]],
                               [mr[0][0]*basal_direction2[1][0]+mr[0][1]*basal_direction2[1][1]+ mr[0][2]*basal_direction2[1][2],
                                mr[1][0]*basal_direction2[1][0]+mr[1][1]*basal_direction2[1][1]+mr[1][2]*basal_direction2[1][2],
                                mr[2][0]*basal_direction2[1][0]+mr[2][1]*basal_direction2[1][1]+mr[2][2]*basal_direction2[1][2]],
                               [mr[0][0]*basal_direction2[2][0]+mr[0][1]*basal_direction2[2][1]+ mr[0][2]*basal_direction2[2][2],
                                mr[1][0]*basal_direction2[2][0]+mr[1][1]*basal_direction2[2][1]+mr[1][2]*basal_direction2[2][2],
                                mr[2][0]*basal_direction2[2][0]+mr[2][1]*basal_direction2[2][1]+mr[2][2]*basal_direction2[2][2]]]
    
    prismatic_directionorientee = [[mr[0][0]*prismatic_direction2[0][0]+mr[0][1]*prismatic_direction2[0][1]+ mr[0][2]*prismatic_direction2[0][2],
                                    mr[1][0]*prismatic_direction2[0][0]+mr[1][1]*prismatic_direction2[0][1]+mr[1][2]*prismatic_direction2[0][2],
                                    mr[2][0]*prismatic_direction2[0][0]+mr[2][1]*prismatic_direction2[0][1]+mr[2][2]*prismatic_direction2[0][2]],
                                   [mr[0][0]*prismatic_direction2[1][0]+mr[0][1]*prismatic_direction2[1][1]+ mr[0][2]*prismatic_direction2[1][2],
                                    mr[1][0]*prismatic_direction2[1][0]+mr[1][1]*prismatic_direction2[1][1]+mr[1][2]*prismatic_direction2[1][2],
                                    mr[2][0]*prismatic_direction2[1][0]+mr[2][1]*prismatic_direction2[1][1]+mr[2][2]*prismatic_direction2[1][2]],
                                   [mr[0][0]*prismatic_direction2[2][0]+mr[0][1]*prismatic_direction2[2][1]+ mr[0][2]*prismatic_direction2[2][2],
                                    mr[1][0]*prismatic_direction2[2][0]+mr[1][1]*prismatic_direction2[2][1]+mr[1][2]*prismatic_direction2[2][2],
                                    mr[2][0]*prismatic_direction2[2][0]+mr[2][1]*prismatic_direction2[2][1]+mr[2][2]*prismatic_direction2[2][2]]]
    
    pyramidal_directionorientee = [[mr[0][0]*pyramidal_direction2[0][0]+mr[0][1]*pyramidal_direction2[0][1]+ mr[0][2]*pyramidal_direction2[0][2],
                                    mr[1][0]*pyramidal_direction2[0][0]+mr[1][1]*pyramidal_direction2[0][1]+mr[1][2]*pyramidal_direction2[0][2],
                                    mr[2][0]*pyramidal_direction2[0][0]+mr[2][1]*pyramidal_direction2[0][1]+mr[2][2]*pyramidal_direction2[0][2]],
                                   [mr[0][0]*pyramidal_direction2[1][0]+mr[0][1]*pyramidal_direction2[1][1]+ mr[0][2]*pyramidal_direction2[1][2],
                                    mr[1][0]*pyramidal_direction2[1][0]+mr[1][1]*pyramidal_direction2[1][1]+mr[1][2]*pyramidal_direction2[1][2],
                                    mr[2][0]*pyramidal_direction2[1][0]+mr[2][1]*pyramidal_direction2[1][1]+mr[2][2]*pyramidal_direction2[1][2]],
                                   [mr[0][0]*pyramidal_direction2[2][0]+mr[0][1]*pyramidal_direction2[2][1]+ mr[0][2]*pyramidal_direction2[2][2],
                                    mr[1][0]*pyramidal_direction2[2][0]+mr[1][1]*pyramidal_direction2[2][1]+mr[1][2]*pyramidal_direction2[2][2],
                                    mr[2][0]*pyramidal_direction2[2][0]+mr[2][1]*pyramidal_direction2[2][1]+mr[2][2]*pyramidal_direction2[2][2]],
                                   [mr[0][0]*pyramidal_direction2[3][0]+mr[0][1]*pyramidal_direction2[3][1]+ mr[0][2]*pyramidal_direction2[3][2],
                                    mr[1][0]*pyramidal_direction2[3][0]+mr[1][1]*pyramidal_direction2[3][1]+mr[1][2]*pyramidal_direction2[3][2],
                                    mr[2][0]*pyramidal_direction2[3][0]+mr[2][1]*pyramidal_direction2[3][1]+mr[2][2]*pyramidal_direction2[3][2]],
                                   [mr[0][0]*pyramidal_direction2[4][0]+mr[0][1]*pyramidal_direction2[4][1]+ mr[0][2]*pyramidal_direction2[4][2],
                                    mr[1][0]*pyramidal_direction2[4][0]+mr[1][1]*pyramidal_direction2[4][1]+mr[1][2]*pyramidal_direction2[4][2],
                                    mr[2][0]*pyramidal_direction2[4][0]+mr[2][1]*pyramidal_direction2[4][1]+mr[2][2]*pyramidal_direction2[4][2]],
                                   [mr[0][0]*pyramidal_direction2[5][0]+mr[0][1]*pyramidal_direction2[5][1]+ mr[0][2]*pyramidal_direction2[5][2],
                                    mr[1][0]*pyramidal_direction2[5][0]+mr[1][1]*pyramidal_direction2[5][1]+mr[1][2]*pyramidal_direction2[5][2],
                                    mr[2][0]*pyramidal_direction2[5][0]+mr[2][1]*pyramidal_direction2[5][1]+mr[2][2]*pyramidal_direction2[5][2]]]
    
    pyramidal_ca_directionorientee = [[mr[0][0]*pyramidal_ca_direction2[0][0]+mr[0][1]*pyramidal_ca_direction2[0][1]+ mr[0][2]*pyramidal_ca_direction2[0][2],
                                       mr[1][0]*pyramidal_ca_direction2[0][0]+mr[1][1]*pyramidal_ca_direction2[0][1]+mr[1][2]*pyramidal_ca_direction2[0][2],
                                       mr[2][0]*pyramidal_ca_direction2[0][0]+mr[2][1]*pyramidal_ca_direction2[0][1]+mr[2][2]*pyramidal_ca_direction2[0][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[1][0]+mr[0][1]*pyramidal_ca_direction2[1][1]+ mr[0][2]*pyramidal_ca_direction2[1][2],
                                       mr[1][0]*pyramidal_ca_direction2[1][0]+mr[1][1]*pyramidal_ca_direction2[1][1]+mr[1][2]*pyramidal_ca_direction2[1][2],
                                       mr[2][0]*pyramidal_ca_direction2[1][0]+mr[2][1]*pyramidal_ca_direction2[1][1]+mr[2][2]*pyramidal_ca_direction2[1][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[2][0]+mr[0][1]*pyramidal_ca_direction2[2][1]+ mr[0][2]*pyramidal_ca_direction2[2][2],
                                       mr[1][0]*pyramidal_ca_direction2[2][0]+mr[1][1]*pyramidal_ca_direction2[2][1]+mr[1][2]*pyramidal_ca_direction2[2][2],
                                       mr[2][0]*pyramidal_ca_direction2[2][0]+mr[2][1]*pyramidal_ca_direction2[2][1]+mr[2][2]*pyramidal_ca_direction2[2][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[3][0]+mr[0][1]*pyramidal_ca_direction2[3][1]+ mr[0][2]*pyramidal_ca_direction2[3][2],
                                       mr[1][0]*pyramidal_ca_direction2[3][0]+mr[1][1]*pyramidal_ca_direction2[3][1]+mr[1][2]*pyramidal_ca_direction2[3][2],
                                       mr[2][0]*pyramidal_ca_direction2[3][0]+mr[2][1]*pyramidal_ca_direction2[3][1]+mr[2][2]*pyramidal_ca_direction2[3][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[4][0]+mr[0][1]*pyramidal_ca_direction2[4][1]+ mr[0][2]*pyramidal_ca_direction2[4][2],
                                       mr[1][0]*pyramidal_ca_direction2[4][0]+mr[1][1]*pyramidal_ca_direction2[4][1]+mr[1][2]*pyramidal_ca_direction2[4][2],
                                       mr[2][0]*pyramidal_ca_direction2[4][0]+mr[2][1]*pyramidal_ca_direction2[4][1]+mr[2][2]*pyramidal_ca_direction2[4][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[5][0]+mr[0][1]*pyramidal_ca_direction2[5][1]+ mr[0][2]*pyramidal_ca_direction2[5][2],
                                       mr[1][0]*pyramidal_ca_direction2[5][0]+mr[1][1]*pyramidal_ca_direction2[5][1]+mr[1][2]*pyramidal_ca_direction2[5][2],
                                       mr[2][0]*pyramidal_ca_direction2[5][0]+mr[2][1]*pyramidal_ca_direction2[5][1]+mr[2][2]*pyramidal_ca_direction2[5][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[6][0]+mr[0][1]*pyramidal_ca_direction2[6][1]+ mr[0][2]*pyramidal_ca_direction2[6][2],
                                       mr[1][0]*pyramidal_ca_direction2[6][0]+mr[1][1]*pyramidal_ca_direction2[6][1]+mr[1][2]*pyramidal_ca_direction2[6][2],
                                       mr[2][0]*pyramidal_ca_direction2[6][0]+mr[2][1]*pyramidal_ca_direction2[6][1]+mr[2][2]*pyramidal_ca_direction2[6][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[7][0]+mr[0][1]*pyramidal_ca_direction2[7][1]+ mr[0][2]*pyramidal_ca_direction2[7][2],
                                       mr[1][0]*pyramidal_ca_direction2[7][0]+mr[1][1]*pyramidal_ca_direction2[7][1]+mr[1][2]*pyramidal_ca_direction2[7][2],
                                       mr[2][0]*pyramidal_ca_direction2[7][0]+mr[2][1]*pyramidal_ca_direction2[7][1]+mr[2][2]*pyramidal_ca_direction2[7][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[8][0]+mr[0][1]*pyramidal_ca_direction2[8][1]+ mr[0][2]*pyramidal_ca_direction2[8][2],
                                       mr[1][0]*pyramidal_ca_direction2[8][0]+mr[1][1]*pyramidal_ca_direction2[8][1]+mr[1][2]*pyramidal_ca_direction2[8][2],
                                       mr[2][0]*pyramidal_ca_direction2[8][0]+mr[2][1]*pyramidal_ca_direction2[8][1]+mr[2][2]*pyramidal_ca_direction2[8][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[9][0]+mr[0][1]*pyramidal_ca_direction2[9][1]+ mr[0][2]*pyramidal_ca_direction2[9][2],
                                       mr[1][0]*pyramidal_ca_direction2[9][0]+mr[1][1]*pyramidal_ca_direction2[9][1]+mr[1][2]*pyramidal_ca_direction2[9][2],
                                       mr[2][0]*pyramidal_ca_direction2[9][0]+mr[2][1]*pyramidal_ca_direction2[9][1]+mr[2][2]*pyramidal_ca_direction2[9][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[10][0]+mr[0][1]*pyramidal_ca_direction2[10][1]+ mr[0][2]*pyramidal_ca_direction2[10][2],
                                       mr[1][0]*pyramidal_ca_direction2[10][0]+mr[1][1]*pyramidal_ca_direction2[10][1]+mr[1][2]*pyramidal_ca_direction2[10][2],
                                       mr[2][0]*pyramidal_ca_direction2[10][0]+mr[2][1]*pyramidal_ca_direction2[10][1]+mr[2][2]*pyramidal_ca_direction2[10][2]],
                                      [mr[0][0]*pyramidal_ca_direction2[11][0]+mr[0][1]*pyramidal_ca_direction2[11][1]+ mr[0][2]*pyramidal_ca_direction2[11][2],
                                       mr[1][0]*pyramidal_ca_direction2[11][0]+mr[1][1]*pyramidal_ca_direction2[11][1]+mr[1][2]*pyramidal_ca_direction2[11][2],
                                       mr[2][0]*pyramidal_ca_direction2[11][0]+mr[2][1]*pyramidal_ca_direction2[11][1]+mr[2][2]*pyramidal_ca_direction2[11][2]]]
    
    loadingdirection = [0,1,0]
    surface = [0,0,1]
    LD = loadingdirection
    bpo = basal_planoriente
    ppo = prismatic_planoriente
    pypo = pyramidal_planoriente
    pcapo = pyramidal_ca_planoriente
    bdo = basal_directionorientee
    pdo = prismatic_directionorientee
    pydo = pyramidal_directionorientee
    pcado = pyramidal_ca_directionorientee
    
    basal_SF = [[(bpo[0][0]*LD[0] + bpo[0][1]*LD[1] + bpo[0][2]*LD[2])*(bdo[0][0]*LD[0] + bdo[0][1]*LD[1] + bdo[0][2]*LD[2])],
                [(bpo[1][0]*LD[0] + bpo[1][1]*LD[1] + bpo[1][2]*LD[2])*(bdo[1][0]*LD[0] + bdo[1][1]*LD[1] + bdo[1][2]*LD[2])],
                [(bpo[2][0]*LD[0] + bpo[2][1]*LD[1] + bpo[2][2]*LD[2])*(bdo[2][0]*LD[0] + bdo[2][1]*LD[1] + bdo[2][2]*LD[2])]]
    
    basal_SF = np.absolute(basal_SF)
    
    prismatic_SF = [[(ppo[0][0]*LD[0] + ppo[0][1]*LD[1] + ppo[0][2]*LD[2])*(pdo[0][0]*LD[0] + pdo[0][1]*LD[1] + pdo[0][2]*LD[2])],
                    [(ppo[1][0]*LD[0] + ppo[1][1]*LD[1] + ppo[1][2]*LD[2])*(pdo[1][0]*LD[0] + pdo[1][1]*LD[1] + pdo[1][2]*LD[2])],
                    [(ppo[2][0]*LD[0] + ppo[2][1]*LD[1] + ppo[2][2]*LD[2])*(pdo[2][0]*LD[0] + pdo[2][1]*LD[1] + pdo[2][2]*LD[2])]]
    
    prismatic_SF = np.absolute(prismatic_SF)
    
    pyramidal_SF = [[(pypo[0][0]*LD[0] + pypo[0][1]*LD[1] + pypo[0][2]*LD[2])*(pydo[0][0]*LD[0] + pydo[0][1]*LD[1] + pydo[0][2]*LD[2])],
                    [(pypo[1][0]*LD[0] + pypo[1][1]*LD[1] + pypo[1][2]*LD[2])*(pydo[1][0]*LD[0] + pydo[1][1]*LD[1] + pydo[1][2]*LD[2])],
                    [(pypo[2][0]*LD[0] + pypo[2][1]*LD[1] + pypo[2][2]*LD[2])*(pydo[2][0]*LD[0] + pydo[2][1]*LD[1] + pydo[2][2]*LD[2])],
                    [(pypo[3][0]*LD[0] + pypo[3][1]*LD[1] + pypo[3][2]*LD[2])*(pydo[3][0]*LD[0] + pydo[3][1]*LD[1] + pydo[3][2]*LD[2])],
                    [(pypo[4][0]*LD[0] + pypo[4][1]*LD[1] + pypo[4][2]*LD[2])*(pydo[4][0]*LD[0] + pydo[4][1]*LD[1] + pydo[4][2]*LD[2])],
                    [(pypo[5][0]*LD[0] + pypo[5][1]*LD[1] + pypo[5][2]*LD[2])*(pydo[5][0]*LD[0] + pydo[5][1]*LD[1] + pydo[5][2]*LD[2])]]
    
    pyramidal_SF = np.absolute(pyramidal_SF)
    
    pyramidal_ca_SF = [[(pcapo[0][0]*LD[0] + pcapo[0][1]*LD[1] + pcapo[0][2]*LD[2])*(pcado[0][0]*LD[0] + pcado[0][1]*LD[1] + pcado[0][2]*LD[2])],
                       [(pcapo[1][0]*LD[0] + pcapo[1][1]*LD[1] + pcapo[1][2]*LD[2])*(pcado[1][0]*LD[0] + pcado[1][1]*LD[1] + pcado[1][2]*LD[2])],
                       [(pcapo[2][0]*LD[0] + pcapo[2][1]*LD[1] + pcapo[2][2]*LD[2])*(pcado[2][0]*LD[0] + pcado[2][1]*LD[1] + pcado[2][2]*LD[2])],
                       [(pcapo[3][0]*LD[0] + pcapo[3][1]*LD[1] + pcapo[3][2]*LD[2])*(pcado[3][0]*LD[0] + pcado[3][1]*LD[1] + pcado[3][2]*LD[2])],
                       [(pcapo[4][0]*LD[0] + pcapo[4][1]*LD[1] + pcapo[4][2]*LD[2])*(pcado[4][0]*LD[0] + pcado[4][1]*LD[1] + pcado[4][2]*LD[2])],
                       [(pcapo[5][0]*LD[0] + pcapo[5][1]*LD[1] + pcapo[5][2]*LD[2])*(pcado[5][0]*LD[0] + pcado[5][1]*LD[1] + pcado[5][2]*LD[2])],
                       [(pcapo[6][0]*LD[0] + pcapo[6][1]*LD[1] + pcapo[6][2]*LD[2])*(pcado[6][0]*LD[0] + pcado[6][1]*LD[1] + pcado[6][2]*LD[2])],
                       [(pcapo[7][0]*LD[0] + pcapo[7][1]*LD[1] + pcapo[7][2]*LD[2])*(pcado[7][0]*LD[0] + pcado[7][1]*LD[1] + pcado[7][2]*LD[2])],
                       [(pcapo[8][0]*LD[0] + pcapo[8][1]*LD[1] + pcapo[8][2]*LD[2])*(pcado[8][0]*LD[0] + pcado[8][1]*LD[1] + pcado[8][2]*LD[2])],
                       [(pcapo[9][0]*LD[0] + pcapo[9][1]*LD[1] + pcapo[9][2]*LD[2])*(pcado[9][0]*LD[0] + pcado[9][1]*LD[1] + pcado[9][2]*LD[2])],
                       [(pcapo[10][0]*LD[0] + pcapo[10][1]*LD[1] + pcapo[10][2]*LD[2])*(pcado[10][0]*LD[0] + pcado[10][1]*LD[1] + pcado[10][2]*LD[2])],
                       [(pcapo[11][0]*LD[0] + pcapo[11][1]*LD[1] + pcapo[11][2]*LD[2])*(pcado[10][0]*LD[0] + pcado[11][1]*LD[1] + pcado[11][2]*LD[2])]]
    
    pyramidal_ca_SF = np.absolute(pyramidal_ca_SF)
    
    basal_vbg1 = [[bpo[0][1]*surface[2] - bpo[0][2]*surface[1]],[bpo[1][1]*surface[2] - bpo[1][2]*surface[1]],[bpo[2][1]*surface[2] - bpo[2][2]*surface[1]]]
    basal_vbg2 = [[bpo[0][2]*surface[0] - bpo[0][0]*surface[2]],[bpo[1][1]*surface[0] - bpo[1][0]*surface[2]],[bpo[2][1]*surface[0] - bpo[2][0]*surface[2]]]
    basal_vbg3 = [[bpo[0][0]*surface[1] - bpo[0][1]*surface[0]],
                         [bpo[1][0]*surface[1] - bpo[1][1]*surface[0]],
                         [bpo[2][0]*surface[1] - bpo[2][1]*surface[0]]]
    
    
    prismatic_vbg1 = [[ppo[0][1]*surface[2] - ppo[0][2]*surface[1]],[ppo[1][1]*surface[2] - ppo[1][2]*surface[1]],[ppo[2][1]*surface[2] - ppo[2][2]*surface[1]]]
    prismatic_vbg2 = [[ppo[0][2]*surface[0] - ppo[0][0]*surface[2]],[ppo[1][1]*surface[0] - ppo[1][0]*surface[2]],[ppo[2][1]*surface[0] - ppo[2][0]*surface[2]]]
    prismatic_vbg3 = [[ppo[0][0]*surface[1] - ppo[0][1]*surface[0]],
                      [ppo[1][0]*surface[1] - ppo[1][1]*surface[0]],
                      [ppo[2][0]*surface[1] - ppo[2][1]*surface[0]]]
    
    
    pyramidal_vbg1 = [[pypo[0][1]*surface[2] - pypo[0][2]*surface[1]],
                      [pypo[1][1]*surface[2] - pypo[1][2]*surface[1]],
                      [pypo[2][1]*surface[2] - pypo[2][2]*surface[1]],
                      [pypo[3][1]*surface[2] - pypo[3][2]*surface[1]],
                      [pypo[4][1]*surface[2] - pypo[4][2]*surface[1]],
                      [pypo[5][1]*surface[2] - pypo[5][2]*surface[1]]]
    
    pyramidal_vbg2 = [[pypo[0][2]*surface[0] - pypo[0][0]*surface[2]],
                      [pypo[1][1]*surface[0] - pypo[1][0]*surface[2]],
                      [pypo[2][1]*surface[0] - pypo[2][0]*surface[2]],
                      [pypo[3][2]*surface[0] - pypo[3][0]*surface[2]],
                      [pypo[4][1]*surface[0] - pypo[4][0]*surface[2]],
                      [pypo[5][1]*surface[0] - pypo[5][0]*surface[2]]]
    
    pyramidal_vbg3 = [[pypo[0][0]*surface[1] - pypo[0][1]*surface[0]],
                      [pypo[1][0]*surface[1] - pypo[1][1]*surface[0]],
                      [pypo[2][0]*surface[1] - pypo[2][1]*surface[0]],
                      [pypo[3][0]*surface[1] - pypo[3][1]*surface[0]],
                      [pypo[4][0]*surface[1] - pypo[4][1]*surface[0]],
                      [pypo[5][0]*surface[1] - pypo[5][1]*surface[0]]]
    
    pyramidal_ca_vbg1 = [[pcapo[0][1]*surface[2] - pcapo[0][2]*surface[1]],
                         [pcapo[1][1]*surface[2] - pcapo[1][2]*surface[1]],
                         [pcapo[2][1]*surface[2] - pcapo[2][2]*surface[1]],
                         [pcapo[3][1]*surface[2] - pcapo[3][2]*surface[1]],
                         [pcapo[4][1]*surface[2] - pcapo[4][2]*surface[1]],
                         [pcapo[5][1]*surface[2] - pcapo[5][2]*surface[1]],
                         [pcapo[6][1]*surface[2] - pcapo[6][2]*surface[1]],
                         [pcapo[7][1]*surface[2] - pcapo[7][2]*surface[1]],
                         [pcapo[8][1]*surface[2] - pcapo[8][2]*surface[1]],
                         [pcapo[9][1]*surface[2] - pcapo[9][2]*surface[1]],
                         [pcapo[10][1]*surface[2] - pcapo[10][2]*surface[1]],
                         [pcapo[11][1]*surface[2] - pcapo[11][2]*surface[1]]]
    
    pyramidal_ca_vbg2 = [[pcapo[0][2]*surface[0] - pcapo[0][0]*surface[2]],
                         [pcapo[1][1]*surface[0] - pcapo[1][0]*surface[2]],
                         [pcapo[2][1]*surface[0] - pcapo[2][0]*surface[2]],
                         [pcapo[3][2]*surface[0] - pcapo[3][0]*surface[2]],
                         [pcapo[4][1]*surface[0] - pcapo[4][0]*surface[2]],
                         [pcapo[5][1]*surface[0] - pcapo[5][0]*surface[2]],
                         [pcapo[6][2]*surface[0] - pcapo[6][0]*surface[2]],
                         [pcapo[7][1]*surface[0] - pcapo[7][0]*surface[2]],
                         [pcapo[8][1]*surface[0] - pcapo[8][0]*surface[2]],
                         [pcapo[9][2]*surface[0] - pcapo[9][0]*surface[2]],
                         [pcapo[10][1]*surface[0] - pcapo[10][0]*surface[2]],
                         [pcapo[11][1]*surface[0] - pcapo[11][0]*surface[2]]]
    
    pyramidal_ca_vbg3 = [[pcapo[0][0]*surface[1] - pcapo[0][1]*surface[0]],
                         [pcapo[1][0]*surface[1] - pcapo[1][1]*surface[0]],
                         [pcapo[2][0]*surface[1] - pcapo[2][1]*surface[0]],
                         [pcapo[3][0]*surface[1] - pcapo[3][1]*surface[0]],
                         [pcapo[4][0]*surface[1] - pcapo[4][1]*surface[0]],
                         [pcapo[5][0]*surface[1] - pcapo[5][1]*surface[0]],
                         [pcapo[6][0]*surface[1] - pcapo[6][1]*surface[0]],
                         [pcapo[7][0]*surface[1] - pcapo[7][1]*surface[0]],
                         [pcapo[8][0]*surface[1] - pcapo[8][1]*surface[0]],
                         [pcapo[9][0]*surface[1] - pcapo[9][1]*surface[0]],
                         [pcapo[10][0]*surface[1] - pcapo[10][1]*surface[0]],
                         [pcapo[11][0]*surface[1] - pcapo[11][1]*surface[0]]]
    
    basal_traction1 = [[(180/math.pi)*math.acos((basal_vbg1[0][0]*LD[0]+basal_vbg2[0][0]*LD[1]+basal_vbg3[0][0]*LD[2])/math.sqrt(((basal_vbg1[0][0])**2)+((basal_vbg2[0][0])**2)+((basal_vbg3[0][0])**2)))],
                       [(180/math.pi)*math.acos((basal_vbg1[1][0]*LD[0]+basal_vbg2[1][0]*LD[1]+basal_vbg3[1][0]*LD[2])/math.sqrt(((basal_vbg1[1][0])**2)+((basal_vbg2[1][0])**2)+((basal_vbg3[1][0])**2)))],
                       [(180/math.pi)*math.acos((basal_vbg1[2][0]*LD[0]+basal_vbg2[2][0]*LD[1]+basal_vbg3[2][0]*LD[2])/math.sqrt(((basal_vbg1[2][0])**2)+((basal_vbg2[2][0])**2)+((basal_vbg3[2][0])**2)))]]
    
    prismatic_traction1 = [[(180/math.pi)*math.acos((prismatic_vbg1[0][0]*LD[0]+prismatic_vbg2[0][0]*LD[1]+prismatic_vbg3[0][0]*LD[2])/math.sqrt(((prismatic_vbg1[0][0])**2)+((prismatic_vbg2[0][0])**2)+((prismatic_vbg3[0][0])**2)))],
                           [(180/math.pi)*math.acos((prismatic_vbg1[1][0]*LD[0]+prismatic_vbg2[1][0]*LD[1]+prismatic_vbg3[1][0]*LD[2])/math.sqrt(((prismatic_vbg1[1][0])**2)+((prismatic_vbg2[1][0])**2)+((prismatic_vbg3[1][0])**2)))],
                           [(180/math.pi)*math.acos((prismatic_vbg1[2][0]*LD[0]+prismatic_vbg2[2][0]*LD[1]+prismatic_vbg3[2][0]*LD[2])/math.sqrt(((prismatic_vbg1[2][0])**2)+((prismatic_vbg2[2][0])**2)+((prismatic_vbg3[2][0])**2)))]]
    
    pyramidal_traction1 = [[(180/math.pi)*math.acos((pyramidal_vbg1[0][0]*LD[0]+pyramidal_vbg2[0][0]*LD[1]+pyramidal_vbg3[0][0]*LD[2])/math.sqrt(((pyramidal_vbg1[0][0])**2)+((pyramidal_vbg2[0][0])**2)+((pyramidal_vbg3[0][0])**2)))],
                           [(180/math.pi)*math.acos((pyramidal_vbg1[1][0]*LD[0]+pyramidal_vbg2[1][0]*LD[1]+pyramidal_vbg3[1][0]*LD[2])/math.sqrt(((pyramidal_vbg1[1][0])**2)+((pyramidal_vbg2[1][0])**2)+((pyramidal_vbg3[1][0])**2)))],
                           [(180/math.pi)*math.acos((pyramidal_vbg1[2][0]*LD[0]+pyramidal_vbg2[2][0]*LD[1]+pyramidal_vbg3[2][0]*LD[2])/math.sqrt(((pyramidal_vbg1[2][0])**2)+((pyramidal_vbg2[2][0])**2)+((pyramidal_vbg3[2][0])**2)))],
                           [(180/math.pi)*math.acos((pyramidal_vbg1[3][0]*LD[0]+pyramidal_vbg2[3][0]*LD[1]+pyramidal_vbg3[3][0]*LD[2])/math.sqrt(((pyramidal_vbg1[3][0])**2)+((pyramidal_vbg2[3][0])**2)+((pyramidal_vbg3[3][0])**2)))],
                           [(180/math.pi)*math.acos((pyramidal_vbg1[4][0]*LD[0]+pyramidal_vbg2[4][0]*LD[1]+pyramidal_vbg3[4][0]*LD[2])/math.sqrt(((pyramidal_vbg1[4][0])**2)+((pyramidal_vbg2[4][0])**2)+((pyramidal_vbg3[4][0])**2)))],
                           [(180/math.pi)*math.acos((pyramidal_vbg1[5][0]*LD[0]+pyramidal_vbg2[5][0]*LD[1]+pyramidal_vbg3[5][0]*LD[2])/math.sqrt(((pyramidal_vbg1[5][0])**2)+((pyramidal_vbg2[5][0])**2)+((pyramidal_vbg3[5][0])**2)))]]
    
    pyramidal_ca_traction1 = [[(180/math.pi)*math.acos((pyramidal_ca_vbg1[0][0]*LD[0]+pyramidal_ca_vbg2[0][0]*LD[1]+pyramidal_ca_vbg3[0][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[0][0])**2)+((pyramidal_ca_vbg2[0][0])**2)+((pyramidal_ca_vbg3[0][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[1][0]*LD[0]+pyramidal_ca_vbg2[1][0]*LD[1]+pyramidal_ca_vbg3[1][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[1][0])**2)+((pyramidal_ca_vbg2[1][0])**2)+((pyramidal_ca_vbg3[1][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[2][0]*LD[0]+pyramidal_ca_vbg2[2][0]*LD[1]+pyramidal_ca_vbg3[2][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[2][0])**2)+((pyramidal_ca_vbg2[2][0])**2)+((pyramidal_ca_vbg3[2][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[3][0]*LD[0]+pyramidal_ca_vbg2[3][0]*LD[1]+pyramidal_ca_vbg3[3][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[3][0])**2)+((pyramidal_ca_vbg2[3][0])**2)+((pyramidal_ca_vbg3[3][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[4][0]*LD[0]+pyramidal_ca_vbg2[4][0]*LD[1]+pyramidal_ca_vbg3[4][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[4][0])**2)+((pyramidal_ca_vbg2[4][0])**2)+((pyramidal_ca_vbg3[4][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[5][0]*LD[0]+pyramidal_ca_vbg2[5][0]*LD[1]+pyramidal_ca_vbg3[5][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[5][0])**2)+((pyramidal_ca_vbg2[5][0])**2)+((pyramidal_ca_vbg3[5][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[6][0]*LD[0]+pyramidal_ca_vbg2[6][0]*LD[1]+pyramidal_ca_vbg3[6][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[6][0])**2)+((pyramidal_ca_vbg2[6][0])**2)+((pyramidal_ca_vbg3[6][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[7][0]*LD[0]+pyramidal_ca_vbg2[7][0]*LD[1]+pyramidal_ca_vbg3[7][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[7][0])**2)+((pyramidal_ca_vbg2[7][0])**2)+((pyramidal_ca_vbg3[7][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[8][0]*LD[0]+pyramidal_ca_vbg2[8][0]*LD[1]+pyramidal_ca_vbg3[8][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[8][0])**2)+((pyramidal_ca_vbg2[8][0])**2)+((pyramidal_ca_vbg3[8][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[9][0]*LD[0]+pyramidal_ca_vbg2[9][0]*LD[1]+pyramidal_ca_vbg3[9][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[9][0])**2)+((pyramidal_ca_vbg2[9][0])**2)+((pyramidal_ca_vbg3[9][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[10][0]*LD[0]+pyramidal_ca_vbg2[10][0]*LD[1]+pyramidal_ca_vbg3[10][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[10][0])**2)+((pyramidal_ca_vbg2[10][0])**2)+((pyramidal_ca_vbg3[10][0])**2)))],
                              [(180/math.pi)*math.acos((pyramidal_ca_vbg1[11][0]*LD[0]+pyramidal_ca_vbg2[11][0]*LD[1]+pyramidal_ca_vbg3[11][0]*LD[2])/math.sqrt(((pyramidal_ca_vbg1[11][0])**2)+((pyramidal_ca_vbg2[11][0])**2)+((pyramidal_ca_vbg3[11][0])**2)))]]
    
    basal_traction2 = [[(180/math.pi)*math.asin(bpo[0][1]/math.sqrt(((bpo[0][0])**2)+((bpo[0][1])**2)))],
                       [(180/math.pi)*math.asin(bpo[1][1]/math.sqrt(((bpo[1][0])**2)+((bpo[1][1])**2)))],
                       [(180/math.pi)*math.asin(bpo[2][1]/math.sqrt(((bpo[2][0])**2)+((bpo[2][1])**2)))]]
    
    prismatic_traction2 = [[(180/math.pi)*math.asin(ppo[0][1]/math.sqrt(((ppo[0][0])**2)+((ppo[0][1])**2)))],
                           [(180/math.pi)*math.asin(ppo[1][1]/math.sqrt(((ppo[1][0])**2)+((ppo[1][1])**2)))],
                           [(180/math.pi)*math.asin(ppo[2][1]/math.sqrt(((ppo[2][0])**2)+((ppo[2][1])**2)))]]
    
    pyramidal_traction2 = [[(180/math.pi)*math.asin(math.sqrt(((pyramidal_vbg2[0][0]*LD[2]-pyramidal_vbg3[0][0]*LD[1])**2)+((pyramidal_vbg3[0][0]*LD[0]-pyramidal_vbg1[0][0]*LD[2])**2)+((pyramidal_vbg1[0][0]*LD[1]-pyramidal_vbg2[0][0]*LD[0])**2))/math.sqrt(pyramidal_vbg1[0][0]**2 + pyramidal_vbg2[0][0]**2 + pyramidal_vbg3[0][0]**2))],
                           [(180/math.pi)*math.asin(math.sqrt(((pyramidal_vbg2[1][0]*LD[2]-pyramidal_vbg3[1][0]*LD[1])**2)+((pyramidal_vbg3[1][0]*LD[0]-pyramidal_vbg1[1][0]*LD[2])**2)+((pyramidal_vbg1[1][0]*LD[1]-pyramidal_vbg2[1][0]*LD[0])**2))/math.sqrt(pyramidal_vbg1[1][0]**2 + pyramidal_vbg2[1][0]**2 + pyramidal_vbg3[1][0]**2))],
                           [(180/math.pi)*math.asin(math.sqrt(((pyramidal_vbg2[2][0]*LD[2]-pyramidal_vbg3[2][0]*LD[1])**2)+((pyramidal_vbg3[2][0]*LD[0]-pyramidal_vbg1[2][0]*LD[2])**2)+((pyramidal_vbg1[2][0]*LD[1]-pyramidal_vbg2[2][0]*LD[0])**2))/math.sqrt(pyramidal_vbg1[2][0]**2 + pyramidal_vbg2[2][0]**2 + pyramidal_vbg3[2][0]**2))],
                           [(180/math.pi)*math.asin(math.sqrt(((pyramidal_vbg2[3][0]*LD[2]-pyramidal_vbg3[3][0]*LD[1])**2)+((pyramidal_vbg3[3][0]*LD[0]-pyramidal_vbg1[3][0]*LD[2])**2)+((pyramidal_vbg1[3][0]*LD[1]-pyramidal_vbg2[3][0]*LD[0])**2))/math.sqrt(pyramidal_vbg1[3][0]**2 + pyramidal_vbg2[3][0]**2 + pyramidal_vbg3[3][0]**2))],
                           [(180/math.pi)*math.asin(math.sqrt(((pyramidal_vbg2[4][0]*LD[2]-pyramidal_vbg3[4][0]*LD[1])**2)+((pyramidal_vbg3[4][0]*LD[0]-pyramidal_vbg1[4][0]*LD[2])**2)+((pyramidal_vbg1[4][0]*LD[1]-pyramidal_vbg2[4][0]*LD[0])**2))/math.sqrt(pyramidal_vbg1[4][0]**2 + pyramidal_vbg2[4][0]**2 + pyramidal_vbg3[4][0]**2))],
                           [(180/math.pi)*math.asin(math.sqrt(((pyramidal_vbg2[5][0]*LD[2]-pyramidal_vbg3[5][0]*LD[1])**2)+((pyramidal_vbg3[5][0]*LD[0]-pyramidal_vbg1[5][0]*LD[2])**2)+((pyramidal_vbg1[5][0]*LD[1]-pyramidal_vbg2[5][0]*LD[0])**2))/math.sqrt(pyramidal_vbg1[5][0]**2 + pyramidal_vbg2[5][0]**2 + pyramidal_vbg3[5][0]**2))]]
    
    pyramidal_ca_traction2 = [[(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[0][0]*LD[2]-pyramidal_ca_vbg3[0][0]*LD[1])**2)+((pyramidal_ca_vbg3[0][0]*LD[0]-pyramidal_ca_vbg1[0][0]*LD[2])**2)+((pyramidal_ca_vbg1[0][0]*LD[1]-pyramidal_ca_vbg2[0][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[0][0]**2 + pyramidal_ca_vbg2[0][0]**2 + pyramidal_ca_vbg3[0][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[1][0]*LD[2]-pyramidal_ca_vbg3[1][0]*LD[1])**2)+((pyramidal_ca_vbg3[1][0]*LD[0]-pyramidal_ca_vbg1[1][0]*LD[2])**2)+((pyramidal_ca_vbg1[1][0]*LD[1]-pyramidal_ca_vbg2[1][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[1][0]**2 + pyramidal_ca_vbg2[1][0]**2 + pyramidal_ca_vbg3[1][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[2][0]*LD[2]-pyramidal_ca_vbg3[2][0]*LD[1])**2)+((pyramidal_ca_vbg3[2][0]*LD[0]-pyramidal_ca_vbg1[2][0]*LD[2])**2)+((pyramidal_ca_vbg1[2][0]*LD[1]-pyramidal_ca_vbg2[2][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[2][0]**2 + pyramidal_ca_vbg2[2][0]**2 + pyramidal_ca_vbg3[2][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[3][0]*LD[2]-pyramidal_ca_vbg3[3][0]*LD[1])**2)+((pyramidal_ca_vbg3[3][0]*LD[0]-pyramidal_ca_vbg1[3][0]*LD[2])**2)+((pyramidal_ca_vbg1[3][0]*LD[1]-pyramidal_ca_vbg2[3][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[3][0]**2 + pyramidal_ca_vbg2[3][0]**2 + pyramidal_ca_vbg3[3][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[4][0]*LD[2]-pyramidal_ca_vbg3[4][0]*LD[1])**2)+((pyramidal_ca_vbg3[4][0]*LD[0]-pyramidal_ca_vbg1[4][0]*LD[2])**2)+((pyramidal_ca_vbg1[4][0]*LD[1]-pyramidal_ca_vbg2[4][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[4][0]**2 + pyramidal_ca_vbg2[4][0]**2 + pyramidal_ca_vbg3[4][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[5][0]*LD[2]-pyramidal_ca_vbg3[5][0]*LD[1])**2)+((pyramidal_ca_vbg3[5][0]*LD[0]-pyramidal_ca_vbg1[5][0]*LD[2])**2)+((pyramidal_ca_vbg1[5][0]*LD[1]-pyramidal_ca_vbg2[5][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[5][0]**2 + pyramidal_ca_vbg2[5][0]**2 + pyramidal_ca_vbg3[5][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[6][0]*LD[2]-pyramidal_ca_vbg3[6][0]*LD[1])**2)+((pyramidal_ca_vbg3[6][0]*LD[0]-pyramidal_ca_vbg1[6][0]*LD[2])**2)+((pyramidal_ca_vbg1[6][0]*LD[1]-pyramidal_ca_vbg2[6][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[6][0]**2 + pyramidal_ca_vbg2[6][0]**2 + pyramidal_ca_vbg3[6][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[7][0]*LD[2]-pyramidal_ca_vbg3[7][0]*LD[1])**2)+((pyramidal_ca_vbg3[7][0]*LD[0]-pyramidal_ca_vbg1[7][0]*LD[2])**2)+((pyramidal_ca_vbg1[7][0]*LD[1]-pyramidal_ca_vbg2[7][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[7][0]**2 + pyramidal_ca_vbg2[7][0]**2 + pyramidal_ca_vbg3[7][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[8][0]*LD[2]-pyramidal_ca_vbg3[8][0]*LD[1])**2)+((pyramidal_ca_vbg3[8][0]*LD[0]-pyramidal_ca_vbg1[8][0]*LD[2])**2)+((pyramidal_ca_vbg1[8][0]*LD[1]-pyramidal_ca_vbg2[8][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[8][0]**2 + pyramidal_ca_vbg2[8][0]**2 + pyramidal_ca_vbg3[8][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[9][0]*LD[2]-pyramidal_ca_vbg3[9][0]*LD[1])**2)+((pyramidal_ca_vbg3[9][0]*LD[0]-pyramidal_ca_vbg1[9][0]*LD[2])**2)+((pyramidal_ca_vbg1[9][0]*LD[1]-pyramidal_ca_vbg2[9][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[9][0]**2 + pyramidal_ca_vbg2[9][0]**2 + pyramidal_ca_vbg3[9][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[10][0]*LD[2]-pyramidal_ca_vbg3[10][0]*LD[1])**2)+((pyramidal_ca_vbg3[10][0]*LD[0]-pyramidal_ca_vbg1[10][0]*LD[2])**2)+((pyramidal_ca_vbg1[10][0]*LD[1]-pyramidal_ca_vbg2[10][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[10][0]**2 + pyramidal_ca_vbg2[10][0]**2 + pyramidal_ca_vbg3[10][0]**2))],
                              [(180/math.pi)*math.asin(math.sqrt(((pyramidal_ca_vbg2[11][0]*LD[2]-pyramidal_ca_vbg3[11][0]*LD[1])**2)+((pyramidal_ca_vbg3[11][0]*LD[0]-pyramidal_ca_vbg1[11][0]*LD[2])**2)+((pyramidal_ca_vbg1[11][0]*LD[1]-pyramidal_ca_vbg2[11][0]*LD[0])**2))/math.sqrt(pyramidal_ca_vbg1[11][0]**2 + pyramidal_ca_vbg2[11][0]**2 + pyramidal_ca_vbg3[11][0]**2))]]
    
    basal_traction3 = [[(180/math.pi)*math.acos(bpo[0][0]/math.sqrt(((bpo[0][0])**2)+((bpo[0][1])**2)))],
                       [(180/math.pi)*math.acos(bpo[1][0]/math.sqrt(((bpo[1][0])**2)+((bpo[1][1])**2)))],
                       [(180/math.pi)*math.acos(bpo[2][0]/math.sqrt(((bpo[2][0])**2)+((bpo[2][1])**2)))]]
    
    prismatic_traction3 = [[(180/math.pi)*math.acos(ppo[0][0]/math.sqrt(((ppo[0][0])**2)+((ppo[0][1])**2)))],
                           [(180/math.pi)*math.acos(ppo[1][0]/math.sqrt(((ppo[1][0])**2)+((ppo[1][1])**2)))],
                           [(180/math.pi)*math.acos(ppo[2][0]/math.sqrt(((ppo[2][0])**2)+((ppo[2][1])**2)))]]
    
    basal_ap = [[-1*np.sign(basal_vbg1[0][0])*np.sign(basal_vbg2[0][0])*np.abs(basal_traction2[0][0])],
                [-1*np.sign(basal_vbg1[1][0])*np.sign(basal_vbg2[1][0])*np.abs(basal_traction2[1][0])],
                [-1*np.sign(basal_vbg1[2][0])*np.sign(basal_vbg2[2][0])*np.abs(basal_traction2[2][0])]]
    
    prismatic_ap = [[-1*np.sign(prismatic_vbg1[0][0])*np.sign(prismatic_vbg2[0][0])*np.abs(prismatic_traction2[0][0])],
                    [-1*np.sign(prismatic_vbg1[1][0])*np.sign(prismatic_vbg2[1][0])*np.abs(prismatic_traction2[1][0])],
                    [-1*np.sign(prismatic_vbg1[2][0])*np.sign(prismatic_vbg2[2][0])*np.abs(prismatic_traction2[2][0])]]
    
    pyramidal_ap = [[-1*np.sign(pyramidal_vbg1[0][0])*np.sign(pyramidal_vbg2[0][0])*np.abs(pyramidal_traction2[0][0])],
                    [-1*np.sign(pyramidal_vbg1[1][0])*np.sign(pyramidal_vbg2[1][0])*np.abs(pyramidal_traction2[1][0])],
                    [-1*np.sign(pyramidal_vbg1[2][0])*np.sign(pyramidal_vbg2[2][0])*np.abs(pyramidal_traction2[2][0])],
                    [-1*np.sign(pyramidal_vbg1[3][0])*np.sign(pyramidal_vbg2[3][0])*np.abs(pyramidal_traction2[3][0])],
                    [-1*np.sign(pyramidal_vbg1[4][0])*np.sign(pyramidal_vbg2[4][0])*np.abs(pyramidal_traction2[4][0])],
                    [-1*np.sign(pyramidal_vbg1[5][0])*np.sign(pyramidal_vbg2[5][0])*np.abs(pyramidal_traction2[5][0])]]
    
    pyramidal_ca_ap = [[-1*np.sign(pyramidal_ca_vbg1[0][0])*np.sign(pyramidal_ca_vbg2[0][0])*np.abs(pyramidal_ca_traction2[0][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[1][0])*np.sign(pyramidal_ca_vbg2[1][0])*np.abs(pyramidal_ca_traction2[1][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[2][0])*np.sign(pyramidal_ca_vbg2[2][0])*np.abs(pyramidal_ca_traction2[2][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[3][0])*np.sign(pyramidal_ca_vbg2[3][0])*np.abs(pyramidal_ca_traction2[3][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[4][0])*np.sign(pyramidal_ca_vbg2[4][0])*np.abs(pyramidal_ca_traction2[4][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[5][0])*np.sign(pyramidal_ca_vbg2[5][0])*np.abs(pyramidal_ca_traction2[5][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[6][0])*np.sign(pyramidal_ca_vbg2[6][0])*np.abs(pyramidal_ca_traction2[6][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[7][0])*np.sign(pyramidal_ca_vbg2[7][0])*np.abs(pyramidal_ca_traction2[7][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[8][0])*np.sign(pyramidal_ca_vbg2[8][0])*np.abs(pyramidal_ca_traction2[8][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[9][0])*np.sign(pyramidal_ca_vbg2[9][0])*np.abs(pyramidal_ca_traction2[9][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[10][0])*np.sign(pyramidal_ca_vbg2[10][0])*np.abs(pyramidal_ca_traction2[10][0])],
                       [-1*np.sign(pyramidal_ca_vbg1[11][0])*np.sign(pyramidal_ca_vbg2[11][0])*np.abs(pyramidal_ca_traction2[11][0])]]
    
    # print(basal_ap, 'basal')
    # print(prismatic_ap, 'prismatic')
    # print(pyramidal_ap, 'pyramidal')
    # print(pyramidal_ca_ap, 'pyramidal_ca')
    
    for i in range(len(basal_ap)):
        if basal_ap[i][0]<0:
            basal_ap[i][0] = basal_ap[i][0] + 180
        else:
            basal_ap[i][0] = basal_ap[i][0] + 0

    for i in range(len(prismatic_ap)):
        if prismatic_ap[i][0]<0:
            prismatic_ap[i][0] = prismatic_ap[i][0] + 180
        else:
            prismatic_ap[i][0] = prismatic_ap[i][0] + 0

    for i in range(len(pyramidal_ca_ap)):
        if pyramidal_ca_ap[i][0]<0:
            pyramidal_ca_ap[i][0] = pyramidal_ca_ap[i][0] + 180
        else:
            pyramidal_ca_ap[i][0] = pyramidal_ca_ap[i][0] + 0

    for i in range(len(pyramidal_ap)):
        if pyramidal_ap[i][0]<0:
            pyramidal_ap[i][0] = pyramidal_ap[i][0] + 180
        else:
            pyramidal_ap[i][0] = pyramidal_ap[i][0] + 0

# print(basal_ap, 'basal')
# print(prismatic_ap, 'prismatic')
# print(pyramidal_ap, 'pyramidal')
# print(pyramidal_ca_ap, 'pyramidal_ca')

    
    basal_min = [[basal_ap[0][0] - 2],
                 [basal_ap[1][0] - 2],
                 [basal_ap[2][0] - 2]] 
    basal_max = [[basal_ap[0][0] + 2],
                 [basal_ap[1][0] + 2],
                 [basal_ap[2][0] + 2]] 
    
    prismatic_min = [[prismatic_ap[0][0] - 2],
                     [prismatic_ap[1][0] - 2],
                     [prismatic_ap[2][0] - 2]] 
    prismatic_max = [[prismatic_ap[0][0] + 2],
                     [prismatic_ap[1][0] + 2],
                     [prismatic_ap[2][0] + 2]] 
    
    pyramidal_min = [[pyramidal_ap[0][0] - 2],
                     [pyramidal_ap[1][0] - 2],
                     [pyramidal_ap[2][0] - 2],
                     [pyramidal_ap[3][0] - 2],
                     [pyramidal_ap[4][0] - 2],
                     [pyramidal_ap[5][0] - 2]] 
    
    pyramidal_max = [[pyramidal_ap[0][0] + 2],
                     [pyramidal_ap[1][0] + 2],
                     [pyramidal_ap[2][0] + 2],
                     [pyramidal_ap[3][0] + 2],
                     [pyramidal_ap[4][0] + 2],
                     [pyramidal_ap[5][0] + 2]]
    
    pyramidal_ca_min = [[pyramidal_ca_ap[0][0] - 2],
                        [pyramidal_ca_ap[1][0] - 2],
                        [pyramidal_ca_ap[2][0] - 2],
                        [pyramidal_ca_ap[3][0] - 2],
                        [pyramidal_ca_ap[4][0] - 2],
                        [pyramidal_ca_ap[5][0] - 2],
                        [pyramidal_ca_ap[6][0] - 2],
                        [pyramidal_ca_ap[7][0] - 2],
                        [pyramidal_ca_ap[8][0] - 2],
                        [pyramidal_ca_ap[9][0] - 2],
                        [pyramidal_ca_ap[10][0] - 2],
                        [pyramidal_ca_ap[11][0] - 2]] 
    
    pyramidal_ca_max = [[pyramidal_ca_ap[0][0] + 2],
                        [pyramidal_ca_ap[1][0] + 2],
                        [pyramidal_ca_ap[2][0] + 2],
                        [pyramidal_ca_ap[3][0] + 2],
                        [pyramidal_ca_ap[4][0] + 2],
                        [pyramidal_ca_ap[5][0] + 2],
                        [pyramidal_ca_ap[6][0] + 2],
                        [pyramidal_ca_ap[7][0] + 2],
                        [pyramidal_ca_ap[8][0] + 2],
                        [pyramidal_ca_ap[9][0] + 2],
                        [pyramidal_ca_ap[10][0] + 2],
                        [pyramidal_ca_ap[11][0] + 2]]
    
# for i in range(len(basal_ap)):
#     if angle_real<= basal_max[i][0] and angle_real>=basal_min[i][0]:
#         print('yes, basal')
#     else:
#         print('no, basal')
        
# for i in range(len(prismatic_ap)):
#     if angle_real<= prismatic_max[i][0] and angle_real>=prismatic_min[i][0]:
#         print('yes, prismatic')
#     else:
#         print('no, prismatic')

# for i in range(len(pyramidal_ap)):
#     if angle_real<= pyramidal_max[i][0] and angle_real>=pyramidal_min[i][0]:
#         print('yes, pyramidal')
#     else:
#         print('no, pyramidal')
        
# for i in range(len(pyramidal_ca_ap)):
#     if angle_real<= pyramidal_ca_max[i][0] and angle_real>=pyramidal_ca_min[i][0]:
#         print('yes, pyramidal_ca')
#     else:
#         print('no, pyramidal_ca')

    arizona = []
    for i in range(len(basal_ap)):
        tuscon = abs(basal_ap[i][0] - angle_real)
        arizona.append(tuscon)
    
    california = []
    for i in range(len(prismatic_ap)):
        berkeley = abs(prismatic_ap[i][0] - angle_real)
        arizona.append(berkeley)
        
    washington = []
    for i in range(len(pyramidal_ap)):
        pullman = abs(pyramidal_ap[i][0] - angle_real)
        arizona.append(pullman)
    
    oregon = []        
    for i in range(len(pyramidal_ca_ap)):
        eugene = abs(pyramidal_ca_ap[i][0] - angle_real)
        arizona.append(eugene)

    michigan =[]
    for i in range(len(basal_SF)):
        annarbor = basal_SF[i][0]
        michigan.append(annarbor)

    for i in range(len(prismatic_SF)):
        lansing = prismatic_SF[i][0]
        michigan.append(lansing)

    for i in range(len(pyramidal_SF)):
        battlecreek = pyramidal_SF[i][0]
        michigan.append(battlecreek)

    for i in range(len(pyramidal_ca_SF)):
        grandrapids = pyramidal_ca_SF[i][0]
        michigan.append(grandrapids)
    
# oregon = [basal_SF, prismatic_SF,pyramidal_SF,pyramidal_ca_SF]    

    for i in range(len(arizona)):
        if arizona[i] == np.min(arizona):
            japan = michigan[i]
            washington.append(japan)
            olympia = np.max(washington)
            
    fridayharbor = np.where(olympia == michigan)
# fridayharbor[0][0]=20
    if fridayharbor[0][0] == 0:
        exeter = 'basal'
    elif fridayharbor[0][0]==1:
        exeter = 'basal'
    elif fridayharbor[0][0]==2:
        exeter = 'basal'
    elif fridayharbor[0][0]==3:
        exeter = 'prismatic'
    elif fridayharbor[0][0]==4:
        exeter =  'prismatic'
    elif fridayharbor[0][0]==5:
        exeter = 'prismatic'
    elif fridayharbor[0][0]==6:
        exeter = 'pyramidal'
    elif fridayharbor[0][0]==7:
        exeter = 'pyramidal'
    elif fridayharbor[0][0]==8:
        exeter = 'pyramidal'
    elif fridayharbor[0][0]==9:
        exeter = 'pyramidal'
    elif fridayharbor[0][0]==10:
        exeter = 'pyramidal'
    elif fridayharbor[0][0]==11:
        exeter = 'pyramidal'
    elif fridayharbor[0][0]==12:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==13:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==14:
        exeter = 'pyramidal_ca'  
    elif fridayharbor[0][0]==15:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==16:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==17:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==18:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==19:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==20:
        exeter = 'pyramidal_ca' 
    elif fridayharbor[0][0]==21:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==22:
        exeter = 'pyramidal_ca'
    elif fridayharbor[0][0]==23:
        exeter = 'pyramidal_ca' 

    print(exeter)
    print(olympia)
    
    thisslipband = kansas[0][j]
    blnintensity = kansas[2][j]
    # thisslipband = j
    final = [[thisslipband],[phasetype],[angle_real],[blnintensity],[exeter],[olympia]]
    final = np.asarray(final).T
    final = pd.DataFrame(final)
    maine.append(final)

rhodes = pd.concat(maine)
#Slip Band, Angle of Slip Band, Slip Plane, Schmid Factor
rhodes.to_csv(directory + '/Deliverables_Ti02_New.csv')
# birch = oak[1][:]
# birch = pd.Series(birch)
# slip_band = pd.Series(slip_band)
# forest = [[slip_band],[birch]]
toc = time.time()
print(toc-tic,'seconds elapsed')