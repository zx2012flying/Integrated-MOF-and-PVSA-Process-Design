# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 17:20:48 2021

@author: zhangx
"""
# %% import package

import numpy as np
import pandas as pd

# %% get movies

Movies = open('Movie_sxb_v1-6c_Fe_2_3_1-FUDQIF_1_2.2.2_300.000000_101325.000000_allcomponents.pdb')
file = Movies.readlines()

AX = []
AY = []
AZ = []

for i in file:
    if 'ATOM' in i:
        X = i[31:38]
        AX.append(X)
        Y = i[39:46]
        AY.append(Y)
        Z = i[47:54]
        AZ.append(Z)

AX = np.asarray(AX, dtype = 'float')
AY = np.asarray(AY, dtype = 'float')
AZ = np.asarray(AZ, dtype = 'float')

del X, Y, Z, Movies, file, i

# # %% convert to X-Y projection

# Right_max = 27.832
# Up_max = 26.928

# Result = np.vstack([AX, AY])
# Result = np.transpose(Result)

# # %% convert to X-Z projection

# Right_max = 27.832
# Up_max = 30.314

# Result = np.vstack([AX, AZ])
# Result = np.transpose(Result)

# %% convert to Y-Z projection

Right_max = 26.928
Up_max = 30.314

Result = np.vstack([AY, AZ])
Result = np.transpose(Result)

# %% calculate intensity

Box = 100
Unit_cell = int(0.5 * Box)

DeltaR = np.linspace(0, Right_max, Box+1)
DeltaU = np.linspace(0, Up_max, Box+1)

Intensity = np.zeros([Box, Box])
Quarter_Intensity = np.zeros([Unit_cell,  Unit_cell])

Count = 0
for i in range(0, int(len(Result) / 3)):
    if 0 < Result[3*i, 0] < Right_max and 0 < Result[3*i, 1] < Up_max:
        Count += 1
        x = len(DeltaR[DeltaR <= Result[3*i, 0]]) 
        y = len(DeltaU[DeltaU <= Result[3*i, 1]])
        Intensity[Box - y, x-1] += 1 
    else:
        if 0 < Result[3*i+1, 0] < Right_max and 0 < Result[3*i+1, 1] < Up_max:
            Count += 1
            x = len(DeltaR[DeltaR <= Result[3*i+1, 0]]) 
            y = len(DeltaU[DeltaU <= Result[3*i+1, 1]])
            Intensity[Box - y, x-1] += 1   
        else:
            if 0 < Result[3*i+2, 0] < Right_max and 0 < Result[3*i+2, 1] < Up_max:
                Count += 1
                x = len(DeltaR[DeltaR <= Result[3*i+2, 0]]) 
                y = len(DeltaU[DeltaU <= Result[3*i+2, 1]])
                Intensity[Box - y, x-1] += 1      


Quarter_Intensity = Intensity[0:Unit_cell, 0:Unit_cell]
Quarter_Intensity += Intensity[0:Unit_cell, Unit_cell:Box]
Quarter_Intensity += Intensity[Unit_cell:Box, 0:Unit_cell]
Quarter_Intensity += Intensity[Unit_cell:Box, Unit_cell:Box]
Quarter_Intensity = Quarter_Intensity / 4

# %% write to excel

Quarter_Intensity = pd.DataFrame(Quarter_Intensity)
Quarter_Intensity.to_excel('Movies_intensity_YZ.xlsx', header = False, index = False)

# Full = np.vstack([AX, AY, AZ])
# Full = np.transpose(Full)
# Full = pd.DataFrame(Full)
# Full.to_excel('Movies_coordinate.xlsx', header = False, index = False)


# # %% convert to X-Z projection

# Result = np.vstack([AX, AZ])
# Result = np.transpose(Result)

# Box = 100
# Unit_cell = int(0.5 * Box)

# Y_max = 26.928
# Z_max = 30.314
# DeltaY = np.linspace(0, Y_max, Box+1)
# DeltaZ = np.linspace(0, Z_max, Box+1)

# Intensity = np.zeros([Box, Box])
# Quarter_Intensity = np.zeros([Unit_cell,  Unit_cell])

# Count = 0
# for i in range(0, int(len(Result) / 3)):
#     if 0 < Result[3*i, 0] < Y_max and 0 < Result[3*i, 1] < Z_max:
#         Count += 1
#         y = len(DeltaY[DeltaY <= Result[3*i, 0]]) 
#         z = len(DeltaZ[DeltaZ <= Result[3*i, 1]])
#         Intensity[y-1, Box - z] += 1        
#     else:
#         if 0 < Result[3*i+1, 0] < Y_max and 0 < Result[3*i+1, 1] < Z_max:
#             Count += 1
#             y = len(DeltaY[DeltaY <= Result[3*i+1, 0]]) 
#             z = len(DeltaZ[DeltaZ <= Result[3*i+1, 1]])
#             Intensity[y-1, Box - z] += 1      
#         else:
#             if 0 < Result[3*i+2, 0] < Y_max and 0 < Result[3*i+2, 1] < Z_max:
#                 Count += 1
#                 y = len(DeltaY[DeltaY <= Result[3*i+2, 0]]) 
#                 z = len(DeltaZ[DeltaZ <= Result[3*i+2, 1]])
#                 Intensity[y-1, Box - z] += 1      

# Quarter_Intensity = Intensity[0:Unit_cell, 0:Unit_cell]
# Quarter_Intensity += Intensity[0:Unit_cell, Unit_cell:Box]
# Quarter_Intensity += Intensity[Unit_cell:Box, 0:Unit_cell]
# Quarter_Intensity += Intensity[Unit_cell:Box, Unit_cell:Box]
# Quarter_Intensity = Quarter_Intensity / 4

# # %% write to excel

# Quarter_Intensity = pd.DataFrame(Quarter_Intensity)
# Quarter_Intensity.to_excel('Movies_intensity_YZ.xlsx', header = False, index = False)

