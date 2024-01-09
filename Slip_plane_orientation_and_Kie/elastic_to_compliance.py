from __future__ import division
import numpy as np
import io
import scipy
from array import *
import os
import math


#W-Re structure 9
C=np.array([[535.7*1e9,246.8*1e9,206.2*1e9,0.0,0.0,0.0],[246.8*1e9,535.7*1e9,206.2*1e9,0.0,0.0,0.0],[206.2*1e9,206.2*1e9,565.4*1e9,0.0,0.0,0.0],[0.0,0.0,0.0,104.9*1e9,0.0,0.0],[0.0,0.0,0.0,0.0,104.9*1e9,0.0],[0.0,0.0,0.0,0.0,0.0,121.7*1e9]])

#gamma-TiAL 
#C=np.array([[195*1e9,107*1e9,113*1e9,0.0,0.0,0.0],[107*1e9,195*1e9,113*1e9,0.0,0.0,0.0],[113*1e9,113*1e9,213*1e9,0.0,0.0,0.0],[0.0,0.0,0.0,92*1e9,0.0,0.0],[0.0,0.0,0.0,0.0,92*1e9,0.0],[0.0,0.0,0.0,0.0,0.0,84*1e9]])

new_cij_1=np.append(C[0,:],C[1,:])
new_cij_2=np.append(new_cij_1,C[2,:])
new_cij_3=np.append(new_cij_2,C[3,:])
new_cij_4=np.append(new_cij_3,C[4,:])
new_cij=np.append(new_cij_4,C[5,:])

np.savetxt ('elastic.dat',new_cij)

S=np.linalg.inv(C)  # Complaince matrix = inverse of elastic constants' matrix

Sij_1=np.append(S[0,:],S[1,:])
Sij_2=np.append(Sij_1,S[2,:])
Sij_3=np.append(Sij_2,S[3,:])
Sij_4=np.append(Sij_3,S[4,:])
Sij=np.append(Sij_4,S[5,:])

np.savetxt ('compliance.dat',Sij)
#print (S)
