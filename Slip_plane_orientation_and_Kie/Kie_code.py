from __future__ import division
import numpy as np
import io
import scipy
from array import *
import os
import math


'''
The code is for PhD research purposes

Author: Praveenkumar Hiremath
Email: praveenkumar.hiremath@mek.lth.se (Email at the University)
       praveenkumar.hiremath2911@gmail.com (Private email)
'''

#crack_prop=np.array([0,-1,0]) #x-axis
perp_crack_plane=np.array([1,1,2]) #y-axis  Enter the vectors such that [1*a,1*a,2*c] i.e indices are multiplied by a,a and c.
u=perp_crack_plane[0]
v=perp_crack_plane[1]
w=perp_crack_plane[2]

crack_front=np.array([1,1,-1]) #z-axis    Enter the vectors such that indices are multiplied by a,a and c.
m=crack_front[0]
n=crack_front[1]
o=crack_front[2]

#Find crack propagation plane
crack_prop=np.cross(perp_crack_plane,crack_front)

Slip_planes=np.array([[1,-1,0],[-1,1,0]])  # Enter the vectors such that indices are multiplied by a,a and c.
Slip_dirs=np.array([[1,1,1],[1,1,1]])      # Enter the vectors such that indices are multiplied by a,a and c.


a=Slip_planes[:,0]  #x-component of slip plane vectors
b=Slip_planes[:,1]  #y-component of slip plane vectors
c=Slip_planes[:,2]  #z-component of slip plane vectors


#burger vectors components
b1=Slip_dirs[:,0]   #x-component of burger vectors
b2=Slip_dirs[:,1]   #y-component of burger vectors
b3=Slip_dirs[:,2]   #z-component of burger vectors


theta_angles=np.zeros(len(Slip_dirs))
theta_angles_radian=np.zeros(len(Slip_dirs))
phi_angles=np.zeros(len(Slip_dirs))

result_text=[]*len(Slip_dirs)  

for i in range(0,len(Slip_dirs),1):
  dot_prod=0.0
  mag_prop=0.0
  mag_slip_plane=0.0

  dot_prod=np.dot(crack_prop,Slip_planes[i,:])                  # Dot product of the two vectors
  mag_prop=np.linalg.norm(crack_prop)                           # Magnitude of the vector
  mag_slip_plane=np.linalg.norm(Slip_planes[i,:])               # Magnitude of the vector

# The following block uses equation (2) in report for calculation of Theta angle.
  cos_theta1=(dot_prod)/(mag_prop*mag_slip_plane)               # Finding angle theta 
  cos_theta1=np.round(cos_theta1,3)
  theta1=(180*np.arccos(cos_theta1))/math.pi                    # Converting to degrees
  theta=90-theta1
  theta_angles[i]=theta
 

#converting into radians
  theta1=((math.pi)/180)*theta1
  theta=((math.pi)/180)*theta
  theta_angles_radian[i]=theta


### Now for phi angle using equations (3), (4), (5) and (6)
  A=np.zeros(3)
  A=np.array([[a[i],b[i],c[i]],[u,v,w],[m,n,o]])    # 3X3 matrix from equation (6) in report 


  C=np.zeros(3)
  C=np.array([0,((u*u)+(v*v)+(w*w)*math.sin(theta)),0])  # 3X1 matrix (vector) on the RHS of '=' in equation (6) of report


  if (theta!=0.0):
    X=np.zeros(3)
    X=np.linalg.solve(A, C)       # Solving for [p,q,r] in equation (6) in report


# Following block uses equation (2) in report
    cos_phi=((X[0]*b1[i])+(X[1]*b2[i])+(X[2]*b3[i]))/((np.sqrt((X[0]*X[0])+(X[1]*X[1])+(X[2]*X[2])))*(np.sqrt((b1[i]*b1[i])+(b1[i]*b1[i])+(b1[i]*b1[i]))))
    cos_phi=np.round(cos_phi,3)
    phi=(180*np.arccos(cos_phi))/math.pi
    phi_angles[i]=phi

  else:
    phi_angles[i]=0.0
 
#  print ('Slip plane: ('+str(Slip_planes[i,0])+str(Slip_planes[i,1])+str(Slip_planes[i,2])+')\t'+'Slip direction(burger vector): ['+ str(Slip_dirs[i,0])+str(Slip_dirs[i,1])+str(Slip_dirs[i,2])+']'+'\t angle theta: '+str(theta_angles[i])+'\t angle phi: '+str(phi_angles[i])+'\n')


################################### KIE calculation using equation (1) in report ###########################################################


K_emit_110=[]*int(len(Slip_dirs)*1)
K_emit_110_test=[]*int(len(Slip_dirs)*1)

print ('theta_angles=',theta_angles)
print ('phi_angles=',phi_angles)

j=0

theta_new=theta_angles
phi_new=phi_angles
G_G=0.0
for j in range(0,int(len(Slip_dirs)*1),1): 
  	angle_theta=theta_angles[j] 
	angle_phi=phi_angles[j]
	if ((angle_theta>90.0)):
	  angle_theta=180-angle_theta
        if (angle_phi>90):
	  angle_phi=180-angle_phi      

        print (angle_theta, angle_phi)
	nu_yx=0.278     # This will be zero for our cases
	gamma_us=1.5775		
	theta=((math.pi)/180)*angle_theta
	phi=((math.pi)/180)*angle_phi
	denom=((1+math.cos(theta))*(pow(math.sin(theta),2)))*gamma_us
	#if (denom==0.0):
    #    G_G=0.0
	if (denom!=0.0):
	    G_G=8*((1+(1-nu_yx)*pow(math.tan(phi),2))/((1+math.cos(theta))*(pow(math.sin(theta),2))))*gamma_us


	C=np.array([[532.55*1e9,204.95*1e9,204.95*1e9,0.0,0.0,0.0],[204.95*1e9,532.55*1e9,204.95*1e9,0.0,0.0,0.0],[204.95*1e9,204.95*1e9,532.55*1e9,0.0,0.0,0.0],[0.0,0.0,0.0,163.13*1e9,0.0,0.0],[0.0,0.0,0.0,0.0,163.13*1e9,0.0],[0.0,0.0,0.0,0.0,0.0,163.13*1e9]])   # Elastic constant matrix in appropriate orientation


	S=np.linalg.inv(C)     # Compliance constant matrix

	s11=S[0,0]
	s12=S[0,1]
	s13=S[0,2]
	s22=S[1,1]
	s33=S[2,2]
	s23=S[1,2]
	s26=S[1,5]
	s66=S[5,5]


##bij are in GPa^-1
	b_11=((s11*s33)-(s13*s13))/s33
	b_22=((s22*s33)-(s23*s23))/s33
	b_12=((s12*s33)-(s13*s23))/s33
	b_66=((s66*s33)-(s26*s26))/s33

#### B from equation (1) in report will be in GPa^-1
	temp_1=(b_11*b_22*0.5)
	temp_2=np.sqrt(b_22/b_11)
	temp_3=((2*b_12)+b_66)
	temp_4=(2*b_11)

	B_1=np.sqrt((temp_1)*(temp_2+(temp_3/temp_4)))
	K_G_N=np.sqrt(G_G/B_1)     # equation (1) in report
	K_G_N_MPa_root_m=K_G_N/1e6     



	print('Slip plane: ('+str(Slip_planes[j,0])+str(Slip_planes[j,1])+str(Slip_planes[j,2])+')\t Slip direction (burger vector): '+'['+ str(Slip_dirs[j,0])+str(Slip_dirs[j,1])+str(Slip_dirs[j,2])+']'+'\t angle theta: '+str(theta_angles[j])+'\t angle phi: '+str(phi_angles[j])+'\t K_IE: '+str(K_G_N_MPa_root_m)+'\n')




