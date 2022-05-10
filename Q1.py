import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#relavent x and y's and edge of square

x_s=[]
y_s=[]
d=0.025
pi=np.pi
e0=8.85*(10**-12)

#Define Circle
R=1
c_point=[0,0]
circ_func=lambda x: x[0]**2+x[1]**2

steps=int(np.ceil(2*R/d))

x_pos=np.linspace(-1,1,steps)
y_pos=np.linspace(-1,1,steps)


for i in range(steps):
    for j in range(steps):
        if circ_func([x_pos[i],y_pos[j]])<R**2:
            x_s.append(x_pos[i])
            y_s.append(y_pos[j])
        

l_matrix=np.zeros([len(x_s),len(x_s)])

for i in range(len(x_s)):
    for j in range(len(x_s)):
        if i != j:
            l_matrix[i][j]=d**2/(4*pi*e0*np.sqrt((x_s[i]-x_s[j])**2+(y_s[i]-y_s[j])**2))
        else:
            l_matrix[i][j]=(d*0.8814)/(pi*e0)

print(l_matrix)            


#inv_l_matrix=np.linalg.inv(l_matrix)


##Building vectors for system
V_potential=np.transpose(np.ones([1,len(x_s)]))

sigma=(d**2)*np.linalg.solve(l_matrix,V_potential)
#sigma=(d**2)*np.matmul(inv_l_matrix,V_potential)
sigma_sum=np.sum(sigma)

print(sigma_sum)


### building pandas for graph of Q(d)
