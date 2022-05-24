import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def calculation(D):
    
    Qs=[]
    
    for d in D:
        
        #relavent x and y's and edge of square
        
        e0=8.85*(10**-12)
        x_s=[]
        y_s=[]
        pi=np.pi
        

        #Define Circle
        R=1
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

                  

        inv_l_matrix=np.linalg.inv(l_matrix)


        ##Building vectors for system
        V_potential=np.transpose(np.ones([1,len(x_s)]))
        ############qi=(d**2)*np.linalg.solve(l_matrix,V_potential)
        qi=(d**2)*np.matmul(inv_l_matrix,V_potential)
        Qs.append(np.sum(qi))
    
    return Qs





e0=8.85*(10**-12)
graphingdf=pd.DataFrame({'d':[0.25,0.15,0.12,0.1,0.075,0.05,0.025,0.02]})
#graphingdf['Q']=[6.0435*(10**-11),6.7425*(10**-11),6.8408*(10**-11),6.8965*(10**-11),6.9397*(10**-11),7.0199*(10**-11),7.0491*(10**-11),7.06004*(10**-11)]
graphingdf['Q']=calculation(graphingdf['d'])
print(graphingdf)
graphingdf['8e0']=8*e0
fig=plt.figure(figsize=(10,10),dpi=200)
plt.xlabel('d[m]')
plt.ylabel('Q Charge [coloumb]')

plt.plot(graphingdf['d'],graphingdf['Q'])
plt.plot(graphingdf['d'],graphingdf['8e0'])
plt.show()