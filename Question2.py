import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def calculate_same_matrix(d):
        
        
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

    return l_matrix              

    
def calculate_relation_matrix(d,D):
    
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
                
                l_matrix[i][j]=d**2/(4*pi*e0*np.sqrt((x_s[i]-x_s[j])**2+(y_s[i]-y_s[j])**2+(D)**2))

    return l_matrix              

def calculate_Q_for_D(d,D):
    
    Qs=[]
    mat_AA,mat_BB=calculate_same_matrix(d),calculate_same_matrix(d)
    
    for ds in D:
        

        mat_AB=calculate_relation_matrix(d,ds)
        mat_BA=np.transpose(mat_AB)
        mat_AAAB=np.concatenate((mat_AA,mat_AB),axis=1)
        mat_BABB=np.concatenate((mat_BA,mat_BB),axis=1)
        big_mat=np.concatenate((mat_AAAB,mat_BABB),axis=0)
        inv_big_mat=np.linalg.inv(big_mat)

        V_potential_positive=np.transpose(np.ones([1,mat_AA.shape[0]]))*(0.5)
        V_potential_negative=np.transpose(np.ones([1,mat_AA.shape[0]]))*(-0.5)
        V_tot=np.concatenate((V_potential_positive,V_potential_negative),axis=0)

        qi=(d**2)*np.matmul(inv_big_mat,V_tot)
        Qs.append(np.sum(qi[0:int(qi.shape[0]/2)]))

    return Qs
        


d=0.025
e0=8.85*(10**-12)

#print(calculate_Q_for_D(d,[0.5]))
#print(calculate_Q_for_D(d,[0.2]))
#print(calculate_Q_for_D(d,[0.1]))
#print(calculate_Q_for_D(d,[0.0001]))

Ds=np.linspace(d/3,1,10)
graphingdf=pd.DataFrame({'D':list(Ds)})
graphingdf['C_theory']=(e0*np.pi)/graphingdf['D']
graphingdf['C_num']=calculate_Q_for_D(d,graphingdf['D'])
graphingdf['Error']=(graphingdf['C_num']/graphingdf['C_theory'])*100

fig=plt.figure(figsize=(10,10),dpi=200)
plt.xlabel('D[m]')
plt.ylabel('C[Farad]')
plt.plot(graphingdf['D'],graphingdf['C_theory'],label="Theory Capacitance")
plt.plot(graphingdf['D'],graphingdf['C_num'],label="Numeric Capacitance")
plt.legend(loc="upper right")
plt.show()