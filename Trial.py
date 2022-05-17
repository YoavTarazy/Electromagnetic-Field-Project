import numpy as np

a=np.transpose(np.array([[1,1,1,1]]))
b=np.transpose(np.array([[2,2,2,2]]))
print(a.shape[0]/2)
print(np.concatenate((a,b),axis=0))