import numpy as np

a=np.transpose(np.array([[1,1,1,1]]))
b=np.transpose(np.array([[2,2,2,2]]))

c=np.concatenate((a,b),axis=0)

v=np.zeros((4))
print(v)
print(np.transpose(v))
