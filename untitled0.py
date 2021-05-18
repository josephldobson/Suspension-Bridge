import numpy as np

A = np.linspace(1,10,20)
B = np.linspace(10,50,20)
C = np.linspace(50,100,20)

D = np.concatenate([A,B,C])
print(D)