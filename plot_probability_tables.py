import numpy as np
import matplotlib.pyplot as plt

UT = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/tables/mortality_UNIVERSAL.csv',skip_header=1,usecols=range(1,24),dtype='float',delimiter=',')

ptiles = [0,0.01]
for p in range(5,100,5):
    ptiles.append(p/100)
ptiles.append(0.99)
ptiles.append(1)

            
