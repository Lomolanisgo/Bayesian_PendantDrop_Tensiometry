from genSingleDrop import *
import numpy as np
#sigma=100 # surface tension [mN/m]
#volume0=32  # prescribed volume in mm^3
rneedle=1 # [mm]

for sigma in np.arange (10,100,5):
    sigma=sigma
    for volume0 in np.arange (5,50,5):
        volume0=volume0
        genSingleDrop(sigma,volume0,rneedle)

print('program finished')