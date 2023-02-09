from fun_genSingleDrop import *
import numpy as np
import warnings
import os.path
warnings.filterwarnings('ignore')


#sigma=100 # surface tension [mN/m]
#volume0=32  # prescribed volume in mm^3

path='../Images_Generated'
path = os.path.abspath(path)

for sigma in np.arange (73,74,1):
    for volume0 in np.arange (20,21,1):
        genSingleDrop(path,sigma,volume0)

print('Program Finished')