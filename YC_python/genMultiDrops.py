from genSingleDrop import *

#sigma=100 # surface tension [mN/m]
#volume0=32  # prescribed volume in mm^3
rneedle=1 # [mm]

for sigma in range (20,22):
    sigma=sigma
    for volume0 in range (7,15):
        volume0=volume0
        genSingleDrop(sigma,volume0,rneedle)