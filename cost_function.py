import os.path
import glob
import cv2
from PIL import Image
import numpy as np

def __init__():
   return

def ResizeImage(filein,width,height):
    fileout=filein.resize((width,height),Image.Resampling.LANCZOS)
    return fileout


def img2bw(img,size):
    img_rs=ResizeImage(img,size,size) #resize 
    img_gray=img_rs.convert('L') # rgb 2 gray
    array_img=np.array(img_gray) # img 2 array
    _,img_BW=cv2.threshold(array_img,0,255,cv2.THRESH_OTSU) # gray 2 bw
    return img_BW

def cost(img_ori,img_gen,size=100):
    '''
    Input: 
    img_ori = Image.open(path)
    img_ori: the picture need be estimate the tension
    img_gen: the picture that function generated
    size: the size of comparation picture
    
    Output
    accurancy:matched pixels/ max drop pixels
    C: the canparation picture in matrix
    '''
    t=0; f=0; drop_ori=0; drop_gen=0 # pixel that same(T) different(F) and original drop (drop)
    A=img2bw(img_ori,size)
    B=img2bw(img_gen,size)
    C=np.zeros((size,size))
    for i in range (size):
        #print(i)
        for j in range(size):
            #print(j)
            if A[i,j]==0:
                drop_ori=drop_ori+1
            if B[i,j]==0:
                drop_gen=drop_gen+1

            if A[i,j]==B[i,j]:
                C[i,j]=255
                t=t+1
            else:
                C[i,j]=0
                f=f+1
    accurancy=1-f/max(drop_gen,drop_ori)
    return accurancy,C
