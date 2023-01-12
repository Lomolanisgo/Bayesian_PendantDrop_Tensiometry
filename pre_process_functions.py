import os.path
import glob
import cv2
from PIL import Image
import numpy as np
from skimage import io

def __init__():
   return

def ResizeImage(filein,width,height):
    fileout=filein.resize((width,height),Image.Resampling.LANCZOS)
    return fileout

def ResizeArray(arrayin,width,height):
    arrayout=arrayin.resize((width,height))
    return arrayout

def img2bw(img):
    img_gray=img.convert('L') # rgb 2 gray
    array_img=np.array(img_gray) # img 2 array
    _,img_BW=cv2.threshold(array_img,0,255,cv2.THRESH_OTSU) # gray 2 bw
    return img_BW

def cost(img_ori,img_gen,size=100,output=0):
    '''
    Input: 
    img_ori = Image.open(path)
    img_ori: the picture need be estimate the tension
    img_gen: the picture that function generated
    size: the size of comparation picture
    
    output = 0 --> lost
    output = 1 --> accurancy
    output = others --> accurancy, C (image matrix)

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
    lost=f/max(drop_gen,drop_ori)
    accurancy=1-lost

    if output == 0:
        return lost
    elif output == 1:
        return accurancy
    else:
        return accurancy,C

def crop_margin(img,path='./images_notupload/cut_edge.png'):
        img2=img.sum(axis=2)
        (row,col)=img2.shape
        row_top=0
        raw_down=0
        col_top=0
        col_down=0
        for r in range(0,row):
                if img2.sum(axis=1)[r]<700*col:
                        row_top=r
                        break
 
        for r in range(row-1,0,-1):
                if img2.sum(axis=1)[r]<700*col:
                        raw_down=r
                        break
 
        for c in range(0,col):
                if img2.sum(axis=0)[c]<700*row:
                        col_top=c
                        break
 
        for c in range(col-1,0,-1):
                if img2.sum(axis=0)[c]<700*row:
                        col_down=c
                        break
 
        new_img=img[row_top:raw_down+1,col_top:col_down+1,0:3]

        io.imsave(path,new_img)

        return new_img,path