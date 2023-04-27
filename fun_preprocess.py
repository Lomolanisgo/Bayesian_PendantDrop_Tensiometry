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

def crop_Margin(syn_array):#(DONE)
        
        (row,col)=syn_array.shape
        img=syn_array

        for i in range(row-1):
                if sum(img[i])!=255*col:
                        cut1=i
                        break               

        for i in range(row-1,0,-1):
                if sum(img[i])!=255*col:
                        cut2=i
                        break 

        _,img1,_=np.split(img,[cut1,cut2])  

        for i in range(col-1):
                if sum(img[:,i])!=255*row:
                        cut3=i
                        break

        for i in range(col-1,0,-1):
                if sum(img[:,i])!=255*row:
                        cut4=i
                        break

        _,cm_arr,_=np.split(img1,[cut3,cut4],axis=1)

        return cm_arr
