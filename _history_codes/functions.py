import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from PIL import Image
from fun_genSingleDrop import *
from matplotlib.backends.backend_agg import FigureCanvasAgg

def removeStandNeedleCentering(exp_path):
    exp_img=Image.open(exp_path).convert('L')
    exp_arr=np.asarray(exp_img)
    h_exp,w_exp=exp_arr.shape

    # initialize the left and right edge array
    left_edge=np.zeros(w_exp)
    right_edge=np.zeros(w_exp)

    # find the left and right position of array
    for i in range(h_exp):
        for j in range (w_exp):
                if exp_arr[i,j-1]>=125 and exp_arr[i,j]<125:
                    left_edge[i]=j
                if exp_arr[i,j-1]<=125 and exp_arr[i,j]>125:
                    right_edge[i]=j

    # calculate the needle array
    needle=np.abs(left_edge-right_edge)

    # find the intersection of droplet and needle
    end=0
    for j in range(needle.shape[0]-1):
        j=j+1
        if end == 0:
            if abs((needle[j]-needle[0])/needle[0])>0.05:
                stopline=j
                end=end+1

    # centering by add edge
    addpad=right_edge[stopline]-needle[stopline]/2-w_exp/2
     # centering
    if addpad>0:
        arr_won=np.hstack( ( exp_arr,np.ones((h_exp,abs(int(addpad*4))))*255 ) )
    elif addpad<0:
        arr_won=np.hstack( ( np.ones((h_exp,abs(int(addpad*4))))*255,exp_arr ) )
    else:
        arr_won=exp_arr


    # use np.split to remove the needle
    _,arr_nn=np.split(arr_won,[stopline])


    # get the pixel of needle
    dneedle=needle[0]

    # reshape the array let the pixel of needle's diameter is 50
    scape=50/dneedle
    arr_50=scipy.ndimage.zoom(arr_nn,scape)
    arr_50[arr_50<125]=0
    arr_50[arr_50>=125]=255

    return arr_50


def genDropletNoMargin(sigma,volume0,rneedle=0.5):
    r,z=genSingleDrop(sigma=sigma,rneedle=rneedle,volume0=volume0,output=1)    
    plt.plot(r,z,'black')
    plt.plot(-r,z,'black')
    plt.axis('equal')
    plt.axis('off')
    figure = plt.gcf().canvas
    ag = figure.switch_backends(FigureCanvasAgg)
    ag.draw()
    plt.close()
    A = np.asarray(ag.buffer_rgba())
    syn_arr = np.rint(A[...,:3] @ [0.2126, 0.7152, 0.0722]).astype(np.uint8)

    h_syn,w_syn=syn_arr.shape

    # fill the droplet to black
    end=0
    left_edge=np.zeros(h_syn)
    right_edge=np.zeros(h_syn)
    p_needle=1
    d_needle=50
    for i in range(h_syn):
        left=0
        right=0 
        for j in range (w_syn):
            if syn_arr[i,j-1]>=125 and syn_arr[i,j]<125 and j < w_syn/2:
                left=j
                end=end+1
        for k in range (w_syn-1,0,-1):
            if syn_arr[i,k]>=125 and syn_arr[i,k-1]<125 and k > w_syn/2:
                right=k
        for l in range (w_syn):
            if l>left and l<right:
                syn_arr[i,l]=0
        left_edge[i]=left
        right_edge[i]=right
        if left!=0 and end==1:
            p_needle=i
            d_needle=right-left


    # cut up and sides margin
    _,syn_re=np.split(syn_arr,[p_needle])
    left_edge[left_edge==0]=w_syn/2
    cut_l=min(left_edge)
    cut_r=max(right_edge)
    _,syn_nm,_=np.split(syn_re,[int(cut_l),int(cut_r)],axis=1)

    # zoom array let rneedle is 50 pixels
    scape=50/d_needle

    syn50=scipy.ndimage.zoom(syn_nm,scape)
    syn50[syn50<125]=0
    syn50[syn50>=125]=255

    # cut bottom margin
    cut_b=h_syn
    h_syn,w_syn=syn50.shape
    for i in range (h_syn):
        if sum(syn50[i,:])==255*w_syn and sum(syn50[i-1])<255*w_syn:
            cut_b=i

    syn_nm,_=np.split(syn50,[int(cut_b)])
    return syn_nm

def lostFunction(exp_arr,syn_arr,output=0):
    h_exp,w_exp=exp_arr.shape
    h_syn,w_syn=syn_arr.shape

    C=np.zeros((h_exp,w_exp))

    droplet=0

    x=int(abs((w_exp-w_syn)/2))

    for i in range(h_exp):
        for j in range(w_exp):
            if exp_arr[i,j]==0:
                droplet=droplet+1


    for i in range(h_exp):
        for j in range(w_exp):
            if j<x or j-x>=w_syn or i>=h_syn:
                C[i,j]=exp_arr[i,j]
            else:
                if exp_arr[i,j]==syn_arr[i,j-x]:
                    C[i,j]=255
                else:
                    C[i,j]=0

    f=0; t=0; 
    for i in range(h_exp):
        for j in range(w_exp):
            if C[i,j]==0:
                f=f+1

    lost=f/droplet
    if output==0:
        return lost
    else:
        return lost,C

def likelihood(sv,path_exp,output=0):    
    sigma=sv[0]; v0=sv[1]

    exp_arr=removeStandNeedleCentering(path_exp)
    syn_arr=genDropletNoMargin(sigma=sigma,volume0=v0)
    lost,C=lostFunction(exp_arr=exp_arr,syn_arr=syn_arr,output=1)

    if output==0:
        return lost
    else:
        return lost,C

def genDropletWithNeedle(sigma,volume0,rneedle=0.5):
    r,z=genSingleDrop(sigma=sigma,rneedle=rneedle,volume0=volume0,output=1)    
    plt.plot(r,z,'black')
    plt.plot(-r,z,'black')
    plt.axis('equal')
    plt.axis('off')
    figure = plt.gcf().canvas
    ag = figure.switch_backends(FigureCanvasAgg)
    ag.draw()
    plt.close()
    A = np.asarray(ag.buffer_rgba())
    syn_arr = np.rint(A[...,:3] @ [0.2126, 0.7152, 0.0722]).astype(np.uint8)

    h_syn,w_syn=syn_arr.shape

    # fill the droplet to black
    end=0
    left_edge=np.zeros(h_syn)
    right_edge=np.zeros(h_syn)
    p_needle=1
    d_needle=50
    for i in range(h_syn):
        left=0
        right=0 
        for j in range (w_syn):
            if syn_arr[i,j-1]>=125 and syn_arr[i,j]<125 and j < w_syn/2:
                left=j
                end=end+1
        for k in range (w_syn-1,0,-1):
            if syn_arr[i,k]>=125 and syn_arr[i,k-1]<125 and k > w_syn/2:
                right=k
        for l in range (w_syn):
            if l>left and l<right:
                syn_arr[i,l]=0
        left_edge[i]=left
        right_edge[i]=right
        if left!=0 and end==1:
            p_needle=i

    for i in range(h_syn):
        if i<=p_needle:
            syn_arr[i,int(left_edge[p_needle]):int(right_edge[p_needle])]=0    

    img_arr=Image.fromarray(syn_arr)
    img_arr.save('./s%.2f_v%.2f_n%.2f.jpg'%(sigma,volume0,rneedle))
    

    return img_arr