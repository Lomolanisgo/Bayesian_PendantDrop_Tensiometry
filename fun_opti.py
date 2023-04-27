import numpy as np
import matplotlib.pyplot as plt
from fun_genSingleDrop import *
from fun_preprocess import *
from matplotlib.backends.backend_agg import FigureCanvasAgg

# Step 1: Remove the needle and centering the experimental (exp) image. 
def remove_Needle_Center(img_ori):#(Need Test)
    '''img_ori=Image.open(path_ori)'''
    w_ori,h_ori=img_ori.size

    A=img2bw(ResizeImage(img_ori,w_ori,h_ori))

    left_edge=np.zeros(w_ori)
    right_edge=np.zeros(w_ori)
    for i in range(h_ori):
        for j in range(w_ori):
            if A[i,j-1]>=125 and A[i,j]<125:
                #print('Left position',i,j)
                left_edge[i]=j
            if A[i,j-1]<=125 and A[i,j]>125:
                #print('Right position',i,j)
                right_edge[i]=j

    needle=np.abs(left_edge-right_edge)

    end=0
    for j in range(needle.shape[0]-1):
        j=j+1
        if end == 0:
            if abs((needle[j]-needle[0])/needle[0])>0.05:
                stopline=j
                end=end+1
    #print(stopline)

    addpad=right_edge[stopline]-needle[stopline]/2-w_ori/2

    img_WON=np.zeros((h_ori-stopline,w_ori))
    for i in range(h_ori):
        if i>=stopline:
            img_WON[i-stopline,:]=A[i,:]
    
    h_won,w_won=img_WON.shape
    
    # centering
    if addpad>0:
        img_won=np.hstack( ( img_WON,np.ones((h_won,abs(int(addpad*4))))*255 ) )
    elif addpad<0:
        img_won=np.hstack( ( np.ones((h_won,abs(int(addpad*4))))*255,img_WON ) )
    else:
        img_won=img_WON

    return img_won,needle

# Step 2: Reshape the exp images.
def reshape_Exp(path_ori,output=0,d_needle=100,threshold=125):
    '''
    d_needle: the amount of pixel of needle's diameter
    output: if not 0, visualize the exp_zoom image 
    threshold: of the binarization
    '''
    # import the exp image
    exp_array0,_=remove_Needle_Center(Image.open(path_ori))
    exp_array=exp_array0
    exp_array[exp_array<threshold]=0
    exp_array[exp_array>=threshold]=255

    # Scaling the syn image to a needle diameter of 100 pixel(r=50 pixel)
    k_exp=d_needle/np.sum(exp_array[0]==0)
    exp_zoom=scipy.ndimage.zoom(exp_array,k_exp)
    # Binarization
    exp_zoom[exp_zoom<threshold]=0
    exp_zoom[exp_zoom>=threshold]=255

    if output!=0:
        plt.imshow(exp_zoom,cmap='gray')
    
    return exp_zoom

# Step 3: Generate the synthetic droplet
def gen_Drop(sigma,volume0,rneedle=0.5):
    # when output=1 gSD return r_a,z_a
    r_a,z_a=genSingleDrop(sigma=sigma,volume0=volume0,rneedle=rneedle,output=1)
    wmax=2*max(abs(r_a))
    r_a0=r_a[0]
    # use plt generate synthetic image and save as arrayd
    plt.fill_between(r_a,z_a,color='black')
    plt.fill_between(-r_a,z_a,color='black')
    plt.axis('equal')
    plt.axis('off')

    figure = plt.gcf().canvas
    ag = figure.switch_backends(FigureCanvasAgg)
    ag.draw()
    plt.close()
    A = np.asarray(ag.buffer_rgba())
    syn_array = np.rint(A[...,:3] @ [0.2126, 0.7152, 0.0722]).astype(np.uint8)
    
    return syn_array,wmax,rneedle,r_a0

# Step 4: Reshape the synthetic (syn) images
def scale_Syn(syn_array,output=0):
    # crop the white margin of the array
    syn_array_cm0=crop_Margin(syn_array)
    # Remove the non-droplet part at the top of the image
    syn_array_cm=syn_array_cm0[2:-1]

    # Scaling the syn image to a needle diameter of 100 pixel(r=50 pixel)

    # The part of the value of 0 is the diameter of the needle
    d_pixel=np.sum(syn_array_cm[0]==0)
    if d_pixel==0:
        d_pixel=1
    k_syn=100/d_pixel
    syn_zoom=scipy.ndimage.zoom(syn_array_cm,k_syn)

    return syn_zoom

# Step 5: Calculate the NOV
def cal_NOV(exp_zoom,syn_zoom,output=0):
    h_exp,w_exp=np.shape(exp_zoom)
    h_syn,w_syn=np.shape(syn_zoom)

    if h_exp<h_syn or w_exp<w_syn:
        NOV=3
        return NOV

    # droplet pixel amount
    droplet=np.sum(exp_zoom==0)

    a0=np.where(exp_zoom[0]==0)[0][0]
    a1=np.where(syn_zoom[0]==0)[0][0]
    x0=a0-a1
    x2=w_exp-x0-w_syn
    y1=h_exp-h_syn
    if x0<0 or x2<0 or y1<0:
        NOV=2
        return NOV

    # syn_arr reshape by add margin
    B1=np.ones((h_syn,x0))*255
    B2=np.ones((h_syn,x2))*255
    B3=np.ones((y1,w_exp))*255
    B_up=np.hstack((B1,syn_zoom,B2))
    B=np.vstack((B_up,B3))

    A=exp_zoom
    A[A==0]=True
    A[A==255]=False

    B[B==0]=True
    B[B==255]=False
    
    C=A+B
    C[C==0]=255
    C[C==2]=255
    C[C==1]=0

    NOV=abs(np.sum(C==0)/droplet)

    # For debugging
    if output!=0:
        return NOV,C
    else:
        return NOV
    

