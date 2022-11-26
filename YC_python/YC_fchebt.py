import numpy as np
import math as m
from YC_fcbt import fcbt
from dif1D import *
pi=m.pi
def fchebt(y,max_order,shift_order):
    '''shift_order=0
    '''
    f=0
    N=len(y)
    
    YC=fcbt('cheb',N,abs(max_order)+1)

    if np.size(y,1)==1:
        y=y.T

    ys0=y

    if max_order == None:
        max_order=N/2

    if shift_order == None:
        shift_order = 0

    revarg = (max_order<0)
    max_order=abs(max_order)
    ds=2/N-1

    d,_,w,s=dif1D('fd',-1+ds,2-2*ds,N-2,5)

    d[0,:]=np.hstack((1,np.zeros(N-3)))
    w2= np.linalg.inv(d)
    w2=w2[-1,:]
    w2[0]=0

    c=np.zeros((max_order+1,1))

    # 1st pass over 6 modes, c are the modes, from y the fitted part is subtracted
    for k in range (0,min(6,max_order)+1):
        y0 = YC [0,k]*y[0]
        dy=(YC[1,k]*y(1)-YC[0,k]*y[0])/ds

        singular_left = (y0+dy)*np.acos(1-ds)-dy*np.sqrt(ds*(2-ds))
        y0 = YC[-1,k]*y[-1]
        dy = (YC[-1,k]*y[-1]-YC[-2,k]*y[-2])/ds
        singular_right = (y0-dy)*np.acos(1-ds)+dy*np.sqrt(ds*(2-ds))

        c[k] = 2/pi*(w2*(y[1:-2]/np.sqrt(1-s**2)*YC[1:-2,k]) + singular_left + singular_right)/(1+(k==0))
        y = y - c[k]*YC[:,k]

    for k in range (0,min(12,max_order)+1):
            y0 = YC [0,k]*y[0]
            dy=(YC[1,k]*y(1)-YC[0,k]*y[0])/ds

            singular_left = (y0+dy)*np.acos(1-ds)-dy*np.sqrt(ds*(2-ds))
            y0 = YC[-1,k]*y[-1]
            dy = (YC[-1,k]*y[-1]-YC[-2,k]*y[-2])/ds
            singular_right = (y0-dy)*np.acos(1-ds)+dy*np.sqrt(ds*(2-ds))

            tmp = 2/pi*(w*(y[1:-2]/np.sqrt(1-s**2)*YC[1:-2,k]) + singular_left + singular_right)/(1+(k==0))
            c[k]=c[k]+tmp
            y = y - c[k]*YC[:,k]

    for k in range (0,max_order+1):
                y0 = YC [0,k]*y[0]
                dy=(YC[1,k]*y(1)-YC[0,k]*y[0])/ds
    
                singular_left = (y0+dy)*np.acos(1-ds)-dy*np.sqrt(ds*(2-ds))
                y0 = YC[-1,k]*y[-1]
                dy = (YC[-1,k]*y[-1]-YC[-2,k]*y[-2])/ds
                singular_right = (y0-dy)*np.acos(1-ds)+dy*np.sqrt(ds*(2-ds))
    
                tmp = 2/pi*(w*(y[1:-2]/np.sqrt(1-s**2)*YC[1:-2,k]) + singular_left + singular_right)/(1+(k==0))
                c[k]=c[k]+tmp
                y = y - c[k]*YC[:,k]

    return c,f