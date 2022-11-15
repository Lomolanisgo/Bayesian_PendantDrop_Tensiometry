# Generated with SMOP  0.41
from libsmop import *
# fchebt.m

    
@function
def fchebt(y=None,max_order=None,shift_order=None,*args,**kwargs):
    varargin = fchebt.varargin
    nargin = fchebt.nargin

    # [coefficient, filtered function] = fchebt(function, max_order, shift_order)
# is a Fast Chebychev Transform, based on first kind Chebyshev Function T.
# The function `y` is approximated up to `max_order`. When `max_order` is
# negative the order is optimized based on the performance of each new
# order. `shift_order` will push the optimal order.
    
    # Chebyshev polynomials have singular end points but are better weighted
# (modes not increasing with wavenumber)
# pass 3 times and make a final error check showed better precision
    
    global g_error
    # if maximum order is negative the response is computed as f
    
    f=0
# fchebt.m:16
    N=length(y)
# fchebt.m:17
    # Matlab Chebyshev is too slow
    YC=fcbt('cheb',N,abs(max_order) + 1)
# fchebt.m:20
    if size(y,1) == 1:
        y=y.T
# fchebt.m:23
    
    # backup function y for final comparison
    ys0=copy(y)
# fchebt.m:27
    # initialize arguments if not passed
    if logical_not(exist('max_order','var')):
        max_order=N / 2
# fchebt.m:31
    
    if logical_not(exist('shift_order','var')):
        shift_order=0
# fchebt.m:35
    
    # revarg detemines if the truncation error is evaluated
    revarg=(max_order < 0)
# fchebt.m:39
    max_order=abs(max_order)
# fchebt.m:40
    ds=2 / (N - 1)
# fchebt.m:41
    
    # create differentiation, integration and ordinate
    d,__,w,s=dif1D('fd',- 1 + ds,2 - dot(2,ds),N - 2,5,nargout=4)
# fchebt.m:44
    
    # w does overall integration
# create w2 for point wise integration
    d[1,arange()]=concat([1,zeros(1,N - 3)])
# fchebt.m:48
    w2=inv(d)
# fchebt.m:49
    w2=w2(end(),arange())
# fchebt.m:50
    w2[1]=0
# fchebt.m:51
    c=zeros(max_order + 1,1)
# fchebt.m:53
    #tic;
#se = [-1; s; 1];
#we = [w(1:10),w(10:11),w(11:end)];
    
    # 1st pass over 6 modes, c are the modes, from y the fitted part is subtracted
    for k in arange(0,min(6,max_order)).reshape(-1):
        # 2nd order singular integration analytically
        y0=dot(YC(1,k + 1),y(1))
# fchebt.m:61
        dy=(dot(YC(2,k + 1),y(2)) - dot(YC(1,k + 1),y(1))) / ds
# fchebt.m:62
        singular_left=dot((y0 + dy),acos(1 - ds)) - dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:63
        y0=dot(YC(end(),k + 1),y(end()))
# fchebt.m:64
        dy=(dot(YC(end(),k + 1),y(end())) - dot(YC(end() - 1,k + 1),y(end() - 1))) / ds
# fchebt.m:65
        singular_right=dot((y0 - dy),acos(1 - ds)) + dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:66
        c[k + 1]=dot(2 / pi,(dot(w2,(multiply(y(arange(2,end() - 1)) / sqrt(1 - s ** 2),YC(arange(2,end() - 1),k + 1)))) + singular_left + singular_right)) / (1 + (k == 0))
# fchebt.m:69
        y=y - dot(c(k + 1),YC(arange(),k + 1))
# fchebt.m:70
    
    # 2nd pass till 12 modes
    for k in arange(0,min(12,max_order)).reshape(-1):
        y0=dot(YC(1,k + 1),y(1))
# fchebt.m:75
        dy=(dot(YC(2,k + 1),y(2)) - dot(YC(1,k + 1),y(1))) / ds
# fchebt.m:76
        singular_left=dot((y0 + dy),acos(1 - ds)) - dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:77
        y0=dot(YC(end(),k + 1),y(end()))
# fchebt.m:78
        dy=(dot(YC(end(),k + 1),y(end())) - dot(YC(end() - 1,k + 1),y(end() - 1))) / ds
# fchebt.m:79
        singular_right=dot((y0 - dy),acos(1 - ds)) + dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:80
        tmp=dot(2 / pi,(dot(w,(multiply(y(arange(2,end() - 1)) / sqrt(1 - s ** 2),YC(arange(2,end() - 1),k + 1)))) + singular_left + singular_right)) / (1 + (k == 0))
# fchebt.m:82
        c[k + 1]=c(k + 1) + tmp
# fchebt.m:83
        y=y - dot(tmp,YC(arange(),k + 1))
# fchebt.m:84
    
    # 3rd pass over all modes
    for k in arange(0,max_order).reshape(-1):
        y0=dot(YC(1,k + 1),y(1))
# fchebt.m:89
        dy=(dot(YC(2,k + 1),y(2)) - dot(YC(1,k + 1),y(1))) / ds
# fchebt.m:90
        singular_left=dot((y0 + dy),acos(1 - ds)) - dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:91
        y0=dot(YC(end(),k + 1),y(end()))
# fchebt.m:92
        dy=(dot(YC(end(),k + 1),y(end())) - dot(YC(end() - 1,k + 1),y(end() - 1))) / ds
# fchebt.m:93
        singular_right=dot((y0 - dy),acos(1 - ds)) + dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:94
        tmp=dot(2 / pi,(dot(w,(multiply(y(arange(2,end() - 1)) / sqrt(1 - s ** 2),YC(arange(2,end() - 1),k + 1)))) + singular_left + singular_right)) / (1 + (k == 0))
# fchebt.m:96
        c[k + 1]=c(k + 1) + tmp
# fchebt.m:97
        y=y - dot(tmp,YC(arange(),k + 1))
# fchebt.m:98
    
    c2=zeros(max_order + 1,1)
# fchebt.m:101
    # 4th pass over all modes (several passes showed lower error)
    for k in arange(0,max_order).reshape(-1):
        y0=dot(YC(1,k + 1),y(1))
# fchebt.m:105
        dy=(dot(YC(2,k + 1),y(2)) - dot(YC(1,k + 1),y(1))) / ds
# fchebt.m:106
        singular_left=dot((y0 + dy),acos(1 - ds)) - dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:107
        y0=dot(YC(end(),k + 1),y(end()))
# fchebt.m:108
        dy=(dot(YC(end(),k + 1),y(end())) - dot(YC(end() - 1,k + 1),y(end() - 1))) / ds
# fchebt.m:109
        singular_right=dot((y0 - dy),acos(1 - ds)) + dot(dy,sqrt(dot(ds,(2 - ds))))
# fchebt.m:110
        c2[k + 1]=dot(2 / pi,(dot(w,(multiply(y(arange(2,end() - 1)) / sqrt(1 - s ** 2),YC(arange(2,end() - 1),k + 1)))) + singular_left + singular_right)) / (1 + (k == 0))
# fchebt.m:112
    
    c=c + c2
# fchebt.m:115
    g_error=- 1
# fchebt.m:117
    errmin=1000000000000.0
# fchebt.m:118
    if revarg:
        errfind=zeros(max_order + 1,1)
# fchebt.m:120
        for k in arange(0,max_order).reshape(-1):
            f=f + dot(c(k + 1),YC(arange(),k + 1))
# fchebt.m:123
            errfind[k + 1]=rms(f - ys0)
# fchebt.m:124
            errmin=min(errfind(k + 1),errmin)
# fchebt.m:125
        # look for the last mode that brought the error within 10# of that
        k=find(errfind < dot(1.5,errmin),1,'first')
# fchebt.m:129
        k=min(k + shift_order,length(c) - 1)
# fchebt.m:132
        c[arange(k + 1,end())]=0
# fchebt.m:135
        g_error=copy(errfind)
# fchebt.m:136
    
    return c,f