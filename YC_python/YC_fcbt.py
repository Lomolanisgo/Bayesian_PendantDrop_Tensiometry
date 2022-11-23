import numpy as np
import os
str=('dif1D.py')
p=os.system(str)
import dif1D

def fcbt(diffmode,N,M,halfmark):
    '''chebyshev_polynomial = fcbt(N,M,halfmark) creates the first 30 chebyshev polynomials
    % for N points and to M order, halfmark is one if only half the drop shape
    % is fitted 0 to 1 instead of -1 to 1.'''
    # create the domain
    if halfmark!= None and var!=None: # there is no input named var
        if halfmark==1 and diffmode=='fd':
            x = np.linspace(0,1,N).reshape(1,N).T
        elif halfmark==1 and diffmode=='cheb':
            [_,_,_,x]=dif1D('cheb',0,1,N)
        else:

    return [Y]