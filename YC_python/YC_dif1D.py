
def dif1D(type,s0,L,N,pts):
    # [d/ds, d/ds^2, integral, s] = dif1D('type',s0,length,N_dof,order), creates
    # the 1D differentiation matrices
    # These functions have been collected from Jerome Hoepffners teaching
    # materials at http://basilisk.fr/sandbox/easystab/README

    if type=='fd': #finite difference
        scale=L/2
        [~,d] = fddif(N,1,pts); 
        [s,dd] = fddif(N,2,pts); 
        s=s*scale; 
        d=(d/scale); dd=dd/scale^2; s=s-s(1)+s0
        d = full(d)
        w=([diff(s'),0]+[0,diff(s')])/2
    elif type=='cheb': # chebychev
        scale=-L/2
    [s,DM] = chebdif(N,2)
    d=DM(:,:,1)
    dd=DM(:,:,2)
    s=s*scale
    d=d/scale 
    dd=dd/scale^2
    s=s-s(1)+s0
    w=L*clencurt(N)/2

    return[d,dd,w,s]

def fddif(N,order,pts):
    # build equispaced grid on [-1,1], and 
    # five points finite diference matrix for N mesh points 

    #pts=5;
    x=linspace(-1,1,N).T
    h=x(2)-x(1)



    return [x,D]