# Generated with SMOP  0.41
from libsmop import *
# dif1D.m

    
@function
def dif1D(type_=None,s0=None,L=None,N=None,pts=None,*args,**kwargs):
    varargin = dif1D.varargin
    nargin = dif1D.nargin

    # [d/ds, d/ds^2, integral, s] = dif1D('type',s0,length,N_dof,order), creates
# the 1D differentiation matrices
# These functions have been collected from Jerome Hoepffners teaching
# materials at http://basilisk.fr/sandbox/easystab/README
    
    if 'fd' == type_:
        scale=L / 2
# dif1D.m:9
        __,d=fddif(N,1,pts,nargout=2)
# dif1D.m:10
        s,dd=fddif(N,2,pts,nargout=2)
# dif1D.m:11
        s=dot(s,scale)
# dif1D.m:12
        d=(d / scale)
# dif1D.m:13
        dd=dd / scale ** 2
# dif1D.m:13
        s=s - s(1) + s0
# dif1D.m:13
        d=full(d)
# dif1D.m:14
        w=(concat([diff(s.T),0]) + concat([0,diff(s.T)])) / 2
# dif1D.m:15
    else:
        if 'cheb' == type_:
            scale=- L / 2
# dif1D.m:17
            s,DM=chebdif(N,2,nargout=2)
# dif1D.m:18
            d=DM(arange(),arange(),1)
# dif1D.m:19
            dd=DM(arange(),arange(),2)
# dif1D.m:20
            s=dot(s,scale)
# dif1D.m:21
            d=d / scale
# dif1D.m:22
            dd=dd / scale ** 2
# dif1D.m:22
            s=s - s(1) + s0
# dif1D.m:22
            w=dot(L,clencurt(N)) / 2
# dif1D.m:23
    
    return d,dd,w,s
    
if __name__ == '__main__':
    pass
    
    
@function
def fddif(N=None,order=None,pts=None,*args,**kwargs):
    varargin = fddif.varargin
    nargin = fddif.nargin

    # build equispaced grid on [-1,1], and 
# five points finite diference matrix for N mesh points
    
    #pts=5;
    x=linspace(- 1,1,N).T
# dif1D.m:33
    h=x(2) - x(1)
# dif1D.m:34
    # subroutine for finite difference weights
    W=ufdwt(h,pts,order)
# dif1D.m:37
    t=(pts + 1) / 2
# dif1D.m:38
    #central difference in the middle
    D=spdiags(dot(ones(N,1),W(t,arange())),arange(- t + 1,t - 1),N,N)
# dif1D.m:41
    for indd in arange(1,t - 1).reshape(-1):
        D[indd,arange(1,pts)]=W(indd,arange())
# dif1D.m:44
        D[N - indd + 1,arange(end() - pts + 1,end())]=W(end() - indd + 1,arange())
# dif1D.m:45
    
    return x,D
    
if __name__ == '__main__':
    pass
    
    # Diffm compute D = differentiation matrix, x = Chebyshev grid
# Order is N, in the interval a to b
#  function [D,x] = dmat(N,a,b)
#  if N==0, D=0; x=1; return, end
#  if a>=b, D=0; x=0; disp('point a smaller than b'); return, end 
#  x = (a-b)/2*cos(pi*(0:N)/N)' + (a+b)/2;
#  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
#  X = repmat(x,1,N+1);
#  dX = X-X';                  
#  D  = (c*(1./c)')./(dX+(eye(N+1)));      # off-diagonal entries
#  D  = D - diag(sum(D,2));                 # diagonal entries
#  end
    
    
@function
def ufdwt(h=None,pts=None,order=None,*args,**kwargs):
    varargin = ufdwt.varargin
    nargin = ufdwt.nargin

    ################################################################################
    
    # ufdwt.m
    
    # Compute Finite Difference Weights for a Uniform Grid
    
    # Input Parameters:
    
    # h     - spacing between FD nodes
# pts   - number of FD points in scheme (3-pt, 5-pt, etc)
# order - order of the derivative operator to compute weights for
#         (note: order<pts-1!)
#         1 computes first derivative differences       
#         2 computes second derivative differences, etc
#  
# Output Parameter:
    
    # W is the weight matrix. Each row contains a different set of weights
# (centered or off). If, for example, the number of finite difference points
# is odd, the centered difference weights will appear in the middle row.
    
    # Written by: Greg von Winckel - 06/16/04
# Contact: gregvw@chtm.unm.edu
    
    ################################################################################
    
    N=dot(2,pts) - 1
# dif1D.m:92
    p1=pts - 1
# dif1D.m:92
    A=repmat((arange(0,p1)).T,1,N)
# dif1D.m:94
    B=repmat(dot((arange(- p1,p1)),h),pts,1)
# dif1D.m:95
    M=(B ** A) / gamma(A + 1)
# dif1D.m:97
    rhs=zeros(pts,1)
# dif1D.m:99
    rhs[order + 1]=1
# dif1D.m:99
    W=zeros(pts,pts)
# dif1D.m:101
    for k in arange(1,pts).reshape(-1):
        W[arange(),k]=numpy.linalg.solve(M(arange(),(arange(0,p1)) + k),rhs)
# dif1D.m:104
    
    W=W.T
# dif1D.m:107
    W[arange(1,pts),arange()]=W(arange(pts,1,- 1),arange())
# dif1D.m:107
    return W
    
if __name__ == '__main__':
    pass
    
    
@function
def poldif(x=None,malpha=None,B=None,*args,**kwargs):
    varargin = poldif.varargin
    nargin = poldif.nargin

    #  The function DM =  poldif(x, maplha, B) computes the
#  differentiation matrices D1, D2, ..., DM on arbitrary nodes.
    
    #  The function is called with either two or three input arguments.
#  If two input arguments are supplied, the weight function is assumed 
#  to be constant.   If three arguments are supplied, the weights should 
#  be defined as the second and third arguments.
    
    #  Input (constant weight):
    
    #  x:        Vector of N distinct nodes.
#  malpha:   M, the number of derivatives required (integer).
#  B:        Omitted.
    
    #  Note:     0 < M < N-1.
    
    #  Input (non-constant weight):
    
    #  x:        Vector of N distinct nodes.
#  malpha:   Vector of weight values alpha(x), evaluated at x = x(k).
#  B:        Matrix of size M x N,  where M is the highest 
#            derivative required.  It should contain the quantities 
#            B(ell,j) = beta(ell,j) = (ell-th derivative
#            of alpha(x))/alpha(x),   evaluated at x = x(j).
    
    #  Output:
#  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
    
    #  J.A.C. Weideman, S.C. Reddy 1998
    
    N=length(x)
# dif1D.m:143
    x=ravel(x)
# dif1D.m:144
    
    if nargin == 2:
        M=copy(malpha)
# dif1D.m:147
        alpha=ones(N,1)
# dif1D.m:148
        B=zeros(M,N)
# dif1D.m:149
    else:
        if nargin == 3:
            alpha=ravel(malpha)
# dif1D.m:151
            M=length(B(arange(),1))
# dif1D.m:152
    
    
    I=eye(N)
# dif1D.m:155
    
    L=logical(I)
# dif1D.m:156
    
    XX=x(arange(),ones(1,N))
# dif1D.m:158
    DX=XX - XX.T
# dif1D.m:159
    
    DX[L]=ones(N,1)
# dif1D.m:161
    
    c=multiply(alpha,prod(DX,2))
# dif1D.m:163
    
    C=c(arange(),ones(1,N))
# dif1D.m:165
    C=C / C.T
# dif1D.m:166
    
    
    Z=1.0 / DX
# dif1D.m:168
    
    Z[L]=zeros(N,1)
# dif1D.m:169
    
    X=Z.T
# dif1D.m:171
    
    X[L]=[]
# dif1D.m:172
    
    X=reshape(X,N - 1,N)
# dif1D.m:173
    Y=ones(N - 1,N)
# dif1D.m:175
    
    D=eye(N)
# dif1D.m:176
    
    DM=zeros(N,N,M)
# dif1D.m:177
    
    for ell in arange(1,M).reshape(-1):
        Y=cumsum(concat([[B(ell,arange())],[multiply(dot(ell,Y(arange(1,N - 1),arange())),X)]]))
# dif1D.m:179
        D=multiply(dot(ell,Z),(multiply(C,repmat(diag(D),1,N)) - D))
# dif1D.m:180
        D[L]=Y(N,arange())
# dif1D.m:181
        DM[arange(),arange(),ell]=D
# dif1D.m:182
    
    return DM
    
if __name__ == '__main__':
    pass
    
    
@function
def clencurt(N=None,*args,**kwargs):
    varargin = clencurt.varargin
    nargin = clencurt.nargin

    
    # Computes the integration weigths for pseudo-chebychev on domain [-1 1]
    
    # INPUTS:
# N  : the number of points
    
    # OUTPUT:
# IW : vector of the integration weigths.
    
    nW=arange(0,N - 1,1)
# dif1D.m:197
    jW=arange(0,N - 1,1)
# dif1D.m:198
    bW=ones(1,N)
# dif1D.m:200
    bW[1]=0.5
# dif1D.m:200
    bW[N]=0.5
# dif1D.m:200
    cW=dot(2,bW)
# dif1D.m:201
    bW=bW / (N - 1)
# dif1D.m:202
    S=cos(dot(dot(nW(arange(3,N)).T,jW),(pi / (N - 1))))
# dif1D.m:204
    IW=multiply(bW,((2 + dot((multiply(cW(arange(3,N)),((1 + (- 1) ** nW(arange(3,N))) / (1 - nW(arange(3,N)) ** 2)))),S))))
# dif1D.m:205
    return IW
    
if __name__ == '__main__':
    pass
    
    
@function
def chebdif(N=None,M=None,*args,**kwargs):
    varargin = chebdif.varargin
    nargin = chebdif.nargin

    #  The function DM =  chebdif(N,M) computes the differentiation 
#  matrices D1, D2, ..., DM on Chebyshev nodes. 
# 
#  Input:
#  N:        Size of differentiation matrix.        
#  M:        Number of derivatives required (integer).
#  Note:     0 < M <= N-1.
    
    #  Output:
#  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
    
    #  The code implements two strategies for enhanced 
#  accuracy suggested by W. Don and S. Solomonoff in 
#  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
#  The two strategies are (a) the use of trigonometric 
#  identities to avoid the computation of differences 
#  x(k)-x(j) and (b) the use of the "flipping trick"
#  which is necessary since sin t can be computed to high
#  relative precision when t is small whereas sin (pi-t) cannot.
    
    #  J.A.C. Weideman, S.C. Reddy 1998.
    
    I=eye(N)
# dif1D.m:233
    
    L=logical(I)
# dif1D.m:234
    
    n1=floor(N / 2)
# dif1D.m:236
    n2=ceil(N / 2)
# dif1D.m:236
    
    k=(arange(0,N - 1)).T
# dif1D.m:238
    
    th=dot(k,pi) / (N - 1)
# dif1D.m:239
    x=sin(dot(pi,(arange(N - 1,1 - N,- 2)).T) / (dot(2,(N - 1))))
# dif1D.m:241
    
    T=repmat(th / 2,1,N)
# dif1D.m:243
    DX=multiply(dot(2,sin(T.T + T)),sin(T.T - T))
# dif1D.m:244
    
    DX=concat([[DX(arange(1,n1),arange())],[- rot90(DX(arange(1,n2),arange()),2)]])
# dif1D.m:245
    
    DX[L]=ones(N,1)
# dif1D.m:246
    
    C=toeplitz((- 1) ** k)
# dif1D.m:248
    
    C[1,arange()]=dot(C(1,arange()),2)
# dif1D.m:249
    C[N,arange()]=dot(C(N,arange()),2)
# dif1D.m:249
    
    C[arange(),1]=C(arange(),1) / 2
# dif1D.m:250
    C[arange(),N]=C(arange(),N) / 2
# dif1D.m:250
    Z=1.0 / DX
# dif1D.m:252
    
    Z[L]=zeros(N,1)
# dif1D.m:253
    
    D=eye(N)
# dif1D.m:255
    
    DM=zeros(N,N,M)
# dif1D.m:256
    for ell in arange(1,M).reshape(-1):
        D=multiply(dot(ell,Z),(multiply(C,repmat(diag(D),1,N)) - D))
# dif1D.m:258
        D[L]=- sum(D,2).T
# dif1D.m:259
        DM[arange(),arange(),ell]=D
# dif1D.m:260
    
    return x,DM
    
if __name__ == '__main__':
    pass
    