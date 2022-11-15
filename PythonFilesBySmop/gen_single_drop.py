# Generated with SMOP  0.41
import sys
sys.path.append('smop\libsmop.py')
from datetime import datetime
import libsmop
# gen_single_drop.m

    #close_('all')
    #clear
addpath('subs/')
    # physical parameters
sigma=100
# gen_single_drop.m:7
    
grav=9807.0
# gen_single_drop.m:8
    
rneedle=1
# gen_single_drop.m:9
    
volume0=32
# gen_single_drop.m:10
    
deltarho=0.001
# gen_single_drop.m:11
    
    # numerical parameters
N=40
# gen_single_drop.m:14
    
Nplot=80
# gen_single_drop.m:15
    
Ncheb=10
# gen_single_drop.m:16
    
alpha=0.1
# gen_single_drop.m:17
    
tic=datetime.now()
# NOTE: the calculation is done in dimensionless form, using the 
# dimensionless surface tension sigma' and volume V'
    
# calculate the dimensionless quantities
sigmaprime=sigma / (dot(dot(deltarho,grav),rneedle ** 2))
# gen_single_drop.m:25
volume0prime=volume0 / rneedle ** 3
# gen_single_drop.m:26
# find the initial guess of the droplet shape
if dot(dot(deltarho,grav),volume0) / (dot(dot(dot(2,pi),sigma),rneedle)) > 0.14:
# predict the maximum length of the interface (empirical Nagel)
smax=dot(sqrt(sigmaprime),2.0) / 0.8701
# gen_single_drop.m:32
D,__,w,s=dif1D('cheb',0,smax,N,5,nargout=4)
# gen_single_drop.m:35
z=dot(dot(- 4 / 3,smax) / pi,(cos(dot(dot(pi,3) / 4,s) / smax)))
# gen_single_drop.m:38
z=z - max(z)
# gen_single_drop.m:39
r=dot(dot(4 / 3,smax) / pi,(sin(dot(dot(pi,3) / 4,s) / smax)))
# gen_single_drop.m:40
psi=dot(dot(pi,3) / 4,s) / smax
# gen_single_drop.m:41
C=1
# gen_single_drop.m:43
p0=dot(sqrt(sigmaprime),1.5)
# gen_single_drop.m:44
else:
# find the roots for the polynomial of a spherical cap
rts=roots(concat([pi / 6,0,pi / 2,- volume0prime]))
# gen_single_drop.m:49
h0=real(rts(3))
# gen_single_drop.m:50
Rguess,xcyc=fit_circle_through_3_points(concat([[1,0],[0,- h0],[- 1,0]]),nargout=2)
# gen_single_drop.m:52
if xcyc(2) < 0:
        theta=acos(1 / Rguess)
# gen_single_drop.m:56
else:
        theta=- acos(1 / Rguess)
# gen_single_drop.m:58
        # predict the maximum length of the interface
smax=dot(Rguess,(dot(2,theta) + pi))
# gen_single_drop.m:62
D,__,w,s=dif1D('fd',0,smax,N,5,nargout=4)
# gen_single_drop.m:65
dtheta=linspace(- pi / 2,theta,N)
# gen_single_drop.m:68
dtheta=dtheta.T
# gen_single_drop.m:69
r=xcyc(1) + dot(Rguess,cos(dtheta))
# gen_single_drop.m:70
z=xcyc(2) + dot(Rguess,sin(dtheta))
# gen_single_drop.m:71
psi=atan2(dot(D,z),dot(D,r))
# gen_single_drop.m:73
C=1
# gen_single_drop.m:75
p0=dot(dot(2,Rguess),sigmaprime)
# gen_single_drop.m:76
# get the differentation/integration matrices and the grid
D,__,w,s=dif1D('cheb',0,smax,N,5,nargout=4)
# gen_single_drop.m:79
    
    
# initialize some variables
Z=zeros(N)
# gen_single_drop.m:84
    
IDL=concat([1,zeros(1,N - 1)])
# gen_single_drop.m:85
    
ZL=zeros(1,N)
# gen_single_drop.m:86
    
u=ones(dot(3,N) + 2,1)
# gen_single_drop.m:87
b=ones(dot(3,N) + 2,1)
# gen_single_drop.m:87
    
iter=0
# gen_single_drop.m:88
crash=0
# gen_single_drop.m:88
while rms(u) > 1e-10:

iter=iter + 1
# gen_single_drop.m:92
if iter > 1200:
        warning('iter > 12000!')
        crash=1
# gen_single_drop.m:96
        break
# determine r from psi
A11=dot(C,D)
# gen_single_drop.m:100
A13=diag(sin(psi))
# gen_single_drop.m:100
A18=dot(D,r)
# gen_single_drop.m:100
b1=- (dot(dot(C,D),r) - cos(psi))
# gen_single_drop.m:100
A22=dot(C,D)
# gen_single_drop.m:103
A23=diag(- cos(psi))
# gen_single_drop.m:103
A28=dot(D,z)
# gen_sngle_drop.m:103
b2=- (dot(dot(C,D),z) - sin(psi))
# gen_sngle_drop.m:103
A31=dot(- sigmaprime,diag(sin(psi) / r ** 2))
# gen_sngle_drop.m:106
A32=diag(ones(N,1))
# gen_sngle_drop.m:107
A33=dot(dot(C,sigmaprime),D) + dot(sigmaprime,diag(cos(psi) / r))
# gen_sngle_drop.m:108
A38=dot(sigmaprime,(dot(D,psi)))
# gen_sngle_drop.m:109
A39=- ones(N,1)
# gen_sngle_drop.m:110
b3=p0 - z - dot(sigmaprime,(dot(dot(C,D),psi) + sin(psi) / r))
# gen_sngle_drop.m:111
# NOTE: the lengths are scaled with the radius, thus its value is one
A81=fliplr(IDL)
# gen_sngle_drop.m:115
b8=(1 - r(end()))
# gen_sngle_drop.m:115
A91=multiply(multiply(dot(2,w),r.T),sin(psi.T))
# gen_sngle_drop.m:118
A93=dot(multiply(w,r.T ** 2.0),cos(psi.T))
# gen_sngle_drop.m:119
A98=- volume0prime / pi
# gen_sngle_drop.m:120
b9=- (dot(w,(dot(r ** 2.0,sin(psi)))) - dot(C,volume0prime) / pi)
# gen_sngle_drop.m:121
A11[1,arange()]=IDL
# gen_sngle_drop.m:124
A13[1,arange()]=ZL
# gen_sngle_drop.m:125
A18[1]=0
# gen_sngle_drop.m:126
b1[1]=- r(1)
# gen_sngle_drop.m:127
A22[1,arange()]=fliplr(IDL)
# gen_sngle_drop.m:130
A23[1,arange()]=ZL
# gen_sngle_drop.m:131
A28[1]=0
# gen_sngle_drop.m:132
b2[1]=- z(end())
# gen_sngle_drop.m:133
A31[1,arange()]=ZL
# gen_sngle_drop.m:136
A32[1,arange()]=ZL
# gen_sngle_drop.m:137
A33[1,arange()]=IDL
# gen_sngle_drop.m:138
A38[1,arange()]=0
# gen_sngle_drop.m:139
A39[1,arange()]=0
# gen_sngle_drop.m:140
b3[1]=- psi(1)
# gen_sngle_drop.m:141
Z1=zeros(N,1)
# gen_sngle_drop.m:144
A=concat([[concat([A11,Z,A13,A18,Z1])],[concat([Z,A22,A23,A28,Z1])],[concat([A31,A32,A33,A38,A39])],[concat([A81,zeros(1,dot(2,N)),- 1,0])],[concat([A91,Z1.T,A93,A98,0])]])
# gen_sngle_drop.m:146
b=concat([[b1],[b2],[b3],[b8],[b9]])
# gen_sngle_drop.m:150
u=numpy.linalg.solve(A,b)
# gen_sngle_drop.m:153
r=r + dot(alpha,u(arange(1,N)))
# gen_sngle_drop.m:156
z=z + dot(alpha,u(arange(N + 1,dot(2,N))))
# gen_sngle_drop.m:157
psi=psi + dot(alpha,u(arange(dot(2,N) + 1,dot(3,N))))
# gen_sngle_drop.m:158
C=C + dot(alpha,u(dot(3,N) + 1))
# gen_sngle_drop.m:159
p0=p0 + dot(alpha,u(dot(3,N) + 2))
# gen_sngle_drop.m:160
if rms(b) > 1000.0:
    crash=1
# gen_sngle_drop.m:163
    break

    
# calculate the Chebyshev coefficients
coefr=fchebt(r,Ncheb,0)
# gen_single_drop.m:169
coefz=fchebt(z,Ncheb,0)
# gen_single_drop.m:170
toc=datetime.now()
print('Elapsed time: %f seconds' % (toc-tic).total_seconds())

# compute volume and area (scaled back to dimensionfull)
disp(concat(['volume = ',num2str(dot(dot(dot(rneedle ** 3,pi),w),(dot(r ** 2.0,sin(psi)))) / C,12),' mm^3']))
disp(concat(['area = ',num2str(dot(dot(dot(dot(rneedle ** 2,pi),2),w),(r)) / C,12),' mm^2']))
disp(concat(['pressure = ',num2str(dot(dot(dot(deltarho,grav),rneedle),p0),12),' Pa']))
# # plot the shape of the drop on the numerical grid
# figure; hold on
# scatter(rneedle*r',rneedle*z','b');
# plot(rneedle*r',rneedle*z','b');
# set(gca,'DataAspectRatio',[1 1 1])
    
# interpolate the numerical solutions on a finer grid. 
# NOTE: the "right" way to interpolate is to fit a higher-orde polynomial 
# though all the points (see book of Trefethen on Spectral Methods in 
# Matlab, page  63). For plotting purposes we use a simpler interpolation
ss=linspace(s(1),s(end()),Nplot).T
# gen_single_drop.m:189
rr=interp1(s,r,ss,'pchip')
# gen_single_drop.m:190
zz=interp1(s,z,ss,'pchip')
# gen_single_drop.m:191
# plot the shape of the drop on the plotting grid
figure
hold('on')
scatter(dot(rneedle,rr.T),dot(rneedle,zz.T),'b')
plot(dot(rneedle,rr.T),dot(rneedle,zz.T),'b')
set(gca,'DataAspectRatio',concat([1,1,1]))