# Generated with SMOP  0.41
from libsmop import *
# fit_circle_through_3_points.m

    
@function
def fit_circle_through_3_points(ABC=None,*args,**kwargs):
    varargin = fit_circle_through_3_points.varargin
    nargin = fit_circle_through_3_points.nargin

    # FIT_CIRCLE_THROUGH_3_POINTS
    # Mathematical background is provided in http://www.regentsprep.org/regents/math/geometry/gcg6/RCir.htm
    
    # Input:
    
    #   ABC is a [3 x 2n] array. Each two columns represent a set of three points which lie on
    #       a circle. Example: [-1 2;2 5;1 1] represents the set of points (-1,2), (2,5) and (1,1) in Cartesian
    #       (x,y) coordinates.
    
    # Outputs:
    
    #   R     is a [1 x n] array of circle radii corresponding to each set of three points.
    #   xcyc  is an [2 x n] array of of the centers of the circles, where each column is [xc_i;yc_i] where i
    #         corresponds to the {A,B,C} set of points in the block [3 x 2i-1:2i] of ABC
    
    # Author: Danylo Malyuta.
    # Version: v1.0 (June 2016)
    # ----------------------------------------------------------------------------------------------------------
    # Each set of points {A,B,C} lies on a circle. Question: what is the circles radius and center?
    # A: point with coordinates (x1,y1)
    # B: point with coordinates (x2,y2)
    # C: point with coordinates (x3,y3)
    # ============= Find the slopes of the chord A<-->B (mr) and of the chord B<-->C (mt)
    #   mt = (y3-y2)/(x3-x2)
    #   mr = (y2-y1)/(x2-x1)
    # /// Begin by generalizing xi and yi to arrays of individual xi and yi for each {A,B,C} set of points provided in ABC array
    x1=ABC(1,arange(1,end(),2))
# fit_circle_through_3_points.m:28
    x2=ABC(2,arange(1,end(),2))
# fit_circle_through_3_points.m:29
    x3=ABC(3,arange(1,end(),2))
# fit_circle_through_3_points.m:30
    y1=ABC(1,arange(2,end(),2))
# fit_circle_through_3_points.m:31
    y2=ABC(2,arange(2,end(),2))
# fit_circle_through_3_points.m:32
    y3=ABC(3,arange(2,end(),2))
# fit_circle_through_3_points.m:33
    
    mr=(y2 - y1) / (x2 - x1)
# fit_circle_through_3_points.m:35
    mt=(y3 - y2) / (x3 - x2)
# fit_circle_through_3_points.m:36
    
    #   (1) First chord is vertical       ==> mr==Inf
    #   (2) Second chord is vertical      ==> mt==Inf
    #   (3) Points are collinear          ==> mt==mr (NB: NaN==NaN here)
    #   (4) Two or more points coincident ==> mr==NaN || mt==NaN
    # Resolve these failure modes case-by-case.
    idf1=isinf(mr)
# fit_circle_through_3_points.m:43
    
    idf2=isinf(mt)
# fit_circle_through_3_points.m:44
    
    idf34=logical_or(logical_or(isequaln(mr,mt),isnan(mr)),isnan(mt))
# fit_circle_through_3_points.m:45
    
    # ============= Compute xc, the circle center x-coordinate
    xcyc=(multiply(multiply(mr,mt),(y3 - y1)) + multiply(mr,(x2 + x3)) - multiply(mt,(x1 + x2))) / (dot(2,(mr - mt)))
# fit_circle_through_3_points.m:47
    xcyc[idf1]=(multiply(mt(idf1),(y3(idf1) - y1(idf1))) + (x2(idf1) + x3(idf1))) / 2
# fit_circle_through_3_points.m:48
    
    xcyc[idf2]=((x1(idf2) + x2(idf2)) - multiply(mr(idf2),(y3(idf2) - y1(idf2)))) / 2
# fit_circle_through_3_points.m:49
    
    xcyc[idf34]=NaN
# fit_circle_through_3_points.m:50
    
    # ============= Compute yc, the circle center y-coordinate
    xcyc[2,arange()]=multiply(- 1.0 / mr,(xcyc - (x1 + x2) / 2)) + (y1 + y2) / 2
# fit_circle_through_3_points.m:52
    idmr0=mr == 0
# fit_circle_through_3_points.m:53
    xcyc[2,idmr0]=multiply(- 1.0 / mt(idmr0),(xcyc(idmr0) - (x2(idmr0) + x3(idmr0)) / 2)) + (y2(idmr0) + y3(idmr0)) / 2
# fit_circle_through_3_points.m:54
    xcyc[2,idf34]=NaN
# fit_circle_through_3_points.m:55
    
    # ============= Compute the circle radius
    R=sqrt((xcyc(1,arange()) - x1) ** 2 + (xcyc(2,arange()) - y1) ** 2)
# fit_circle_through_3_points.m:57
    R[idf34]=Inf
# fit_circle_through_3_points.m:58
    
    return R,xcyc
    
if __name__ == '__main__':
    pass
    