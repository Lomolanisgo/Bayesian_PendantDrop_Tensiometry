# Generated with SMOP  0.41
from libsmop import *
# fcbt.m

    
@function
def fcbt(diffmode=None,N=None,M=None,halfmark=None,*args,**kwargs):
    varargin = fcbt.varargin
    nargin = fcbt.nargin

    # chebyshev_polynomial = fcbt(N,M,halfmark) creates the first 30 chebyshev polynomials
# for N points and to M order, halfmark is one if only half the drop shape
# is fitted 0 to 1 instead of -1 to 1.
    
    # create the domain
    if exist('halfmark','var'):
        if (halfmark == 1) and strcmp(diffmode,'fd'):
            x=linspace(0,1,N).T
# fcbt.m:9
        else:
            if (halfmark == 1) and strcmp(diffmode,'cheb'):
                __,__,__,x=dif1D('cheb',0,1,N,nargout=4)
# fcbt.m:11
            else:
                __,__,__,x=dif1D('cheb',- 1,2,N,nargout=4)
# fcbt.m:13
                #x = linspace(-1,1,N)';
    else:
        x=linspace(- 1,1,N).T
# fcbt.m:17
    
    # reserve memory and create powers of x
    Y=zeros(N,31)
# fcbt.m:21
    x2=multiply(x,x)
# fcbt.m:23
    x3=multiply(x,x2)
# fcbt.m:24
    x4=multiply(x2,x2)
# fcbt.m:25
    x5=multiply(x2,x3)
# fcbt.m:26
    x6=multiply(x3,x3)
# fcbt.m:27
    x7=multiply(x3,x4)
# fcbt.m:28
    x8=multiply(x4,x4)
# fcbt.m:29
    x9=multiply(x4,x5)
# fcbt.m:30
    x10=multiply(x5,x5)
# fcbt.m:31
    x11=multiply(x6,x5)
# fcbt.m:32
    x12=multiply(x6,x6)
# fcbt.m:33
    x13=multiply(x6,x7)
# fcbt.m:34
    x14=multiply(x7,x7)
# fcbt.m:35
    x15=multiply(x7,x8)
# fcbt.m:36
    x16=multiply(x8,x8)
# fcbt.m:37
    x17=multiply(x8,x9)
# fcbt.m:38
    x18=multiply(x9,x9)
# fcbt.m:39
    x19=multiply(x9,x10)
# fcbt.m:40
    x20=multiply(x10,x10)
# fcbt.m:41
    x21=multiply(x10,x11)
# fcbt.m:42
    x22=multiply(x11,x11)
# fcbt.m:43
    x23=multiply(x11,x12)
# fcbt.m:44
    x24=multiply(x12,x12)
# fcbt.m:45
    x25=multiply(x12,x13)
# fcbt.m:46
    x26=multiply(x13,x13)
# fcbt.m:47
    x27=multiply(x13,x14)
# fcbt.m:48
    x28=multiply(x14,x14)
# fcbt.m:49
    x29=multiply(x14,x15)
# fcbt.m:50
    x30=multiply(x15,x15)
# fcbt.m:51
    # First 20 polynoms
    Y[arange(),arange(1,5)]=concat([ones(N,1),x,dot(2,x2) - 1,dot(4,x3) - dot(3,x),dot(8,(x4 - x2)) + 1])
# fcbt.m:54
    Y[arange(),arange(6,9)]=concat([dot(16,x5) - dot(20,x3) + dot(5,x),dot(32,x6) - dot(48,x4) + dot(18,x2) - 1,dot(64,x7) - dot(112,x5) + dot(56,x3) - dot(7,x),dot(128,x8) - dot(256,x6) + dot(160,x4) - dot(32,x2) + 1])
# fcbt.m:55
    Y[arange(),arange(10,11)]=concat([dot(9,x) - dot(120,x3) + dot(432,x5) - dot(576,x7) + dot(256,x9),- 1 + dot(50,x2) - dot(400,x4) + dot(1120,x6) - dot(1280,x8) + dot(512,x10)])
# fcbt.m:56
    Y[arange(),12]=dot(- 11,x) + dot(220,x3) - dot(1232,x5) + dot(2816,x7) - dot(2816,x9) + dot(1024,x11)
# fcbt.m:57
    Y[arange(),13]=1 - dot(72,x2) + dot(840,x4) - dot(3584,x6) + dot(6912,x8) - dot(6144,x10) + dot(2048,x12)
# fcbt.m:58
    Y[arange(),14]=dot(13,x) - dot(364,x3) + dot(2912,x5) - dot(9984,x7) + dot(16640,x9) - dot(13312,x11) + dot(4096,x13)
# fcbt.m:59
    Y[arange(),15]=- 1 + dot(98,x2) - dot(1568,x4) + dot(9408,x6) - dot(26880,x8) + dot(39424,x10) - dot(28672,x12) + dot(8192,x14)
# fcbt.m:60
    Y[arange(),16]=dot(- 15,x) + dot(560,x3) - dot(6048,x5) + dot(28800,x7) - dot(70400,x9) + dot(92160,x11) - dot(61440,x13) + dot(16384,x15)
# fcbt.m:61
    Y[arange(),17]=1 - dot(128,x2) + dot(2688,x4) - dot(21504,x6) + dot(84480,x8) - dot(180224,x10) + dot(212992,x12) - dot(131072,x14) + dot(32768,x16)
# fcbt.m:62
    Y[arange(),18]=dot(17,x) - dot(816,x3) + dot(11424,x5) - dot(71808,x7) + dot(239360,x9) - dot(452608,x11) + dot(487424,x13) - dot(278528,x15) + dot(65536,x17)
# fcbt.m:63
    Y[arange(),19]=- 1 + dot(162,x2) - dot(4320,x4) + dot(44352,x6) - dot(228096,x8) + dot(658944,x10) - dot(1118208,x12) + dot(1105920,x14) - dot(589824,x16) + dot(131072,x18)
# fcbt.m:64
    Y[arange(),20]=dot(- 19,x) + dot(1140,x3) - dot(20064,x5) + dot(160512,x7) - dot(695552,x9) + dot(1770496,x11) - dot(2723840,x13) + dot(2490368,x15) - dot(1245184,x17) + dot(262144,x19)
# fcbt.m:65
    Y[arange(),21]=1 - dot(200,x2) + dot(6600,x4) - dot(84480,x6) + dot(549120,x8) - dot(2050048,x10) + dot(4659200,x12) - dot(6553600,x14) + dot(5570560,x16) - dot(2621440,x18) + dot(524288,x20)
# fcbt.m:66
    # polynoms till 30
    Y[arange(),22]=dot(21,x) - dot(1540,x3) + dot(33264,x5) - dot(329472,x7) + dot(1793792,x9) - dot(5870592,x11) + dot(12042240,x13) - dot(15597568,x15) + dot(12386304,x17) - dot(5505024,x19) + dot(1048576,x21)
# fcbt.m:69
    Y[arange(),23]=- 1 + dot(242,x2) - dot(9680,x4) + dot(151008,x6) - dot(1208064,x8) + dot(5637632,x10) - dot(16400384,x12) + dot(30638080,x14) - dot(36765696,x16) + dot(27394048,x18) - dot(11534336,x20) + dot(2097152,x22)
# fcbt.m:70
    Y[arange(),24]=dot(- 23,x) + dot(2024,x3) - dot(52624,x5) + dot(631488,x7) - dot(4209920,x9) + dot(17145856,x11) - dot(44843008,x13) + dot(76873728,x15) - dot(85917696,x17) + dot(60293120,x19) - dot(24117248,x21) + dot(4194304,x23)
# fcbt.m:71
    Y[arange(),25]=1 - dot(288,x2) + dot(13728,x4) - dot(256256,x6) + dot(2471040,x8) - dot(14057472,x10) + dot(50692096,x12) - dot(120324096,x14) + dot(190513152,x16) - dot(199229440,x18) + dot(132120576,x20) - dot(50331648,x22) + dot(8388608,x24)
# fcbt.m:72
    Y[arange(),26]=dot(25,x) - dot(2600,x3) + dot(80080,x5) - dot(1144000,x7) + dot(9152000,x9) - dot(45260800,x11) + dot(146227200,x13) - dot(317521920,x15) + dot(466944000,x17) - dot(458752000,x19) + dot(288358400,x21) - dot(104857600,x23) + dot(16777216,x25)
# fcbt.m:73
    Y[arange(),27]=- 1 + dot(338,x2) - dot(18928,x4) + dot(416416,x6) - dot(4759040,x8) + dot(32361472,x10) - dot(141213696,x12) + dot(412778496,x14) - dot(825556992,x16) + dot(1133117440,x18) - dot(1049624576,x20) + dot(627048448,x22) - dot(218103808,x24) + dot(33554432,x26)
# fcbt.m:74
    Y[arange(),28]=dot(- 27,x) + dot(3276,x3) - dot(117936,x5) + dot(1976832,x7) - dot(18670080,x9) + dot(109983744,x11) - dot(428654592,x13) + dot(1143078912,x15) - dot(2118057984,x17) + dot(2724986880,x19) - dot(2387607552,x21) + dot(1358954496,x23) - dot(452984832,x25) + dot(67108864,x27)
# fcbt.m:75
    Y[arange(),29]=1 - dot(392,x2) + dot(25480,x4) - dot(652288,x6) + dot(8712704,x8) - dot(69701632,x10) + dot(361181184,x12) - dot(1270087680,x14) + dot(3111714816,x16) - dot(5369233408,x18) + dot(6499598336,x20) - dot(5402263552,x22) + dot(2936012800,x24) - dot(939524096,x26) + dot(134217728,x28)
# fcbt.m:76
    Y[arange(),30]=dot(29,x) - dot(4060,x3) + dot(168896,x5) - dot(3281408,x7) + dot(36095488,x9) - dot(249387008,x11) + dot(1151016960,x13) - dot(3683254272,x15) + dot(8341487616,x17) - dot(13463453696,x19) + dot(15386804224,x21) - dot(12163481600,x23) + dot(6325010432,x25) - dot(1946157056,x27) + dot(268435456,x29)
# fcbt.m:77
    Y[arange(),31]=- 1 + dot(450,x2) - dot(33600,x4) + dot(990080,x6) - dot(15275520,x8) + dot(141892608,x10) - dot(859955200,x12) + dot(3572121600,x14) - dot(10478223360,x16) + dot(22052208640,x18) - dot(33426505728,x20) + dot(36175872000,x22) - dot(27262976000,x24) + dot(13589544960,x26) - dot(4026531840,x28) + dot(536870912,x30)
# fcbt.m:78
    if M <= 31:
        Y=Y(arange(),arange(1,M))
# fcbt.m:81
    
    # For larger polynoms use MATLAB (slow)
    if M > 31:
        for k in arange(31,M).reshape(-1):
            Y[arange(),k + 1]=chebyshevT(k,x)
# fcbt.m:87
    
    return Y
    
if __name__ == '__main__':
    pass
    