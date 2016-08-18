# -*- coding: utf-8 -*-
"""
classes
"""


class CartVector():
    x = None
    y = None
    z = None
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z
    
    def x(self, x):
        self.x = x

    
    def y(self, y):
        self.y = y

    
    def z(self, z):
        self.z = z
    
class CylVector():
    rho = None
    zeta = None
    h = None
    def __init__(self, rho=None, zeta=None, h=None):
        self.rho = rho
        self.zeta = zeta
        self.h = h
    
    def rho(self, rho):
        self.rho = rho

    
    def zeta(self, zeta):
        self.zeta = zeta

    
    def h(self, h):
        self.h = h

class SphereVector():
    r = None
    theta = None
    phi = None
    def __init__(self, r=None, theta=None, phi=None):
        self.r = r
        self.theta = theta
        self.phi = phi

    
    def r(self, r):
        self.r = r

    
    def theta(self, theta):
        self.theta = theta

    
    def phi(self, phi):
        self.phi = phi

class filenames():
    parfile = None
    syspars = None
    lightcurves = None
    def __init__(self, parfile=None, syspars=None, lightcurves=None):
        self.parfile = parfile
        self.syspars = syspars
        self.lightcurves = lightcurves

    
    def parfile(self, parfile):
        self.parfile = parfile

    
    def syspars(self, syspars):
        self.syspars = syspars

    
    def lightcurves(self, lightcurves):
        self.lightcurves = lightcurves

class flowcontrol():
    star1 = None
    star2 = None
    star2spots = None
    disk = None
    diskrim = None
    disktorus = None
    innerdisk = None
    diskspots = None
    adc = None
    thirdlight = None
    irradiation = None
    diagnostics = None
    diagnosephase = None
    diagnoseband = None
    diagnoseindex = None    
    def __init__(self, star1=None, star2=None, star2spots=None, disk=None, diskrim=None, disktorus=None, innerdisk=None,
            diskspots=None, adc=None, thirdlight=None, irradiation=None, diagnostics=None, diagnosephase=None, diagnoseband=None, diagnoseindex=None):
                      self.star1 = star1
                      self.star2 = star2
                      self.star2spots = star2spots
                      self.disk = disk
                      self.diskrim = diskrim
                      self.disktorus = disktorus
                      self.innerdisk = innerdisk
                      self.diskspots = diskspots
                      self.adc = adc
                      self.thirdlight = thirdlight
                      self.irradiation = irradiation
                      self.diagnostics = diagnostics
                      self.diagnosephase = diagnosephase
                      self.diagnoseband = diagnoseband
                      self.diagnoseindex = diagnoseindex

    
    def star1(self, star1):
        self.star1 = star1
    
    
    def star2(self, star2):
        self.star2 = star2

    
    def star2spots(self, star2spots):
        self.star2spots = star2spots

    
    def disk(self, disk):
        self.disk = disk

    
    def diskrim(self, diskrim):
        self.diskrim = diskrim

    
    def disktorus(self, disktorus):
        self.disktorus = disktorus

    
    def innerdisk(self, innerdisk):
        self.innerdisk = innerdisk

    
    def diskspots(self, diskspots):
        self.diskspots = diskspots

    
    def adc(self, adc):
        self.adc = adc

    
    def thirdlight(self, thirdlight):
        self.thirdlight = thirdlight

    
    def irradiation(self, irradiation):
        self.irradiation = irradiation

    
    def diagnostics(self, diagnostics):
        self.diagnostics = diagnostics

    
    def diagnosephase(self, diagnosephase):
        self.diagnosephase = diagnosephase

    
    def diagnoseband(self, diagnoseband):
        self.diagnoseband = diagnoseband

    
    def diagnoseindex(self, diagnoseindex):
        self.diagnoseindex = diagnoseindex

class orbitparams():
    phasemin = None
    deltaphase = None
    maxpindex = None
    phaseoffset = None
    nbands = None
    filtermax = None
    minlambda = [0 for i in range(22)]
    maxlambda = [0 for i in range(22)]
    normalize = [0 for i in range(21)]
    normfilter = [0 for i in range(21)]
    normMinlambda = None
    normMaxlambda = None
    normvalue = None
    def __init__(self, phasemin=None, phasemax=None, deltaphase=None, maxpindex=None,
                 phaseoffset=None, nbands=None, filtermax = None, minlambda = None,
                 maxlambda=None, normalize=None, normfilter=None, normMinlambda=None,
                 normMaxlambda=None, normvalue=None):
                     self.phasemin = phasemin
                     self.phasemax = phasemax
                     self.deltaphase = deltaphase
                     self.maxpindex = maxpindex
                     self.phaseoffset = phaseoffset
                     self.nbands = nbands
                     self.filtermax = filtermax
                     self.minlambda = minlambda
                     self.maxlambda = maxlambda
                     self.normalize = normalize
                     self.normfilter = normfilter
                     self.normMinlambda = normMinlambda
                     self.normMaxlambda = normMaxlambda
                     self.normvalue = normvalue
                
    
    def phasemin(self, phasemin):
        self.phasemin = phasemin

    
    def phasemax(self, phasemax):
        self.phasemax = phasemax

    
    def deltaphase(self, deltaphase):
        self.deltaphase = deltaphase

    
    def maxpindex(self, maxpindex):
        self.maxpindex = maxpindex

    
    def phaseoffset(self, phaseoffset):
        self.phaseoffset = phaseoffset

    
    def nbands(self, nbands):
        self.nbands = nbands

    
    def filtermax(self, filtermax):
        self.filtermax = filtermax

    
    def minlambda(self, minlambda):
        self.minlambda = minlambda

    
    def maxlambda(self, maxlambda):
        self.maxlambda = maxlambda

    
    def normalize(self, normalize):
        self.normalize = normalize

    
    def normfilter(self, normfilter):
        self.normfilter = normfilter

    
    def normMinlambda(self, normMinlambda):
        self.normMinlambda = normMinlambda

    
    def normMaxlambda(self, normMaxlambda):
        self.normMaxlambda = normMaxlambda

    
    def normvalue(self, normvalue):
        self.normvalue = normvalue

class systemparams():
    p = None
    omega = None
    k2 = None
    q = None
    i = None
    a = None
    zcm = None
    M1 = None
    M2 = None
    VL1 = None
    rL1 = None
    MeanLobe1Radius = None
    MeanLobe2Radius = None    
    def __init__(self, p=None, omega=None, K2=None, q=None, i=None, a=None, zcm=None,
                 M1=None, M2=None, rL1=None, VL1=None, MeanLobe1Radius=None, MeanLobe2Radius=None):
                     self.p = p
                     self.omega = omega
                     self.K2 = K2
                     self.q = q
                     self.i = i
                     self.a = a
                     self.zcm = zcm
                     self.M1 = M1
                     self.M2 = M2
                     self.rL1 = rL1
                     self.VL1 = VL1
                     self.MeanLobe1Radius = MeanLobe1Radius
                     self.MeanLobe2Radius = MeanLobe2Radius

    
    def p(self, p):
        self.p = p

    
    def omega(self, omega):
        self.omega = omega

    
    def K2(self, K2):
        self.K2 = K2
        
    
    def q(self, q):
        self.q = q

    
    def i(self, i):
        self.i = i

    
    def a(self, a):
        self.a = a

    
    def zcm(self, zcm):
        self.zcm = zcm

    
    def M1(self, M1):
        self.M1 = M1

    
    def M2(self, M2):
        self.M2 = M2

    
    def rL1(self, rL1):
        self.rL1 = rL1
        
    
    def VL1(self, VL1):
        self.VL1 = VL1

    
    def MeanLobe1Radius(self, r):
        self.MeanLobe1Radius = r

    
    def MeanLobe2Radius(self, r):
        self.MeanLobe2Radius = r

class star2spotparams():
    nspots = None
    theta = [0 for i in range(21)]
    phi = [0 for i in range(21)]
    radius = [0 for i in range(21)]
    SpotoverStarT = [0 for i in range(21)]
    def __init__(self, nspots=None, theta=None, phi=None, radius=None, SpotToverStarT=None):
        self.nspots = nspots
        self.theta = theta
        self.phi = phi
        self.radius = radius
        self.SpotToverStarT = SpotToverStarT
    
    def nspots(self, nspots):
        self.nspots = nspots

    
    def theta(self, theta):
        self.theta = theta

    
    def phi(self, phi):
        self.phi = phi

    
    def radius(self, radius):
        self.radius = radius

    
    def SpotToverStarT(self, s):
        self.SpotToverStarT = s

class wholediskpars():
    targetNtiles = None
    e = None
    zetazero = None
    albedo = None
    L = None
    TopTmax = None
    TopTmin = None
    def __init__(self, targetNtiles=None, Ntiles=None, e=None, zetazero=None, albedo=None,
                 L=None, TopTmax=None, TopTmin=None):
                     self.targetNtiles = targetNtiles
                     self.Ntiles = Ntiles
                     self.e = e
                     self.zetazero = zetazero
                     self.albedo = albedo
                     self.L = L
                     self.TopTmax = TopTmax
                     self.TopTmin = TopTmin

    
    def targetNtiles(self, t):
        self.targetNtiles = t

    
    def Ntiles(self, n):
        self.Ntiles = n

    
    def e(self, e):
        self.e = e

    
    def zetazero(self, z):
        self.zetazero = z

    
    def albedo(self, albedo):
        self.albedo = albedo

    
    def L(self, L):
        self.L = L

    
    def TopTmax(self, t):
        self.TopTmax = t

     
    def TopTmin(self, t):
        self.Toptmin = t

class diskedgepars():
    T = None
    Tspot = None
    ZetaMid = None
    ZetaWidth = None
    def __init__(self, T=None, Tspot=None, ZetaMid=None, ZetaWidth=None):
        self.T = T
        self.Tspot = Tspot
        self.ZetaMid = ZetaMid
        self.ZetaWidth = ZetaWidth
 
    
    def T(self, T):
        self.T = T

    
    def Tspot(self, Tspot):
        self.Tspot = Tspot

    
    def ZetaMid(self, z):
        self.ZetaMid = z

    
    def ZetaWidth(self, ZetaWidth):
        self.ZetaWidth = ZetaWidth

class diskrimpars():
    Type = None
    awidth = None
    Hmax = None
    Hmin = None
    ZetaHmax = None
    Tmax = None
    Tmin = None
    ZetaTmax = None
    points = None
    PointZeta = [0 for i in range(102)]
    PointH = [0 for i in range(102)]
    PointT = [0 for i in range(102)]
    def __init__(self, Type=None, awidth=None, Hmax=None, Hmin=None, ZetaHmax=None,
                 Tmax=None, Tmin=None, ZetaTmax=None, points=None, PointZeta=None,
                 PointH=None, PointT=None):
                     self.Type = Type
                     self.awidth = awidth
                     self.Hmax = Hmax
                     self.Hmin = Hmin
                     self.ZetaHmax = ZetaHmax
                     self.Tmax = Tmax
                     self.Tmin = Tmin
                     self.ZetaTmax = ZetaTmax
                     self.points = points
                     self.PointZeta = PointZeta
                     self.PointH = PointH
                     self.PointT = PointT
    
    def Type(self, Type):
        self.Type = Type

    
    def awidth(self, awidth):
        self.awidth = awidth

    
    def Hmax(self, Hmax):
        self.Hmax = Hmax

    
    def Hmin(self, Hmin):
        self.Hmin = Hmin

    
    def ZetaHmax(self, z):
        self.ZetaHmax = z

    
    def Tmax(self, t):
        self.Tmax = t

    
    def Tmin(self, t):
        self.Tmin = t

    
    def ZetaTmax(self, z):
        self.ZetaTmax = z

    
    def points(self, p):
        self.points = p

    
    def PointZeta(self, p):
        self.PointZeta = p

    
    def PointH(self, p):
        self.PointH = p

    
    def PointT(self, p):
        self.PointT = p

class disktorusparams():
    Type = None
    awidth = None
    Hmax = None
    Hmin = None
    ZetaHmax = None
    Tmax = None
    Tmin = None
    ZetaTmax = None
    points = None
    PointZeta = [0 for i in range(102)]
    PointH = [0 for i in range(102)]
    PointT = [0 for i in range(102)]
    azero = None
    def __init(self, Type=None, awidth=None, Hmax=None, Hmin=None, ZetaHmax=None,
                 Tmax=None, Tmin=None, ZetaTmax=None, points=None, PointZeta=None,
                 PointH=None, PointT=None, azero=None):
                     self.Type = Type
                     self.awidth = awidth
                     self.Hmax = Hmax
                     self.Hmin = Hmin
                     self.ZetaHmax = ZetaHmax
                     self.Tmax = Tmax
                     self.Tmin = Tmin
                     self.ZetaTmax = ZetaTmax
                     self.points = points
                     self.PointZeta = PointZeta
                     self.PointH = PointH
                     self.PointT = PointT
                     self.azero = azero
    
    def Type(self, Type):
        self.Type = Type

    
    def awidth(self, awidth):
        self.awidth = awidth

    
    def Hmax(self, Hmax):
        self.Hmax = Hmax

    
    def Hmin(self, Hmin):
        self.Hmin = Hmin

    
    def ZetaHmax(self, z):
        self.ZetaHmax = z

    
    def Tmax(self, t):
        self.Tmax = t

    
    def Tmin(self, t):
        self.Tmin = t

    
    def ZetaTmax(self, z):
        self.ZetaTmax = z

    
    def points(self, p):
        self.points = p

    
    def PointZeta(self, p):
        self.PointZeta = p

    
    def PointH(self, p):
        self.PointH = p

    
    def PointT(self, p):
        self.PointT = p

    
    def azero(self, azero):
        self.azero = azero

class diskspotpars():
    nspots = None
    zetamin = [0 for i in range(21)]
    zetamax = [0 for i in range(21)]
    amin = [0 for i in range(21)]
    amax = [0 for i in range(21)]
    spotToverT = [0 for i in range(21)]
    def __init__(self, nspots=None, zetamin=None, zetamax=None, amin=None, amax=None,
                 spotToverT=None):
                     self.nspots = nspots
                     self.zetamin = zetamin
                     self.zetamax = zetamax
                     self.amin = amin
                     self.amax = amax
                     self.spotToverT = spotToverT

    
    def nspots(self, n):
        self.nspots = n

    
    def zetamin(self, z):
        self.zetamin = z

    
    def zetamax(self, z):
        self.zetamax = z

    
    def amin(self, amin):
        self.amin = amin

    
    def amax(self, amax):
        self.amax = amax

    
    def spotToverT(self, s):
        self.spotToverT = s

class innerdiskpars():
    T = None
    L = None
    sigmaT4 = None
    radius = None
    def __init__(self, T=None, L=None, sigmaT4=None, radius=None):
        self.T = T
        self.L = L
        self.sigmaT4 = sigmaT4
        self.radius = radius

    
    def T(self, T):
        self.T = T

    
    def L(self, L):
        self.L = L

    
    def sigmaT4(self, s):
        self.sigmaT4 = s

    
    def radius(self, r):
        self.radius = r

class adcpars():
    L = None
    height = None
    def __init__(self, L=None, height=None):
        self.L = L
        self.height = height

    
    def L(self, L):
        self.L = L

    
    def height(self, height):
        self.height = height

class thirdlightparams():
    orbphase = None
    nbands = None
    Filter = None
    minlambda = [0 for i in range(22)]
    maxlambda = [0 for i in range(22)]
    fraction = [0 for i in range(22)]
    addFlux = [0 for i in range(22)]
    def __init__(self, orbphase=None, nbands=None, Filter=None, minlambda=None,
                 maxlambda=None, fraction=None, addFlux=None):
                     self.orbphase = orbphase
                     self.nbands = nbands
                     self.Filter = Filter
                     self.minlambda = minlambda
                     self.maxlambda = maxlambda
                     self.fraction = fraction
                     self.addFlux = addFlux

    
    def orbphase(self, o):
        self.orbphase = o

    
    def nbands(self, n):
        self.nbands = n

    
    def Filter(self, f):
        self.Filter = f

    
    def minlambda(self, m):
        self.minlambda = m

    
    def maxlambda(self, m):
        self.maxlambda = m

    
    def fraction(self, fraction):
        self.fraction = fraction

    
    def addFlux(self, addFlux):
        self.addFlux = addFlux

class XYGrid():
    Nxtiles = None
    Nztiles = None
    deltax = None
    deltaz = None
    deltal = None
    xmin = None
    xmax = None
    ymin = None
    ymax = None
    zmin = None
    zmax = None
    Topy = [[0 for i in range(402)] for j in range(402)]
    Bottomy = [[0 for i in range(402)] for j in range(402)]
    def __init__(self, Nxtiles=None, Nztiles=None, deltax=None, deltaz=None, deltal=None,
                 xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None,
                 Topy=[], Bottomy=[]):
                     self.Nxtiles = Nxtiles
                     self.Nztiles = Nztiles
                     self.deltax = deltax
                     self.deltaz = deltaz
                     self.deltal = deltal
                     self.xmin = xmin
                     self.xmax = xmax
                     self.ymin = ymin
                     self.ymax = ymax
                     self.zmin = zmin
                     self.zmax = zmax
                     self.Topy = Topy
                     self.Bottomy = Bottomy
    
    
    def Nxtiles(self, x):
        self.Nxtiles = x

    
    def Nztiles(self, z):
        self.Nztiles = z

    
    def deltax(self, x):
        self.deltax = x

    
    def deltaz(self, z):
        self.deltaz = z

    
    def deltal(self, l):
        self.deltal = l

    
    def xmin(self, x):
        self.xmin = x

    
    def xmax(self, x):
        self.xmax = x

    
    def ymin(self, y):
        self.ymin = y

    
    def ymax(self, y):
        self.ymax = y

    
    def zmin(self, z):
        self.zmin = z

    def zmax(self, z):
        self.zmax = z

    
    def Topy(self, y):
        self.Topy = y

    
    def Bottomy(self, y):
        self.Bottomy = y

class dataparams():
    nbands = None
    filename = None
    Filter = None
    minlambda = [0 for i in range(22)]
    maxlambda = [0 for i in range(22)]
    npoints = [0 for i in range(22)]
    phase = [[0 for i in range(22)] for j in range(1002)]
    flux = [[0 for i in range(22)] for j in range(1002)]
    standdev = [[0 for i in range(22)] for j in range(1002)]
    chisquare = None
    def __init__(self, nbands=None, filename=None, Filter=None, minlambda=None,
                 maxlambda=None, npoints=None, phase=None, flux=None, standdev=None, chisquare=None):
                     self.nbands = nbands
                     self.filename = filename
                     self.Filter = Filter
                     self.minlambda = minlambda
                     self.maxlambda = maxlambda
                     self.npoints = npoints
                     self.phase = phase
                     self.flux = flux
                     self.standdev = standdev
                     self.chisquare = chisquare
    
    def nbands(self, b):
        self.nbands = b

    
    def filename(self, filename):
        self.filename = filename

    
    def Filter(self, Filter):
        self.Filter = Filter

    
    def minlambda(self, minlambda):
        self.minlambda = minlambda

    
    def maxlambda(self, maxlambda):
        self.maxlambda = maxlambda

    
    def npointa(self, npoints):
        self.npoints = npoints

    
    def phase(self, phase):
        self.phase = phase

    
    def flux(self, flux):
        self.flux = flux

    
    def standdev(self, standdev):
        self.standdev = standdev

    
    def chisquare(self, chisquare):
        self.chisquare = chisquare

class globalvar:
    verbose = None
    maxGDindex = None
    GDT = [0 for i in range(31)]
    fourbeta = [0 for i in range(31)]
    LDfilterName = None
    maxLDgindex = None
    maxLDTindex = None
    maxLDfilterindex = None
    LDlogg = [0 for i in range(21)]
    LDT = [0 for i in range(101)]
    LDtable = [[[[0 for i in range(21)] for j in range(101)] for k in range(13)] for l in range(6)]
    IperpfilterName = None
    Iperplogg = [0 for i in range(21)]
    IperpT = [0 for i in range(101)]
    maxIperpgindex = None
    maxIperpTindex = None
    maxIperpfilterindex = None
    Iperptable = [[[0 for i in range(21)] for j in range(101)] for k in range(13)]
    IBBfilterName = None
    maxIBBTindex = None
    maxIBBfilterindex = None
    IBBTmin  = None
    IBBTmax = None
    IBBdeltaT = None
    IBBT = [0 for i in range(15002)]
    IBBtable = [[0 for i in range(15002)] for j in range(13)]
    maxBBzetaindex = None
    deltaBBzeta = None
    BBzetamax = None
    ZBBzeta = [0 for i in range(40001)]
    T1I = [0 for i in range(22)]
    a = CartVector()
    b = CylVector()
    c = SphereVector()
    T2normCart = [a for i in range(40507)]
    T2gradV = [c for i in range(40507)]
    T2normSphere = [c for i in range(40507)]
    T2r = [0 for i in range(40507)]
    T2theta = [0 for i in range(40507)]
    T2phi = [0 for i in range(40507)] 
    T2x = [0 for i in range(40507)]
    T2y = [0 for i in range(40507)]
    T2z = [0 for i in range(40507)]
    T2g = [0 for i in range(40507)]
    T2logg = [0 for i in range(40507)]
    T2dS = [0 for i in range(40507)]
    T2T = [0 for i in range(40507)]
    T2I = [[0 for i in range(22)] for j in range(40507)]
    TDisknormCyl = [b for i in range(40507)]
    TDisknormCart = [a for i in range(40507)]
    TDiskRho = [0 for i in range(40507)]
    TDiskZeta = [0 for i in range(40507)]
    TDiskH = [0 for i in range(40507)] 
    TDiskx = [0 for i in range(40507)]
    TDisky = [0 for i in range(40507)]
    TDiskz = [0 for i in range(40507)] 
    TDiskdS = [0 for i in range(40507)]
    TDiskT = [0 for i in range(40507)]
    TDiskT4 = [0 for i in range(40507)]
    TDiska = [0 for i in range(40507)]
    TDiskI = [[0 for i in range(22)] for j in range(40507)]
    LCphase = [0 for i in range(502)]
    LCflux = [[0 for i in range(22)] for j in range(502)]
    NormLC = [[0 for i in range(22)] for j in range(502)]
