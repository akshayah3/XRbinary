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
    @set
    def x(self, x):
        self.x = x

    @set
    def y(self, y):
        self.y = y

    @set
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
    @set
    def rho(self, rho):
        self.rho = rho

    @set
    def zeta(self, zeta):
        self.zeta = zeta

    @set
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

    @set
    def r(self, r):
        self.r = r

    @set
    def theta(self, theta):
        self.theta = theta

    @set
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

    @set
    def parfile(self, parfile):
        self.parfile = parfile

    @set
    def syspars(self, syspars):
        self.syspars = syspars

    @set
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

    @set
    def star1(self, star1):
        self.star1 = star1
    
    @set
    def star2(self, star2):
        self.star2 = star2

    @set
    def star2spots(self, star2spots):
        self.star2spots = star2spots

    @set
    def disk(self, disk):
        self.disk = disk

    @set
    def diskrim(self, diskrim):
        self.diskrim = diskrim

    @set
    def disktorus(self, disktorus):
        self.disktorus = disktorus

    @set
    def innerdisk(self, innerdisk):
        self.innerdisk = innerdisk

    @set
    def diskspots(self, diskspots):
        self.diskspots = diskspots

    @set
    def adc(self, adc):
        self.adc = adc

    @set
    def thirdlight(self, thirdlight):
        self.thirdlight = thirdlight

    @set
    def irradiation(self, irradiation):
        self.irradiation = irradiation

    @set
    def diagnostics(self, diagnostics):
        self.diagnostics = diagnostics

    @set
    def diagnosephase(self, diagnosephase):
        self.diagnosephase = diagnosephase

    @set
    def diagnoseband(self, diagnoseband):
        self.diagnoseband = diagnoseband

    @set
    def diagnoseindex(self, diagnoseindex):
        self.diagnoseindex = diagnoseindex

class orbitparams():
    phasemin = None
    deltaphase = None
    maxpindex = None
    phaseoffset = None
    nbands = None
    filtermax = None
    minlambda = []
    maxlambda = []
    normalize = None
    normfilter = None
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
                
    @set
    def phasemin(self, phasemin):
        self.phasemin = phasemin

    @set
    def phasemax(self, phasemax):
        self.phasemax = phasemax

    @set
    def deltaphase(self, deltaphase):
        self.deltaphase = deltaphase

    @set
    def maxpindex(self, maxpindex):
        self.maxpindex = maxpindex

    @set
    def phaseoffset(self, phaseoffset):
        self.phaseoffset = phaseoffset

    @set
    def nbands(self, nbands):
        self.nbands = nbands

    @set
    def filtermax(self, filtermax):
        self.filtermax = filtermax

    @set
    def minlambda(self, minlambda):
        self.minlambda = minlambda

    @set
    def maxlambda(self, maxlambda):
        self.maxlambda = maxlambda

    @set
    def normalize(self, normalize):
        self.normalize = normalize

    @set
    def normfilter(self, normfilter):
        self.normfilter = normfilter

    @set
    def normMinlambda(self, normMinlambda):
        self.normMinlambda = normMinlambda

    @set
    def normMaxlambda(self, normMaxlambda):
        self.normMaxlambda = normMaxlambda

    @set
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

    @set
    def p(self, p):
        self.p = p

    @set
    def omega(self, omega):
        self.omega = omega

    @set
    def K2(self, K2):
        self.K2 = K2
        
    @set
    def q(self, q):
        self.q = q

    @set
    def i(self, i):
        self.i = i

    @set
    def a(self, a):
        self.a = a

    @set
    def zcm(self, zcm):
        self.zcm = zcm

    @set
    def M1(self, M1):
        self.M1 = M1

    @set
    def M2(self, M2):
        self.M2 = M2

    @set
    def rL1(self, rL1):
        self.rL1 = rL1
        
    @set
    def VL1(self, VL1):
        self.VL1 = VL1

    @set
    def MeanLobe1Radius(self, r):
        self.MeanLobe1Radius = r

    @set
    def MeanLobe2Radius(self, r):
        self.MeanLobe2Radius = r

class star2spotparams():
    nspots = None
    theta = []
    phi = []
    radius = []
    SpotoverStarT = None
    def __init__(self, nspots=None, theta=None, phi=None, radius=None, SpotToverStarT=None):
        self.nspots = nspots
        self.theta = theta
        self.phi = phi
        self.radius = radius
        self.SpotToverStarT = SpotToverStarT
    @set
    def nspots(self, nspots):
        self.nspots = nspots

    @set
    def theta(self, theta):
        self.theta = theta

    @set
    def phi(self, phi):
        self.phi = phi

    @set
    def radius(self, radius):
        self.radius = radius

    @set
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

    @set
    def targetNtiles(self, t):
        self.targetNtiles = t

    @set
    def Ntiles(self, n):
        self.Ntiles = n

    @set
    def e(self, e):
        self.e = e

    @set
    def zetazero(self, z):
        self.zetazero = z

    @set
    def albedo(self, albedo):
        self.albedo = albedo

    @set
    def L(self, L):
        self.L = L

    @set
    def TopTmax(self, t):
        self.TopTmax = t

    @set 
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
 
    @set
    def T(self, T):
        self.T = T

    @set
    def Tspot(self, Tspot):
        self.Tspot = Tspot

    @set
    def ZetaMid(self, z):
        self.ZetaMid = z

    @set
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
    PointZeta = []
    PointH = []
    PointT = []
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
    @set
    def Type(self, Type):
        self.Type = Type

    @set
    def awidth(self, awidth):
        self.awidth = awidth

    @set
    def Hmax(self, Hmax):
        self.Hmax = Hmax

    @set
    def Hmin(self, Hmin):
        self.Hmin = Hmin

    @set
    def ZetaHmax(self, z):
        self.ZetaHmax = z

    @set
    def Tmax(self, t):
        self.Tmax = t

    @set
    def Tmin(self, t):
        self.Tmin = t

    @set
    def ZetaTmax(self, z):
        self.ZetaTmax = z

    @set
    def points(self, p):
        self.points = p

    @set
    def PointZeta(self, p):
        self.PointZeta = p

    @set
    def PointH(self, p):
        self.PointH = p

    @set
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
    PointZeta = []
    PointH = []
    PointT = []
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
    @set
    def Type(self, Type):
        self.Type = Type

    @set
    def awidth(self, awidth):
        self.awidth = awidth

    @set
    def Hmax(self, Hmax):
        self.Hmax = Hmax

    @set
    def Hmin(self, Hmin):
        self.Hmin = Hmin

    @set
    def ZetaHmax(self, z):
        self.ZetaHmax = z

    @set
    def Tmax(self, t):
        self.Tmax = t

    @set
    def Tmin(self, t):
        self.Tmin = t

    @set
    def ZetaTmax(self, z):
        self.ZetaTmax = z

    @set
    def points(self, p):
        self.points = p

    @set
    def PointZeta(self, p):
        self.PointZeta = p

    @set
    def PointH(self, p):
        self.PointH = p

    @set
    def PointT(self, p):
        self.PointT = p

    @set
    def azero(self, azero):
        self.azero = azero

class diskspotpars():
    nspots = None
    zetamin = []
    zetamax = []
    amin = []
    amax = []
    spotToverT = []
    def __init__(self, nspots=None, zetamin=None, zetamax=None, amin=None, amax=None,
                 spotToverT):
                     self.nspots = nspots
                     self.zetamin = zetamin
                     self.zetamax = zetamax
                     self.amin = amin
                     self.amax = amax
                     self.spotToverT = spotToverT

    @set
    def nspots(self, n):
        self.nspots = n

    @set
    def zetamin(self, z):
        self.zetamin = z

    @set
    def zetamax(self, z):
        self.zetamax = z

    @set
    def amin(self, amin):
        self.amin = amin

    @set
    def amax(self, amax):
        self.amax = amax

    @set
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

    @set
    def T(self, T):
        self.T = T

    @set
    def L(self, L):
        self.L = L

    @set
    def sigmaT4(self, s):
        self.sigmaT4 = s

    @set
    def radius(self, r):
        self.radius = r

class adcpars():
    L = None
    height = None
    def __init__(self, L=None, height=None):
        self.L = L
        self.height = height

    @set
    def L(self, L):
        self.L = L

    @set
    def height(self, height):
        self.height = height

class thirdlightparams():
    orbphase = None
    nbands = None
    Filter = []
    minlambda = []
    maxlambda = []
    fraction = []
    addFlux = None
    def __init__(self, orbphase=None, nbands=None, Filter=None, minlambda=None,
                 maxlambda=None, fraction=None, addFlux=None):
                     self.orbphase = orbphase
                     self.nbands = nbands
                     self.Filter = Filter
                     self.minlambda = minlambda
                     self.maxlambda = maxlambda
                     self.fraction = fraction
                     self.addFlux = addFlux

    @set
    def orbphase(self, o):
        self.orbphase = o

    @set
    def nbands(self, n):
        self.nbands = n

    @set
    def Filter(self, f):
        self.Filter = f

    @set
    def minlambda(self, m):
        self.minlambda = m

    @set
    def maxlambda(self, m):
        self.maxlambda = m

    @set
    def fraction(self, fraction):
        self.fraction = fraction

    @set
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
    Topy = []
    Bottomy = []
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
    
    @set
    def Nxtiles(self, x):
        self.Nxtiles = x

    @set
    def Nztiles(self, z):
        self.Nztiles = z

    @set
    def deltax(self, x):
        self.deltax = x

    @set
    def deltaz(self, z):
        self.deltaz = z

    @set
    def deltal(self, l):
        self.deltal = l

    @set
    def xmin(self, x):
        self.xmin = x

    @set
    def xmax(self, x):
        self.xmax = x

    @set
    def ymin(self, y):
        self.ymin = y

    @set
    def ymax(self, y):
        self.ymax = y

    @set
    def zmin(self, z):
        self.zmin = z

    @set
    def zmax(self, z):
        self.zmax = z

    @set
    def Topy(self, y):
        self.Topy = y

    @set
    def Bottomy(self, y):
        self.Bottomy = y

class dataparams():
    nbands = None
    filename = []
    Filter = []
    minlambda = []
    maxlambda = []
    npoints = None
    phase = []
    flux = []
    standdev = []
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
    @set
    def nbands(self, b):
        self.nbands = b

    @set
    def filename(self, filename):
        self.filename = filename

    @set
    def Filter(self, Filter):
        self.Filter = Filter

    @set
    def minlambda(self, minlambda):
        self.minlambda = minlambda

    @set
    def maxlambda(self, maxlambda):
        self.maxlambda = maxlambda

    @set
    def npointa(self, npoints):
        self.npoints = npoints

    @set
    def phase(self, phase):
        self.phase = phase

    @set
    def flux(self, flux):
        self.flux = flux

    @set
    def standdev(self, standdev):
        self.standdev = standdev

    @set
    def chisquare(self, chisquare):
        self.chisquare = chisquare