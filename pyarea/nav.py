__author__ = 'jano'
import gc
import math
import struct
from multiprocessing import cpu_count
import logging
import datetime
import numexpr as ne
import numpy as np
from numpy import nan

from pyarea import utils

ne.set_num_threads(cpu_count())
logger = logging.getLogger('nav')

class AreaNav():


    # Code value in AREA files used to designate DMSP navigation
    DMSP = 'DSMP'

    # Code value in AREA files used to designate GMSX (GMS) navigation
    GMSX = 'GMSX'

    # Code value in AREA files used to designate GOES (GOES D-H) navigation
    GOES = 'GOES'

    # Code value in AREA files used to designate GEOS navigation
    GEOS = 'GEOS'

    # Code value in AREA files used to designate GVAR (GOES I-M) navigation
    GVAR = 'GVAR'

    # Code value in AREA files used to designate MOLL (Mollweide) navigation
    MOLL = 'MOLL'

    # Code value in AREA files used to designate MSAT (Meteosat) navigation
    MSAT = 'MSAT'

    # Code value in AREA files used to designate MSGT  navigation
    MSGT = 'MSGT'

    # Code value in AREA files used to designate MSG navigation
    MSG  = 'MSG'

    # Code value in AREA files used to designate POES navigation
    POES =  'POES'

    # Code value in AREA files used to designate RADR (radar) navigation
    RADR = 'RADR'

    # Code value in AREA files used to designate RECT (rectilinear) navigation
    RECT =  'RECT'

    # Code value in AREA files used to designate PS (polar stereographic) navigation
    PS = 'PS'

    # Code value in AREA files used to designate MERC (mercator) navigation
    MERC = 'MERC'

    # Code value in AREA files used to designate TANC (tangent cone) navigation
    TANC =  'TANC'

    # Code value in AREA files used to designate SIN (sinusoidal cone) navigation
    SIN =  'SIN'

    # Code value in AREA files used to designate LAMB (lambert conformal)navigation
    LAMB = 'LAMB'

    # Code value in AREA files used to designate Lat/Lon
    LALO = 'LALO'

    # Code value in AREA files used to designate KALP
    KALP = 'KALP'

    # Code value in AREA files used to designate ABIS
    ABIS = 'ABIS'

    #ABI GOES R
    ABIN = 'ABIN'

    # Code value for specifying Latitude/Longitude transformations
    LL = 123

    # Code value for specifying Cartesian (X/Y) transformations
    XY = 234

    # "Line" index in line/element array
    indexLine=1
    # "Element" index in line/element array
    indexEle=0
    # "Latitude" index in latitude/longitude array
    indexLat=0
    # "Longitude" index in latitude/longitude array
    indexLon=1

    isLineFlipped = False
    lineOffset = 0.0

    # the following are ancillary info set by the set/get
    # public methods Res, Mag, and Start
    resLine = 1.
    resElement = 1
    magLine = 1.
    magElement = 1.
    startLine = 0
    startElement = 0.
    startImageLine = 0.
    startImageElement = 0.
    ##COMMON attrs added by jano purposefully
    PI = math.pi
    DEG = 180. / PI
    RAD = PI / 180.  # degrees to radians conversion pi/180
    STTYPE = 0  # word position of satellite type

    #java project uses makeAreaNav factory static method. I prefer different alternative,
    #to use __new__()
    def __new__(cls, *args, **kwargs):
        if not 'navblock' in kwargs:
            raise RuntimeError('Can not create {} instance without a navigation block'.format( cls.__name__))
        if not 'byte_order' in kwargs:
            raise RuntimeError('Can not create {} instance without knowing byte_order. Please supply "byte_order" argument.'.format( cls.__name__))
        if not 'int_format' in kwargs:
            raise RuntimeError(
                'Can not create {} instance without knowing format. Please supply "int_format" argument.'.format(cls.__name__))

        navblock= kwargs.get('navblock')


        byte_order= kwargs.get('byte_order')
        int_format= kwargs.get('int_format')
        #navtype = struct.unpack('<4s', struct.pack('<i', navblock[0]))[0]
        #navtype = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=byte_order, size=4,format_str='s'), navblock[:4])[0]
        navtype = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=byte_order, size=4, format_str='s'),
                      struct.pack('{byte_order:s}{format_str:s}'.format(byte_order=byte_order, format_str=int_format),
                                  navblock[0]))[0]
        #navtype = .from_bytes(navblock[:4])
        navtype = navtype.decode()

        #evaluate dynamically, does have the advantage of being self updateable
        #nav_classes = [klass.code_string() for name, klass in globals().items() if 'nav' in name.lower() and issubclass(klass, AreaNav)]
        nav_classes = cls.get_subclasses()



        # dynamically create the new class
        for klass in nav_classes:
            if klass.__name__.startswith(navtype):
                inst = object.__new__(klass)
                inst.navtype = navtype
                inst.byte_order = byte_order
                inst.int_format = int_format
                return  inst



        #raise RuntimeError(f'Navigation for {navtype} data is not implemented. Supported implementations are: {str(nav_classes)}')


    def toLatLon(self, linele=None):
        raise RuntimeError('Abstract method. Should be implemented by subclasses')

    def toLinEle(self, latlon=None):
        raise RuntimeError('Abstract method. Should be implemented by subclasses')

    def areaCoordtoImageCoord_np(self, linele=None):
        c, l = linele
        l1 =  self.lineOffset - l if self.isLineFlipped else l
        l1 = self.startImageLine + (self.resLine * (l1 - self.startLine)) / self.magLine
        c1 = self.startImageElement + (self.resElement * (c - self.startElement)) / self.magElement
        return [l1, c1]

    def areaCoordtoImageCoord(self, linele=None, newvals = None):
        l = len(linele[0])
        if newvals is None:
            newvals = [[nan for i in range(l)], [nan for i in range(l)]]

        for i in  range(l):
            if linele[self.indexLine][i] == linele[self.indexLine][i]:
                # account for flipped coordinates
                line  =  self.lineOffset - linele[self.indexLine][i] if self.isLineFlipped else linele[self.indexLine][i]
                #print 'line %s startImageLine %s startLine %s resLine %s magLine %s' % (line, self.startImageLine, self.startLine, self.resLine, self.magLine)
                newvals[self.indexLine][i] = self.startImageLine + (self.resLine * (line - self.startLine)) / self.magLine

            if linele[self.indexEle][i] == linele[self.indexEle][i]:
                #print 'elem %s startImageElement %s startElement %s resElement %s magElement %s' % (linele[self.indexEle][i], self.startImageElement, self.startElement, self.resElement, self.magElement)
                newvals[self.indexEle][i] = self.startImageElement + (self.resElement * (linele[self.indexEle][i] - self.startElement)) / self.magElement


        return newvals


    def imageCoordToAreaCoord(self, linele=None, newvals = None):
        l = len(linele[0])
        if newvals is None:
            newvals = [[nan for i in range(l)], [nan for i in range(l)]]

        for i in range(l):
            if linele[self.indexLine][i] == linele[self.indexLine][i]:
                newvals[self.indexLine][i] = self.startLine + ( self.magLine * (linele[self.indexLine][i] - self.startImageLine)) / self.resLine
            # account for flipped coordinates
                if self.isLineFlipped:
                    newvals[self.indexLine][i] =  self.lineOffset - newvals[self.indexLine][i]

            if linele[self.indexEle][i] == linele[self.indexEle][i]:
                newvals[self.indexEle][i] = self.startElement + ( self.magElement * (linele[self.indexEle][i] - self.startImageElement)) / self.resElement
        return newvals


    def imageCoordToAreaCoord_np(self, linele=None):
        l, c = linele
        l1 = self.startLine + ( self.magLine * (l - self.startImageLine)) / self.resLine
        if self.isLineFlipped:
            l1 -= self.lineOffset
        c1 = self.startElement + ( self.magElement * (c - self.startImageElement)) / self.resElement

        return [l1,c1]

    def canApproximateWithSpline(self):
        return True

    def __str__(self):
        return self.navtype

    def setImageStart(self, start_line, start_elem):
        self.startImageLine = float(start_line)
        self.startImageElement = float(start_elem)

    def setStart(self, start_line, start_elem):
        self.startLine = start_line
        self.startElement = start_elem
    def setRes(self, resLine, resElement):
        self.resLine = resLine
        self.resElement = resElement
    def setMag(self, magLine, magElement):
        self.magLine = magLine
        self.magElement = magElement

    @classmethod
    def get_subclasses(cls):

        subs = set(cls.__subclasses__())
        return subs.union(*(i.get_subclasses() for i in subs))



class GVARNav(AreaNav):

    isEastPositive = True


    NOMORB=42164.365 # nominal radial distance of instrument (km)
    AE=6378.137 # earth equatorial radius (km)
    FER=1.-(6356.7533/AE) # earth flattening coeff = 1-(be/ae)
    AEBE2= 1./(1-FER)**2
    AEBE3= AEBE2-1.
    AEBE4=((1.-FER)**4.)-1.

    xs = [nan]  * 3 # normalized s/c position in ecef coordinates
    #bt= [[nan for i in range(3)], [nan for i in range(3)], [nan for i in range(3)] ]
    bt = np.empty((3,3))
    bt[:] = np.nan
    #print bt
    # ecef to instrument coordinates transformation
    q3 = nan # used in subrtn lpoint
    pitch = roll = yaw = nan # pitch,roll,yaw angles of instrument (rad)
    pma = nan # pitch misalignment of instrument (rad)
    rma = nan # roll misalignment of instrument (rad)

    incmax = [6136, 2805] # number of increments per cycle

    elvmax = [0.220896, 0.22089375] # bounds in elevation (radians)
    scnmax= [0.24544,0.2454375] # bounds in scan angle (radians)
    elvinc = [8.e-6, 17.5e-6] # change in elevation angle per increment (rad)
    scninc = [16.e-6, 35.e-6] # change in scan angle per increment (radians)
    elvln = [28.e-6, 280.e-6] # elevation angle per detector line (radians)
    scnpx = [16.e-6, 280.e-6 ]# scan angle per pixel (radians)
    nsnom = [.220896, .22089375] # north-south center of instrument (=4.5 x incmax x elvinc)

    ewnom = [.24544, .2454375] # east-west center of instrument (=2.5 x incmax x sncinc)

    STTYPE = 0 # position of instrument type
    IDNTFR = 1
    IMCACT = 2 # position of imc active flag
    IYFLIP = 3 # position of yaw flip enabled flag
    REFLON = 5 # position of reference longitude
    REFDIS = 6 # position of reference distance from nominal
    REFLAT = 7 # position of reference latitude
    REFYAW = 8 # position of reference yaw
    RATROL = 9 # position of reference attitude roll
    RATPTC = 10 # position of reference attitude pitch
    RATYAW = 11 # position of reference attitude yaw
    ETIME  = 12 # position of epoch time
    EDTIME = 14 # location of delta from epoch time
    IMCROL = 15 # location of image motion compensation roll
    IMCPTC = 16 # location of image motion compensation pitch
    IMCYAW = 17 # location of image motion compensation yaw

    # ldr1-13: location of longitude delta from reference parameters
    LDR1   = 18
    LDR2   = 19
    LDR3   = 20
    LDR4   = 21
    LDR5   = 22
    LDR6   = 23
    LDR7   = 24
    LDR8   = 25
    LDR9   = 26
    LDR10  = 27
    LDR11  = 28
    LDR12  = 29
    LDR13  = 30

    # rddr1-11: location of radial distance delta from reference parameters
    RDDR1  = 31
    RDDR2  = 32
    RDDR3  = 33
    RDDR4  = 34
    RDDR5  = 35
    RDDR6  = 36
    RDDR7  = 37
    RDDR8  = 38
    RDDR9  = 39
    RDDR10 = 40
    RDDR11 = 41

    # dgl1-9: location of geocentric latitude delta parameters
    DGL1   = 42
    DGL2   = 43
    DGL3   = 44
    DGL4   = 45
    DGL5   = 46
    DGL6   = 47
    DGL7   = 48
    DGL8   = 49
    DGL9   = 50

    # doy1-9: location of orbit yaw delta parameters
    DOY1   = 51
    DOY2   = 52
    DOY3   = 53
    DOY4   = 54
    DOY5   = 55
    DOY6   = 56
    DOY7   = 57
    DOY8   = 58
    DOY9   = 59

    EXPTIM = 61 # exponential start time from epoch
    RAAWDS = 62 # location of start of roll attitude angle words
    PAAWDS = 129 # location of start of pitch attitude angle words
    YAAWDS = 184 # location of start of yaw attitude angle words
    RMAWDS = 257 # location of start of roll misalignment angle words
    PMAWDS = 312 # location of start of pitch misalignment angle words

    IMGDAY = 367 # position of image day value  (yyddd)
    IMGTM  = 368 # position of image time value (hhmmss)
    IMGSND = 369 # location of imager/sounder instrument flag
    ISTIM = 226

    # the following four words were added 5-26-94 to comply w/ the new elug...
    # numbering started at 380 because these same parameters are used
    # in the nav message sent from the ingestor to evx, and we had to
    # start somewhere after the 378 nav parameters

    IOFNC = 379
    IOFEC = 380
    IOFNI = 381
    IOFEI = 382

    MXCDSZ=5*128 # maximum directory block entry size

    OASIZE = 336 # size of gvar orbit and attitude set
    PCOEFS = 117 # start of pitch coefficients
    RMACFS = 227 # start of rma coefficients
    CUTOF1 = 115 # first dividing poin o&a set (128 word sets)
    CUTOF2 = 225 # second dividing poin o&a set (128 word sets)

    IMCFLG = 7 #  bit # in sscan for imc active flag
    FLPFLG = 15 #  bit # in sscan for yaw flip enabled flat
    iflip = None # value of FLPFLG

    # now some variables that are shared among the methods...

    aec = ts = dr = lam = dlat = dyaw = phi = psi = imc = nan
    aebe2c = aebe3c = aebe4c = ferc = nan
    instr = itype = None
    sublat = sublon = nan
    rlat = rlon = gam = alf = nan
    subpoint = [nan] * 2

    RELLST = (  ( 4,  10), ( 13,  63), ( 65,  94),
                 ( 98, 100), (103, 105), (108, 110), (113, 115),
                 (116, 118), (120, 149), (153, 155), (158, 160),
                 (163, 165), (168, 170), (171, 173), (175, 204),
                 (208, 210), (213, 215), (218, 220), (223, 225),
                 (226, 228), (230, 259), (263, 265), (268, 270),
                 (273, 275), (278, 283), (285, 314), (318, 320),
                 (323, 325), (328, 330), (333, 335))

    def __init__(self, navblock=None, auxblock=None, **kwargs):
        """
            doc
        """

        #python is not java but i declare this so I do not lose the track of instance vars

        te  = u = sinu = cosu = nan
        sinoi = cosoi = slat = asc = sinasc = cosasc = nan
        syaw = w = sw = cw = s2w = c2w = nan
        sw1 = cw1 = sw3 = cw3 = nan
        secs = wa = nan
        imgtim = epoch = timex = time50 = r = nan
        #b = [[nan for i in range(3)], [nan for i in range(3)], [nan for i in range(3)] ]
        b = np.empty((3,3))
        b[:] = np.nan
        lit = year = day = nan
        iftok = hour = minuts = nan
        count = offset = loop = nan
        time = [nan] * 2
        rparms = [nan] * self.MXCDSZ
        # rellst two dimensional array of location of float in o&a set

        # initialize the earth radius constants
        self.aec = self.AE
        self.ferc = self.FER
        self.aebe2c = self.AEBE2
        self.aebe3c = float(self.AEBE3)
        self.aebe4c = float(self.AEBE4)


        #navtype = struct.unpack( '<4s', struct.pack('<i', navblock[self.STTYPE]))[0]

        navtype = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=4, format_str='s'), struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=1, format_str=self.int_format), navblock[self.STTYPE]))[0]
        navtype = navtype.decode()

        if  navtype != self.GVAR:
            raise Exception(f'Invalid navigation type {navtype}')


        self.itype = 1
        #is this used somewhere???
        self.flywheel = None
        #convert to float, no Gould bullshit like in C, just int to float, same order of course

        for i in range(self.MXCDSZ):
            #rparms[i] = struct.unpack('<f', struct.pack('<i', navblock[i]))[0]
            rparms[i] = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=1, format_str='f'), struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=1, format_str=self.int_format), navblock[i]))[0]



        #divide floats as per OGE
        for i, e in enumerate(self.RELLST):
            li, ui = e
            offset = 1
            if li > self.CUTOF1: offset = 13
            if li > self.CUTOF2: offset = 31
            for j in range(li, ui):

                if (j ==  13 or j == 60 or (j - 7) % 55 == 0) and j != 7:
                    rparms[j+offset] = navblock[j+offset] / 100.
                else:
                    rparms[j + offset] = navblock[j + offset] / 10000000.

        self.instr = navblock[self.IMGSND]




        # ----------------------------------------------------------
        # new code from kath kelly for sounder nav - 10/27/??
        # because change to elug -- nadir postion is avabile in
        # signal and should be used to compute values instead
        # of making them hardwired...


        self.nadnsc = navblock[self.IOFNC]
        self.nadnsi = navblock[self.IOFNI]
        self.nadewc = navblock[self.IOFEC]
        self.nadewi = navblock[self.IOFEI]

        if self.nadnsc != 0 and self.nadnsi != 0 and self.nadewc != 0 and self.nadewi != 0:
            if self.instr == 1:
                self.elvmax[0] = (self.nadnsc* self.incmax[0]+self.nadnsi)*self.elvinc[0]
            else:
                self.elvmax[1]=( (9-self.nadnsc)*self.incmax[1]-self.nadnsi)*self.elvinc[1]
            self.scnmax[self.instr-1] = (self.nadewc*self.incmax[self.instr-1]+self.nadewi)*self.scninc[self.instr-1]


        # end of new code from kathy kelly
        # ----------------------------------------------------------
        #compute image time and epoch

        #java and c code uses this computation to get to epoch time in minutes
        #this could nevertheless be a source of rounding error
        #but hey, python has all this figured out through timedelta

        year = 1900 + navblock[self.IMGDAY] // 1000
        day = navblock[self.IMGDAY] - navblock[self.IMGDAY] / 1000 * 1000
        hour = int(rparms[self.IMGTM] / 10000)
        minuts = int(rparms[self.IMGTM] / 100 - hour * 100)
        secs = rparms[self.IMGTM] - float(100 * minuts) - float(10000 * hour)





        j = day + 1461 * (year + 4799) / 4 - 3 * ((year + 4899) / 100) / 4 - 2465022

        imgtim = float( j * 1440. + hour * 60. + minuts + (secs / 60.))


        #sdt = '%4d-%02d-%02d %02d:%02d:%0.2f' % (year, month, mday, hour, minuts, secs)
        #self.scan_datetime = datetime.datetime.strptime(sdt, '%Y-%m-%d %H:%M:%S.%f')
        self.scan_datetime = None

        jan_1_1950 = datetime.datetime(1950, 1, 1)
        '''
        my version of imgtime in minutres is not quite the same like the java code so i use that one

        imgtime = utils.yyyddd_hhmmss2datetime(navblock[self.IMGDAY], navblock[self.IMGTM])

        diff = imgtime - jan_1_1950
        imgtim1 = diff.total_seconds()/60.
        '''
        bcdtime = self.bcd2datetime(navblock[self.ETIME], navblock[self.ETIME+1])



        diff1 = bcdtime - jan_1_1950
        epoch1 = diff1.total_seconds()/60.
        '''
        # convert the BCD to binary integer
        #heer the my code and java is identical so i use mine
        #gain this is a port from java, but because the results are identical to my own code
        t0 = 0
        t1 = 0
        power10 = 1
        e0 = navblock[self.ETIME]
        e1 = navblock[self.ETIME + 1]
        for i in range(8):
            t0 = t0 + ( e0 & int('0xf', 16) ) * power10
            t1 = t1 + ( e1 & int('0xf', 16) ) * power10
            e0 = e0 >>  4
            e1 = e1 >>  4
            power10 = power10 * 10
        year = t0 / 10000
        day = int((t0 - (year * 10000)) * 0.1)
        iaa = t0 - (year * 10000)
        iab = (iaa - (day * 10)) * 10
        nbc = t1 / 10000000
        iac = t1 - (nbc * 10000000)
        def1 = t1 - iac
        hour = iab + nbc
        minuts = int(iac * 0.00001)
        s = (t1 - (def1 + (minuts * 100000))) * 0.001;
        j = day + 1461 * (year + 4799) / 4 - 3 * ((year + 4899) / 100) / 4 - 2465022
        epoch = float(j * 1440. + hour * 60. + float(minuts) + (s / 60.))
        print epoch, epoch1
        '''


        self.imc = 0 if ( (navblock[self.IMCACT] & (1 << self.IMCFLG)) != 0 ) else  1
        self.iflip = -1 if ( (navblock[self.IYFLIP] & (1 << self.FLPFLG)) != 0) else 1
        #print self.imc, self.iflip, imgtim
        # print self.imc
        # print self.iflip

        # assign reference values to the subsatellite longitude and
        # latitude, the radial distance and the orbit yaw.
        self.lam = rparms[self.REFLON]

        self.dr = rparms[self.REFDIS]

        self.phi = rparms[self.REFLAT]

        self.psi = rparms[self.REFYAW]


        self.subpoint[0] = rparms[self.REFLAT] / self.RAD
        self.subpoint[1] = rparms[self.REFLON] / self.RAD


        # assign reference values to the attitudes and misalignments
        self.roll = rparms[self.RATROL]
        self.pitch = rparms[self.RATPTC]
        self.yaw = rparms[self.RATYAW]


        self.rma = 0.
        self.pma = 0.
        # if imc is off, compute changes in the instrument orbit
        if self.imc != 0:
            # set reference radial distance, latitude and orbit yaw to zero
            self.dr = 0.
            self.phi = 0.
            self.psi = 0.

            # compute time since epoch (in minutes)
            self.ts = imgtim - epoch
            # computes orbit angle and the related trigonometric functions.
            # earth rotational rate=.729115e-4 (RAD/s)
            w = 0.729115e-4 * 60.0 * self.ts
            sw = math.sin(w)
            cw = math.cos(w)
            sw1 = math.sin(0.927 * w)
            cw1 = math.cos(0.927 * w)
            s2w = math.sin(2. * w)
            c2w = math.cos(2. * w)
            sw3 = math.sin(1.9268 * w)
            cw3 = math.cos(1.9268 * w)

            # computes change in the imc longitude from the reference
            self.lam = self.lam + rparms[self.LDR1] + (rparms[self.LDR2] + rparms[self.LDR3] * w) * w \
                        + (rparms[self.LDR10] * sw1 + rparms[self.LDR11] * cw1 + rparms[self.LDR4] * sw \
                        + rparms[self.LDR5] * cw + rparms[self.LDR6] * s2w + rparms[self.LDR7] * c2w \
                        + rparms[self.LDR8] * sw3 + rparms[self.LDR9] * cw3 + w * (rparms[self.LDR12] * sw + rparms[self.LDR13] * cw)) * 2.

            # computes change in radial distance from the reference (km)
            self.dr = self.dr + rparms[self.RDDR1] + rparms[self.RDDR2] * cw + rparms[self.RDDR3] * sw \
                        + rparms[self.RDDR4] * c2w + rparms[self.RDDR5] * s2w + rparms[self.RDDR6]*cw3 \
                        + rparms[self.RDDR7] * sw3 + rparms[self.RDDR8] * cw1 \
                        + rparms[self.RDDR9] * sw1 + w * (rparms[self.RDDR10] * cw + rparms[self.RDDR11] * sw)
            # computes the sine of the change in the geocentric latitude
            self.dlat = rparms[self.DGL1] + rparms[self.DGL2] * cw + rparms[self.DGL3] * sw \
                        + rparms[self.DGL4] * c2w + rparms[self.DGL5] * s2w + w * (rparms[self.DGL6] * cw \
                        + rparms[self.DGL7] * sw) + rparms[self.DGL8] * cw1 + rparms[self.DGL9] * sw1

            # computes geocentric latitude by using an expansion for arcsine
            self.phi = self.phi + self.dlat * (1. + self.dlat * self.dlat / 6.)

            # computes sine of the change in the orbit yaw
            self.dyaw = rparms[self.DOY1] + rparms[self.DOY2] * sw + rparms[self.DOY3] * cw \
                        + rparms[self.DOY4] * s2w + rparms[self.DOY5] * c2w \
                        + w * (rparms[self.DOY6] * sw + rparms[self.DOY7] * cw) \
                        + rparms[self.DOY8] * sw1 + rparms[self.DOY9] * cw1
            # computes the orbit yaw by using an expansion for arcsine.
            self.psi = self.psi + self.dyaw * (1. + self.dyaw * self.dyaw / 6.)

            # calculation of changes in the instrument orbit ends here

        # conversion of the imc longitude and orbit yaw to the subinstrument
        # longitude and the orbit inclination (ref: goes-pcc-tm-2473, inputs
        # required for earth location and gridding by sps, june 6, 1988)
        slat = math.sin(self.phi)
        syaw = math.sin(self.psi)
        sinoi = slat * slat + syaw * syaw
        cosoi = math.sqrt(1. - sinoi)
        sinoi = math.sqrt(sinoi)

        if slat == 0.0 and syaw == 0.0:
            u = 0.0
        else:
            u = math.atan2(slat, syaw)
        sinu = math.sin(u)
        cosu = math.cos(u)

        # computes longitude of the ascending node
        asc = self.lam - u
        sinasc = math.sin(asc)
        cosasc = math.cos(asc)

        # computes the subinstrument geographic latitude
        self.sublat = math.atan(self.aebe2c * math.tan(self.phi))

        # computes the subinstrument longitude
        self.sublon = asc + math.atan2(cosoi * sinu, cosu)

        # computes the spacecraft to earth fixed coordinates transformation
        # matrix:
        #     (vector in ecef coordinates) = b * (vector in s/c coordinates)

        b[0][1] = -sinasc * sinoi
        b[1][1] = cosasc * sinoi
        b[2][1] = -cosoi
        b[0][2] = -cosasc * cosu + sinasc * sinu * cosoi
        b[1][2] = -sinasc * cosu - cosasc * sinu * cosoi
        b[2][2] = -slat
        b[0][0] = -cosasc * sinu - sinasc * cosu * cosoi
        b[1][0] = -sinasc * sinu + cosasc * cosu * cosoi
        b[2][0] = cosu * sinoi

        # computes the normalized spacecraft position vector in earth fixed
        # coordinates - xs.
        r = (self.NOMORB + self.dr) / self.aec
        self.xs[0] = -b[0][2] * r
        self.xs[1] = -b[1][2] * r
        self.xs[2] = -b[2][2] * r

        # precomputes q3 (used in lpoint funciton (now in navToLatLon() )

        self.q3 = self.xs[0] * self.xs[0] + self.xs[1] * self.xs[1] + self.aebe2c * self.xs[2] * self.xs[2] - 1.0

        # computes the attitudes and misalignments if imc is off
        if self.imc != 0: # computes the solar orbit angle
            wa = rparms[61-1] * self.ts  # computes the difference between current time, ts, and the  # exponential time, iparms(62). note that both times are since epoch.
            te = self.ts - rparms[self.EXPTIM]  # computes roll + roll misalignment
            self.roll = self.roll + self.gatt(self.RAAWDS, rparms, navblock, wa, te)  # computes pitch + pitch misalignment
            self.pitch = self.pitch + self.gatt(self.PAAWDS, rparms, navblock, wa, te)  # computes yaw
            self.yaw = self.yaw + self.gatt(self.YAAWDS, rparms, navblock, wa, te)  # computes roll misalignment
            self.rma = float(self.gatt(self.RMAWDS, rparms, navblock, wa, te))# computes pitch misalignment
            self.pma = float(self.gatt(self.PMAWDS, rparms, navblock, wa, te)) # apply the earth sensor compensation if needed
            self.roll = self.roll + rparms[self.IMCROL]
            self.pitch = self.pitch + rparms[self.IMCPTC]
            self.yaw = self.yaw + rparms[self.IMCYAW]
            # end if (imc...)
        #this needs to be fixed, those params are passed for no reason as they are instance
        self.inst2e(self.roll, self.pitch, self.yaw, b, self.bt)




    def bcd2datetime(self, w1, w2):
        """
            converts a 2word 4 bytes binary string representing a datetime in  BCD format
            into a datetime object

        """
        names = 'year_1000', 'year_100', 'year_10', 'year_1', 'day_100', 'day_10', 'day_1', 'hour_10', 'hour_1', 'min_10', 'min_1', 'sec_10', 'sec_1', 'msec_100', 'msec_10', 'msec_1'
        #w1 and 2 should be  little endian, the code handles big endian so convert conditionally
        if self.byte_order == '<':
            #chars = struct.pack('>i', w1) + struct.pack('>i', w2)
            chars = struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='>', size=1, format_str=self.int_format), w1) + struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='>', size=1, format_str=self.int_format), w2)
        else:
            chars = struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='>', size=1, format_str=self.int_format),w1) + struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='>', size=1, format_str=self.int_format), w2)
        values = list()
        k = 0
        for i, byte in enumerate(chars):
            #char = ord(byte)  # int values
            for val in (byte >> 4, byte & 15):  # split in 2 nibbles, use shift by 4 and bitewise & 15 to get the values
                if i == 3:  # the most significant bit of the first nibble day_100 is the flywheel. GOD knows what that means
                    #set flywheel
                    mask = 1 << 0
                    self.flywheel = True if (val & mask) else False


                values.append(val)
        #create a dict from names and nibbles
        d = dict(zip(names, values))

        #compose year
        year = d['year_1000'] * 1000 + d['year_100'] * 100 + d['year_10'] * 10 + d['year_1']
        yday = d['day_100'] * 100 + d['day_10'] * 10 + d['day_1']
        month, day = utils.get_month_day(year, yday, True)
        hour = d['hour_10'] * 10 + d['hour_1']
        minute = d['min_10'] * 10 + d['min_1']
        second = d['sec_10'] * 10 + d['sec_1']
        milisecond = d['msec_100'] * 100 + d['msec_10'] * 10 + d['msec_1']
        return datetime.datetime(year, month, day, hour, minute, second, milisecond*1000)


    def gatt(self, k, floatblock, intblock, wa, te):
        """
            k =  starting position of a parameter subset in the real o&a set
            floatblock(mxcdxz) = input o&a parameter set
            intblock(mxcdxz) = input o&a parameter set
            double wa = input solar orbit angle in radians
            double te = input exponential time delay from epoch (minutes)
        """
        # local variables

        #ir, jr = mr = att

        # constant component
        att = floatblock[k + 2]

        #  computes the exponential term
        if te >= 0:
            att = att + floatblock[k] * math.exp(-te / floatblock[k + 1])


        # extracts the number of sinusoids
        ir = float(intblock[k + 3])
        i = int(ir)

        # calculation of sinusoids

        for j in range(1, i + 1):
            att = att + floatblock[k + 2 * j + 2] * math.cos(wa * float(j) + floatblock[k + 2 * j + 3])


        # pointer to the number of monomial sinusoids
        k = k + 34

        # extarct number of monomial sinusoids
        ir = float(intblock[k])
        kkk = int(intblock[k])


        # computes monomial sinusoids

        for l in range(1, kkk+1):

            ll = k + 5 * l  # order of sinusoid
            jr = float(intblock[ll - 4])

            # order of monomial sinusoid
            mr = float(intblock[ll - 3])
            att = att + floatblock[ll - 2] * math.pow((wa - floatblock[ll]), mr) * math.cos(jr * wa + floatblock[ll - 1])


        return att

    def inst2e(self, r, p, y, a, at):
        """
         r =  roll angle in radians
         p = pitch angle in radians
         y = yaw angle in radians
         a(3,3) = spacecraft to ecef coordinates transformation matrix
         at(3,3)= instrument to ecef coordinates transformation matrix
        """
        #rpy = [[nan for i in range(3)], [nan for i in range(3)], [nan for i in range(3)]]
        rpy = np.empty((3,3))
        rpy[:] = np.nan

        # we compute instrument to body coordinates transformation
        # matrix by using a small angle approximation of trigonometric
        # funktions of the roll, pitch and yaw.
        rpy[0][0] = 1. - 0.5 * (p * p + y * y)
        rpy[0][1] = -y
        rpy[0][2] = p
        rpy[1][0] = y + p * r
        rpy[1][1] = 1. - 0.5 * (y * y + r * r)
        rpy[1][2] = -r
        rpy[2][0] = -p + r * y
        rpy[2][1] = r + p * y
        rpy[2][2] = 1. - 0.5 * (p * p + r * r)

        # multiplication of matrices a and rpy


        for i in range(3):
            for j in range(3):
                at[i][j] = a[i][0] * rpy[0][j] + a[i][1] * rpy[1][j] +a[i][2] * rpy[2][j]


    def toLatLon2(self, lines=None, cols=None):
        """

        :param lines:
        :param cols:
        :return:
        """
        # try:
        #     # broadcasting is faster here
        #     lines = lines[:, np.newaxis]
        #     cols = cols[np.newaxis, :]
        # except Exception as e:  # scalar version, r
        #     print e, 1
        #     pass
        #dude, make sure the inout pyarea pyarea coordinates not GVAR image coordinates
        rl, rp = self.areaCoordtoImageCoord_np(linele=[cols, lines])
        #print rl, rp
        #rl, rp = lines, cols
        cshape = 3, rl.shape[0], rp.shape[1]
        g1 = np.empty(cshape)
        g1[:] = nan
        g = np.empty(cshape)
        g[:] = nan

        u = np.empty(cshape)
        u[:] = nan

        #  compute elevation and scan angles (e,s) related to input
        #  line and pixel numbers
        alpha0 = self.elvmax[0] - (rl - 4.5) * self.elvln[0]
        zeta0 = (rp - 1.0) * self.scnpx[self.instr - 1] - self.scnmax[self.instr - 1]

        # compute sign of misalignment corrections and origin offset
        ff = float(-self.iflip) if self.instr == 2 else float(self.iflip)
        doff = self.scnmax[self.instr - 1] - self.ewnom[self.instr - 1]

        # add new second order origin offset correction

        alpha = alpha0 - alpha0 * zeta0 * doff
        zeta = zeta0 + 0.5 * alpha0 * alpha0 * doff

        #  transform elevation and scan angles to geographic coordinates
        #  (this is the old 'lpoint' routine...

        # computes trigonometric funktions of the scan and elevation
        # angles corrected for the roll and pitch misalignments
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        cz = np.cos(zeta)
        da = alpha - self.pma * sa * (ff / cz + np.tan(zeta)) - self.rma * (1.0 - ca / cz)
        dz = zeta + ff * self.rma * sa

        # corrected scan angle
        cz = np.cos(dz)


        # computes pointing vector in instrument coordinates
        g[0] = np.sin(dz)
        g[1] = -cz * np.sin(da)
        g[2] = cz * np.cos(da)

        # transforms the pointing vector to earth fixed coordinates
        g1[0] = self.bt[0][0] * g[0] + self.bt[0][1] * g[1] + self.bt[0][2] * g[2]
        g1[1] = self.bt[1][0] * g[0] + self.bt[1][1] * g[1] + self.bt[1][2] * g[2]
        g1[2] = self.bt[2][0] * g[0] + self.bt[2][1] * g[1] + self.bt[2][2] * g[2]

        # computes coefficients and solves a quadratic equation to
        # find the intersect of the pointing vector with the earth
        # surface
        q1 = g1[0] * g1[0] + g1[1] * g1[1] + self.aebe2c * g1[2] * g1[2]
        q2 = self.xs[0] * g1[0] + self.xs[1] * g1[1] + self.aebe2c * self.xs[2] * g1[2]
        d = q2 * q2 - q1 * self.q3

        # if the discriminant of the equation, d, is negative, the
        # instrument points off the earth
        dm = np.abs(d) < 1. - 9
        d = np.where(dm, 0, d)
        dm1 = d >= 0


        d = np.sqrt(d)
        #print 'd2 is %s ' % d

        # slant distance from the instrument to the earth point
        h = -(q2 + d) / q1

        # cartesian coordinates of the earth point
        u[0] = self.xs[0] + h * g1[0]
        u[1] = self.xs[1] + h * g1[1]
        u[2] = self.xs[2] + h * g1[2]

        # sinus of geocentric latitude
        d1 = u[2] / np.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])


        rlat = np.where(dm1, np.arctan(self.aebe2c * d1 / np.sqrt(1. - d1 * d1)), nan)
        rlon = np.where(dm1, np.arctan2(u[1], u[0]), nan)

        rlat = rlat * self.DEG
        rlon = rlon * self.DEG

        #  put longitude into mcidas form
        if not self.isEastPositive:
            rlon = -rlon
        return rlat, rlon


    def toLatLon(self, linele=None):
        """
            Converts form satelite coordinates to Lat/Lon

        """

        rl = rp = rlat = rlon = nan
        number = len(linele[0])

        # alpha = elevation angle (rad)
        # zeta = scan angle (rad)
        q1 = q2 = d = h = alpha = zeta = alpha0 = zeta0 = ff = doff = nan
        g1 = [nan for i in range(3)]
        g = [nan for i in range(3)]
        u = [nan for i in range(3)]
        sa = ca =  da = dz = d1 = cz = nan

        latlon = [[nan for i in range(number)], [nan for i in range(number)]]

        #  transform line/pixel to geographic coordinates:
        imglinele = self.areaCoordtoImageCoord(linele=linele)

        for point in range(number):

            #  set input line/pixel numbers
            rl = imglinele[self.indexLine][point]
            rp = imglinele[self.indexEle][point]

            #  if doing sounder nav, have to trick routines into thinking image is
            #  at res 1, because nav routines take sounder res into account
            if self.instr == 2:
                rl = (rl + 9.) / 10.
                rp = (rp + 9.) / 10.


            #  compute elevation and scan angles (e,s) related to input
            #  line and pixel numbers

            if self.instr == 1:
                alpha0 = self.elvmax[0] - (rl - 4.5) * self.elvln[0]
            else:
                alpha0 = self.elvmax[1] - (rl - 2.5) * self.elvln[1]

            zeta0 = (rp - 1.0) * self.scnpx[self.instr - 1] - self.scnmax[self.instr - 1]
            #print alpha0, zeta0
            # compute sign of misalignment corrections and origin offset
            ff = float(-self.iflip) if self.instr == 2 else float(self.iflip)
            doff = self.scnmax[self.instr - 1] - self.ewnom[self.instr - 1]
            #print doff

            # add new second order origin offset correction
            alpha = alpha0 - alpha0 * zeta0 * doff
            zeta = zeta0 + 0.5 * alpha0 * alpha0 * doff
            #print alpha, zeta
            #  transform elevation and scan angles to geographic coordinates
            #  (this is the old 'lpoint' routine...

            # computes trigonometric funktions of the scan and elevation
            # angles corrected for the roll and pitch misalignments
            ca = math.cos(alpha)
            sa = math.sin(alpha)
            cz = math.cos(zeta)
            da = alpha - self.pma * sa * (ff / cz + math.tan(zeta)) - self.rma * (1.0 - ca / cz)
            dz = zeta + ff * self.rma * sa

            # corrected scan angle
            cz = math.cos(dz)


            # computes pointing vector in instrument coordinates
            g[0] = math.sin(dz)
            g[1] = -cz * math.sin(da)
            g[2] = cz * math.cos(da)

            # transforms the pointing vector to earth fixed coordinates
            g1[0] = self.bt[0][0] * g[0] + self.bt[0][1] * g[1] + self.bt[0][2] * g[2]
            g1[1] = self.bt[1][0] * g[0] + self.bt[1][1] * g[1] + self.bt[1][2] * g[2]
            g1[2] = self.bt[2][0] * g[0] + self.bt[2][1] * g[1] + self.bt[2][2] * g[2]

            # computes coefficients and solves a quadratic equation to
            # find the intersect of the pointing vector with the earth
            # surface
            q1 = g1[0] * g1[0] + g1[1] * g1[1] + self.aebe2c * g1[2] * g1[2]
            q2 = self.xs[0] * g1[0] + self.xs[1] * g1[1] + self.aebe2c * self.xs[2] * g1[2]
            d = q2 * q2 - q1 * self.q3

            if abs(d) < 1. - 9:
                d=0.
            # if the discriminant of the equation, d, is negative, the
            # instrument points off the earth


            if d >= 0.0:
                d = math.sqrt(d)

                # slant distance from the instrument to the earth point
                h = -(q2 + d) / q1

                # cartesian coordinates of the earth point
                u[0] = self.xs[0] + h * g1[0]
                u[1] = self.xs[1] + h * g1[1]
                u[2] = self.xs[2] + h * g1[2]

                # sinus of geocentric latitude
                d1 = u[2] / math.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])

                # geographic (geodetic) coordinates of the point
                rlat = math.atan(self.aebe2c * d1 / math.sqrt(1. - d1 * d1))
                rlon = math.atan2(u[1], u[0])
            else:
                latlon[self.indexLat][point] = nan
                latlon[self.indexLon][point] = nan
                continue

            rlat = rlat * self.DEG
            rlon = rlon * self.DEG

            #  put longitude into mcidas form
            if not self.isEastPositive:
                rlon = -rlon

            #  see if we have to convert to x y z coordinates
            if self.itype == 2: #what the heck is this? the func is commented
                pass
                # llcart(ylat,ylon,xlat,xlon,z)
            else:
                latlon[self.indexLat][point] = rlat
                latlon[self.indexLon][point] = rlon


        # end point for loop

        return latlon

    def toLinEle2(self, lats=None, lons=None):
        # print '#'*50,
        # print 'AREA',
        # print '#' * 50

        #print self.resLine, self.resElement, self.magLine, self.magElement, self.startLine, self.startElement, self.startImageLine, self.startImageElement, self.isLineFlipped, self.lineOffset
        #print self.nadnsc, self.nadnsi, self.nadewc, self.nadewi, self.elvmax[0], self.scnmax[0]
        h_, w_ = lats.shape

        latsm = lats != nan
        lonsm = lons != nan


        hw = lats[latsm].size
        # f = np.empty((3,h_,w_))
        f = np.empty((3,hw))
        # ft = np.empty((3,h_,w_))
        ft = np.empty((3,hw))
        # u = np.empty((3,h_,w_))
        u = np.empty((3,hw))
        f[:] = nan
        ft[:] = nan
        u[:] = nan

        #lats and lons can contain nans
        ff = float(-self.iflip) if self.instr == 2 else float(self.iflip)
        doff = self.scnmax[self.instr - 1] - self.ewnom[self.instr - 1]
        rlat = lats[latsm] * self.RAD
        rlon = lons[lonsm] * self.RAD
        if not self.isEastPositive:
            rlon = -rlon

        # transform lat/lon to elevation and scan angles
        # (used to be the gpoint routine...)

        # computes sinus of geographic (geodetic) latitude
        sing = np.sin(rlat)
        w1 = self.aebe4c * sing * sing
        # sinus of the geocentric latitude
        slat = ((0.375 * w1 - 0.5) * w1 + 1.) * sing / self.aebe2c

        # computes local earth radius at specified point
        w2 = slat * slat
        w1 = self.aebe3c * w2
        w1 = (0.375 * w1 - 0.5) * w1 + 1.

        # computes cartesian coordinates of the point
        u[2] = slat * w1
        w2 = w1 * np.sqrt(1. - w2)
        u[0] = w2 * np.cos(rlon)
        u[1] = w2 * np.sin(rlon)
        # pointing vector from instrument to the earth point
        f[0] = u[0] - self.xs[0]
        f[1] = u[1] - self.xs[1]
        f[2] = u[2] - self.xs[2]
        w2 = u[0] * f[0] + u[1] * f[1] + u[2] * f[2] * self.aebe2c

        #wmask = w2 <= 0.
        wmask = ~np.isnan(w2)


        # converts pointing vector to instrument coordinates
        ft[0] = self.bt[0][0] * f[0] + self.bt[1][0] * f[1] + self.bt[2][0] * f[2]
        ft[1] = self.bt[0][1] * f[0] + self.bt[1][1] * f[1] + self.bt[2][1] * f[2]
        ft[2] = self.bt[0][2] * f[0] + self.bt[1][2] * f[1] + self.bt[2][2] * f[2]


        # converts pointing vector to scan and elevation angles and
        # corrects for the roll and pitch misalignments
        gam = np.arctan(ft[0] / np.sqrt(ft[1] * ft[1] + ft[2] * ft[2]))
        alf = -np.arctan(ft[1] / ft[2])
        w1 = np.sin(alf)
        w2 = np.cos(gam)
        alpha1 = alf + self.rma * (1. - np.cos(alf) / w2) + self.pma * w1 * (doff / w2 + np.arctan(gam))
        gam = gam - ff * self.rma * w1
        alf = alpha1 + alpha1 * gam * doff
        gam = gam - 0.5 * alpha1 * alpha1 * doff

        # convert elevation and scan angles to line/pixel coordinates

        # compute fractional line number

        tmplin = (self.elvmax[self.instr - 1] - alf) / self.elvln[self.instr - 1]
        tmplin = tmplin + 4.5 if self.instr == 1 else tmplin + 2.5

        # compute fractional pixel number
        tmpele = (self.scnmax[self.instr - 1] + gam) / self.scnpx[self.instr - 1] + 1.

        # if doing sounder nav, change lin & ele returned to res 10 values
        if self.instr == 2:
            tmplin *= 10. - 9.
            tmpele *= 10. - 9.

        lin = np.where(wmask, tmplin, 0)
        ele = np.where(wmask, tmpele, 0)
        # print self.startLine, self.startElement
        # print self.magLine, self.magElement, self.resLine
        # print self.startImageLine, self.startImageElement
        alin, aele= self.imageCoordToAreaCoord_np(linele=[lin, ele])
        #print alin, aele


        olin = np.empty_like(lats)

        olin[:] = 0
        oele = np.empty_like(lons)
        oele[:] = 0
        olin[latsm] = alin
        oele[lonsm] = aele
        return olin, oele



    def toLinEle(self, latlon=None):

        tmplin = tmpele = sing = slat = w1 = w2 = ff = doff = alpha1 = rlat = rlon = gam = alf = nan
        f = [nan] * 3
        ft = [nan] * 3
        u = [nan] * 3
        number = len(latlon[0])
        linele = [[nan for i in range(number)], [nan for i in range(number)]]
        ff = float(-self.iflip) if self.instr == 2 else float(self.iflip)
        doff = self.scnmax[self.instr - 1] - self.ewnom[self.instr - 1]

        for point in range(number):

            if abs(latlon[self.indexLat][point]) > 90.:
                linele[self.indexLine][point] = nan
                linele[self.indexEle][point] = nan
                continue

            rlat = float(latlon[self.indexLat][point] * self.RAD)
            rlon = float(latlon[self.indexLon][point] * self.RAD)
            if not self.isEastPositive:
                rlon = -rlon

            # transform lat/lon to elevation and scan angles
            # (used to be the gpoint routine...)

            # computes sinus of geographic (geodetic) latitude
            sing = math.sin(rlat)
            w1 = self.aebe4c * sing * sing

            # sinus of the geocentric latitude
            slat = ((0.375 * w1 - 0.5) * w1 + 1.) * sing / self.aebe2c

            # computes local earth radius at specified point
            w2 = slat * slat
            w1 = self.aebe3c * w2
            w1 = (0.375 * w1 - 0.5) * w1 + 1.

            # computes cartesian coordinates of the point
            u[2] = slat * w1
            w2 = w1 * math.sqrt(1. - w2)
            u[0] = w2 * math.cos(rlon)
            u[1] = w2 * math.sin(rlon)

            # pointing vector from instrument to the earth point
            f[0] = u[0] - self.xs[0]
            f[1] = u[1] - self.xs[1]
            f[2] = u[2] - self.xs[2]
            w2 = u[0] * f[0] + u[1] * f[1] + u[2] * f[2] * self.aebe2c

            # verifies visibility of the point
            if w2 <= 0.0:
                # converts pointing vector to instrument coordinates
                ft[0] = self.bt[0][0] * f[0] + self.bt[1][0] * f[1] + self.bt[2][0] * f[2]
                ft[1] = self.bt[0][1] * f[0] + self.bt[1][1] * f[1] + self.bt[2][1] * f[2]
                ft[2] = self.bt[0][2] * f[0] + self.bt[1][2] * f[1] + self.bt[2][2] * f[2]

                # converts pointing vector to scan and elevation angles and
                # corrects for the roll and pitch misalignments
                gam = math.atan(ft[0] / math.sqrt(ft[1] * ft[1] + ft[2] * ft[2]))
                alf = -math.atan(ft[1] / ft[2])
                w1 = math.sin(alf)
                w2 = math.cos(gam)
                alpha1 = alf + self.rma * (1. - math.cos(alf) / w2) + self.pma * w1 * (doff / w2 + math.tan(gam))
                gam = gam - ff * self.rma * w1
                alf = alpha1 + alpha1 * gam * doff
                gam = gam - 0.5 * alpha1 * alpha1 * doff

            else:
                # not visible...
                linele[self.indexLine][point] = nan
                linele[self.indexEle][point] = nan
                continue
            # convert elevation and scan angles to line/pixel coordinates

            # compute fractional line number

            tmplin = (self.elvmax[self.instr - 1] - alf) / self.elvln[self.instr - 1]
            tmplin = tmplin + 4.5 if self.instr == 1 else tmplin + 2.5


            # compute fractional pixel number
            tmpele = (self.scnmax[self.instr - 1] + gam) / self.scnpx[self.instr - 1] + 1.

            # convert internal 8 byte values to 4 bytes
            linele[self.indexLine][point] = float(tmplin)
            linele[self.indexEle][point] = float(tmpele)

            # if doing sounder nav, change lin & ele returned to res 10 values
            if self.instr == 2:
                linele[self.indexLine][point] = linele[self.indexLine][point] * 10. - 9.
                linele[self.indexEle][point] = linele[self.indexEle][point] * 10. - 9.


            # end for loop on points
        # Return in 'File' coordinates

        return self.imageCoordToAreaCoord(linele, linele)

    def getSubpoint(self):

        return self.sublat * self.DEG, self.sublon * self.DEG



class LALONav(AreaNav):
    def __init__(self, navblock=None, auxblock=None, **kwargs):
        pass


class GOESNav(AreaNav):
    isEastPositive = True

    # NAVCOM variables
    navday = nan
    lintot = nan
    deglin = nan
    ieltot = nan
    degele = nan
    spinra = nan
    ietimy = nan
    ietimh = nan
    semima = nan
    oeccen = nan
    orbinc = nan
    perhel = nan
    asnode = nan
    nopcln = nan
    declin = nan
    rascen = nan
    piclin = nan
    prerat = nan
    predir = nan
    pitch = nan
    yaw = nan
    roll = nan
    skew = nan

    # BETCOM variables
    iajust = nan
    ibtcon = nan
    negbet = nan
    iseang = nan

    # VASCOM variables
    scan1 = nan
    time1 = nan
    scan2 = nan
    time2 = nan

    # NAVINI variables
    emega = nan
    ab = nan
    asq = nan
    bsq = nan
    r = nan
    rsq = nan
    rdpdg = nan
    numsen = nan
    totlin = nan
    radlin = nan
    totele = nan
    radele = nan
    picele = nan
    cpitch = nan
    cyaw = nan
    croll = nan
    pskew = nan
    rfact = nan
    roasin = nan
    tmpscl = nan
    b11 = nan
    b12 = nan
    b13 = nan
    b21 = nan
    b22 = nan
    b23 = nan
    b31 = nan
    b32 = nan
    b33 = nan
    gamma = nan
    gamdot = nan
    rotm11 = nan
    rotm13 = nan
    rotm21 = nan
    rotm23 = nan
    rotm31 = nan
    rotm33 = nan
    pictim = nan
    xref = nan

    iold = 0

    # variables needed for satvec
    tdife = nan
    xmmc = nan
    epsiln = nan
    srome2 = nan
    pz = nan
    py = nan
    px = nan
    qz = nan
    qy = nan
    qx = nan

    def __init__(self, navblock=None, auxblock=None, **kwargs):

        navtype = \
        struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=4, format_str='s'),
                      struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=1,
                                                                                format_str=self.int_format),
                                  navblock[self.STTYPE]))[0]
        if navtype != self.GOES:
            raise Exception(f'Invalid navigation type {navtype}')
        jday = navblock[1]
        jtime = navblock[2]

        # INTIALIZE NAVCOM
        self.navday = jday % 100000
        # print jday, jtime, self.navday

        # for i in range(6,12):
        #     if navblock[i] <=0:
        #         raise RuntimeError('Invalid orbital parameter {}:{}. Should be greater than 0'.format(i, navblock[i]) )

        self.ietimy = self.icon1(navblock[4])
        self.ietimh = 100 * (navblock[5] / 100) + round(.6 * (navblock[5] % 100))
        self.semima = navblock[6] / 100.
        self.oeccen = navblock[7] / 1000000.
        self.orbinc = navblock[8] / 1000.
        self.xmeana = navblock[9] / 1000.
        self.perhel = navblock[10] / 1000.
        self.asnode = navblock[11] / 1000.
        # print self.ietimy
        # print self.ietimh
        # call epoch
        self.epoch()
        # print self.ietimy
        # print self.ietimh
        self.declin = utils.mcPI2F(navblock[12])
        self.rascen = utils.mcPI2F(navblock[13])
        self.piclin = navblock[14]
        if navblock[14] > 1000000:
            self.piclin /= 10000
        if (navblock[12] == 0 and navblock[13] == 0 and navblock[14] == 0):
            raise ValueError('Invalid ascension/declination parameters')
        if navblock[15] == 0:
            raise ValueError('Invalid spin period')
        self.spinra = navblock[15] / 1000.
        if navblock[15] != 0 and self.spinra < 300:
            self.spinra = 60000. / self.spinra
        self.deglin = utils.mcPI2F(navblock[16])
        self.lintot = navblock[17]
        self.degele = utils.mcPI2F(navblock[18])
        self.ieltot = navblock[19]
        self.pitch = utils.mcPI2F(navblock[20])
        self.yaw = utils.mcPI2F(navblock[21])
        self.roll = utils.mcPI2F(navblock[22])
        self.skew = navblock[28] / 100000.
        if navblock[28] == utils.MCMISSING:
            self.skew = 0.
        ## BETCOM
        self.iajust = navblock[24]
        self.iseang = navblock[27]
        self.ibtcon = 6289920
        self.negbet = 3144960

        ## ----- NAVINI
        self.emega = .26251617
        self.ab = 40546851.22
        self.asq = 40683833.48
        self.bsq = 40410330.18
        self.r = 6371.22
        self.rsq = self.r ** 2
        self.rdpdg = 1.745329252E-02
        self.numsen = (self.lintot / 100000) % 100
        if self.numsen < 1:
            self.numsen = 1
        self.totlin = self.numsen * (self.lintot % 100000)
        self.radlin = self.rdpdg * self.deglin / (self.totlin - 1.)
        self.totele = self.ieltot
        self.radele = self.rdpdg * self.degele / (self.totele - 1.)
        self.picele = (1. + self.totele) / 2.
        self.cpitch = self.rdpdg * self.pitch
        self.cyaw = self.rdpdg * self.yaw
        self.croll = self.rdpdg * self.roll

        self.pskew = math.atan2(self.skew, self.radlin / self.radele)

        stp = math.sin(self.cpitch)
        ctp = math.cos(self.cpitch)

        sty = math.sin(self.cyaw - self.pskew)
        cty = math.cos(self.cyaw - self.pskew)

        _str = math.sin(self.croll)
        ctr = math.cos(self.croll)

        self.rotm11 = ctr * ctp
        self.rotm13 = sty * _str * ctp + cty * stp
        self.rotm21 = -_str
        self.rotm23 = sty * ctr
        self.rotm31 = -ctr * stp
        self.rotm33 = cty * ctp - sty * _str * stp
        self.rfact = self.rotm31 ** 2 + self.rotm33 ** 2
        self.roasin = math.atan2(self.rotm31, self.rotm33)
        self.tmpscl = self.spinra / 3600000.
        dec = self.declin * self.rdpdg
        sindec = math.sin(dec)
        cosdec = math.cos(dec)
        ras = self.rascen * self.rdpdg
        sinras = math.sin(ras)
        cosras = math.cos(ras)

        self.b11 = -sinras
        self.b12 = cosras
        self.b13 = 0.
        self.b21 = -sindec * cosras
        self.b22 = -sindec * sinras
        self.b23 = cosdec
        self.b31 = cosdec * cosras
        self.b32 = cosdec * sinras
        self.b33 = sindec

        jan_1_1974 = datetime.datetime(year=1974, month=1, day=1)
        nav_daytime = utils.yyyddd_hhmmss2datetime(self.navday, 0)
        diff = nav_daytime - jan_1_1974
        time_diff_mins = diff.total_seconds() / 60.
        raha = time_diff_mins * 1.00273791 / 4.0 + 100.26467

        rac = raha % 360.
        if rac < 0:
            rac += 360.
        self.xref = rac * self.rdpdg
        # TIME SPECIFIC
        self.pictim = utils.mcPI2F(jtime)
        self.gamma = navblock[38] / 100.
        self.gamdot = navblock[39] / 100.

        # VASCOM
        iss = jday / 100000
        if (iss > 25 or iss == 12) and navblock[30] > 0:

            #       THIS SECTION DOES VAS BIRDS AND GMS
            #       IT USES TIMES AND SCAN LINE FROM BETA RECORDS

            self.scan1 = float(navblock[30])
            self.time1 = utils.mcPI2F(navblock[31])
            self.scan2 = float(navblock[34])
            self.time2 = utils.mcPI2F(navblock[35])
        else:
            # THIS SECTION DOES THE OLD GOES BIRDS
            self.scan1 = 1.
            self.time1 - utils.mcPI2F(jtime)
            self.scan2 = float(self.lintot % 100000)
            self.time2 = self.time1 + self.scan2 * self.tmpscl

        self.iold = 0

    def toLatLon(self, linele=None):

        number = len(linele[0])
        latlon = [[nan for i in range(number)], [nan for i in range(number)]]

        imglinele = self.areaCoordtoImageCoord(linele)
        for point in range(number):
            xlin = imglinele[self.indexLine][point]
            xele = imglinele[self.indexEle][point]
            ilin = round(xlin)
            parlin = (ilin - 1) / self.numsen + 1
            framet = self.tmpscl * parlin
            samtim = framet + self.pictim
            xyz = self.satvec(samtim)
            # print xyz
            ylin = (xlin - self.piclin) * self.radlin
            yele = (xele - self.picele + self.gamma + self.gamdot * samtim) * self.radele

            xcor = self.b11 * xyz[0] + self.b12 * xyz[1] + self.b13 * xyz[2]
            ycor = self.b21 * xyz[0] + self.b22 * xyz[1] + self.b23 * xyz[2]
            rot = math.atan2(ycor, xcor) + math.pi
            yele = yele - rot
            coslin = math.cos(ylin)
            sinlin = math.sin(ylin)
            sinele = math.sin(yele)
            cosele = math.cos(yele)
            eli = self.rotm11 * coslin - self.rotm13 * sinlin
            emi = self.rotm21 * coslin - self.rotm23 * sinlin
            eni = self.rotm31 * coslin - self.rotm33 * sinlin
            temp = eli
            eli = cosele * eli + sinele * emi
            emi = -sinele * temp + cosele * emi
            elo = self.b11 * eli + self.b21 * emi + self.b31 * eni
            emo = self.b12 * eli + self.b22 * emi + self.b32 * eni
            eno = self.b13 * eli + self.b23 * emi + self.b33 * eni
            basq = self.bsq / self.asq
            onemsq = 1.0 - basq
            aq = basq + onemsq * math.pow(eno, 2)
            bq = 2.0 * ((elo * xyz[0] + emo * xyz[1]) * basq + eno * xyz[2])
            cq = (math.pow(xyz[0], 2) + math.pow(xyz[1], 2)) * basq + math.pow(xyz[2], 2) - self.bsq
            rad = math.pow(bq, 2) - 4.0 * aq * cq
            if rad < 1.0:
                latlon[self.indexLat][point] = nan
                latlon[self.indexLon][point] = nan
            else:
                s = -(bq + math.sqrt(rad)) / (2.0 * aq)
                x = xyz[0] + elo * s
                y = xyz[1] + emo * s
                z = xyz[2] + eno * s

                ct = math.cos(self.emega * samtim + self.xref)
                st = math.sin(self.emega * samtim + self.xref)
                x1 = ct * x + st * y
                y1 = -st * x + ct * y
                ll = self.nxyzll(x1, y1, z)
                latlon[self.indexLat][point] = float(ll[0])
                # put longitude into East Positive (form)
                if self.isEastPositive:
                    rlat = -ll[1]
                    latlon[self.indexLon][point] = float(rlat)
                else:
                    latlon[self.indexLon][point] = float(ll[1])

        return latlon

    def nxyzll(self, x, y, z):
        '''
        SUBROUTINE NXYZLL(X,Y,Z,XLAT,XLON)
        CONVERT EARTH-CENTERED X,Y,Z TO LAT & LON
        X,Y,Z ARE GIVEN IN KM. THEY ARE THE COORDINATES IN A RECTANGULAR
        COORDINATE SYSTEM WITH ORIGIN AT THE EARTH CENTER, WHOSE POS.
        X-AXIS PIERCES THE EQUATOR AT LON 0 DEG, WHOSE POSITIVE Y-AXIS
        PIERCES THE EQUATOR AT LON 90 DEG, AND WHOSE POSITIVE Z-AXIS
        INTERSECTS THE NORTH POLE.
        XLAT,XLON ARE IN DEGREES, WITH NORTH AND WEST POSITIVE
        '''
        xlat = nan
        xlon = nan
        if x != 0 and y != 0 and z != 0:
            a = math.atan(z / math.sqrt(x * x + y * y))
            xlat = math.atan2(self.asq * math.sin(a), self.bsq * math.cos(a)) / self.rdpdg
            xlon = -math.atan2(y, x) / self.rdpdg

        return xlat, xlon

    def satvec(self, samtim):

        if self.iold != 1:
            self.iold = 1

            PI = 3.14159265
            rdpdg = PI / 180.0
            re = 6378.388
            gracon = .07436574
            sha = 100.26467
            sha = rdpdg * sha
            irayd = 74001
            irahms = 0
            o = rdpdg * self.orbinc
            p = rdpdg * self.perhel
            a = rdpdg * self.asnode

            so = math.sin(o)
            co = math.cos(o)
            sp = math.sin(p) * self.semima
            cp = math.cos(p) * self.semima
            sa = math.sin(a)
            ca = math.cos(a)
            self.px = cp * ca - sp * sa * co
            self.py = cp * sa + sp * ca * co
            self.pz = sp * so
            self.qx = -sp * ca - cp * sa * co
            self.qy = -sp * sa + cp * ca * co
            self.qz = cp * so
            self.srome2 = ((1. - self.oeccen) ** .5) * ((1 + self.oeccen) ** .5)
            self.xmmc = gracon * re * ((re / self.semima) ** .5) / self.semima
            iey = (self.ietimy / 1000) % 100
            ied = self.ietimy % 1000
            iefac = (iey - 1) / 4 + 1
            de = 365 * (iey - 1) + iefac + ied - 1

            te = 1440. * de + 60. * utils.mcPI2F(self.ietimh)

            iray = irayd / 1000
            irad = irayd % 1000
            irafac = (iray - 1) / 4 + 1
            dra = 365 * (iray - 1) + irafac + irad - 1
            tra = 1440. * dra + 60. * utils.mcPI2F(irahms)

            inavy = (self.navday / 1000) % 100
            inavd = self.navday % 1000
            infac = (inavy - 1) / 4 + 1
            dnav = 365 * (inavy - 1) + infac + inavd - 1
            self.tdife = float(dnav) * 1440. - te  # tdifra = dnav * 1440. - tra
        self.epsiln = 1.0E-8  # guess this is a BUG in java code because this line is below the upper if a then it is used below

        timsam = samtim * 60.
        diftim = self.tdife + timsam
        xmanom = self.xmmc * diftim
        ecanm1 = xmanom
        ecanom = 0

        for i in range(20):
            ecanom = xmanom + self.oeccen * np.sin(ecanm1)
            if np.all(ecanom - ecanm1 < self.epsiln):
                break
            ecanm1 = ecanom

        xomega = np.cos(ecanom) - self.oeccen
        yomega = self.srome2 * np.sin(ecanom)

        z = xomega * self.pz + yomega * self.qz
        y = xomega * self.py + yomega * self.qy
        x = xomega * self.px + yomega * self.qx
        return x, y, z

    def toLinEle(self, latlon=None):

        number = len(latlon[0])
        linele = [[nan for i in range(number)], [nan for i in range(number)]]

        for point in range(number):
            xpar = latlon[self.indexLat][point]
            # expects positive West Longitude.
            ypar = latlon[self.indexLon][point]
            if self.isEastPositive:
                ypar = -ypar
            xlin = nan
            xele = nan

            # if abs(xpar <= 90.:) original
            _xpar = np.abs(xpar)
            if np.any(_xpar <= 90.):
                oldlin = 910.
                orbtim = -99999.
                xsat = ysat = zsat = 0.0
                x = y = z = 0.0
                xht = znorm = 0.0
                xyz = self.nllxyz(xpar, ypar)
                x1 = xyz[0]
                y1 = xyz[1]
                z = xyz[2]
                samtim = self.time1

                for i in range(2):
                    if abs(samtim - orbtim) >= 0.0005:
                        xyzsat = self.satvec(samtim)
                        xsat = xyzsat[0]
                        ysat = xyzsat[1]
                        zsat = xyzsat[2]
                        orbtim = samtim
                        xht = math.sqrt(math.pow(xyzsat[0], 2) + math.pow(xyzsat[1], 2) + math.pow(xyzsat[2], 2))

                    ct = math.cos(self.emega * samtim + self.xref)
                    st = math.sin(self.emega * samtim + self.xref)
                    x = ct * x1 - st * y1
                    y = st * x1 + ct * y1
                    vcste1 = x - xsat
                    vcste2 = y - ysat
                    vcste3 = z - zsat
                    vcses3 = self.b31 * vcste1 + self.b32 * vcste2 + self.b33 * vcste3
                    znorm = np.sqrt(np.power(vcste1, 2) + np.power(vcste2, 2) + np.power(vcste3, 2))
                    x3 = vcses3 / znorm
                    umv = np.arctan2(x3, np.sqrt(self.rfact - np.power(x3, 2))) - self.roasin
                    xlin = self.piclin - umv / self.radlin

                    if i == 0:
                        samtim = self.time2
                        oldlin = xlin

                scnnum = ((oldlin + xlin) / 2.0 - 1.0) / self.numsen
                scnfrc = (scnnum - self.scan1) / (self.scan2 - self.scan1)
                xlin = oldlin + scnfrc * (xlin - oldlin)
                samtim = self.time1 + self.tmpscl * (scnnum - self.scan1)
                xyzsat = self.satvec(samtim)
                xsat = xyzsat[0]
                ysat = xyzsat[1]
                zsat = xyzsat[2]
                cosa = x * xsat + y * ysat + z * zsat
                ctst = 0.0001 * self.r * xht + self.rsq
                if np.any(cosa >= ctst):
                    xsats1 = self.b11 * xsat + self.b12 * ysat + self.b13 * zsat
                    ysats2 = self.b21 * xsat + self.b22 * ysat + self.b23 * zsat
                    ct = np.cos(self.emega * samtim + self.xref)
                    st = np.sin(self.emega * samtim + self.xref)
                    x = ct * x1 - st * y1
                    y = st * x1 + ct * y1
                    vcste1 = x - xsat
                    vcste2 = y - ysat
                    vcste3 = z - zsat
                    vcses1 = self.b11 * vcste1 + self.b12 * vcste2 + self.b13 * vcste3
                    vcses2 = self.b21 * vcste1 + self.b22 * vcste2 + self.b23 * vcste3
                    vcses3 = self.b31 * vcste1 + self.b32 * vcste2 + self.b33 * vcste3
                    xnorm = np.sqrt(np.power(znorm, 2) - np.power(vcses3, 2))
                    ynorm = np.sqrt(np.power(xsats1, 2) + np.power(ysats2, 2))
                    znorm = np.sqrt(np.power(vcste1, 2) + np.power(vcste2, 2) + np.power(vcste3, 2))
                    x3 = vcses3 / znorm
                    umv = np.arctan2(x3, np.sqrt(self.rfact - np.power(x3, 2))) - self.roasin
                    slin = np.sin(umv)
                    clin = np.cos(umv)
                    u = self.rotm11 * clin + self.rotm13 * slin
                    v = self.rotm21 * clin + self.rotm23 * slin
                    xele = self.picele + np.arcsin((xsats1 * vcses2 - ysats2 * vcses1) / (xnorm * ynorm)) / self.radele
                    xele = xele + np.arctan2(v, u) / self.radele
                    xele = xele - self.gamma - self.gamdot * samtim

            linele[self.indexLine][point] = xlin
            linele[self.indexEle][point] = xele

        # Return in 'File' coordinates
        return self.imageCoordToAreaCoord(linele, linele)

    def toLinEle2(self, latlon=None):

        ne_e = ne.evaluate

        point = 0

        linele = [[nan] for i in range(len(latlon))]

        xpar = latlon[self.indexLat]
        # expects positive West Longitude.
        ypar = latlon[self.indexLon]
        if self.isEastPositive:
            ypar = -ypar
        xlin = 0
        xele = 0

        # _xpar=ne_e("abs(xpar)")
        if np.any(abs(xpar) <= 90.):
            # del _xpar
            oldlin = 910.
            orbtim = -99999.
            xsat = ysat = zsat = 0.0, 0.0, 0.0
            x = y = z = 0.0, 0.0, 0.0
            xht = znorm = 0.0, 0.0
            xyz = self.nllxyz_ne(xpar, ypar)

            x1 = xyz[0]
            y1 = xyz[1]
            z = xyz[2]
            del xyz
            samtim = self.time2

            if abs(samtim - orbtim) >= 0.0005:
                xyzsat = self.satvec(samtim)
                xsat = xyzsat[0]
                ysat = xyzsat[1]
                zsat = xyzsat[2]
                orbtim = samtim  # xht = ne_e("sqrt( ((xsat)**2) +((ysat)**2) + ((zsat)**2))")

            emega = self.emega
            xref = self.xref
            ct = ne_e("cos(emega*samtim + xref)")
            st = ne_e("sin(emega*samtim + xref)")
            # x = ne_e("ct*x1 - st*y1")
            # y = ne_e("st*x1 + ct*y1")
            # vcste1 = ne_e("(ct*x1 - st*y1) - xsat")
            # vcste2 = ne_e("(st*x1 + ct*y1) - ysat")
            # vcste3 = ne_e("z - zsat")
            b31 = self.b31
            b32 = self.b32
            b33 = self.b33
            vcses3 = ne_e("b31*((ct*x1 - st*y1) - xsat) + b32*((st*x1 + ct*y1) - ysat) + b33*(z - zsat)")
            znorm = ne_e("sqrt( (((ct*x1 - st*y1) - xsat)**2) + (((st*x1 + ct*y1) - ysat)**2) + ((z - zsat)**2))")
            # x3 = ne_e("vcses3/znorm")

            rfact = self.rfact
            roasin = self.roasin
            piclin = self.piclin
            radlin = self.radlin
            scan1 = self.scan1
            scan2 = self.scan2
            numsen = self.numsen
            # umv = ne_e("arctan2((vcses3/znorm),sqrt(rfact - ((vcses3/znorm)**2))) - roasin")
            xlin = ne_e("piclin - (arctan2((vcses3/znorm),sqrt(rfact - ((vcses3/znorm)**2))) - roasin)/radlin")
            # xlin = ne_e("piclin - umv/radlin")
            # oldlin = xlin
            del vcses3

            # scnnum = ne_e("((xlin + xlin)/2.0 - 1.0)/numsen")
            # scnfrc = ne_e("((((xlin + xlin)/2.0 - 1.0)/numsen) - scan1)/(scan2 - scan1)")
            xlin = ne_e("xlin + (((((xlin + xlin)/2.0 - 1.0)/numsen) - scan1)/(scan2 - scan1))*(xlin - xlin)")
            time1 = self.time1
            tmpscl = self.tmpscl
            scan1 = self.scan1
            samtim = ne.evaluate("time1 + tmpscl*((((xlin + xlin)/2.0 - 1.0)/numsen) - scan1)")
            xyzsat = self.satvec(samtim)
            xsat = xyzsat[0]
            ysat = xyzsat[1]
            zsat = xyzsat[2]
            cosa = ne_e("(ct*x1 - st*y1)*xsat + (st*x1 + ct*y1)*ysat + z*zsat")
            r = self.r
            rsq = self.rsq
            ctst = ne_e("0.0001*r*(sqrt( ((xsat)**2) +((ysat)**2) + ((zsat)**2))) + rsq")

            del xyzsat

            if np.any(cosa >= ctst):
                del cosa
                del ctst
                b11 = self.b11
                b12 = self.b12
                b13 = self.b13
                b21 = self.b21
                b22 = self.b22
                b23 = self.b23
                b31 = self.b31
                b32 = self.b32
                b33 = self.b33
                # xsats1 = ne_e("b11*xsat + b12*ysat + b13*zsat")
                # ysats2 = ne_e("b21*xsat + b22*ysat + b23*zsat")
                # ct = ne_e("cos(emega*samtim + xref)")
                # st = ne_e("sin(emega*samtim + xref)")
                x = ne_e("(cos(emega*samtim + xref))*x1 - (sin(emega*samtim + xref))*y1")
                y = ne_e("(sin(emega*samtim + xref))*x1 + (cos(emega*samtim + xref))*y1")
                # vcste1 = ne_e("x - xsat")
                # vcste2 = ne_e("y - ysat")
                # vcste3 = ne_e("z - zsat")

                del x1
                del y1

                vcses1 = ne_e("b11*(x - xsat) + b12*(y - ysat) + b13*(z - zsat)")
                vcses2 = ne_e("b21*(x - xsat) + b22*(y - ysat) + b23*(z - zsat)")
                vcses3 = ne_e("b31*(x - xsat) + b32*(y - ysat) + b33*(z - zsat)")

                xnorm = ne_e("sqrt((znorm**2) - (vcses3**2))")
                ynorm = ne_e("sqrt(((b11*xsat + b12*ysat + b13*zsat)**2) + ((b21*xsat + b22*ysat + b23*zsat)**2))")

                znorm = ne_e("sqrt( ((x - xsat)**2) + ((y - ysat)**2) + ((z - zsat)**2))")
                # x3 = vcses3/znorm

                rfact = self.rfact
                roasin = self.roasin
                umv = ne_e("arctan2((vcses3/znorm),sqrt(rfact - ((vcses3/znorm)**2))) - roasin")
                # slin = ne_e("sin(umv)")
                # clin = ne_e("cos(umv)")

                rotm11 = self.rotm11
                rotm13 = self.rotm13
                rotm21 = self.rotm21
                rotm23 = self.rotm23
                # u = ne_e("rotm11*cos(umv) + rotm13*sin(umv)")
                # v = ne_e("rotm21*cos(umv) + rotm23*sin(umv)")

                picele = self.picele
                radele = self.radele
                gamma = self.gamma
                gamdot = self.gamdot
                xele = ne_e(
                    "picele + arcsin(((b11*xsat + b12*ysat + b13*zsat)*vcses2 - (b21*xsat + b22*ysat + b23*zsat)*vcses1)/(xnorm*ynorm))/radele")
                xele = ne_e(
                    "xele + arctan2((rotm21*cos(umv) + rotm23*sin(umv)),(rotm11*cos(umv) + rotm13*sin(umv)))/radele")
                xele = ne_e("xele-gamma-gamdot*samtim")

                del xsat
                del ysat
                del zsat
                del znorm
                del samtim
                del vcses1
                del vcses2
                del vcses3
                del xnorm
                del ynorm
                del umv
                del z

        linele[self.indexLine][point] = xlin
        linele[self.indexEle][point] = xele

        del xlin
        del xele
        # Return in 'File' coordinates
        gc.collect()

        return self.imageCoordToAreaCoord(linele, linele)

    def nllxyz(self, xlat=None, xlon=None):
        '''
        SUBROUTINE NLLXYZ(XLAT,XLON,X,Y,Z)
        CONVERT LAT,LON TO EARTH CENTERED X,Y,Z
        (DALY, 1978)
        XLAT,XLON ARE IN DEGREES, WITH NORTH AND WEST POSITIVE
        X,Y,Z ARE GIVEN IN KM. THEY ARE THE COORDINATES IN A RECTANGULAR
        FRAME WITH ORIGIN AT THE EARTH CENTER, WHOSE POSITIVE
        X-AXIS PIERCES THE EQUATOR AT LON 0 DEG, WHOSE POSITIVE Y-AXIS
        PIERCES THE EQUATOR AT LON 90 DEG, AND WHOSE POSITIVE Z-AXIS
        INTERSECTS THE NORTH POLE.
        '''

        ylat = self.rdpdg * xlat
        ylat1 = np.arctan2(self.bsq * np.sin(ylat), self.asq * np.cos(ylat))
        ylon = -self.rdpdg * xlon
        snlt = np.sin(ylat1)
        cslt = np.cos(ylat1)
        csln = np.cos(ylon)
        snln = np.sin(ylon)
        tnlt = np.power((snlt / cslt), 2)
        r = self.ab * np.sqrt((1.0 + tnlt) / (self.bsq + self.asq * tnlt))
        x = r * cslt * csln
        y = r * cslt * snln
        z = r * snlt
        return x, y, z

    def nllxyz_ne(self, xlat=None, xlon=None):
        '''
        SUBROUTINE NLLXYZ(XLAT,XLON,X,Y,Z)
        CONVERT LAT,LON TO EARTH CENTERED X,Y,Z
        (DALY, 1978)
        XLAT,XLON ARE IN DEGREES, WITH NORTH AND WEST POSITIVE
        X,Y,Z ARE GIVEN IN KM. THEY ARE THE COORDINATES IN A RECTANGULAR
        FRAME WITH ORIGIN AT THE EARTH CENTER, WHOSE POSITIVE
        X-AXIS PIERCES THE EQUATOR AT LON 0 DEG, WHOSE POSITIVE Y-AXIS
        PIERCES THE EQUATOR AT LON 90 DEG, AND WHOSE POSITIVE Z-AXIS
        INTERSECTS THE NORTH POLE.
        '''

        bsq = self.bsq
        asq = self.asq
        ab = self.ab
        rdpdg = self.rdpdg
        # ylat = self.rdpdg*xlat
        ylat1 = ne.evaluate("arctan2(bsq*sin(rdpdg*xlat), asq*cos(rdpdg*xlat))")
        # ylon = -self.rdpdg*xlon
        snlt = ne.evaluate("sin(arctan2(bsq*sin(rdpdg*xlat), asq*cos(rdpdg*xlat)))")
        cslt = ne.evaluate("cos(arctan2(bsq*sin(rdpdg*xlat), asq*cos(rdpdg*xlat)))")
        csln = ne.evaluate("cos(-rdpdg*xlon)")
        snln = ne.evaluate("sin(-rdpdg*xlon)")
        tnlt = ne.evaluate("(snlt/cslt)**2")
        r = ne.evaluate("ab*sqrt((1.0+tnlt)/(bsq+asq*tnlt))")
        x = ne.evaluate("r*cslt*csln")
        y = ne.evaluate("r*cslt*snln")
        z = ne.evaluate("r*snlt")
        return x, y, z

    def getSubpoint(self):

        samtim = self.time1
        xyzsat = self.satvec(samtim)

        ct = math.cos(self.emega * samtim + self.xref)
        st = math.sin(self.emega * samtim + self.xref)
        x = xyzsat[0]
        y = xyzsat[1]
        z = xyzsat[2]
        x1 = ct * x + st * y
        y1 = -st * x + ct * y

        ll = self.nxyzll(x1, y1, z)

        ssp_lat = ll[0]
        ssp_lon = ll[1]
        if self.isEastPositive:
            ssp_lon = -ssp_lon

        return ssp_lat, ssp_lon

    def epoch(self):
        """
        Epoch does time computations. Totally unintuitive I choose to implement it as instance method in spite of the
        fact it does change some class attributes.
        This is because this code is ported form java where it was ported from Fortran which really sucks
        :return:
        """
        # print self.ietimy
        # print self.ietimh
        # print self.semima
        # print self.oeccen
        # print self.xmeana
        PI = math.pi
        RDPDG = PI / 180.
        RE = 6378.388
        GRACON = 0.07436574
        axmmc = GRACON * ((RE / self.semima) ** .5) ** 3
        xmanom = RDPDG * self.xmeana
        time = (xmanom - self.oeccen * math.sin(xmanom)) / (60. * axmmc)
        time1 = utils.mcPI2F(self.ietimh)
        time = time1 - time

        iday = 0
        if time > 48:
            time = - 48
            iday = 2
        elif time > 24:
            time -= 24
            iday = 1
        elif time < -24:
            time += 48
            iday = -2
        elif time < 0:

            time += 24
            iday = -1
        self.ietimh = utils.mcF2PI(time)
        if iday != 0:
            jyear = (self.ietimy / 1000) % 100

            jyear += 1000

            jday = self.ietimy % 1000

            jday += iday
            if jday < 1:
                jyear -= 1
                jday = utils.leapyear(jyear) + jday
            else:
                jtot = utils.leapyear(jyear)

                if jday > jtot:
                    jyear += 1
                    jday += jtot

            jyear = jyear % 100
            self.ietimy = 1000 * jyear + jday

    def icon1(self, yymmdd=None):
        """

        :param yymmdd:
        :return:
        """
        num = 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
        year = (yymmdd / 10000) % 100
        month = (yymmdd / 100) % 100
        day = yymmdd % 100;
        if (month < 0 or month > 12):
            month = 1
        julday = day + num[month - 1]
        if year % 4 == 0 and month > 2:
            julday = julday + 1
        return 1000 * year + julday



class RECTNav(AreaNav):
    labels =  { 	1:'test',
                    2: 'a particular image row number',
                    3: 'latitude corresponding to word 2, degrees x10000',
                    4: 'a particular image column number',
                    5: 'longitude corresponding to word 4, degrees x10000',
                    6:'latitude degrees/image line, degrees x10000',
                    7:'longitude degrees/image line, degrees x10000',
                    8:'radius of the planet, meters',
                    9:'eccentricity of the planet, x1000000',
                    10:'coordinate type, >= 0 planetodetic, < 0 planetocentric',
                    11:'longitude convention, >= 0 west positive, < 0 west negative'
                }

    def __init__(self, navblock=None, auxblock=None, **kwargs):
        inavblock = navblock
        #for n, v in self.labels.items():
        #	print n, v, inavblock[n-1]
        self.itype = 2
        self.iwest = 1 if inavblock[10] >= 0 else -1

        self.xrow = inavblock[1]
        self.xcol = inavblock[3]
        ipowlat = inavblock[11]
        if ipowlat == 0:
            ipowlat = 4
        self.zslat = inavblock[2] / (10. ** ipowlat)
        ipowlon = inavblock[12]
        if ipowlon == 0:
            ipowlon == 4
        self.zslon = inavblock[4] / (10. ** ipowlon)

        ipowdlin = inavblock[13]
        if ipowdlin == 0:
            ipowdlin == 4
        self.zdlat = float(inavblock[5]) / (10. ** ipowdlin)
        ipowdele = inavblock[14]
        if ipowdele == 0:
            ipowdele == 4
        self.zdlon = float(inavblock[6]) / (10 ** ipowdele)
        if self.xcol == 1:
            self.zslon = self.zslon - 180.0 * self.iwest

    def toLinEle(self, latlon=None):

        number = len(latlon[0])
        linele = [[nan for i in range(number)], [nan for i in range(number)]]
        for point in range(number):
            xlat = latlon[self.indexLat][point]
            xlon = -latlon[self.indexLon][point] if self.iwest == 1 else latlon[self.indexLon][point]

            if xlon > (self.zslon+180):
                xlon = xlon - 360
            elif xlon < (self.zslon-180):
                xlon = xlon + 360

            xlin = self.xrow - (xlat - self.zslat)/self.zdlat
            if self.xcol == 1:
                xele = self.xcol - (xlon - self.zslon-180*self.iwest) / (self.zdlon*self.iwest)
            else:
                xele = self.xcol - (xlon - self.zslon)/(self.zdlon*self.iwest)
            linele[self.indexLine][point] = float(xlin)
            linele[self.indexEle][point] = float(xele)
        #return linele

        return self.imageCoordToAreaCoord(linele=linele)

    def toLatLon(self, linele=None, lon_circular=True ):

        number = len(linele[0])
        latlon = [[nan for i in range(number)], [nan for i in range(number)]]


        # Convert array to Image coordinatesfor computations

        imglinele = self.areaCoordtoImageCoord(linele=linele)

        for point in range(number):
            xlin = imglinele[self.indexLine][point]
            xele = imglinele[self.indexEle][point]

            xldif = self.xrow - xlin
            if self.xcol == 1:
                xedif = self.iwest * (xele-self.xcol)
                xlon = self.zslon + 180 * self.iwest - xedif*self.zdlon
            else:
                xedif = self.iwest * (self.xcol - xele)
                xlon = self.zslon + xedif * self.zdlon

            xlat = self.zslat + xldif * self.zdlat
            if  xlat > 90 or xlat < -90:
                xlat = nan

            if xlon > (self.zslon + 180) or xlon <  (self.zslon -180):
                xlon = nan
            if not math.isnan(xlon):
                if lon_circular:
                    if xlon < -180:
                        xlon = xlon + 360
                    if xlon > 180:
                        xlon = xlon - 360

            if math.isnan(xlon) or math.isnan(xlat):
                latlon[self.indexLat][point] = nan
                latlon[self.indexLon][point] = nan
            else:
                latlon[self.indexLat][point] = float(xlat)
                latlon[self.indexLon][point] = -xlon if self.iwest == 1 else xlon

        return latlon








    def __str__(self):
        s = ''
        s+= 'itype %i \n' % self.itype
        s+= 'xrow %i \n' % self.xrow
        s+= 'xcol %i \n' % self.xcol
        s+= 'iwest %i \n' % self.iwest
        s+= 'zslat %.6f \n' % self.zslat
        s+= 'zslon %.6f \n' % self.zslon
        s+= 'zdlat %.6f \n' % self.zdlat
        s+= 'zdlon %.6f \n' % self.zdlon
        return s



class GMSXNav(AreaNav):



    def __init__(self, navblock=None, auxblock=None, **kwargs):
        self.bParms = np.zeros((3200), dtype=np.int8)
        self.subLon = nan
        self.subLat = nan
        self.resLin = np.zeros((4,), dtype=np.float32)
        self.resEle = np.zeros((4,), dtype=np.float32)
        self.rlic = np.zeros((4,), dtype=np.float32)
        self.relmfc = np.zeros((4,), dtype=np.float32)
        self.senssu = np.zeros((4,), dtype=np.float32)
        self.rline = np.zeros((4,), dtype=np.float32)
        self.relem = np.zeros((4,), dtype=np.float32)
        self.vmis = np.zeros((3,), dtype=np.float32)
        self.elmis = np.zeros((3, 3), dtype=np.float32)
        self.dtims = 0.
        self.dspin = 0.
        self.sitagt = 0.
        self.sunalp = 0.
        self.sundel = 0.
        self.sat = [nan, nan, nan]
        self.sp = [nan, nan, nan]
        self.ss = [nan, nan, nan]
        self.orbt1 = np.zeros((35, 8), dtype=np.float64)
        self.atit = np.zeros((10, 10), dtype=np.float64)

        self.cdr = math.pi / 180.
        self.crd = 180. / math.pi
        self.hpai = math.pi / 2.
        self. dpai = math.pi / 2.
        self.ea = 6378136.0
        self.ef = 1 / 298.257




        navtype = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=4, format_str='s'),
                      struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=1, format_str=self.int_format),navblock[self.STTYPE]))[0]
        if navtype != self.GMSX:
            raise Exception(f'Invalid navigation type {navtype}')




        #I will use a completely different approach, numpy based
        #convert the string MORE to a integer
        #excluded_value = struct.unpack('>1i',struct.pack('>4s', 'MORE'))[0]
        excluded_value = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=1, format_str=self.int_format),struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=4, format_str='s'), 'MORE'))[0]

        #convert initial nav block to numpy array
        navblock_array = np.array(navblock, dtype=np.int64)
        #compute mask by excluding elements whose value == 'MORE'
        navblock_mask = navblock_array!=excluded_value

        #filter the initial block end exclude first value (GMSX), then convert to list
        cleaned_navblock = navblock_array[navblock_mask][1:].tolist()



        #convert back to byte string
        fmt = '{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=len(cleaned_navblock), format_str=self.int_format)
        self.strbParms = struct.pack(fmt, *cleaned_navblock)
        fmt1 = '{byte_order:s}{size:d}{format_str:s}'.format(byte_order=self.byte_order, size=len(cleaned_navblock)*4, format_str='b')
        #convert to  bytes
        self.bParms = struct.unpack(fmt1, self.strbParms)

        self.decOABlock(self.bParms, self.strbParms, 0)
        #int to 4 bytes intToBytes
        #print struct.unpack('>4b', struct.pack('>1i', 123277))

    @staticmethod
    def sv01000(byte_str_array=None, mult_pow=None, byte_order=None):
        """
        Converts a byte string array of a certain length to a float number as per
        http://www.ssec.wisc.edu/mcidas/doc/prog_man/current/formats-11.html#25566
        I have not followed a close to Java implementation because it would have been slow
        Rather than that i used numpy and struct
        The original Java code can deal with 4 or 6 bytes, this code can deal with arbitrary lenghts
        The idea ois that given a byte string of a certain length it is converted to a float number
        by anding &  the most significant byte with 127 '01111111' which is ignoring the most significand bit for BIG endian (it seems the GMSX data is in BIG endian)
        na then multypling with  2  at the power corresponding to the index of each byte expressed as the maximum number of bits (see pows var)
        and the summbin this multiplication
        This implementation is byte order specific which is bad but as the data are just 1 byte  it holds!
        "The Type column in the tables below shows scaled integers in the format R*M.N. The R indicates real numbers, M is the number of bytes and N is the exponent."


        :param byte_str_array: a
        :param mult_pow:
        :param byte_order:
        :return:
        """
        size = len(byte_str_array)

        #get unsigned and signed versions

        unsigned_byte_array = np.array(struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=byte_order, size=size,format_str='B'), byte_str_array))
        signed_byte_array = np.array(struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order=byte_order, size=size,format_str='b'), byte_str_array))

        #patch the MSB of the unsigned
        unsigned_byte_array[0] = unsigned_byte_array[0]&127
        #print signed_byte_array
        negative = True if signed_byte_array[0] < 0 else False
        #the power
        pows = 2. ** (np.arange(size)[::-1] * 8)

        sm = sum(unsigned_byte_array*pows)/(10**mult_pow)
        if negative:
            sm= -sm
        return np.asscalar(sm)

    @staticmethod
    def normalize_lon(lon=None, lon_range=180):
        rtn = lon % lon_range
        if rtn <=0:
            rtn+=lon_range
        if rtn> lon_range/2.:
            rtn-=lon_range
        return rtn

    def decOABlock(self, byte_array, byte_string_array, form):
        """
        decode Orbit and Attitude data block
        :param byte_array:
        :param form:
        :return:
        """

        self.dtims = self.sv01000(byte_string_array[:6], 8, self.byte_order)
        #print self.dtims
        self.dspin = self.sv01000(byte_string_array[240:246], 8, self.byte_order)
        #print self.dspin
        self.resLin[0] = self.sv01000(byte_string_array[6:10],8, self.byte_order)
        self.resLin[1] = self.sv01000(byte_string_array[10:14],8, self.byte_order)
        self.resLin[2] = self.sv01000(byte_string_array[10:14],8, self.byte_order)
        self.resLin[3] = self.sv01000(byte_string_array[10:14],8, self.byte_order)
        #print self.resLin
        self.resEle[0] = self.sv01000(byte_string_array[14:18], 10, self.byte_order)
        self.resEle[1] = self.sv01000(byte_string_array[18:22], 10, self.byte_order)
        self.resEle[2] = self.sv01000(byte_string_array[18:22], 10, self.byte_order)
        self.resEle[3] = self.sv01000(byte_string_array[18:22], 10, self.byte_order)
        #print self.resEle
        self.rlic[0] = self.sv01000(byte_string_array[22:26], 4, self.byte_order)
        self.rlic[1] = self.sv01000(byte_string_array[26:30], 4, self.byte_order)
        self.rlic[2] = self.sv01000(byte_string_array[110:114], 4, self.byte_order)
        self.rlic[3] = self.sv01000(byte_string_array[114:118], 4, self.byte_order)
        #print self.rlic
        self.relmfc[0] = self.sv01000(byte_string_array[30:34], 4, self.byte_order)
        self.relmfc[1] = self.sv01000(byte_string_array[34:38], 4, self.byte_order)
        self.relmfc[2] = self.sv01000(byte_string_array[118:122], 4, self.byte_order)
        self.relmfc[3] = self.sv01000(byte_string_array[122:126], 4, self.byte_order)
        #print self.relmfc
        self.senssu[0] = self.sv01000(byte_string_array[38:42], 0, self.byte_order)
        self.senssu[1] = self.sv01000(byte_string_array[42:46], 0, self.byte_order)
        self.senssu[2] = self.sv01000(byte_string_array[42:46], 0, self.byte_order)
        self.senssu[3] = self.sv01000(byte_string_array[42:46], 0, self.byte_order)
        #print self.senssu
        self.rline[0] = self.sv01000(byte_string_array[46:50], 0, self.byte_order)
        self.rline[1] = self.sv01000(byte_string_array[50:54], 0, self.byte_order)
        self.rline[2] = self.sv01000(byte_string_array[50:54], 0, self.byte_order)
        self.rline[3] = self.sv01000(byte_string_array[50:54], 0, self.byte_order)
        #print self.rline
        self.relem[0] = self.sv01000(byte_string_array[54:58], 0, self.byte_order)
        self.relem[1] = self.sv01000(byte_string_array[58:62], 0, self.byte_order)
        self.relem[2] = self.sv01000(byte_string_array[58:62], 0, self.byte_order)
        self.relem[3] = self.sv01000(byte_string_array[58:62], 0, self.byte_order)
        #print self.relem
        self.vmis[0] = self.sv01000(byte_string_array[62:66], 10, self.byte_order)
        self.vmis[1] = self.sv01000(byte_string_array[66:70], 10, self.byte_order)
        self.vmis[2] = self.sv01000(byte_string_array[70:74], 10, self.byte_order)
        #print self.vmis
        self.elmis[0][0] = self.sv01000(byte_string_array[74:78], 7, self.byte_order)
        self.elmis[1][0] = self.sv01000(byte_string_array[78:82], 10, self.byte_order)
        self.elmis[2][0] = self.sv01000(byte_string_array[82:86], 10, self.byte_order)
        self.elmis[0][1] = self.sv01000(byte_string_array[86:90], 10, self.byte_order)
        self.elmis[1][1] = self.sv01000(byte_string_array[90:94], 7, self.byte_order)
        self.elmis[2][1] = self.sv01000(byte_string_array[94:98], 10, self.byte_order)
        self.elmis[0][2] = self.sv01000(byte_string_array[98:102], 10, self.byte_order)
        self.elmis[1][2] = self.sv01000(byte_string_array[102:106], 10, self.byte_order)
        self.elmis[2][2] = self.sv01000(byte_string_array[106:110], 7, self.byte_order)
        #for i in range(3):
        #    print self.elmis[i]
        self.subLon = self.sv01000(byte_string_array[198:204], 6, self.byte_order)
        #print self.subLon
        self.subLat = self.sv01000(byte_string_array[204:210], 6, self.byte_order)
        #print self.subLat
        for i in range(10):
            if form == 1:
                j = (i - 1) * 64 + 256
            if form == 0:
                j = (i - 1) * 48 + 256
            self.atit[0][i] = self.sv01000(byte_string_array[j:j+6], 8, self.byte_order)
            self.atit[2][i] = self.sv01000(byte_string_array[j+12:j+12+6], 8, self.byte_order)
            self.atit[3][i] = self.sv01000(byte_string_array[j+18:j+18+6], 11, self.byte_order)
            self.atit[4][i] = self.sv01000(byte_string_array[j+24:j+24+6], 8, self.byte_order)
            self.atit[5][i] = self.sv01000(byte_string_array[j+30:j+30+6], 8, self.byte_order)

            # if i == 0:
            #     print self.atit[0][i]
            #     print self.atit[2][i]
            #     print self.atit[3][i]
            #     print self.atit[4][i]
            #     print self.atit[5][i]
            #
        for i in range(8):
            if form == 1:
                j = (i-1)*256 + 896
            if form == 0:
                j = (i-1) * 200 + 736
            self.orbt1[0][i] = self.sv01000(byte_string_array[j:j+6], 8, self.byte_order)
            self.orbt1[8][i] = self.sv01000(byte_string_array[j+48:j+48+6], 6, self.byte_order)
            self.orbt1[9][i] = self.sv01000(byte_string_array[j+54:j+54+6], 6, self.byte_order)
            self.orbt1[10][i] = self.sv01000(byte_string_array[j+60:j+60+6], 6, self.byte_order)
            self.orbt1[14][i] = self.sv01000(byte_string_array[j+84:j+84+6], 8, self.byte_order)
            self.orbt1[17][i] = self.sv01000(byte_string_array[j+102:j+102+6], 8, self.byte_order)
            self.orbt1[18][i] = self.sv01000(byte_string_array[j+108:j+108+6], 8, self.byte_order)
            self.orbt1[19][i] = self.sv01000(byte_string_array[j+128:j+128+6], 12, self.byte_order)
            self.orbt1[20][i] = self.sv01000(byte_string_array[j+134:j+134+6], 14, self.byte_order)
            self.orbt1[21][i] = self.sv01000(byte_string_array[j+140:j+140+6], 14, self.byte_order)
            self.orbt1[22][i] = self.sv01000(byte_string_array[j+146:j+146+6], 14, self.byte_order)
            self.orbt1[23][i] = self.sv01000(byte_string_array[j+152:j+152+6], 12, self.byte_order)
            self.orbt1[24][i] = self.sv01000(byte_string_array[j+158:j+158+6], 16, self.byte_order)
            self.orbt1[25][i] = self.sv01000(byte_string_array[j+164:j+164+6], 12, self.byte_order)
            self.orbt1[26][i] = self.sv01000(byte_string_array[j+170:j+170+6], 16, self.byte_order)
            self.orbt1[27][i] = self.sv01000(byte_string_array[j+176:j+176+6], 12, self.byte_order)


            # print '-'*50
            # for mmm in (0, 8, 9, 10, 14, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27):
            #     print 'orbt1[%d][%d] %.5f %s ' % (mmm, i, self.orbt1[mmm][i], type(self.orbt1[mmm][i]))
            #

    def linele_scantime(self, linele=None):
        mode = -1
        imglinele = self.areaCoordtoImageCoord_np(linele=linele)

        rPix = imglinele[1]
        rLin = imglinele[0]

        shp = rPix.shape
        # cmask = np.ones(lat.shape).astype(np.bool)
        rPix = rPix.ravel()
        rLin = rLin.ravel()
        # vissr frame
        lMode = abs(mode) - 1

        rsamp = self.resEle[lMode]

        sens = self.senssu[lMode]

        rtim = ((self.rint_np(rLin - 1) / sens) + (rPix * rsamp) / self.dpai) / (self.dspin * 1440.) + self.dtims
        return rtim.reshape(shp)




    def latlon_scantime(self, latlon=None):
        mode = 1
        lon = ((latlon[self.indexLon] + 180.0) % 360.0) - 180.0
        lat = latlon[self.indexLat]
        sz = lat.size
        shp = lon.shape


        lon = lon.ravel()
        lat = lat.ravel()

        # now allocate matrices for the computataion dependant attributes cecessary in array mode
        # these attributes are scalars in point mode but they become array in array mode
        self.sat = np.zeros((3, sz))
        self.sp = np.zeros((3, sz))
        self.ss = np.zeros((3, sz))
        self.sitagt = np.zeros((sz))
        self.sunalp = np.zeros((sz))
        self.sundel = np.zeros((sz))
        tri = np.zeros((sz))
        trj = np.zeros((sz))

        # cmask = np.ones(lat.shape).astype(np.bool)

        stn1 = [nan, nan, nan]


        m = (lat > 90) | (lat < -90)
        if m[m].size > 0:
            return
        # vissr frame
        lMode = abs(mode) - 1
        rstep = self.resLin[lMode]
        rfcl = self.rlic[lMode]
        sens = self.senssu[lMode]
        dLat = lat * self.cdr
        dLon = lon * self.cdr
        ee = 2. * self.ef - (self.ef * self.ef)

        sin_lat = np.sin(dLat)
        cos_lat = np.cos(dLat)
        sin_lon = np.sin(dLon)
        cos_lon = np.cos(dLon)
        en = (self.ea / np.sqrt(1. - ee * sin_lat * sin_lat))

        stn1[0] = en * cos_lat * cos_lon
        stn1[1] = en * cos_lat * sin_lon
        stn1[2] = en * (1. - ee) * sin_lat
        rio = rfcl - np.arctan(sin_lat / (6.610689 - cos_lat)) / rstep  # i guess 6.1 ms is the scan ytime for a swath, it makes sens if I look to this line of code

        rtim = self.dtims + (rio / sens / 1440.) / self.dspin
        return rtim.reshape(shp)

    def compute_orbital_zones(self, rscantime=None):
        _shp = rscantime.shape
        rscantime = rscantime.ravel()
        orbt_values = self.orbt1[0]
        orbt_zones = np.digitize(rscantime, orbt_values) - 1  # ?? should i extract one
        return orbt_zones.reshape(_shp)
        # unique_orb_zones = np.unique(orbt_zones)
        # for orbt_zone in unique_orb_zones:
        #     m = orbt_zones == orbt_zone

    def toLinEle2(self, latlon=None):
        #print 'indexLat', self.indexLat
        #print 'indexLoin', self.indexLon

        mode = 1
        lon = ((latlon[self.indexLon] + 180.0) % 360.0) - 180.0
        lat = latlon[self.indexLat]
        sz = lat.size
        shp = lon.shape



        lon = lon.ravel()
        lat = lat.ravel()


        #now allocate matrices for the computataion dependant attributes cecessary in array mode
        #these attributes are scalars in point mode but they become array in array mode
        self.sat = np.zeros((3, sz))
        self.sp = np.zeros((3, sz))
        self.ss = np.zeros((3, sz))
        self.sitagt = np.zeros((sz))
        self.sunalp = np.zeros((sz))
        self.sundel = np.zeros((sz))
        tri = np.zeros((sz))
        trj = np.zeros((sz))


        #cmask = np.ones(lat.shape).astype(np.bool)

        stn1 = np.zeros((3, sz))

        slv = np.zeros((3, sz))
        sw3 = np.zeros((3, sz))
        m = (lat > 90) | (lat <-90)
        if m[m].size> 0:
            return None, None



        # vissr frame
        lMode = abs(mode) - 1

        rstep = self.resLin[lMode]
        rsamp = self.resEle[lMode]
        rfcl = self.rlic[lMode]

        rfcp = self.relmfc[lMode]

        sens = self.senssu[lMode]
        rftl = self.rline[lMode] + .5
        rftp = self.relem[lMode] + .5

        dLat = lat * self.cdr
        dLon = lon * self.cdr
        ee = 2. * self.ef - (self.ef * self.ef)

        sin_lat = np.sin(dLat)
        cos_lat = np.cos(dLat)
        sin_lon = np.sin(dLon)
        cos_lon = np.cos(dLon)
        en = (self.ea / np.sqrt(1. - ee * sin_lat * sin_lat))

        stn1[0] = en * cos_lat * cos_lon
        stn1[1] = en * cos_lat * sin_lon
        stn1[2] = en * (1. - ee) * sin_lat
        rio = rfcl - np.arctan(sin_lat / (6.610689 - cos_lat)) / rstep  # i guess 6.1 ms is the scan ytime for a swath, it makes sens if I look to this line of code

        rtim = self.dtims + (rio / sens / 1440.) / self.dspin



        gg = 0
        while True:
            #print 'iter', gg, 'rio size', rio.size
            # if gg == 0:
			#
            #     print 'sp', self.sp
            #     print 'ss', self.ss
            #     print 'sat', self.sat
            beta = self.mg1100_np(rtim)





            # if gg == 0:
            #     print 'beta', beta, beta.dtype
            #     print 'sp', self.sp
            #     print 'ss', self.ss
            #     print 'sat', self.sat

            sw1 = self.mg1220_np(self.sp, self.ss)

            sw2 = self.mg1220_np(sw1, self.sp)

            bc = np.cos(beta)
            bs = np.sin(beta)
            sw3[0] = (sw1[0] * bs) + (sw2[0] * bc)
            sw3[1] = (sw1[1] * bs) + (sw2[1] * bc)
            sw3[2] = (sw1[2] * bs) + (sw2[2] * bc)
            sx = self.mg1200_np(sw3)
            sy = self.mg1220_np(self.sp, sx)

            # try:
			#
            #     slv[0] = stn1[0][cmask] - self.sat[0][cmask]
            #     slv[1] = stn1[1][cmask] - self.sat[1][cmask]
            #     slv[2] = stn1[2][cmask] - self.sat[2][cmask]
            # except UnboundLocalError:
            #     slv[0] = stn1[0] - self.sat[0]
            #     slv[1] = stn1[1] - self.sat[1]
            #     slv[2] = stn1[2] - self.sat[2]
            slv[0] = stn1[0] - self.sat[0]
            slv[1] = stn1[1] - self.sat[1]
            slv[2] = stn1[2] - self.sat[2]
            sl = self.mg1200_np(slv)

            sw2 = self.mg1210(self.sp, sl)
            sw3 = self.mg1210(sy, sw2)
            tp = self.mg1230_np(sy, sw2)
            tf = (self.sp[0] * sw3[0]) + (self.sp[1] * sw3[1]) + (self.sp[2] * sw3[2])

            if np.all(tf < 0):
                tp = -tp

            tl = self.mg1230_np(self.sp, sl)




            ri = (self.hpai - tl) / rstep + rfcl - (self.vmis[1] / rstep)


            rj = tp / rsamp + rfcp + (self.vmis[2] / rsamp) - (self.hpai - tl) * np.tan(self.vmis[0]) / rsamp


            #now this is a bit tricky as I ned to update the cmask in order to correctly emulate the original while loop
            gg += 1
            diff = np.abs(ri-rio)


            cmask = diff >=1
            # from pylab import show, imshow, title, colorbar
            # imshow(cmask.reshape((shp)))
            # title('cmask')
            # colorbar()
            # show()
            # imshow(rio.reshape((shp)))
            # title('rio')
            # colorbar()
            # show()
            #print 'diff', diff.mean()


            if cmask[cmask].size > 0:
                print ('masking %s for %s in percent %s' % (cmask[cmask].size, lat.size, cmask[cmask].size/float(lat.size) * 100))

                # from pylab import show, imshow, title, colorbar
                # imshow(cmask.reshape((3900, 2700)))
                # colorbar()
                # title('cmask%s' % gg)
                # show()
                #if abs(ri - rio) >= 1.:
                #mask the rtim adn rio, i am not sure if this is enough thoufh
                #rtim = ((self.rint_np((ri - 1) / sens) + (rj * rsamp) / self.dpai) / (self.dspin * 1440.0) + self.dtims)[cmask]
                rtim1 = (self.rint_np((ri - 1) / sens) + rj * rsamp / self.dpai) / (self.dspin * 1440.0) + self.dtims

                rtim[cmask] = rtim1[cmask]

                rio[cmask] = ri[cmask]
                #rio = ri[cmask]
                continue

            break

        #print 'end rio size', rio.size
        #sanity check for bounds
        if np.any((ri < 0) | (ri > rftl)):
            rc = 4
            print ('inv line')
            return None, None

        if np.any((rj < 0) | (rj > rftp)):
            rc = 5
            print ('inv ele')
            return None, None


        rri = ri.reshape(shp)
        rrj = rj.reshape(shp)

        return self.imageCoordToAreaCoord_np(linele=[rri, rrj])

    def toLinEle(self, latlon=None):
        mode = 1
        count = len(latlon[0])
        linele = [[nan for i in range(count)], [nan for i in range(count)]]

        for point in range(count):
            #init as not navigable
            linele[self.indexLine][point] = nan
            linele[self.indexEle][point] = nan
            if latlon[self.indexLat][point] > 90:
                continue
            #normalize
            lon =  ((latlon[self.indexLon][point] + 180.0) % 360.0) - 180.0
            lat = latlon[self.indexLat][point]

            rtnPoint = self.msivsr(mode,0,0,lat, lon)

            linele[self.indexLine][point] = rtnPoint[0]
            linele[self.indexEle][point] = rtnPoint[1]

        return self.imageCoordToAreaCoord(linele=linele)

    def toLatLon(self, linele=None):
        mode = -1
        count = len(linele[0])
        latlon = [[nan ], [nan]]
        #transform file to image coordinates
        imgLinEle = self.areaCoordtoImageCoord(linele=linele)

        #print imgLinEle
        for point in range(count):
            rtnPoint = self.msivsr(mode, imgLinEle[self.indexEle][point], imgLinEle[self.indexLine][point], 0,0,)
            #print rtnPoint
            lon = rtnPoint[1]
            lon = ((lon + 180.0) % 360.0) - 180.0
            if abs(lon) > 180 or (-lon > 90-self.subLon and -lon < 270 - self.subLon):
                return [nan, nan]

            latlon[self.indexLat][point] = rtnPoint[0]
            latlon[self.indexLon][point] = lon
        return latlon

    def toLatLon2(self, linele=None, ):

        mode = -1
        imglinele = self.areaCoordtoImageCoord_np(linele=linele)

        rPix = imglinele[1]
        rLin = imglinele[0]


        shp = rPix.shape
        # cmask = np.ones(lat.shape).astype(np.bool)
        rPix = rPix.ravel()
        rLin = rLin.ravel()

        stn1 = [nan, nan, nan]
        sw3 = [nan, nan, nan]


        # vissr frame
        lMode = abs(mode) - 1
        rstep = self.resLin[lMode]
        rsamp = self.resEle[lMode]
        rfcl = self.rlic[lMode]
        rfcp = self.relmfc[lMode]
        sens = self.senssu[lMode]
        #print rLin, sens, rPix, rsamp, self.dpai, self.dspin, self.dtims

        rtim = ((self.rint_np(rLin - 1) / sens) + (rPix * rsamp) / self.dpai) / (self.dspin * 1440.) + self.dtims
        #print 'rtim', rtim
        beta = self.mg1100_np(rtim)
        # print 'beta: ', beta
        # print 'sp:', self.sp
        # print 'ss:', self.ss
        sw1 = self.mg1220_np(self.sp, self.ss)
        sw2 = self.mg1220_np(sw1, self.sp)

        bc = np.cos(beta)
        bs = np.sin(beta)
        sw3[0] = (sw1[0] * bs) + (sw2[0] * bc)
        sw3[1] = (sw1[1] * bs) + (sw2[1] * bc)
        sw3[2] = (sw1[2] * bs) + (sw2[2] * bc)
        # print 'sw1:', sw1
        # print 'sw2:', sw2
        # print 'sw3', sw3
        sx = self.mg1200_np(sw3)
        sy = self.mg1220_np(self.sp, sx)
        # print 'sx:', sx
        # print 'sy:', sy
        pc = np.cos(rstep * (rLin - rfcl))
        ps = np.sin(rstep * (rLin - rfcl))
        qc = np.cos(rsamp * (rPix - rfcp))
        qs = np.sin(rsamp * (rPix - rfcp))
        sw1[0] = (self.elmis[0][0] * pc) + (self.elmis[0][2] * ps)
        sw1[1] = (self.elmis[1][0] * pc) + (self.elmis[1][2] * ps)
        sw1[2] = (self.elmis[2][0] * pc) + (self.elmis[2][2] * ps)
        # print sw1
        sw2[0] = (qc * sw1[0]) - (qs * sw1[1])
        sw2[1] = (qs * sw1[0]) + (qc * sw1[1])
        sw2[2] = sw1[2]
        # print sw2
        sw3[0] = (sx[0] * sw2[0]) + (sy[0] * sw2[1]) + (self.sp[0] * sw2[2])
        sw3[1] = (sx[1] * sw2[0]) + (sy[1] * sw2[1]) + (self.sp[1] * sw2[2])
        sw3[2] = (sx[2] * sw2[0]) + (sy[2] * sw2[1]) + (self.sp[2] * sw2[2])

        # print 'sw1:', sw1
        # print 'sw2:', sw2
        # print 'sw3', sw3

        # print sw3
        sl = self.mg1200_np(sw3)
        # print 'sl:', sl
        # print 'sat %.5f, %.5f, %.5f:' % (self.sat[0], self.sat[1], self.sat[2])
        def1 = (1. - self.ef) * (1 - self.ef)
        # print "def %.9f" % def1
        dda = def1 * ((sl[0] * sl[0]) + (sl[1] * sl[1])) + (sl[2] * sl[2])
        ddb = def1 * ((self.sat[0] * sl[0]) + (self.sat[1] * sl[1])) + (self.sat[2] * sl[2])
        # print def1, self.sat[0], self.sat[1], self.sat[2], self.ea
        ddc = def1 * ((self.sat[0] * self.sat[0]) + (self.sat[1] * self.sat[1]) - (self.ea * self.ea)) + (self.sat[2] * self.sat[2])
        dd = (ddb * ddb) - (dda * ddc)

        # print 'dda %.5f' % dda
        # print 'ddb %.5f' % ddb
        # print 'ddc %.5f' % ddc
        # print 'dd %.5f' % dd

        # if dd >= 0 and dda != 0:
        #     dk1 = (-ddb + np.sqrt(dd)) / dda
        #     dk2 = (- ddb - np.sqrt(dd)) / dda
        # else:
        #     rc = 6
        #     return (None, None)

        #this code is equivalent to the above  commented python code
        ddm = (dd >= 0) & (dda != 0)
        # now we need to be careful here as we will return invalid data
        dk1 = np.where(ddm, (-ddb + np.sqrt(dd)) / dda, nan)
        dk2 = np.where(ddm, (-ddb - np.sqrt(dd)) / dda, nan)


        # if abs(dk1) < abs(dk2):
        #     dk = dk1
        # else:
        #     dk = dk2
        # print 'dk %.5f' % dk
        abscondm = np.abs(dk1) < np.abs(dk2)
        dk = np.where(abscondm, dk1, dk2)
        stn1[0] = self.sat[0] + (dk * sl[0])
        stn1[1] = self.sat[1] + (dk * sl[1])
        stn1[2] = self.sat[2] + (dk * sl[2])
        # print 'stn1', stn1
        dLat = np.arctan(stn1[2] / (def1 * np.sqrt((stn1[0] * stn1[0]) + (stn1[1] * stn1[1]))))
        # print dLat




        # if stn1[0] != 0:
        #     dLon = np.arctan(stn1[1] / stn1[0])
        #
        #     if stn1[0] < 0 and stn1[1] >= 0:
        #         dLon = dLon + self.PI
        #     if stn1[0] < 0 and stn1[1] < 0:
        #         dLon = dLon - self.PI
        # else:
        #     if stn1[1] > 0:
        #         dLon = self.hpai
        #     else:
        #         dLon = -self.hpai

        #this code is equivalent to above python code
        stn1_0_m = stn1[0] != 0
        stn1_1_m = stn1[1] > 0
        dLon = np.arctan(stn1[1] / stn1[0])
        acon = (stn1[0] < 0) & (stn1[1] >= 0)[stn1_0_m]
        bcon = (stn1[0] < 0) & (stn1[1] < 0)[stn1_0_m]
        if np.any(acon):
            dLon[acon] += self.PI
        if np.any(bcon):
            dLon[bcon] -= self.PI
        if np.any(~stn1_0_m):
            dLon[~stn1_0_m] = np.where(stn1_1_m, self.hpai, -self.hpai)




        lat = dLat * self.crd
        lon = dLon * self.crd

        #normalize
        lon = ((lon + 180.0) % 360.0) - 180.0
        return [lat.reshape(shp), lon.reshape(shp)]


    def msivsr(self, iMode, rPix, rLin, rLat, rLon):
        """
        actual conversion to/from lat/lon or line/elem
        @param iMode - conversion mode, to lat/lon or to line/elem
        @param rPix - float pixel or element value
        @param rLin - float line value
        @param rLat - float latitude value
        @param rLon - float longitude value

        @return array of two floating point values, lat/lon or line/elem pair
        """
        name = 'msivsr'
        point = [nan, nan]
        stn1 = [nan, nan, nan]
        slv = [nan, nan, nan]
        sw3  = [nan, nan, nan]
        #sanity check
        if iMode > 4:
            return point
        if np.any(np.abs(rLat)>90) and iMode > 0:
            return point
        #vissr frame
        lMode = abs(iMode) -1
        rstep = self.resLin[lMode]
        rsamp = self.resEle[lMode]
        rfcl = self.rlic[lMode]
        rfcp = self.relmfc[lMode]
        sens = self.senssu[lMode]
        rftl = self.rline[lMode] + .5
        rftp = self.relem[lMode] + .5



        if iMode > 0: #transformation, geographical -> VISSR
            dLat = rLat * self.cdr
            dLon = rLon * self.cdr
            ee = 2.*self.ef - (self.ef * self.ef)
            sin_lat = np.sin(dLat)
            cos_lat = np.cos(dLat)
            sin_lon = np.sin(dLon)
            cos_lon = np.cos(dLon)
            en = self.ea / np.sqrt(1. - ee*sin_lat*sin_lat)
            stn1[0] = en * cos_lat * cos_lon
            stn1[1] = en * cos_lat * sin_lon
            stn1[2] = en * (1. - ee) * sin_lat
            rio = rfcl-np.arctan(sin_lat / (6.610689 - cos_lat)) / rstep # i guess 6.1 ms is the scan ytime for a swath, it makes sens if I look to this line of code
            rtim = self.dtims + (rio/sens/1440.)/ self.dspin

            gg = 0
            while True:
                #print 'mg1100'
                #print 'classic iter ', gg
                # if gg == 0:
                #     print 'sp', self.sp
                #     print 'ss', self.ss
                #     print 'sat', self.sat

                beta = self.mg1100(rtim)
                # if gg == 0:

                #     print 'sp', self.sp
                #     print 'ss', self.ss
                #     print 'sat', self.sat

                sw1 = self.mg1220(self.sp, self.ss)

                sw2 = self.mg1220(sw1, self.sp)

                bc = np.cos(beta)
                bs = np.sin(beta)
                sw3[0] = (sw1[0] * bs) + (sw2[0] * bc)
                sw3[1] = (sw1[1] * bs) + (sw2[1] * bc)
                sw3[2] = (sw1[2] * bs) + (sw2[2] * bc)
                sx = self.mg1200(sw3)
                sy = self.mg1220(self.sp, sx)

                slv[0] = stn1[0] - self.sat[0]
                slv[1] = stn1[1] - self.sat[1]
                slv[2] = stn1[2] - self.sat[2]

                sl = self.mg1200(slv)

                sw2 = self.mg1210(self.sp, sl)
                sw3 = self.mg1210(sy, sw2)
                tp = self.mg1230(sy, sw2)
                tf = (self.sp[0] * sw3[0]) + (self.sp[1] * sw3[1]) + (self.sp[2] * sw3[2])

                if tf < 0:
                    tp = -tp


                tl = self.mg1230(self.sp, sl)

                ri = float(self.hpai - tl) / rstep + rfcl - (self.vmis[1] / rstep)
                rj = float(tp / rsamp + rfcp + (self.vmis[2] / rsamp) - (self.hpai - tl) * np.tan(self.vmis[0]) / rsamp)

                gg+=1
                if abs(ri - rio) >= 1.:
                    rtim = (self.rint((ri - 1) / sens) + (rj * rsamp) / self.dpai) / (self.dspin * 1440.0) + self.dtims
                    rio = ri
                    continue
                break

            rLin = ri
            rPix = rj

            if rLin > 0 and rLin < rftl:
                if rPix > 0 and rPix < rftp:
                    point[0] = rLin
                    point[1] = rPix

                else:
                    rc = 5
            else:
                rc = 4

        if iMode < 0: #transformation, VISSR -> geographical
            #print rLin, sens, rPix, rsamp, self.dpai, self.dspin, self.dtims
            rtim = ((self.rint(rLin-1) / sens) + (rPix*rsamp) / self.dpai) / (self.dspin*1440.) + self.dtims
            #print 'rtim: ', rtim
            beta = self.mg1100(rtim)
            #print 'beta: ', beta
            #print 'sp:', self.sp
            #print 'ss:', self.ss
            sw1 = self.mg1220(self.sp, self.ss)
            sw2 = self.mg1220(sw1, self.sp)

            bc = np.cos(beta)
            bs = np.sin(beta)
            sw3[0] = (sw1[0] * bs) + (sw2[0] * bc)
            sw3[1] = (sw1[1] * bs) + (sw2[1] * bc)
            sw3[2] = (sw1[2] * bs) + (sw2[2] * bc)
            #print 'sw1:', sw1
            #print 'sw2:', sw2
            #print 'sw3', sw3
            sx = self.mg1200(sw3)
            sy = self.mg1220(self.sp, sx)
            #print 'sx:', sx
            #print 'sy:', sy
            pc = np.cos(rstep * (rLin-rfcl))
            ps = np.sin(rstep * (rLin-rfcl))
            qc = np.cos(rsamp * (rPix-rfcp))
            qs = np.sin(rsamp * (rPix-rfcp))
            sw1[0] = (self.elmis[0][0] * pc) + (self.elmis[0][2] * ps)
            sw1[1] = (self.elmis[1][0] * pc) + (self.elmis[1][2] * ps)
            sw1[2] = (self.elmis[2][0] * pc) + (self.elmis[2][2] * ps)
            #print sw1
            sw2[0] = (qc * sw1[0]) - (qs * sw1[1])
            sw2[1] = (qs * sw1[0]) + (qc * sw1[1])
            sw2[2] = sw1[2]
            #print sw2
            sw3[0] = (sx[0] * sw2[0]) + (sy[0] * sw2[1]) + (self.sp[0] * sw2[2])
            sw3[1] = (sx[1] * sw2[0]) + (sy[1] * sw2[1]) + (self.sp[1] * sw2[2])
            sw3[2] = (sx[2] * sw2[0]) + (sy[2] * sw2[1]) + (self.sp[2] * sw2[2])

            # print 'sw1:', sw1
            # print 'sw2:', sw2
            # print 'sw3', sw3

            #print sw3
            sl = self.mg1200(sw3)
            #print 'sl:', sl
            #print 'sat %.5f, %.5f, %.5f:' % (self.sat[0], self.sat[1], self.sat[2])
            def1 = (1.- self.ef) * (1-self.ef)
            #print "def %.9f" % def1
            dda = def1 * ((sl[0] * sl[0]) + (sl[1] * sl[1])) + (sl[2] * sl[2])
            ddb = def1 * ((self.sat[0] * sl[0]) + (self.sat[1] * sl[1])) + (self.sat[2] * sl[2])
            #print def1, self.sat[0], self.sat[1], self.sat[2], self.ea
            ddc = def1 * ((self.sat[0] * self.sat[0]) + (self.sat[1] * self.sat[1]) - (self.ea * self.ea)) + (self.sat[2] * self.sat[2])
            dd = (ddb * ddb) - (dda*ddc)

            # print 'dda %.5f' % dda
            # print 'ddb %.5f' % ddb
            # print 'ddc %.5f' % ddc
            #print 'dd %.5f' % dd

            if dd >= 0 and dda != 0:
                dk1 = (-ddb + np.sqrt(dd))/ dda
                dk2 =  (- ddb - np.sqrt(dd)) / dda
            else:
                rc= 6
                return point

            if abs(dk1) < abs(dk2):
                dk = dk1
            else:
                dk = dk2
            #print 'dk %.5f' % dk
            stn1[0] = self.sat[0] + (dk * sl[0])
            stn1[1] = self.sat[1] + (dk * sl[1])
            stn1[2] = self.sat[2] + (dk * sl[2])
            #print 'stn1', stn1
            dLat = np.arctan(stn1[2] / (def1 * np.sqrt((stn1[0]* stn1[0]) + (stn1[1]*stn1[1])) ))
            #print dLat
            if stn1[0] != 0:
                dLon = math.atan(stn1[1]/stn1[0])
                if stn1[0] < 0 and stn1[1] >= 0:
                    dLon = dLon + self.PI
                if stn1[0] < 0 and stn1[1] < 0:
                    dLon = dLon - self.PI
            else:
                if stn1[1] > 0:
                    dLon = self.hpai
                else:
                    dLon = -self.hpai

            point[0]= dLat * self.crd
            point[1]= dLon * self.crd

        return point

    @staticmethod
    def rint(num):
        """Rounds toward the even number if equidistant"""
        return round(num + (num % 2 - 1 if (num % 1 == 0.5) else 0))

    @staticmethod
    def rint_np(num):
        """Rounds toward the even number if equidistant"""

        t = np.where(num % 1 == 0.5, num % 2 - 1 , 0)
        return np.round(num + t)

    def mg1110_np(self, i, rtim, orbta, mask):

        npa = np.zeros((3, 3))

        #print i, rtim, (rtim-orbta[0][0:6]) / ( orbta[0][1:7]-orbta[0][0:6]).tolist()

        if i != 7:

            # print 'i sat %d' % i
            # print 'rtim %.5f' % rtim
            # print 'orbta[0][%d] %.8f , orbta[0][%d+1] %.8f' % ( i, orbta[0][i], i, orbta[0][i+1])

            delt = (rtim - orbta[0][i]) / (orbta[0][i + 1] - orbta[0][i])


            # delt = np.where(mask, delt, 0)

            # print 'delt', delt

            self.sat[0] = orbta[8][i] + (orbta[8][i + 1] - orbta[8][i]) * delt
            self.sat[1] = orbta[9][i] + (orbta[9][i + 1] - orbta[9][i]) * delt
            self.sat[2] = orbta[10][i] + (orbta[10][i + 1] - orbta[10][i]) * delt


            # print 'setting sat %.5f, %.5f, %.5f:' %  (self.sat[0], self.sat[1], self.sat[2])
            self.sitagt = (orbta[14][i] + (orbta[14][i + 1] - orbta[14][i]) * delt) * self.cdr
            if orbta[14][i + 1] - orbta[14][i] < 0:
                self.sitagt = (orbta[14][i] + (orbta[14][i + 1] - orbta[14][i] + 360.0) * delt) * self.cdr

            self.sunalp = (orbta[17][i] + (orbta[17][i + 1] - orbta[17][i]) * delt) * self.cdr

            if orbta[17][i + 1] - orbta[17][i] > 0:
                self.sunalp = (orbta[17][i] + (orbta[17][i + 1] - orbta[17][i] - 360.0) * delt) * self.cdr

            self.sundel = (orbta[18][i] + (orbta[18][i + 1] - orbta[18][i]) * delt) * self.cdr

            # compute npa for every point/coodinate?

            npa[0][0] = orbta[19][i]
            npa[1][0] = orbta[20][i]
            npa[2][0] = orbta[21][i]
            npa[0][1] = orbta[22][i]
            npa[1][1] = orbta[23][i]
            npa[2][1] = orbta[24][i]
            npa[0][2] = orbta[25][i]
            npa[1][2] = orbta[26][i]
            npa[2][2] = orbta[27][i]
        return npa

    def mg1100_np(self, rtim):
        """
        I have no clue , but it seem to do some conversion. Im mcidas check the file gms5_nav.for
        for further details.
        :param rtim:
        :return:
        """
        sz = rtim.size

        #beta = delta = attalp = attdel = wkcos = wksin = 0.
        beta = np.zeros((sz))
        delta_orbt = np.zeros((sz))
        delta_atit = np.zeros((sz))
        attalp = np.zeros((sz))
        attdel = np.zeros((sz))
        #wkcos = np.zeros((sz))
        #wksin = np.zeros((sz))

        #att1 = [nan, nan, nan]
        att1 = np.zeros((3, sz))
        #att2 = [nan, nan, nan]
        att2 = np.zeros((3, sz))
        #att3 = [nan, nan, nan]
        att3 = np.zeros((3,sz))
        # npa = [[nan for _i in range(3)] for _ii in range(3)]
        #npa = np.zeros((3, 3))


        """
            This is where the array based code is pretty much different from original code. The backbone of the navigation are the orbital
            and attitude coefficients. One has to be aware that as the scan pattern (swaths) plays a central role in determining the correct navigation

            While in single point mode the computation happens in a loop that exits when the orbit prediction sub-block for the specific pixel is determined.

            There are 8 orbit prediction sub blocks.
        """



        # step 2 determine and calcultae attitude corrections
        atit_values = self.atit[0]
        atit_zones = np.digitize(rtim, atit_values) - 1  # ?? should i extract one
        unique_atit_zones = np.unique(atit_zones)

        for atit_zone in unique_atit_zones:
            m_ = atit_zones == atit_zone
            delta_atit[m_] = (rtim[m_] - atit_values[atit_zone]) / (atit_values[atit_zone + 1] - atit_values[atit_zone])

            attalp[m_] = self.atit[2][atit_zone] + (self.atit[2][atit_zone + 1] - self.atit[2][atit_zone]) * delta_atit[m_]
            attdel[m_] = self.atit[3][atit_zone] + (self.atit[3][atit_zone + 1] - self.atit[3][atit_zone]) * delta_atit[m_]
            beta[m_] = self.atit[4][atit_zone] + (self.atit[4][atit_zone + 1] - self.atit[4][atit_zone]) * delta_atit[m_]
            if self.atit[4][atit_zone + 1] - self.atit[4][atit_zone] > 0:
                beta[m_] = self.atit[4][atit_zone] + (self.atit[4][atit_zone + 1] - self.atit[4][atit_zone] - 360.0 * self.cdr) * delta_atit[m_]

        wkcos = np.cos(attdel)
        att1[0] = np.sin(attdel)
        att1[1] = wkcos * -np.sin(attalp)
        att1[2] = wkcos * np.cos(attalp)
        # from pylab import show, imshow, title, colorbar
        # imshow(attdel.reshape((4800, 2400)))
        # colorbar()
        # title('attdel')
        # show()
        # print wkcos, att1
        '''
		for i in range(9):
			# print 'atit', i, rtim, self.atit[0][i], self.atit[0][i + 1]
			m__ = (rtim > self.atit[0][i]) & (rtim < self.atit[0][i + 1])
			mrtim = rtim[m__]

			if np.all(m__) :

				if mrtim.size != rtim.size:  # what happens here?? the npa seems to be som
					print 'Inconsistent1 results are possible! mrtim.size %s rtim.size %s' % (mrtim.size, rtim.size)

				delt_orbt = (rtim - self.atit[0][i]) / (self.atit[0][i + 1] - self.atit[0][i])
				# delt = np.where(~m__, delt,0)
				attalp = self.atit[2][i] + (self.atit[2][i + 1] - self.atit[2][i]) * delt_orbt
				attdel = self.atit[3][i] + (self.atit[3][i + 1] - self.atit[3][i]) * delt_orbt
				beta = self.atit[4][i] + (self.atit[4][i + 1] - self.atit[4][i]) * delt_orbt
				if self.atit[4][i + 1] - self.atit[4][i] > 0:
					beta = self.atit[4][i] + (self.atit[4][i + 1] - self.atit[4][i] - 360.0 * self.cdr) * delt_orbt
				break
		'''

        # step 1 determine the orbit subblock fir each element in he rtims
        orbt_values = self.orbt1[0]
        orbt_zones = np.digitize(rtim, orbt_values) - 1  # minus one because digitize gives you zones +1
        unique_orb_zones = np.unique(orbt_zones)

        for orbt_zone in unique_orb_zones:
            m = orbt_zones == orbt_zone
            delta_orbt[m] = (rtim[m] - orbt_values[orbt_zone]) / (orbt_values[orbt_zone + 1] - orbt_values[orbt_zone])

            # sate positions
            self.sat[0][m] = self.orbt1[8][orbt_zone] + (self.orbt1[8][orbt_zone + 1] - self.orbt1[8][orbt_zone]) * delta_orbt[m]
            self.sat[1][m] = self.orbt1[9][orbt_zone] + (self.orbt1[9][orbt_zone + 1] - self.orbt1[9][orbt_zone]) * delta_orbt[m]
            self.sat[2][m] = self.orbt1[10][orbt_zone] + (self.orbt1[10][orbt_zone + 1] - self.orbt1[10][orbt_zone]) * delta_orbt[m]

            self.sitagt[m] = (self.orbt1[14][orbt_zone] + (self.orbt1[14][orbt_zone + 1] - self.orbt1[14][orbt_zone]) * delta_orbt[m]) * self.cdr
            if self.orbt1[14][orbt_zone + 1] - self.orbt1[14][orbt_zone] < 0:
                self.sitagt[m] = (self.orbt1[14][orbt_zone] + (self.orbt1[14][orbt_zone + 1] - self.orbt1[14][orbt_zone] + 360.0) * delta_orbt[m]) * self.cdr

            self.sunalp[m] = (self.orbt1[17][orbt_zone] + (self.orbt1[17][orbt_zone + 1] - self.orbt1[17][orbt_zone]) * delta_orbt[m]) * self.cdr

            if self.orbt1[17][orbt_zone + 1] - self.orbt1[17][orbt_zone] > 0:
                self.sunalp[m] = (self.orbt1[17][orbt_zone] + (self.orbt1[17][orbt_zone + 1] - self.orbt1[17][orbt_zone] - 360.0) * delta_orbt[m]) * self.cdr

            self.sundel[m] = (self.orbt1[18][orbt_zone] + (self.orbt1[18][orbt_zone + 1] - self.orbt1[18][orbt_zone]) *delta_orbt[m]) * self.cdr
            '''
			npa[0][0] = self.orbt1[19][orbt_zone]
			npa[1][0] = self.orbt1[20][orbt_zone]
			npa[2][0] = self.orbt1[21][orbt_zone]
			npa[0][1] = self.orbt1[22][orbt_zone]
			npa[1][1] = self.orbt1[23][orbt_zone]
			npa[2][1] = self.orbt1[24][orbt_zone]
			npa[0][2] = self.orbt1[25][orbt_zone]
			npa[1][2] = self.orbt1[26][orbt_zone]
			npa[2][2] = self.orbt1[27][orbt_zone]
			'''
            # just replace directly, not reason to use npa. for this the atit was moved before the orbt computation, as opposed to the java code to avoid looping again through orbtila zones
            att2[0][m] = (self.orbt1[19][orbt_zone] * att1[0][m]) + (self.orbt1[22][orbt_zone] * att1[1][m]) + (self.orbt1[25][orbt_zone] * att1[2][m])
            att2[1][m] = (self.orbt1[20][orbt_zone] * att1[0][m]) + (self.orbt1[23][orbt_zone] * att1[1][m]) + (self.orbt1[26][orbt_zone] * att1[2][m])
            att2[2][m] = (self.orbt1[21][orbt_zone] * att1[0][m]) + (self.orbt1[24][orbt_zone] * att1[1][m]) + (self.orbt1[27][orbt_zone] * att1[2][m])




            # from pylab import show, imshow, title, colorbar
        # imshow(rtim.reshape((4800, 4200)))
        # colorbar()
        # title('sitagt')
        # show()


        '''
        for i in range(7):
            m_ = (rtim > self.orbt1[0][i]) & (rtim < self.orbt1[0][i + 1])
            # print i, rtim,  self.orbt1[0][i],  self.orbt1[0][i+1], m_
            # from pylab import show, imshow, title, colorbar
            # imshow(m_.reshape((3900, 4200)))
            # colorbar()
            # title('m_%s' % i)
            # show()
            mrtim = rtim[m_]
            npa = self.mg1110_np(i, rtim, self.orbt1, m_)
            print i, mrtim.size>0, npa.mean()
            if np.all(m_):
                if mrtim.size != rtim.size:  # what happens here?? the npa seems to be som
                    print 'Inconsistent results are possible! mrtim.size %s rtim.size %s, orbt1[0] size %s for i %s' % (mrtim.size, rtim.size, self.orbt1[0].size, i)

                npa = self.mg1110_np(i, rtim, self.orbt1, m_)
                print i,  npa
                # print m_npa.shape
                break
        '''




        # print att1
        # print npa[0][0], npa[0][1], npa[0][2]
        wksin = np.sin(self.sitagt)
        wkcos = np.cos(self.sitagt)

        att3[0] = (wkcos * att2[0]) + (wksin * att2[1])
        att3[1] = (-wksin * att2[0]) + (wkcos * att2[1])
        att3[2] = att2[2]

        self.sp = self.mg1200_np(att3)

        wkcos = np.cos(self.sundel)
        self.ss[0] = wkcos * np.cos(self.sunalp)
        self.ss[1] = wkcos * np.sin(self.sunalp)
        self.ss[2] = np.sin(self.sundel)

        return beta


    def mg1200_np(self, vec):

        vecu = [nan, nan, nan]

        rv1 = (vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2])
        if np.all(rv1 == 0):
            print ('issue')
            return vecu
        rv2 = np.sqrt(rv1)
        vecu[0] = vec[0] / rv2
        vecu[1] = vec[1] / rv2
        vecu[2] = vec[2] / rv2
        return vecu

    def mg1220_np(self, va, vb):

        vc = [nan, nan, nan]

        vc[0] = (va[1] * vb[2]) - (va[2] * vb[1])
        vc[1] = (va[2] * vb[0]) - (va[0] * vb[2])
        vc[2] = (va[0] * vb[1]) - (va[1] * vb[0])

        vd = self.mg1200_np(vc)

        return vd

    def mg1230_np(self, va, vb):

        as1 = (va[0] * vb[0]) + (va[1] * vb[1]) + (va[2] * vb[2])
        as2 = ((va[0] * va[0]) + (va[1] * va[1]) + (va[2] * va[2])) * (
            (vb[0] * vb[0]) + (vb[1] * vb[1]) + (vb[2] * vb[2]))

        if np.all(as2 == 0):
            return as2
        temp = as1 / np.sqrt(as2)
        tmask = (temp > 1) & (temp - 1. < 0.000001)
        temp = np.where(tmask, 1, temp)

        return np.arccos(temp)

    def mg1110(self, i, rtim, orbta):
        """

        :param i:
        :param rtim:
        :param orbta:
        :return: 3 by 3 array
        """

        npa  = [[nan for _i in range(3)] for _i in range(3)]
        aa = []
        for j in range(6):
            aa.append( (rtim - orbta[0][j]) / (orbta[0][j+1] - orbta[0][j]))
        if i != 7:
            #print 'i sat %d' % i
            #print 'rtim %.5f' % rtim
            #print 'orbta[0][%d] %.8f , orbta[0][%d+1] %.8f' % ( i, orbta[0][i], i, orbta[0][i+1])
            delt = (rtim - orbta[0][i]) / (orbta[0][i+1] - orbta[0][i])
            #print 'classic delt ', delt
            self.sat[0] = orbta[8][i] + (orbta[8][i+1] - orbta[8][i]) * delt
            self.sat[1] = orbta[9][i] + (orbta[9][i+1] - orbta[9][i]) * delt
            self.sat[2] = orbta[10][i] + (orbta[10][i+1] - orbta[10][i]) * delt
            #print 'setting sat %.5f, %.5f, %.5f:' %  (self.sat[0], self.sat[1], self.sat[2])
            self.sitagt = (orbta[14][i] + (orbta[14][i+1] - orbta[14][i]) * delt) * self.cdr
            if orbta[14][i+1] - orbta[14][i] < 0:
                self.sitagt = (orbta[14][i] + (orbta[14][i + 1] - orbta[14][i] + 360.0) * delt) * self.cdr

            self.sunalp = (orbta[17][i] + (orbta[17][i + 1] - orbta[17][i]) * delt) * self.cdr

            if orbta[17][i + 1] - orbta[17][i] > 0:
                self.sunalp = (orbta[17][i] + (orbta[17][i + 1] - orbta[17][i] - 360.0) * delt) * self.cdr

            self.sundel = (orbta[18][i] + (orbta[18][i+1] - orbta[18][i]) * delt) * self.cdr

            npa[0][0] = orbta[19][i]
            npa[1][0] = orbta[20][i]
            npa[2][0] = orbta[21][i]
            npa[0][1] = orbta[22][i]
            npa[1][1] = orbta[23][i]
            npa[2][1] = orbta[24][i]
            npa[0][2] = orbta[25][i]
            npa[1][2] = orbta[26][i]
            npa[2][2] = orbta[27][i]

        return npa

    def mg1200(self, vec):

        vecu = [nan, nan, nan]
        rv1 = (vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2])
        if rv1 == 0:
            return vecu
        rv2 = np.sqrt(rv1)
        vecu[0] = vec[0] / rv2
        vecu[1] = vec[1] / rv2
        vecu[2] = vec[2] / rv2
        return vecu

    def mg1100(self, rtim):
        """
        I have no clue , but it seem to do some conversion. Im mcidas check the file gms5_nav.for
        for further details.
        :param rtim:
        :return:
        """
        beta = delta = attalp = attdel = wkcos = wksin = 0.
        att1 = [nan, nan, nan]
        att2 = [nan, nan, nan]
        att3 = [nan, nan, nan]
        npa  = [[nan for _i in range(3)] for _ii in range(3)]

        #print rtim
        for i in range(7):
            #print i, rtim,  self.orbt1[0][i],  self.orbt1[0][i+1], rtim > self.orbt1[0][i] and rtim < self.orbt1[0][i+1]

            if rtim > self.orbt1[0][i] and rtim < self.orbt1[0][i+1]:
                #print 'computing npa', i
                npa = self.mg1110(i, rtim, self.orbt1)
                break
        #print 'npa classic', npa
        for i in range(9):
            #print 'atit', i, rtim, self.atit[0][i], self.atit[0][i + 1]
            if rtim >= self.atit[0][i] and rtim < self.atit[0][i+1]:
                delt = (rtim - self.atit[0][i]) / (self.atit[0][i + 1] - self.atit[0][i])
                attalp = self.atit[2][i] + (self.atit[2][i + 1] - self.atit[2][i]) * delt
                attdel = self.atit[3][i] + (self.atit[3][i + 1] - self.atit[3][i]) * delt
                beta = self.atit[4][i] + (self.atit[4][i + 1] - self.atit[4][i]) * delt
                if self.atit[4][i + 1] - self.atit[4][i] > 0:
                    beta = self.atit[4][i] + (self.atit[4][i + 1] - self.atit[4][i] - 360.0 * self.cdr) * delt
                break
        #print beta
        wkcos = np.cos(attdel)
        att1[0] = np.sin(attdel)
        att1[1] = wkcos * -np.sin(attalp)
        att1[2] =  wkcos * np.cos(attalp)
        #print wkcos, att1

        att2[0] = (npa[0][0] * att1[0]) + (npa[0][1] * att1[1]) + (npa[0][2] * att1[2])
        att2[1] = (npa[1][0] * att1[0]) + (npa[1][1] * att1[1]) + (npa[1][2] * att1[2])
        att2[2] = (npa[2][0] * att1[0]) + (npa[2][1] * att1[1]) + (npa[2][2] * att1[2])
        #print att1
        #print npa[0][0], npa[0][1], npa[0][2]
        wksin = np.sin(self.sitagt)
        wkcos = np.cos(self.sitagt)

        att3[0] = (wkcos * att2[0]) + (wksin * att2[1])
        att3[1] = (-wksin * att2[0]) + (wkcos * att2[1])
        att3[2] = att2[2]
        self.sp = self.mg1200(att3)

        wkcos = np.cos(self.sundel)
        self.ss[0] = wkcos * np.cos(self.sunalp)
        self.ss[1] = wkcos * np.sin(self.sunalp)
        self.ss[2] = np.sin(self.sundel)

        return beta


    def mg1220(self, va, vb):
        vc = [nan, nan, nan]

        vc[0] = (va[1] * vb[2]) - (va[2] * vb[1])
        vc[1] = (va[2] * vb[0]) - (va[0] * vb[2])
        vc[2] = (va[0] * vb[1]) - (va[1] * vb[0])

        vd = self.mg1200(vc)

        return vd


    def mg1210(self, va, vb):

        vc = [nan, nan, nan]
        vc[0] = (va[1] * vb[2]) - (va[2] * vb[1])
        vc[1] = (va[2] * vb[0]) - (va[0] * vb[2])
        vc[2] = (va[0] * vb[1]) - (va[1] * vb[0])
        return vc

    def mg1230(self, va, vb):

        as1 = (va[0] * vb[0]) + (va[1] * vb[1]) + (va[2] * vb[2])
        as2 = ((va[0] * va[0]) + (va[1] * va[1]) + (va[2] * va[2])) * ((vb[0] * vb[0]) + (vb[1] * vb[1]) + (vb[2] * vb[2]))

        if as2 ==0:
            return as2
        temp = as1/np.sqrt(as2)
        if temp > 1 and temp-1. < 0.000001:
            temp = 1
        return np.arccos(temp)



class ABINNav(AreaNav):

    def __init__(self, navblock=None, auxblock=None, **kwargs):
        self.loff = navblock[1] / 100000000.0;
        self.coff = navblock[2] / 100000000.0;
        self.lfac = navblock[3] / 100000000.0;
        self.cfac = navblock[4] / 100000000.0;
        self.plon = navblock[5] / 10.0;
        self.bres = navblock[6];




