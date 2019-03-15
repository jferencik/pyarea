__author__ = 'jano'
import math

from numpy import nan

from pyarea import utils


class Calibrator(object):

    CAL_NONE = -1
    CAL_MIN = 1
    CAL_RAW = 1
    CAL_RAD = 2
    CAL_ALB = 3
    CAL_TEMP = 4
    CAL_BRIT = 5
    CAL_MAX = 5
    # FY-2D
    SENSOR_FY2D = 36
    # FY-2E
    SENSOR_FY2E = 37
    # FY-2F
    SENSOR_FY2F = 38
    # FY-2G
    SENSOR_FY2G = 39
    # FY-2H
    SENSOR_FY2H = 40
    # Meteosat Second Generation imager.
    SENSOR_MSG_IMGR = 51
    # GOES 8 imager.
    SENSOR_GOES8_IMGR = 70
    # GOES 8 sounder.
    SENSOR_GOES8_SNDR = 71
    # GOES 9 imager.
    SENSOR_GOES9_IMGR = 72
    # GOES 9 sounder.
    SENSOR_GOES9_SNDR = 73
    # GOES 10 imager.
    SENSOR_GOES10_IMGR = 74
    # GOES 10 sounder.
    SENSOR_GOES10_SNDR = 75
    # GOES 12 imager.
    SENSOR_GOES12_IMGR = 78
    # GOES 12 sounder.
    SENSOR_GOES12_SNDR = 79
    # GOES 13 imager.
    SENSOR_GOES13_IMGR = 180
    # GOES 13 sounder.
    SENSOR_GOES13_SNDR = 181

    def setCalType(self, calType):
        raise Exception( 'Function selCalType is abstract')

    def calibrate_pixel(self, inputPixel, band, calTypeOut):
        raise Exception('Function calibrate_pixel is abstract')

    def calibrate_array(self, input_array, band, calTypeOut):
        raise Exception( 'Function calibrate_array is abstract')


class DefaultCalibrator(Calibrator):
    curCalType = 0

    def setCalType(self, calType):
        if calType < Calibrator.CAL_MIN or calType > Calibrator.CAL_MAX:
            return -1
        self.curCalType = calType

    def calibrate_pixel(self, inputPixel, band, calTypeOut ):
        """
            Calibrate one pixel
        """
        if self.curCalType in [self.CAL_RAW, self.CAL_RAD, self.CAL_ALB, self.CAL_TEMP, self.CAL_BRIT]:
            return inputPixel
        else:
            return 0.


    def calibrate_array(self, input_array, band, calTypeOut):

        out = list()
        for e in input_array:
            out.append(self.calibrate_pixel(e,band, calTypeOut))
        return  out


class GVARCalibrator(Calibrator):
    NUM_BANDS_IMAGER = 5
    NUM_BANDS_SOUNDER = 18

    NUM_VIS_DETECTORS = 8
    NUM_IR_DETECTORS = 2
    NUM_IR_BANDS = 4

    LOOKUP_TABLE_SZ_IMGR = 1024
    LOOKUP_TABLE_SZ_SNDR = 32768

    # var to store current cal type

    curCalType = 0
    index = 0

    visBiasCoef = [None] * NUM_VIS_DETECTORS
    visGain1Coef = [None] * NUM_VIS_DETECTORS
    visGain2Coef = [None] * NUM_VIS_DETECTORS
    visRadToAlb = 0.
    irBiasCoef = [[None] * NUM_IR_DETECTORS] * NUM_IR_BANDS
    irGainCoef = [[None] * NUM_IR_DETECTORS] * NUM_IR_BANDS
    sBiasCoef = [None] * NUM_BANDS_SOUNDER
    sGainCoef = [None] * NUM_BANDS_SOUNDER
    lookupTable = None

    #used in calibrator method
    gain = 0.
    bias = 0.
    scale = 1
    bandNum = 0
    sid = 0


    def __init__(self, sensorId=None, calBlock=None):
        """
            Default constructor
        """
        calIndex = 0
        self.sid = sensorId
        #correction for G12 onwards
        irOffset = 0 if self.sid > 77 else 2

        #intililize the lookpuptable
        self.lookupTable = [[nan] * self.NUM_BANDS_IMAGER] * self.LOOKUP_TABLE_SZ_IMGR
        if self.sid % 2 == 0:

            #read in an imager format cal block
            for i in range(self.NUM_VIS_DETECTORS):
                self.visBiasCoef[i] = utils.gould2native(calBlock[i])
                calIndex += 1
            for i in range(self.NUM_VIS_DETECTORS):
                self.visGain1Coef[i] = utils.gould2native(calBlock[i])
                calIndex += 1
            for i in range(self.NUM_VIS_DETECTORS):
                self.visGain2Coef[i] = utils.gould2native(calBlock[i])
                calIndex+=1

            self.visRadToAlb = utils.gould2native(calBlock[calIndex])
            calIndex+=1

            for i in range(self.NUM_IR_BANDS):
                self.irBiasCoef[0][(i+irOffset) % self.NUM_IR_BANDS] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            for i in range(self.NUM_IR_BANDS):
                self.irBiasCoef[1][(i + irOffset) % self.NUM_IR_BANDS] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            for i in range(self.NUM_IR_BANDS):
                self.irGainCoef[0][(i + irOffset) % self.NUM_IR_BANDS] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            for i in range(self.NUM_IR_BANDS):
                self.irGainCoef[1][(i + irOffset) % self.NUM_IR_BANDS] = utils.gould2native(calBlock[calIndex])
                calIndex += 1
        else:
            # read in a sounder format cal block
            for i in range(self.NUM_VIS_DETECTORS/2):
                self.visBiasCoef[i] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            for i in range(self.NUM_VIS_DETECTORS / 2):
                self.visGain1Coef[i] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            for i in range(self.NUM_VIS_DETECTORS / 2):
                self.visGain2Coef[i] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            self.visRadToAlb = utils.gould2native(calBlock[calIndex])
            calIndex += 1

            for i in range(self.NUM_BANDS_SOUNDER):
                self.sBiasCoef[i] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

            for i in range(self.NUM_BANDS_SOUNDER):
                self.sGainCoefCoef[i] = utils.gould2native(calBlock[calIndex])
                calIndex += 1

    def setCalType(self, calType):
        if calType < Calibrator.CAL_MIN or calType > Calibrator.CAL_MAX:
            return -1
        self.curCalType = calType
        return 0

    def radToTemp(self, inval=None, band=None, sId=None):
        raise Exception('Method radToTemp is abstract')

    def calibrate_array(self, input_array, band, calTypeOut):
        out = list()
        for e in input_array:
            out.append(self.calibrate_pixel(e, band, calTypeOut))
        return out

    def calibrate_pixel(self, inputPixel, band, calTypeOut):
        out = 0.
        if band != self.bandNum:
            self.bandNum = band
            if self.sid % 2 == 0:
                if band == 1:
                    self.gain = self.visGain1Coef[0]
                    self.bias = self.visBiasCoef[0]
                else:
                    self.gain = self.irGainCoef[0][band-2]
                    self.bias = self.irBiasCoef[0][band-2]
                self.scale = 32
            else:
                if band == 19:
                    self.gain = self.visGain1Coef[0]
                    self.bias = self.visBiasCoef[0]
                else:
                    self.gain = self.sGainCoef[band-1]
                    self.bias = self.sBiasCoef[band-1]

                self.scale = 2
        if self.curCalType == self.CAL_BRIT:
            self.index = int(inputPixel + 128)
        else:
            self.index = int(inputPixel / self.scale)

        if self.lookupTable[band -1][self.index] is not nan:
            return  self.lookupTable[band -1][self.index]

        if self.curCalType == self.CAL_RAW:
            out = inputPixel
            if calTypeOut == self.CAL_RAW:
                pass
            elif calTypeOut == self.CAL_RAD:
                if self.sid % 2 == 0:
                    out = inputPixel / self.scale
                out = (out - self.bias) / self.gain
            elif calTypeOut == self.CAL_BRIT:
                out = max((660 - int(out * 2)), 0) if out >= 242. else min((418 - int(out), 255))

        elif self.curCalType in (self.CAL_RAD, self.CAL_ALB, self.CAL_TEMP, self.CAL_BRIT):
            out = inputPixel
        self.lookupTable[band-1][self.index] = out
        return out


class G8GVARCalibrator(GVARCalibrator):

    imager8FK1 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder8FK1 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager8FK2 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder8FK2 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager8TC1 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder8TC1 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager8TC2 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder8TC2 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER

    imager8FK1[0] = 0.0000000E0
    imager8FK1[1] = 0.1999862E6
    imager8FK1[2] = 0.3879239E5
    imager8FK1[3] = 0.9737930E4
    imager8FK1[4] = 0.6944640E4

    sounder8FK1[0] = 0.3756810E4
    sounder8FK1[1] = 0.4011100E4
    sounder8FK1[2] = 0.4296870E4
    sounder8FK1[3] = 0.4681130E4
    sounder8FK1[4] = 0.4975250E4
    sounder8FK1[5] = 0.5881410E4
    sounder8FK1[6] = 0.6787440E4
    sounder8FK1[7] = 0.8873710E4
    sounder8FK1[8] = 0.1299794E5
    sounder8FK1[9] = 0.2862932E5
    sounder8FK1[10] = 0.3424830E5
    sounder8FK1[11] = 0.4311430E5
    sounder8FK1[12] = 0.1242353E6
    sounder8FK1[13] = 0.1281235E6
    sounder8FK1[14] = 0.1351482E6
    sounder8FK1[15] = 0.1691671E6
    sounder8FK1[16] = 0.1882350E6
    sounder8FK1[17] = 0.2257944E6

    imager8FK2[0] = 0.0000000E0
    imager8FK2[1] = 0.3684270E4
    imager8FK2[2] = 0.2132720E4
    imager8FK2[3] = 0.1345370E4
    imager8FK2[4] = 0.1201990E4

    sounder8FK2[0] = 0.3765120E4
    sounder8FK2[1] = 0.3981160E4
    sounder8FK2[2] = 0.4281880E4
    sounder8FK2[3] = 0.4678910E4
    sounder8FK2[4] = 0.4962590E4
    sounder8FK2[5] = 0.5860420E4
    sounder8FK2[6] = 0.6770320E4
    sounder8FK2[7] = 0.8958910E4
    sounder8FK2[8] = 0.1296593E5
    sounder8FK2[9] = 0.2839828E5
    sounder8FK2[10] = 0.3420134E5
    sounder8FK2[11] = 0.4252514E5
    sounder8FK2[12] = 0.1240574E6
    sounder8FK2[13] = 0.1280114E6
    sounder8FK2[14] = 0.1348497E6
    sounder8FK2[15] = 0.1678142E6
    sounder8FK2[16] = 0.1888012E6
    sounder8FK2[17] = 0.2258565E6

    imager8TC1[0] = 0.0000000E0
    imager8TC1[1] = 0.6357000E0
    imager8TC1[2] = 0.6060000E0
    imager8TC1[3] = 0.3735000E0
    imager8TC1[4] = 0.2217000E0

    sounder8TC1[0] = 0.1230000E-1
    sounder8TC1[1] = 0.1330000E-1
    sounder8TC1[2] = 0.1860000E-1
    sounder8TC1[3] = 0.1500000E-1
    sounder8TC1[4] = 0.1650000E-1
    sounder8TC1[5] = 0.4740000E-1
    sounder8TC1[6] = 0.1318000E0
    sounder8TC1[7] = 0.1200000E0
    sounder8TC1[8] = 0.4260000E-1
    sounder8TC1[9] = 0.1505000E0
    sounder8TC1[10] = 0.2743000E0
    sounder8TC1[11] = 0.1447000E0
    sounder8TC1[12] = 0.2240000E-1
    sounder8TC1[13] = 0.2200000E-1
    sounder8TC1[14] = 0.2170000E-1
    sounder8TC1[15] = 0.5790000E-1
    sounder8TC1[16] = 0.6230000E-1
    sounder8TC1[17] = 0.3675000E0

    imager8TC2[0] = 0.0000000E0
    imager8TC2[1] = 0.9991000E0
    imager8TC2[2] = 0.9986000E0
    imager8TC2[3] = 0.9987000E0
    imager8TC2[4] = 0.9992000E0

    sounder8TC2[0] = 0.9999000E0
    sounder8TC2[1] = 0.9999000E0
    sounder8TC2[2] = 0.9999000E0
    sounder8TC2[3] = 0.9999000E0
    sounder8TC2[4] = 0.9999000E0
    sounder8TC2[5] = 0.9998000E0
    sounder8TC2[6] = 0.9995000E0
    sounder8TC2[7] = 0.9996000E0
    sounder8TC2[8] = 0.9999000E0
    sounder8TC2[9] = 0.9996000E0
    sounder8TC2[10] = 0.9993000E0
    sounder8TC2[11] = 0.9997000E0
    sounder8TC2[12] = 0.1000000E1
    sounder8TC2[13] = 0.1000000E1
    sounder8TC2[14] = 0.1000000E1
    sounder8TC2[15] = 0.9999000E0
    sounder8TC2[16] = 0.9999000E0
    sounder8TC2[17] = 0.9995000E0

    def radToTemp(self, inval=None, band=None, sId=None):
        outval = None
        if sId % 2 == 0:
            expn = (self.imager8FK1[band-1] / inval) + 1.
            temp = self.imager8FK2[band-1] / math.log(expn)
            outval = float(temp - (self.imager8TC1[band-1] / self.imager8TC2[band-1]))
        else:
            expn = (self.sounder8FK1[band - 1] / inval) + 1.
            temp = self.sounder8FK2[band - 1] / math.log(expn)
            outval = float(temp - (self.sounder8TC1[band - 1] / self.sounder8TC2[band - 1]))

        return outval

class G9GVARCalibrator(GVARCalibrator):

    imager9FK1 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder9FK1 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager9FK2 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder9FK2 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager9TC1 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder9TC1 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager9TC2 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder9TC2 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER

    imager9FK1[0] = 0.0000000E0
    imager9FK1[1] = 0.1988078E6
    imager9FK1[2] = 0.3873241E5
    imager9FK1[3] = 0.9717210E4
    imager9FK1[4] = 0.6899470E4

    sounder9FK1[0] = 0.3765120E4
    sounder9FK1[1] = 0.3981160E4
    sounder9FK1[2] = 0.4281880E4
    sounder9FK1[3] = 0.4678910E4
    sounder9FK1[4] = 0.4962590E4
    sounder9FK1[5] = 0.5860420E4
    sounder9FK1[6] = 0.6770320E4
    sounder9FK1[7] = 0.8958910E4
    sounder9FK1[8] = 0.1296593E5
    sounder9FK1[9] = 0.2839828E5
    sounder9FK1[10] = 0.3420134E5
    sounder9FK1[11] = 0.4252514E5
    sounder9FK1[12] = 0.1240574E6
    sounder9FK1[13] = 0.1280114E6
    sounder9FK1[14] = 0.1348497E6
    sounder9FK1[15] = 0.1678142E6
    sounder9FK1[16] = 0.1888012E6
    sounder9FK1[17] = 0.2258565E6

    imager9FK2[0] = 0.0000000E0
    imager9FK2[1] = 0.3677020E4
    imager9FK2[2] = 0.2131620E4
    imager9FK2[3] = 0.1344410E4
    imager9FK2[4] = 0.1199380E4

    sounder9FK2[0] = 0.9801200E3
    sounder9FK2[1] = 0.9985200E3
    sounder9FK2[2] = 0.1023050E4
    sounder9FK2[3] = 0.1053740E4
    sounder9FK2[4] = 0.1074620E4
    sounder9FK2[5] = 0.1135870E4
    sounder9FK2[6] = 0.1191850E4
    sounder9FK2[7] = 0.1308490E4
    sounder9FK2[8] = 0.1480080E4
    sounder9FK2[9] = 0.1922130E4
    sounder9FK2[10] = 0.2045030E4
    sounder9FK2[11] = 0.2199040E4
    sounder9FK2[12] = 0.3142140E4
    sounder9FK2[13] = 0.3175180E4
    sounder9FK2[14] = 0.3230740E4
    sounder9FK2[15] = 0.3475050E4
    sounder9FK2[16] = 0.3614260E4
    sounder9FK2[17] = 0.3836740E4

    imager9TC1[0] = 0.0000000E0
    imager9TC1[1] = 0.5864000E0
    imager9TC1[2] = 0.4841000E0
    imager9TC1[3] = 0.3622000E0
    imager9TC1[4] = 0.2014000E0

    sounder9TC1[0] = 0.9900000E-2
    sounder9TC1[1] = 0.1190000E-1
    sounder9TC1[2] = 0.1220000E-1
    sounder9TC1[3] = 0.1190000E-1
    sounder9TC1[4] = 0.1350000E-1
    sounder9TC1[5] = 0.4400000E-1
    sounder9TC1[6] = 0.1345000E0
    sounder9TC1[7] = 0.1193000E0
    sounder9TC1[8] = 0.4070000E-1
    sounder9TC1[9] = 0.1438000E0
    sounder9TC1[10] = 0.2762000E0
    sounder9TC1[11] = 0.1370000E0
    sounder9TC1[12] = 0.1890000E-1
    sounder9TC1[13] = 0.1980000E-1
    sounder9TC1[14] = 0.1910000E-1
    sounder9TC1[15] = 0.5310000E-1
    sounder9TC1[16] = 0.6120000E-1
    sounder9TC1[17] = 0.3042000E0

    imager9TC2[0] = 0.0000000E0
    imager9TC2[1] = 0.9992000E0
    imager9TC2[2] = 0.9989000E0
    imager9TC2[3] = 0.9988000E0
    imager9TC2[4] = 0.9992000E0

    sounder9TC2[0] = 0.1000000E1
    sounder9TC2[1] = 0.9999000E0
    sounder9TC2[2] = 0.9999000E0
    sounder9TC2[3] = 0.9999000E0
    sounder9TC2[4] = 0.9999000E0
    sounder9TC2[5] = 0.9998000E0
    sounder9TC2[6] = 0.9995000E0
    sounder9TC2[7] = 0.9996000E0
    sounder9TC2[8] = 0.9999000E0
    sounder9TC2[9] = 0.9996000E0
    sounder9TC2[10] = 0.9993000E0
    sounder9TC2[11] = 0.9997000E0
    sounder9TC2[12] = 0.1000000E1
    sounder9TC2[13] = 0.1000000E1
    sounder9TC2[14] = 0.1000000E1
    sounder9TC2[15] = 0.9999000E0
    sounder9TC2[16] = 0.9999000E0
    sounder9TC2[17] = 0.9996000E0

    def radToTemp(self, inval=None, band=None, sId=None):
        outval = None
        if sId % 2 == 0:
            expn = (self.imager9FK1[band - 1] / inval) + 1.
            temp = self.imager9FK2[band - 1] / math.log(expn)
            outval = float(temp - (self.imager9TC1[band - 1] / self.imager9TC2[band - 1]))
        else:
            expn = (self.sounder9FK1[band - 1] / inval) + 1.
            temp = self.sounder9FK2[band - 1] / math.log(expn)
            outval = float(temp - (self.sounder9TC1[band - 1] / self.sounder9TC2[band - 1]))

        return outval


class G10GVARCalibrator(GVARCalibrator):

    imager10FK1 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder10FK1 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager10FK2 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder10FK2 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager10TC1 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder10TC1 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER
    imager10TC2 = [None] * GVARCalibrator.NUM_BANDS_IMAGER
    sounder10TC2 = [None] * GVARCalibrator.NUM_BANDS_SOUNDER

    imager10FK1[0] = 0.00000E0
    imager10FK1[1] = 0.19841E6
    imager10FK1[2] = 0.39086E5
    imager10FK1[3] = 0.97744E4
    imager10FK1[4] = 0.68286E4

    sounder10FK1[0] = 0.37305E4
    sounder10FK1[1] = 0.40039E4
    sounder10FK1[2] = 0.43124E4
    sounder10FK1[3] = 0.46616E4
    sounder10FK1[4] = 0.49734E4
    sounder10FK1[5] = 0.58698E4
    sounder10FK1[6] = 0.68161E4
    sounder10FK1[7] = 0.89404E4
    sounder10FK1[8] = 0.12973E5
    sounder10FK1[9] = 0.28708E5
    sounder10FK1[10] = 0.34401E5
    sounder10FK1[11] = 0.43086E5
    sounder10FK1[12] = 0.12468E6
    sounder10FK1[13] = 0.12882E6
    sounder10FK1[14] = 0.13532E6
    sounder10FK1[15] = 0.16853E6
    sounder10FK1[16] = 0.18862E6
    sounder10FK1[17] = 0.22487E6

    imager10FK2[0] = 0.00000E0
    imager10FK2[1] = 0.36745E4
    imager10FK2[2] = 0.21381E4
    imager10FK2[3] = 0.13470E4
    imager10FK2[4] = 0.11953E4

    sounder10FK2[0] = 0.97710E3
    sounder10FK2[1] = 0.10004E4
    sounder10FK2[2] = 0.10255E4
    sounder10FK2[3] = 0.10524E4
    sounder10FK2[4] = 0.10754E4
    sounder10FK2[5] = 0.11365E4
    sounder10FK2[6] = 0.11945E4
    sounder10FK2[7] = 0.13076E4
    sounder10FK2[8] = 0.14804E4
    sounder10FK2[9] = 0.19291E4
    sounder10FK2[10] = 0.20490E4
    sounder10FK2[11] = 0.22087E4
    sounder10FK2[12] = 0.31474E4
    sounder10FK2[13] = 0.31818E4
    sounder10FK2[14] = 0.32345E4
    sounder10FK2[15] = 0.34800E4
    sounder10FK2[16] = 0.36131E4
    sounder10FK2[17] = 0.38311E4

    imager10TC1[0] = 0.00000
    imager10TC1[1] = 0.62226
    imager10TC1[2] = 0.61438
    imager10TC1[3] = 0.27791
    imager10TC1[4] = 0.21145

    sounder10TC1[0] = 0.00988
    sounder10TC1[1] = 0.01196
    sounder10TC1[2] = 0.01245
    sounder10TC1[3] = 0.01245
    sounder10TC1[4] = 0.01366
    sounder10TC1[5] = 0.04311
    sounder10TC1[6] = 0.13973
    sounder10TC1[7] = 0.11707
    sounder10TC1[8] = 0.03979
    sounder10TC1[9] = 0.14968
    sounder10TC1[10] = 0.27603
    sounder10TC1[11] = 0.13049
    sounder10TC1[12] = 0.02008
    sounder10TC1[13] = 0.01834
    sounder10TC1[14] = 0.02017
    sounder10TC1[15] = 0.05292
    sounder10TC1[16] = 0.05330
    sounder10TC1[17] = 0.28683

    imager10TC2[0] = 0.00000
    imager10TC2[1] = 0.99912
    imager10TC2[2] = 0.99857
    imager10TC2[3] = 0.99905
    imager10TC2[4] = 0.99919

    sounder10TC2[0] = 0.99995
    sounder10TC2[1] = 0.99994
    sounder10TC2[2] = 0.99994
    sounder10TC2[3] = 0.99995
    sounder10TC2[4] = 0.99994
    sounder10TC2[5] = 0.99983
    sounder10TC2[6] = 0.99947
    sounder10TC2[7] = 0.99959
    sounder10TC2[8] = 0.99988
    sounder10TC2[9] = 0.99962
    sounder10TC2[10] = 0.99933
    sounder10TC2[11] = 0.99970
    sounder10TC2[12] = 0.99997
    sounder10TC2[13] = 0.99997
    sounder10TC2[14] = 0.99997
    sounder10TC2[15] = 0.99992
    sounder10TC2[16] = 0.99992
    sounder10TC2[17] = 0.99961


    def radToTemp(self, inval=None, band=None, sId=None):

        outval = None
        if sId % 2 == 0:
            expn = (self.imager10FK1[band - 1] / inval) + 1.
            temp = self.imager10FK2[band - 1] / math.log(expn)
            outval = float(temp - (self.imager10TC1[band - 1] / self.imager10TC2[band - 1]))
        else:
            expn = (self.sounder10FK1[band - 1] / inval) + 1.
            temp = self.sounder10FK2[band - 1] / math.log(expn)
            outval = float(temp - (self.sounder10TC1[band - 1] / self.sounder10TC2[band - 1]))

        return outval