#######################################################################################################################
# Copyright (C) 2016  Regents of the University of California
#
# This is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License (GNU GPL) as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE.TXT in the root of the source code repository.
# Additionally, it can be found at http://www.gnu.org/licenses/.
#
# NOTES: Per GNU GPLv3 terms:
#   * This notice must be kept in this source file
#   * Changes to the source must be clearly noted with date & time of change
#
# If you use this software in a product, an explicit acknowledgment in the product documentation of the contribution
# by Project IDA, Institute of Geophysics and Planetary Physics, UCSD would be appreciated but is not required.
#######################################################################################################################

from collections import namedtuple
from numpy import sqrt

"""Instrument properties and related application constants and structures.

    Adding support for a new seismometer model:
        1) Add a SEISTYPE_ constant for model
        2) Add SEISTYPE_ constant into SEISMOMETER_MODELS list
        3) If Triaxial model:
            a) add to TRIAXIAL_SEIS_MODELS list
            b) add transformations of type XFRM_TYPE_XYZ2UVW and XFRM_TYPE_UVW2ENZ_ABS
                named SEISTYPE_ + 'XYZ2UVW' and SEISTYPE_ + 'UVW2ENZ_ABS'
            c) add transformations to TRIAXIAL_TRANSFORMS dict with SEISTYPE_ as key
        4) if model is CTBTO model, add to CTBTO_SEIS_MODELS list
        5) Add model to SEIS_INVERT_ lists, as appropriate
        6) Add entry to SEISMOMETER_RESPONSES dict
            a) Set embedded nominal full response file
            b) Set Fitting poles/zeros indices into full response paz
            c) Set (default) perturbing poles/zeros indices into full response paz
        7) Add model entry to INSRTUMENT_NOMINAL_GAINS dict
        8) Add model entry to Q330_GCALIB_FOR_SEIS (for Q330 <=> sensor impedance adjustment)
"""

# making these match existing instrument abbreviations in DataScope ABBREV table
# only include seismometer models deployed by IDA on 2016-05-26, plus TR360
SEISTYPE_STS1E3 = 'STS1E3'
SEISTYPE_STS1 = 'STS1'
SEISTYPE_STS2 = 'STS2'
SEISTYPE_STS25 = 'STS2_5'
SEISTYPE_STS25F = 'STS2_5_FAST'
SEISTYPE_STS5A = 'STS-5A'
SEISTYPE_STS6 = 'STS6'
SEISTYPE_TR12PA = 'TR12PB'
SEISTYPE_TR12PB = 'TR12PA'
SEISTYPE_TR240 = 'TR240'
SEISTYPE_TR360 = 'TR360'
SEISTYPE_GS13 = 'GS13'
SEISTYPE_3ESPC = '3ESPC'
SEISTYPE_KS54A = 'KS54A'
SEISTYPE_KS54B = 'KS54B'
SEISTYPE_KS54BEFI = 'KS54BEFI'
SEISTYPE_CMG3T = 'CMG3T'
SEISTYPE_FBAEST = 'FBAEST'
SEISTYPE_FBA23 = 'FBA23'
SEISTYPE_M2166 = 'M2166'

DIGITYPE_Q330 = 'Q330HR'

SEISMOMETER_MODELS = [
    SEISTYPE_STS1E3,
    # SEISTYPE_STS1HB,
    # SEISTYPE_STS1VB,
    SEISTYPE_STS1,
    SEISTYPE_STS2,
    SEISTYPE_STS25,
    # SEISTYPE_STS25F,
    SEISTYPE_STS5A,
    SEISTYPE_STS6,
    SEISTYPE_TR12PA,
    SEISTYPE_TR12PB,
    SEISTYPE_TR240,
    SEISTYPE_TR360,
    # SEISTYPE_GS13,
    # SEISTYPE_3ESPC,
    SEISTYPE_KS54A,
    SEISTYPE_KS54B,
    SEISTYPE_KS54BEFI,
    SEISTYPE_3ESPC,
    SEISTYPE_CMG3T,
    # SEISTYPE_FBAEST,
    # SEISTYPE_FBA23,
    SEISTYPE_M2166,
]

# only include currently supported triaxial list
# There are others. Need to define transforms
# as done with STS 2.5 below before adding addl seismometers
TRIAXIAL_SEIS_MODELS = [
    SEISTYPE_STS2,
    SEISTYPE_STS25,
    SEISTYPE_STS25F,
    SEISTYPE_STS5A,
    SEISTYPE_STS6,
    SEISTYPE_TR12PA,
    SEISTYPE_TR12PB,
    SEISTYPE_TR240,
    SEISTYPE_TR360
]

# need set of these transform for each supported triaxial sensor
STS2_XYZ2UVW = [
    [-sqrt(2/3),  0,         sqrt(1/3)],
    [ sqrt(1/6),  sqrt(1/2), sqrt(1/3)],
    [ sqrt(1/6),  -sqrt(1/2), sqrt(1/3)]
]
# STS2_XYZ2UVW = [
#     [-sqrt(6)/3,  sqrt(6)/6,   sqrt(6)/6],
#     [         0,  sqrt(2)/2,  -sqrt(2)/2],
#     [ sqrt(3)/3,  sqrt(3)/3,  sqrt(3)/3]
# ]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
STS2_UVW2ENZ_ABS = [
    [ sqrt(2/3), sqrt(1/6), sqrt(1/6)],
    [         0, sqrt(1/2), sqrt(1/2)],
    [ sqrt(1/3), sqrt(1/3), sqrt(1/3)]
]
STS2_5_XYZ2UVW = [
    [0,             -sqrt(6)/6,  sqrt(3)/6],
    [-sqrt(2)/4,  sqrt(6)/12, sqrt(3)/6],
    [ sqrt(2)/4,  sqrt(6)/12, sqrt(3)/6]
]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
STS2_5_UVW2ENZ_ABS = [
    [0,           sqrt(2),     sqrt(2)],
    [2*sqrt(6)/3, sqrt(6)/3,   sqrt(6)/3],
    [2*sqrt(3)/3, 2*sqrt(3)/3, 2*sqrt(3)/3]
]

STS5A_XYZ2UVW = [
    [0,             -sqrt(6)/6,  sqrt(3)/6],
    [-sqrt(2)/4,  sqrt(6)/12, sqrt(3)/6],
    [ sqrt(2)/4,  sqrt(6)/12, sqrt(3)/6]
]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
STS5A_UVW2ENZ_ABS = [
    [0,           sqrt(2),     sqrt(2)],
    [2*sqrt(6)/3, sqrt(6)/3,   sqrt(6)/3],
    [2*sqrt(3)/3, 2*sqrt(3)/3, 2*sqrt(3)/3]
]

STS6_XYZ2UVW = [
    [0,             -sqrt(6)/6,  sqrt(3)/6],
    [-sqrt(2)/4,  sqrt(6)/12, sqrt(3)/6],
    [ sqrt(2)/4,  sqrt(6)/12, sqrt(3)/6]
]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
STS6_UVW2ENZ_ABS = [
    [0,           sqrt(2),     sqrt(2)],
    [2*sqrt(6)/3, sqrt(6)/3,   sqrt(6)/3],
    [2*sqrt(3)/3, 2*sqrt(3)/3, 2*sqrt(3)/3]
]

# TRILLIUM TRIAXIAL TRANSFORMS
# for both A nad B models
TR120PH_XYZ2UVW = [
    [ sqrt(6) / 3,  0,           sqrt(3) / 3],
    [-sqrt(6) / 6,  sqrt(2) / 2, sqrt(3) / 3],
    [-sqrt(6) / 6, -sqrt(2) / 2, sqrt(3) / 3]
]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
TR120PH_UVW2ENZ_ABS = [
    [sqrt(6) / 3, sqrt(6) / 6, sqrt(6) / 6],
    [          0, sqrt(2) / 2, sqrt(2) / 2],
    [sqrt(3) / 3, sqrt(3) / 3, sqrt(3) / 3]
]

TR240_XYZ2UVW = [
    [ sqrt(6) / 3,  0,           sqrt(3) / 3],
    [-sqrt(6) / 6,  sqrt(2) / 2, sqrt(3) / 3],
    [-sqrt(6) / 6, -sqrt(2) / 2, sqrt(3) / 3]
]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
TR240_UVW2ENZ_ABS = [
    [sqrt(6) / 3, sqrt(6) / 6, sqrt(6) / 6],
    [          0, sqrt(2) / 2, sqrt(2) / 2],
    [sqrt(3) / 3, sqrt(3) / 3, sqrt(3) / 3]
]

TR360_XYZ2UVW = [
    [ sqrt(6) / 3,  0,           sqrt(3) / 3],
    [-sqrt(6) / 6,  sqrt(2) / 2, sqrt(3) / 3],
    [-sqrt(6) / 6, -sqrt(2) / 2, sqrt(3) / 3]
]
# from UVW back to ENZ, but ABS values, so all going in same direction at same time to maximize signal
TR360_UVW2ENZ_ABS = [
    [sqrt(6) / 3, sqrt(6) / 6, sqrt(6) / 6],
    [          0, sqrt(2) / 2, sqrt(2) / 2],
    [sqrt(3) / 3, sqrt(3) / 3, sqrt(3) / 3]
]


XFRM_TYPE_XYZ2UVW = 'XYZ2UVW'
XFRM_TYPE_UVW2ENZ_ABS = 'UVW2ENZ_ABS'

TRIAXIAL_TRANSFORMS = {
    SEISTYPE_STS2: {
        XFRM_TYPE_XYZ2UVW: STS2_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS : STS2_UVW2ENZ_ABS
    },
    SEISTYPE_STS25: {
        XFRM_TYPE_XYZ2UVW: STS2_5_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS : STS2_5_UVW2ENZ_ABS
    },
    SEISTYPE_STS25F: {
        XFRM_TYPE_XYZ2UVW: STS2_5_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS : STS2_5_UVW2ENZ_ABS
    },
    SEISTYPE_STS5A: {
        XFRM_TYPE_XYZ2UVW: STS5A_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS: STS5A_UVW2ENZ_ABS
    },
    SEISTYPE_STS6: {
        XFRM_TYPE_XYZ2UVW: STS6_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS: STS6_UVW2ENZ_ABS
    },
    SEISTYPE_TR12PA: {
        XFRM_TYPE_XYZ2UVW: TR120PH_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS: TR120PH_UVW2ENZ_ABS
    },
    SEISTYPE_TR12PB: {
        XFRM_TYPE_XYZ2UVW: TR120PH_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS: TR120PH_UVW2ENZ_ABS
    },
    SEISTYPE_TR240: {
        XFRM_TYPE_XYZ2UVW: TR240_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS: TR240_UVW2ENZ_ABS
    },
    SEISTYPE_TR360: {
        XFRM_TYPE_XYZ2UVW: TR360_XYZ2UVW,
        XFRM_TYPE_UVW2ENZ_ABS: TR360_UVW2ENZ_ABS
    }

}

CTBTO_SEIS_MODELS = [
    SEISTYPE_STS25,
    SEISTYPE_STS25F
]

SEIS_INVERT_CAL_CHAN = [SEISTYPE_GS13, SEISTYPE_TR12PA, SEISTYPE_TR12PB, SEISTYPE_TR240, SEISTYPE_TR360]
SEIS_INVERT_NORTH_CHAN = [SEISTYPE_STS1E3, SEISTYPE_STS1]
SEIS_INVERT_EAST_CHAN = [SEISTYPE_STS1E3, SEISTYPE_STS1, SEISTYPE_KS54BEFI]

CALTYPE_RBHF = 'rbhf'
CALTYPE_RBLF = 'rblf'
CALIBRATION_TYPES = [CALTYPE_RBHF, CALTYPE_RBLF]
CALTYPE_RBHF_SAMPLERATE = 100.0
CALTYPE_RBLF_SAMPLERATE = 1.0

# obtained from DB stage 1 rec
INSTRUMENT_NOMINAL_GAINS = {
    SEISTYPE_STS1E3: 2400,
    SEISTYPE_STS1: 2400,
    SEISTYPE_STS2: 1500,
    SEISTYPE_STS25: 1500,
    SEISTYPE_STS25F: 1500,
    SEISTYPE_STS5A: 1500,
    SEISTYPE_STS6: 1200,

    SEISTYPE_TR12PA: 1200,
    SEISTYPE_TR12PB: 1200,
    SEISTYPE_TR240: 1200,
    SEISTYPE_TR360: 1200,

    SEISTYPE_KS54A: 2400,
    SEISTYPE_KS54B: 2400,
    SEISTYPE_KS54BEFI: 2400,

    SEISTYPE_CMG3T: 1500,
    SEISTYPE_3ESPC: 20000,

    SEISTYPE_M2166: 2400,
}

# obtained from IDA filter file q330.40/20 May/Aug 2016
Q330_40_FIR_FILTER_DELAY = 17.218
Q330_FIR_FILTER_DELAY = {
    40: 17.218,
    20: 32.6090
}
Q330_40_FIR_FILTER_COEFFS = [ 4.18952E-13,  3.30318E-04,  1.02921E-03, -3.14123E-03,  2.05709E-04,  1.52521E-03, -6.23193E-03,  1.04801E-02,
 -1.31202E-02,  1.07821E-02, -1.44455E-03, -1.58729E-02,  3.95074E-02, -6.51036E-02,  8.53716E-02, -8.91913E-02,
 5.00619E-02,  8.37233E-01,  2.66723E-01, -1.66693E-01,  9.52840E-02, -5.09218E-02,  1.61458E-02,  7.06362E-03,
 -1.83877E-02,  1.99414E-02, -1.54895E-02,  8.52735E-03, -2.55789E-03, -1.81103E-03,  2.42649E-03, -3.75769E-03,
 4.67293E-04,  6.33072E-04, -1.56874E-06, -1.25480E-05,  3.21041E-07, -2.63324E-08, -5.09997E-08]

Q330_20_FIR_FILTER_COEFFS = [ -3.65342E-17,  3.67488E-08, -4.27060E-07,  1.14502E-06, -1.87594E-07, -3.37274E-07,  2.78747E-06, -3.74403E-06,
5.41172E-06,  7.47336E-06, -5.17759E-04,  2.10677E-04,  4.63258E-05, -6.08222E-04,  1.44175E-03, -2.40627E-03,
3.22534E-03, -3.50639E-03,  2.81441E-03, -7.71971E-04, -2.80512E-03,  7.77805E-03, -1.35815E-02,  1.91765E-02,
-2.29704E-02,  2.40398E-02, -2.20986E-02,  8.60734E-03,  1.17525E-02, -4.47787E-02,  9.64923E-02, -1.91755E-01,
5.27652E-01,  7.24167E-01, -1.56905E-01,  4.42574E-02,  3.14168E-03, -2.66714E-02,  3.61532E-02, -3.85687E-02,
3.10842E-02, -2.35259E-02,  1.53211E-02, -7.40398E-03,  1.09645E-03,  3.09797E-03, -5.19320E-03,  5.56131E-03,
-4.76110E-03,  3.38213E-03, -1.92052E-03,  7.15218E-04,  7.67719E-05, -4.51897E-04,  5.02700E-04, -5.65037E-04,
-5.56800E-05,  1.57736E-05, -1.41985E-06,  8.14909E-07,  6.80795E-07, -1.25273E-06,  1.52435E-06, -2.83336E-07,
-1.06384E-08,  1.25712E-09, -5.42954E-11]
Q330_FIR_COEFFS = {
    40: Q330_40_FIR_FILTER_COEFFS,
    20: Q330_20_FIR_FILTER_COEFFS
}

# obtained from IDA database. This value is specific to each Q330-Seis_Model combination
# and is in the IDA database as gcalib from stage 3 rec
Q330_GCALIB_FOR_SEIS = {
    SEISTYPE_KS54A: 0.9953,
    SEISTYPE_KS54B: 0.9953,
    SEISTYPE_KS54BEFI : 0.9953,
    SEISTYPE_STS1E3: 0.9481,
    SEISTYPE_STS1: 0.9481,
    SEISTYPE_STS2: 0.9911,  # value for Q330/STS2 combination
    SEISTYPE_STS25: 0.9911,  # value for Q330/STS2.5 combination
    SEISTYPE_STS25F: 0.9911,  # value for Q330/STS2.5 combination
    SEISTYPE_STS5A: 0.9911,  # value for Q330/STS5 combination
    SEISTYPE_STS6: 0.9911,  # value for Q330/STS6 combination
    SEISTYPE_TR12PA: 1.0011,
    SEISTYPE_TR12PB: 1.0011,
    SEISTYPE_TR240: 1.0011,
    SEISTYPE_TR360: 1.0011,
    SEISTYPE_CMG3T: 1.0011,
    SEISTYPE_3ESPC: 1.0011,
    SEISTYPE_M2166: 0.9481
}


# compute_response_fir(Q330_40_FIR_COEFFS,...)
# normalize_response() at 0Hz

# can be obtained by computing response of FIR filter in file q330.40
#   1) taking FFT of coefficients ==> Freq Resp
#   2) normalize on bin 0 (0 hz) value
#   3) Below is response amplitude value at 1Hz
#   OR
#   4) Run evalresp for stage 4 only
Q330_40HZ_NOMINAL_FIR_GAIN_1HZ = 1.00666915769
Q330_20HZ_NOMINAL_FIR_GAIN_1HZ = 1.006939

Q330_FIR_GAIN_1HZ = {
    40: Q330_40HZ_NOMINAL_FIR_GAIN_1HZ,
    20: Q330_20HZ_NOMINAL_FIR_GAIN_1HZ
}

# get from stages rec 3 gnom
Q330_NOMINAL_GAIN = {
    'A': 1.67e6,
    'B': 4.175e5
}
ComponentsTpl = namedtuple('Components', ['north', 'east', 'vertical'])
