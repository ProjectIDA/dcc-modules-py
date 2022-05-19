#######################################################################################################################
# Copyright (C) 2020  Regents of the University of California
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
import os.path
import sys
import math

from obspy.signal.invsim import paz_2_amplitude_value_of_freq_resp


def response_file_type(filename: str) -> str:
    """
    retrieve file type information from first line in header

    Args:
        filename: fully qualified path to response file

    Returns:
        str: 2-byte hex file type asc an ASCII string in format '0xABCD'
    """

    header = 'NO HEADER INFO FOUND'
    ftype = None
    with open(filename, 'rt') as rfl:
        header = rfl.readline()
        ftype = header.split('#')[0].strip()

    return ftype

class CoefficientsFile(object):
    """
    Generic (and simple) container that reads IDA filter coefficient files
    """

    # CONSTANTS
    FILTER_COMB = "0x4001"
    FILTER_FIR_SYM = "0x4002"
    FILTER_FIR_ASYM = "0x4004"


    def __init__(self, filename: str):
        """
        Coefficients initializer
        """

        self._filename = filename
        self.type = None
        self.group_delay = 0.0
        self._coeffs = []

        if not isinstance(self._filename, str):
            raise TypeError("Str type expected for filename: '{}'".format(self._filename))
        else:
            self._filename = os.path.abspath(os.path.expanduser(self._filename))

        if not os.path.exists(self._filename):
            raise Exception("File not found: '{}'".format(self._filename))

        self.load_file()

    def load_file(self):
        """
        Read response file
        """

        with open(self._filename, 'rt') as pzfl:
            lines = pzfl.readlines()
            self._parse_ida_coeff(lines)

    def filter_gain(self, nom_freq, in_srate):

        r, i = self.complex_filter_gain(nom_freq, in_srate) 
        return 1.0/math.sqrt(r**2 + i**2)

    def complex_filter_gain(self, nom_freq, in_srate):
  
        if not self.type in [self.FILTER_COMB, self.FILTER_FIR_SYM, self.FILTER_FIR_ASYM]:
            raise Exception('Transform type must either FILTER_COMB, FILTER_FIR_SYM or FILTER_FIR_ASYM')

        w = math.pi * 2.0 * float(nom_freq)
        wsint = w/float(in_srate)

        if self.type == self.FILTER_COMB:
            r, i = self._gain_complex_comb(wsint)

        elif self.type == self.FILTER_FIR_SYM:
            r, i = self._gain_complex_fir_sym(wsint)

        elif self.type == self.FILTER_FIR_ASYM:
            r, i = self._gain_complex_fir_asym(wsint)

        return r, i

    def _gain_complex_comb(self, wsint):
        """Ported from NRTS src/lib/filter/response.c func: RunMean(...)
        """

        r = math.sin((wsint/2.0)*self.num_coeffs) / math.sin(wsint / 2.0)
        i = 0.0
        return r, i

    def _gain_complex_fir_sym(self, wsint):
        """Ported from NRTS src/lib/filter/response.c func: FIRSymTrans(...)
        """

        R = 0.0
        n0 = math.floor( (self.num_coeffs - 1) / 2)
        for ndx in range(1, self.num_coeffs - int(n0)):
            R += self._coeffs[ndx + n0] * math.cos(wsint * ndx)

        r = self._coeffs[n0] + 2 * R
        i = 0.0

        return r, i

    def _gain_complex_fir_asym(self, wsint):
        """Ported from NRTS src/lib/filter/response.c func: FIRASymTrans(...)
        """

        R = 0.0
        I = 0.0
        ndx = 0
        break_ndx = 0
        for ndx in range(1, self.num_coeffs):
            break_ndx = ndx
            if self._coeffs[ndx] != self._coeffs[0]: 
                break

        if ((break_ndx + 1) == self.num_coeffs):

            if (wsint == 0):
                r = 1.0
            else:
                r = math.sin((wsint/2.0) * self.num_coeffs) / math.sin(wsint / 2.0) * self._coeffs[0]
            
            i = 0.0
        
        else:

            for ndx in range(self.num_coeffs):
                y = wsint * ndx
                R += self._coeffs[ndx] * -math.cos(y)
                I += self._coeffs[ndx] * -math.sin(y)

            mod = math.sqrt(R**2 + I**2)
            pha = math.atan2(I, R)
            r = mod * math.cos(pha)
            i = mod * math.sin(pha)

        return r, i

    def _parse_ida_coeff(self, lines: [str]):
        """

        Args:
            lines: Lines read from ida filter coefficient response file

        Returns:
            nothing
        """

        if ' # type ' in lines[0]:

            self.type = lines[0].split('#')[0].strip()

            # get num of coefficients
            parts = lines[1].split('#')
            coeff_cnt = int(parts[0])

            # get group delay of filter
            parts = lines[2].split('#')
            self.group_delay = float(parts[0])

            # skip over blank lines
            ndx = 3
            while lines[ndx].strip() == '':
                ndx += 1

            for _ in range(0, coeff_cnt):
                self.add_coeff(float(lines[ndx].strip()))
                ndx += 1

            # dbl check we got the num of values we expected
            if (coeff_cnt != self.num_coeffs):
                msg = "Error reading correct number of coefficients in file '{}'".format(self._filename)
                print(msg, file=sys.stderr)
                raise Exception(msg)

        else:
            raise Exception('Format error reading "ipaz" type paz file')

    @property
    def num_coeffs(self):
        return len(self._coeffs)

    @property
    def coeffs(self):
        return self._coeffs

    def add_coeff(self, coeff: float):
        self._coeffs.append(coeff)

    def __str__(self):
        txt = ''
        txt = f'{txt}File: {self._filename}\n'
        txt = f'{txt}Type: {self.type}\n'
        txt = f'{txt}Delay: {self.group_delay}\n'
        txt = f'{txt}Count: {self.num_coeffs}\n'
        txt = f'{txt}Coefficients:\n'
        for coeff in self._coeffs:
            txt = f'{txt}   {coeff:>16.6E}\n'

        return txt


class PolesZerosFile(object):

    """ Generic (and simple) container that reads IDA ipaz
    Poles and Zeros files. For something more fully functional, 
    see ida.signals.paz.PAZ."""

    # CONSTANTS
    TRANSFORM_ANALOG = "0x8001"
    TRANSFORM_IIR_DIGITAL = "0x8002"
    TRANSFORM_LAPLACE = "0x8004"

    def __init__(self, filename: str):
        """
        PolesZeros initializer
        """

        self._filename = filename
        self.type = None
        self._poles = []
        self._zeros = []

        if not isinstance(self._filename, str):
            raise TypeError("Str type expected for filename: '{}'".format(self._filename))
        else:
            self._filename = os.path.abspath(os.path.expanduser(self._filename))

        if not os.path.exists(self._filename):
            raise Exception("File not found: '{}'".format(self._filename))

        self.load_file()

    def load_file(self):
        """
        Read response file
        """

        with open(self._filename, 'rt') as pzfl:
            pz_lines = pzfl.readlines()
            self._parse_ida_paz(pz_lines)

    def _parse_ida_paz(self, pz_lines : [str]):
        """

        Args:
            pz_lines: Lines read from ida paz response file

        Returns:
            nothing
        """

        if ' # type ' in pz_lines[0]:

            self.type = pz_lines[0].split('#')[0].strip()

            # get num zeros and if any to skip when fitting
            parts = pz_lines[1].split('#')
            z_num = int(parts[0])

            # get num poles; perturbed indices (lf:hf); cnts to exclude when fitting (lf:hf)
            parts = pz_lines[2].split('#')
            p_num = int(parts[0])

            # skip over blank lines
            ndx = 3
            while 'zeros' not in pz_lines[ndx]:
                ndx += 1

            ndx += 1
            for _ in range(0, z_num):  # skip first zero, so one less to process
                vals = [val.strip() for val in pz_lines[ndx].split(',')]
                self.add_zero(complex(float(vals[0]), float(vals[1])))
                ndx += 1

            # skip over blank lines
            while 'poles' not in pz_lines[ndx]:
                ndx += 1

            ndx += 1
            for _ in range(0, p_num):
                vals = [val.strip() for val in pz_lines[ndx].split(',')]
                self.add_pole(complex(float(vals[0]), float(vals[1])))
                ndx += 1

            # dbl check we got the num of values we expected
            if (z_num != self.num_zeros) or (p_num != self.num_poles):
                msg = "Error reading correct number of poles and zeros in file '{}'".format(self._filename)
                print(msg, file=sys.stderr)
                raise Exception(msg)

        else:
            raise Exception('Format error reading "ipaz" type paz file')

    def a0(self, nomfreq, in_srate):

        r, i = self.a0_inv_complex(nomfreq, in_srate)
        return 1.0/math.sqrt(r**2 + i**2)

    def a0_inv_complex(self, nomfreq, in_srate):

        if not self.type in [self.TRANSFORM_ANALOG, self.TRANSFORM_LAPLACE, self.TRANSFORM_IIR_DIGITAL]:
            raise Exception('Transform type must either TRANSFORM_ANALOG, TRANSFORM_LAPLACE or TRANSFORM_IIR_DIGITAL')

        if self.type in [self.TRANSFORM_ANALOG, self.TRANSFORM_LAPLACE]:
            r, i = self._a0inv_complex_analog(nomfreq)
            # return self._a0_analog(nomfreq)
        elif self.type == self.TRANSFORM_IIR_DIGITAL:
            r, i = self._a0inv_complex_iir_digital(nomfreq, in_srate)

        return r, i

    def _a0inv_complex_iir_digital(self, freq, in_srate):
        """Ported from NRTS src/lib/filter/response.c `
            funcs: filterA0(...) and IIRtrans(...)
        """

        w = math.pi * 2 * freq  # convert to radians
        wsint = w / float(in_srate)        # sample interval in radians at freq 

        wcos = math.cos(wsint)
        wsin = math.sin(wsint)

        paz = {'poles': self.poles, 'zeros': self.zeros, 'gain': 1.0}

        mod = 1.0
        pha = 0.0

        for zero in paz['zeros']:
            R = wcos + zero.real
            I = wsin + zero.imag
            mod *= math.sqrt(R**2 + I**2)
            if R == 0.0 and I == 0.0:
                pha += 0.0  # keeping just to match c code in response.c
            else:
                pha += math.atan2(I, R)

        for pole in paz['poles']:
            R = wcos + pole.real
            I = wsin + pole.imag
            mod /= math.sqrt(R**2 + I**2)
            if R == 0.0 and I == 0.0:
                pha -= 0.0  # keeping just to match c code in response.c
            else:
                pha -= math.atan2(I, R)

        recipa0real = mod * math.cos(pha)
        recipa0imag = mod * math.sin(pha)

        return recipa0real, recipa0imag


    def _a0inv_complex_analog(self, freq):
        """Ported from NRTS src/lib/filter/response.c `
            funcs: filterA0(...) and AnalogTrans(...)
        """


        paz = {'poles': self.poles, 'zeros': self.zeros, 'gain': 1.0}

        mod = 1.0
        pha = 0.0

        for zero in paz['zeros']:
            R = -zero.real
            I = -zero.imag
            if R == 0.0:
                mod *= freq
                pha += math.pi/2.0
            else:
                mod *= math.sqrt((freq+I)**2 + R**2)
                pha += math.atan2(freq+I, R)

        for pole in paz['poles']:
            R = -pole.real
            I = -pole.imag
            if R == 0.0:
                mod /= freq
                pha -= math.pi/2.0
            else:
                mod /= math.sqrt((freq+I)**2 + R**2)
                pha == math.atan2(freq+I, R)

        recipa0real = mod * math.cos(pha)
        recipa0imag = mod * math.sin(pha)

        return recipa0real, recipa0imag

    @property
    def num_poles(self):
        return len(self._poles)

    @property
    def num_zeros(self):
        return len(self._zeros)

    @property
    def poles(self):
        return self._poles

    @property
    def zeros(self):
        return self._zeros

    def add_pole(self, pole: complex):
        self._poles.append(pole)

    def add_zero(self, zero: complex):
        self._zeros.append(zero)

    def __str__(self):
        txt = ''
        txt = f'{txt}File: {self._filename}\n'
        txt = f'{txt}Type: {self.type}\n'
        txt = f'{txt}Poles ({self.num_poles} cnt):\n'
        for pole in self._poles:
            txt = f'{txt}   {pole}\n'
        txt = f'{txt}Zeros ({self.num_zeros} cnt):\n'
        for zero in self._zeros:
            txt = f'{txt}   {zero}\n'

        return txt
