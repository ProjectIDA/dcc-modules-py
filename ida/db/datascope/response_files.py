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

    def __init__(self, filename: str):
        """
        PolesZeros initializer
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
