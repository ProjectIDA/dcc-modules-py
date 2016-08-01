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
import os.path
import logging
from copy import copy, deepcopy

from numpy import zeros as npzeros
from numpy import array, pi, complex128, concatenate

import ida.signals.utils


class PAZ(object):
    """ Generic object represnting frequency response in
    Poles and Zeros format"""

    PARTIAL_ALL = 1
    PARTIAL_FITTING_LF = 2
    PARTIAL_FITTING_HF = 3
    PARTIAL_PERTURBING_LF = 4
    PARTIAL_PERTURBING_HF = 4

    FILE_FORMATS = ['ida']
    MODE_ZEROS = dict(acc=0, vel=1, disp=2)
    UNITS = ['hz', 'rad']
    PAZ_HEADER_IDA = '0x8001  # type 0x8001 = analog'

    def __init__(self, mode, units, pzfilename=None, fileformat=None):
        """

        :param mode:
        :type mode:
        :param units:
        :type units:
        :param pzfilename:
        :type pzfilename:
        :param fileformat:
        :type fileformat:
        """
        self._filename = pzfilename
        self.fileformat = fileformat
        self._mode = mode
        self._units = units
        self._h0 = 1
        self._poles = npzeros(0, dtype=complex128)
        self._zeros = npzeros(0, dtype=complex128)
        self._poles_no_fitting_count = (0, 0)   # (lf, hf)
        self._zeros_no_fitting_count = (0, 0)
        self._poles_default_pert_ndxs = ([], [])   # (lf, hf)
        self._zeros_default_pert_ndxs = ([], [])

        if mode not in PAZ.MODE_ZEROS.keys():
            raise ValueError("Invalid MODE requested: '{}'. Valid values: {}".format(mode, PAZ.MODE_ZEROS.keys()))

        if units not in PAZ.UNITS:
            raise ValueError("Invalid units specified: '{}'. Valid values: {}".format(units, PAZ.UNITS))

        if pzfilename and fileformat:

            if not isinstance(self._filename, str):
                raise TypeError("Str type expected for filename: '{}'".format(self._filename))
            else:
                self._filename = os.path.abspath(os.path.expanduser(self._filename))

            if not os.path.exists(self._filename):
                raise Exception("File not found: '{}'".format(self._filename))

            if fileformat not in PAZ.FILE_FORMATS:
                raise ValueError("Unsupported fileformat requested: '{}'.  Supported values: {}".format(
                    fileformat,
                    PAZ.FILE_FORMATS))

            self._load_paz_file()

    def _load_paz_file(self):
        """
        Load me a file...

        :return: none
        :rtype: none
        """

        with open(self._filename, 'rt') as pzfl:
            pzlines = pzfl.readlines()
            if self.fileformat == 'ida':
                self._parse_ida_paz(pzlines)

    def _parse_ida_paz(self, pzlines):
        """

        :param pzlines:
        :type pzlines:
        :return:
        :rtype:
        """

        #todo:  should read non-fitting counts before perturb indices to
        #       check valid indices within fitting set of values
        if ' # type ' in pzlines[0]:

            # get num zeros and if any to skip when fitting
            parts = pzlines[1].split('#')
            z_num = int(parts[0])

            if len(parts) > 3:  # has extra non-fitting value counts
                nofit_parts = parts[3].split(':')  # split for lf:hf values
                if nofit_parts[0].strip() != '':
                    val = int(nofit_parts[0])
                    if (val < 0) or (val > z_num):
                        raise ValueError('LF non-fitting count out of range on zeros line":'+ pzlines[1])
                    else:
                        self._zeros_no_fitting_count = (val, self._zeros_no_fitting_count[1])
                else:
                    self._zeros_no_fitting_count = (0, self._zeros_no_fitting_count[1])
                if nofit_parts[1].strip() != '':
                    val = int(nofit_parts[1])
                    if (val < 0) or (val > z_num):
                        raise ValueError('HF non-fitting count out of range on zeros line":'+ pzlines[1])
                    else:
                        self._zeros_no_fitting_count = (self._zeros_no_fitting_count[0], val)
                else:
                    self._zeros_no_fitting_count = (self._zeros_no_fitting_count[0], 0)
            else:
                self._zeros_no_fitting_count = (0, 0)

            if len(parts) > 2:  # has perturbed indice info
                pert_parts = parts[2].split(':')  # split for lf:hf values
                if pert_parts[0].strip() != '':  # process lf piece
                    ndxs = pert_parts[0].strip().split(',')
                    lst = [int(ndx)-1 for ndx in ndxs]  # make sure valid index
                    if (min(lst) < 0) or (max(lst) > (z_num-self._zeros_no_fitting_count[0]-1)):
                        raise ValueError('LF Pert index out of range on "# of zeros line":'+ pzlines[1])
                    else:
                        self._zeros_default_pert_ndxs = (lst, self._zeros_default_pert_ndxs[1])
                else:
                    self._zeros_default_pert_ndxs = ([], self._zeros_default_pert_ndxs[1])
                if pert_parts[1].strip() != '':  # process hf piece
                    ndxs = pert_parts[1].strip().split(',')
                    lst = [int(ndx)-1 for ndx in ndxs]  # make sure valid index
                    if (min(lst) < 0) or (max(lst) > (z_num-self._zeros_no_fitting_count[1]-1)):
                        raise ValueError('HF Pert index out of range on "# of zeros line":'+ pzlines[1])
                    else:
                        self._zeros_default_pert_ndxs = (self._zeros_default_pert_ndxs[0], lst)
                else:
                    self._zeros_default_pert_ndxs = (self._zeros_default_pert_ndxs[0], [])
            else:
                self._zeros_default_pert_ndxs = ([], [])
                self._zeros_no_fitting_count = (0, 0)

            # get num poles; perturbed indices (lf:hf); cnts to exclude when fitting (lf:hf)
            parts = pzlines[2].split('#')
            p_num = int(parts[0])

            if len(parts) > 3:  # has extra non-fitting value counts
                nofit_parts = parts[3].split(':')  # split for lf:hf values
                if nofit_parts[0].strip() != '':
                    val = int(nofit_parts[0])
                    if (val < 0) or (val > p_num):
                        raise ValueError('LF non-fitting count out of range on poles line":'+ pzlines[2])
                    else:
                        self._poles_no_fitting_count = (val, self._poles_no_fitting_count[1])
                else:
                    self._poles_no_fitting_count = (0, self._poles_no_fitting_count[1])
                if nofit_parts[1].strip() != '':
                    val = int(nofit_parts[1])
                    if (val < 0) or (val > p_num):
                        raise ValueError('HF non-fitting count out of range on poles line":'+ pzlines[2])
                    else:
                        self._poles_no_fitting_count = (self._poles_no_fitting_count[0], val)
                else:
                    self._poles_no_fitting_count = (self._poles_no_fitting_count[0], 0)
            else:
                self._poles_no_fitting_count = (0, 0)

            if len(parts) > 2:  # has perturbed indice info
                pert_parts = parts[2].split(':')  # split for lf:hf values
                if pert_parts[0].strip() != '':  # process lf piece
                    ndxs = pert_parts[0].strip().split(',')
                    lst = [int(ndx)-1 for ndx in ndxs]
                    if (min(lst) < 0) or (max(lst) > (p_num-self._poles_no_fitting_count[0]-1)):
                        raise ValueError('LF Pert index out of range on poles line":'+ pzlines[2])
                    else:
                        self._poles_default_pert_ndxs = (lst, self._poles_default_pert_ndxs[1])
                else:
                    self._poles_default_pert_ndxs = ([], self._poles_default_pert_ndxs[1])
                if pert_parts[1].strip() != '':  # process hf piece
                    ndxs = pert_parts[1].strip().split(',')
                    lst = [int(ndx)-1 for ndx in ndxs]
                    if (min(lst) < 0) or (max(lst) > (p_num-self._poles_no_fitting_count[1]-1)):
                        raise ValueError('HF Pert index out of range on poles line":'+ pzlines[2])
                    else:
                        self._poles_default_pert_ndxs = (self._poles_default_pert_ndxs[0], lst)
                else:
                    self._poles_default_pert_ndxs = (self._poles_default_pert_ndxs[0], [])
            else:
                self._zeros_default_pert_ndxs = ([], [])
                self._zeros_no_fitting_count = (0, 0)

            ndx = 3
            while 'zeros' not in pzlines[ndx]:
                ndx += 1

            ndx += 1
            for z_ndx in range(0, z_num):  # skip first zero, so one less to process
                vals = [val.strip() for val in pzlines[ndx].split(',')]
                self.add_zero(complex(float(vals[0]), float(vals[1])))
                ndx += 1

            while 'poles' not in pzlines[ndx]:
                ndx += 1

            ndx += 1
            for p_ndx in range(0, p_num):
                vals = [val.strip() for val in pzlines[ndx].split(',')]
                self.add_pole(complex(float(vals[0]), float(vals[1])))
                ndx += 1

            # dbl check we got the num of values we expected
            if (z_num != self.num_zeros) or (p_num != self.num_poles):
                msg = "Error reading correct number of poles and zeros in file '{}'".format(self._filename)
                logging.error(msg)
                raise Exception(msg)

        else:
            raise Exception('Format error reading "ipaz" type paz file')

    def save(self, filename):
        """
        Save me a pizza!

        :param filename: fully qualified location of file in which to save PAZ
        :type filename: str
        """

        with open(filename, 'wt') as ofl:
            ofl.write(self.PAZ_HEADER_IDA + '\n')
                                                                          # add 1 to indices to make 1-based for humans.
            ofl.write('{:<3} # number of zeros # {}:{} # {}:{} \n'.format(self.num_zeros,
                                                                          ','.join(
                                                                              [str(v+1) for v in self._zeros_default_pert_ndxs[0]]
                                                                          ),
                                                                          ','.join(
                                                                              [str(v+1) for v in self._zeros_default_pert_ndxs[1]]
                                                                          ),
                                                                          self._zeros_no_fitting_count[0],
                                                                          self._zeros_no_fitting_count[1]
                                                                          ))
            ofl.write('{:<3} # number of poles # {}:{} # {}:{} \n'.format(self.num_poles,
                                                                          ','.join(
                                                                              [str(v+1) for v in self._poles_default_pert_ndxs[0]]
                                                                          ),
                                                                          ','.join(
                                                                              [str(v+1) for v in self._poles_default_pert_ndxs[1]]
                                                                          ),
                                                                          self._poles_no_fitting_count[0],
                                                                          self._poles_no_fitting_count[1]
                                                                ))
            ofl.write('\n')
            ofl.write('# zeros\n')
            for zero in self._zeros:
                ofl.write('{:>12.5E}, {:>12.5E}\n'.format(zero.real, zero.imag))
            ofl.write('\n')
            ofl.write('# poles\n')
            for pole in self._poles:
                ofl.write('{:>12.5E}, {:>12.5E}\n'.format(pole.real, pole.imag))

        return

    @property
    def h0(self):
        return self._h0

    @h0.setter
    def h0(self, h0):
        self._h0 = h0

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode):
        self._mode = mode

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        self._units = units

    @property
    def num_poles(self):
        return self._poles.size

    @property
    def num_zeros(self):
        return self._zeros.size

    def add_pole(self, pole):
        ndx = len(self._poles)
        self._poles.resize((ndx + 1,))
        self._poles[ndx] = complex128(pole)

    def add_zero(self, zero):
        ndx = len(self._zeros)
        self._zeros.resize((ndx + 1,))
        self._zeros[ndx] = complex128(zero)

    def zeros(self, mode=None, units=None):

        if mode:
            zero_cnt_dif = PAZ.MODE_ZEROS[mode] - PAZ.MODE_ZEROS[self.mode]
            if zero_cnt_dif > 0:
                zeros = concatenate((npzeros(zero_cnt_dif, dtype=complex128), self._zeros))
            elif zero_cnt_dif < 0:
                if self._zeros.size >= abs(zero_cnt_dif):
                    zeros = self._zeros[abs(zero_cnt_dif):self._zeros.size].copy()
                else:
                    raise Exception("Can't convert PAZ with mode '{}' to mode '{}'".format(self.mode, mode))
            else:
                zeros = self._zeros.copy()
        else:
            zeros = self._zeros.copy()

        if (units == 'hz') and (self.units == 'rad'):
            zeros /= 2 * pi
        elif (units == 'rad') and (self.units == 'hz'):
            zeros *= 2 * pi

        return zeros

    def poles(self, mode=None, units=None):

        poles = self._poles.copy()
        if (units == 'hz') and (self.units == 'rad'):
            poles /= 2 * pi
        elif (units == 'rad') and (self.units == 'hz'):
            poles *= 2 * pi

        return poles

    def perturb_defaults(self):

        poleperdef = deepcopy(self._poles_default_pert_ndxs)
        zeroperdef = deepcopy(self._zeros_default_pert_ndxs)

        return poleperdef, zeroperdef

    def merge_paz_partial(self, paz_partial, paz_map, norm_freq):

        if len(paz_map[0]) > 0:
            self._poles[paz_map[0]] = paz_partial._poles

        if len(paz_map[1]) > 0:
            self._zeros[paz_map[1]] = paz_partial._zeros

        resp = ida.signals.utils.compute_response(array([norm_freq]), self, mode=self.mode)
        self.h0 = 1.0 / abs(resp)

    def make_partial(self, paz_map, norm_freq):

        newpaz = PAZ(mode=self.mode, units=self.units, fileformat=self.fileformat)
        newpaz.h0 = self.h0
        newpaz._poles = self._poles[paz_map[0]]
        newpaz._zeros = self._zeros[paz_map[1]]

        resp = ida.signals.utils.compute_response(array([norm_freq]), newpaz, mode=self.mode)
        newpaz.h0 = 1.0 / abs(resp)

        return newpaz

    def make_partial2(self, norm_freq, partial_mode=PARTIAL_ALL):

        newpaz = PAZ(mode=self.mode, units=self.units, fileformat=self.fileformat)

        if partial_mode == self.PARTIAL_ALL:
            newpaz._poles = copy(self._poles)
            newpaz._zeros = copy(self._zeros)
            newpaz._poles_no_fitting_count = self._poles_no_fitting_count
            newpaz._zeros_no_fitting_count = self._zeros_no_fitting_count
            newpaz._poles_default_pert_ndxs = deepcopy(self._poles_default_pert_ndxs)
            newpaz._zeros_default_pert_ndxs = deepcopy(self._zeros_default_pert_ndxs)
            newpaz.h0 = self.h0
        elif partial_mode == self.PARTIAL_FITTING_LF:
            newpaz._poles = copy(self._poles)[:self.num_poles - self._poles_no_fitting_count[0]]
            newpaz._zeros = copy(self._zeros)[:self.num_zeros - self._zeros_no_fitting_count[0]]
            newpaz._poles_no_fitting_count = (0, 0)
            newpaz._zeros_no_fitting_count = (0, 0)
            newpaz._poles_default_pert_ndxs = deepcopy(self._poles_default_pert_ndxs)
            newpaz._zeros_default_pert_ndxs = deepcopy(self._zeros_default_pert_ndxs)
            resp = ida.signals.utils.compute_response(array([norm_freq]), newpaz, mode=self.mode)
            newpaz.h0 = 1.0 / abs(resp)
        elif partial_mode == self.PARTIAL_FITTING_HF:
            newpaz._poles = copy(self._poles)[:self.num_poles - self._poles_no_fitting_count[1]]
            newpaz._zeros = copy(self._zeros)[:self.num_zeros - self._zeros_no_fitting_count[1]]
            newpaz._poles_no_fitting_count = (0, 0)
            newpaz._zeros_no_fitting_count = (0, 0)
            newpaz._poles_default_pert_ndxs = deepcopy(self._poles_default_pert_ndxs)
            newpaz._zeros_default_pert_ndxs = deepcopy(self._zeros_default_pert_ndxs)
            resp = ida.signals.utils.compute_response(array([norm_freq]), newpaz, mode=self.mode)
            newpaz.h0 = 1.0 / abs(resp)
        elif partial_mode == self.PARTIAL_PERTURBING_LF:
            newpaz._poles = self._poles[self._poles_default_pert_ndxs[0]]
            newpaz._zeros = self._zeros[self._zeros_default_pert_ndxs[0]]
            newpaz._poles_no_fitting_count = (0, 0)
            newpaz._zeros_no_fitting_count = (0, 0)
            newpaz._poles_default_pert_ndxs = ([], [])
            newpaz._zeros_default_pert_ndxs = ([], [])
            resp = ida.signals.utils.compute_response(array([norm_freq]), newpaz, mode=self.mode)
            newpaz.h0 = 1.0 / abs(resp)
        elif partial_mode == self.PARTIAL_PERTURBING_HF:
            newpaz._poles = self._poles[self._poles_default_pert_ndxs[1]]
            newpaz._zeros = self._zeros[self._zeros_default_pert_ndxs[1]]
            newpaz._poles_no_fitting_count = (0, 0)
            newpaz._zeros_no_fitting_count = (0, 0)
            newpaz._poles_default_pert_ndxs = ([], [])
            newpaz._zeros_default_pert_ndxs = ([], [])
            resp = ida.signals.utils.compute_response(array([norm_freq]), newpaz, mode=self.mode)
            newpaz.h0 = 1.0 / abs(resp)

        return newpaz

    def copy(self):

        newpaz = PAZ(mode=self.mode, units=self.units, fileformat=self.fileformat)

        newpaz._poles = self._poles.copy()
        newpaz._zeros = self._zeros.copy()
        newpaz._poles_no_fitting_count =  copy(self._poles_no_fitting_count)
        newpaz._zeros_no_fitting_count =  copy(self._zeros_no_fitting_count)
        newpaz._poles_default_pert_ndxs = copy(self._poles_default_pert_ndxs)
        newpaz._zeros_default_pert_ndxs = copy(self._zeros_default_pert_ndxs)
        newpaz.h0 = self.h0

        return newpaz
