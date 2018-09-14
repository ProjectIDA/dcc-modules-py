#!/usr/bin/env python3
#######################################################################################################################
# Copyright (C) 2018-  Regents of the University of California
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

import argparse
import datetime
import calendar
from os import environ, path, makedirs

################################################################################
""" custom argparse timestamp, check for both month/day date and ordinal date formats """
def datetimeType(stamp):
    try:
        return datetime.datetime.strptime(stamp, "%Y-%m-%dT%H-%M-%S")
    except ValueError:
        try:
            return datetime.datetime.strptime(stamp, "%Y-%m-%d-%H-%M-%S")
        except ValueError:
            try:
                return datetime.datetime.strptime(stamp, "%Y-%jT%H-%M-%S")
            except ValueError:
                try:
                    return datetime.datetime.strptime(stamp, "%Y-%j-%H-%M-%S")
                except ValueError:
                    raise argparse.ArgumentTypeError("".join("Invalid date: \"", stamp, "\""))


################################################################################
""" custom argparse date type """
def dateType(date_str):
    try:
        return datetime.datetime.strptime(date_str, "%Y-%m-%d")
    except ValueError:
        try:
            return datetime.datetime.strptime(date_str, "%Y-%j")
        except ValueError:
            msg = "Date string ({0}) not valid. Expected format, YYYY-MM-DD or YYYY-JJJ".format(date_str)
            raise argparse.ArgumentTypeError(msg)


################################################################################
def makeStationList(stations):
    try:
        return stations.split(',')
    except ValueError:
        raise argparse.ArgumentTypeError("".join("Invalid stations list: \"", stations, "\""))


################################################################################
def makeLocList(locations):
    try:
        return locations.split(',')
    except ValueError:
        raise argparse.ArgumentTypeError("".join("Invalid locations list: \"", locations, "\""))


################################################################################
def makeChanList(channels):
    try:
        return channels.split(',')
    except ValueError:
        raise argparse.ArgumentTypeError("".join("Invalid channels list: \"", channels, "\""))


################################################################################
def mseedOutputFilename(sta, chanloc, day):
    # environ.get second paramter is default if first doesn't exist
    outdir = environ.get('IDA_PPSDPLOT_DATA_DIR', './waveforms')
    if not path.exists(outdir):
        makedirs(outdir)
        
    return "{}/{}.{}.{}.ms".format(outdir, sta, chanloc, day)


################################################################################
def mseedInputFilename(filename):
    # environ.get second paramter is default if first doesn't exist
    outdir = environ.get('IDA_ARCHIVE_MS_DIR', 'mseed')
    if path.exists(outdir):
        return outdir + "/" + filename
    else:
        raise ValueError("Directory missing: {}".format(outdir))
        

################################################################################
def imageOutputFile(sta, loc, chan, start, end, outdirname):
    if not outdirname:
        outdirname = "./images"
    outdir = environ.get('IDA_PPSDPLOT_IMAGE_DIR', outdirname)

    if not path.exists(outdir):
        makedirs(outdir)
        
    return "{}/{}.{}.{}_{}-{}.png".format(outdir, sta, loc, chan, 
                                          start.strftime("%Y%m%dT%H%M%S"), 
                                          end.strftime("%Y%m%dT%H%M%S"))

################################################################################
def makePageImagePath(start, end, outDir):
    if start == end:
        # doing a day, dir = base/year/month/day
        pathName = outDir + "/generated_files/" + str(start.year) + "/" + str(start.month).rjust(2, '0') + "/" + str(start.day).rjust(2, '0')
    elif (end.toordinal() - start.toordinal()) > 1:
        # dir = base/startyear/startmonth
        pathName = outDir + "/generated_files/" + str(start.year) + "/" + str(start.month).rjust(2, '0')

    return pathName

################################################################################
def makeDaysList(month):
    daylist = []
    for day in range(1, calendar.monthrange(month.year, month.month)[1]+1):
        daylist.append(day)
    return daylist


################################################################################
def createDateString(month, daynum):
    return str(datetime.datetime(month.year, month.month, daynum).date())


################################################################################
def lastDayOfMonth(singleDay):
    nextMonth = singleDay.replace(day=28) + datetime.timedelta(days=4)
    return nextMonth - datetime.timedelta(days=nextMonth.day)


