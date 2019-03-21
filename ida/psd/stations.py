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

import os
import glob
import datetime

CHANNELLIST = ["bh1", "bh2", "bhz", "bhn", "bhe"]

LOCATIONLIST = ["00", "10"]

def isValidStationList(stations, start, end):
    for station in stations:
        if not isValidStation(station, start, end):
            return False,station
    return True,None

def isValidStation(station_id, start, end):
    stationList = getStationList(start, end)
    return True if station_id in stationList else False

def isValidChanList(chans):
    for chan in chans:
       if not isValidChan(chan):
           return False, chan
    return True, None

def isValidChan(channel):
    return True if channel in CHANNELLIST else False

def isValidLocList(locs):
    for loc in locs:
       if not isValidLoc(loc):
           return False, loc
    return True, None

def isValidLoc(location):
    return True if location in LOCATIONLIST else False

def getMetadataFilename(station):
    webRoot = os.environ.get('IDA_WEB_ROOT')
    filename = "{}/pickup/dataless/II_{}.*".format(webRoot, station.upper())
    filelist = glob.glob(filename)
    if len(filelist) > 0:
        newest = max(filelist, key=os.path.getctime)
        return(newest)
    else:
        raise OSError("Missing Metadata file for station {}".format(station.upper()))



################################################################################
def getStationList(monthStart, monthEnd):
 
    dbDir = os.environ.get('IDA_DATASCOPEDB_DIR')
    dbSite = dbDir + '/' + "IDA.site"

    try:
        file = open(dbSite, 'r')
    except FileNotFoundError:
        print("File does not exist: {}".format(dbSite))
        return False

    siteList = []
    for line in file:
        siteName = line.split()[0].lower()
        siteActiveStart = datetime.datetime.utcfromtimestamp(float(line.split()[1]))
        siteActiveEnd = datetime.datetime.utcfromtimestamp(float(line.split()[2]))

        if siteActiveStart < monthEnd and siteActiveEnd > monthStart:
            siteList.append(siteName)
        siteList.sort()

    file.close()
    return siteList

################################################################################
