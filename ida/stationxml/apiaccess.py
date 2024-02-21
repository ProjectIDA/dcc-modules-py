#!/usr/bin/env python3
#######################################################################################################################
# Copyright (C) 2020-  Regents of the University of California
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

""" system imports """
import argparse
import urllib, json
import pprint

""" local imports """

""" other imports """

def main():
    response = getJSONFromAPI("stations")
    print(response.json())

def hostURL():
    #return 'http://bargahack:8000/'
    return 'https://idastatus-dev.ucsd.edu:8000/'

def getJSONFromAPI(apiType):
    import requests

    validTypes = ["networks", "stations", "channelepochs", "stages", "instypes", "units"]
    if apiType not in validTypes:
        return False
    
    response = requests.get(url=hostURL() + 'api/' + apiType + '/?format=json')
    return response
    
def getJSONFromFile(apiType):
    validTypes = ["networks", "stations", "channelepochs", "stages", "instypes", "units"]
    if apiType not in validTypes:
        return False
    
    with open(apiType + ".json", 'rt') as jsonFile:
        data = jsonFile.read()
    
    jsonObj = json.loads(data)
    
    return jsonObj 

if __name__ == '__main__':
    main()
