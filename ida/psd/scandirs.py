#!/usr/bin/env python3

import argparse
import os
import glob
import json
import pprint
from calendar import month_abbr, monthrange

def process_dir(rootdir):
    """ must start in root of tree you're processing. it should be the directory containing YYYY subdirs
        Constructs dir datatree that can be used by jinga templates to construct the nav menu

        Does not use recursion due to somewhat specific logic needed at year/month and day levels.
    """
    savedir = os.getcwd()
    os.chdir(rootdir)

    datetree = []

    for ayear in sorted(glob.glob('????')):

        if os.path.isdir(ayear):

            yearinfo = { 'year': ayear,
                         'months': []
                       }
            os.chdir(ayear)

            for amonth in sorted(glob.glob('??')):

                monthinfo = { 'month': amonth,
                              'monthname': month_abbr[int(amonth)],
                              'days': [],
                              'complete': False
                            }

                os.chdir(amonth)
                for aday in sorted(glob.glob('??')):
                    monthinfo['days'].append(aday)

                # if we have all days, mark month as "complete" used in template for partial month
                if len(monthinfo['days']) == monthrange(int(ayear), int(monthinfo['month']))[1]:
                    monthinfo['complete'] = True

                yearinfo['months'].append(monthinfo)
                os.chdir("..")

            datetree.append(yearinfo)
            os.chdir("..")

    os.chdir(savedir)

    return datetree



#def process_dir(idroot):
#    # must start in root of tree you're processing
#    entries = []
#    separator = '/'
#
#    for entry in sorted(glob.glob('*')):
#        if entry in ['.', '..']:
#            continue
#
#        if os.path.isdir(entry):
#            items = [idroot, entry]
#            nextroot = separator.join( str(x) for x in items)
#
#            new_entry = {}
#            new_entry["id"] = nextroot.lstrip(separator) # trim leading '/'
#            dateList = new_entry["id"].split(separator)
#            dateStr = ""
#            # need to do this more pythonically
#            if len(dateList) is 1:
#                dateStr = dateList[0]
#            if len(dateList) is 2:
#                dateStr = dateList[1] + separator + dateList[0]
#            if len(dateList) is 3:
#                dateStr = dateList[1] + separator + dateList[2] + separator + dateList[0]
#            
#            new_entry["text"] = dateStr
#
#            os.chdir(entry)
#            new_entry['children'] = process_dir(nextroot)
#            os.chdir("..")
#            
#            entries.append(new_entry)
#
#    pp = pprint.PrettyPrinter(indent=4)
#    print("\n\nENTRIES: \n\n")
#    pp.pprint(entries)
#    print("\n\n")
#
#    return entries
#
#################################################################################
