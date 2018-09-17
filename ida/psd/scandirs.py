#!/usr/bin/env python3

import argparse
import os
import glob
import json

def main():
    parser = argparse.ArgumentParser(description="Walk the given directory and gather dict of subdirectories")
    parser.add_argument("-d", "--datadir",
                        help="the directory of PSD data", 
                        action="store",
                        type = str,
                        required=True)
    args = parser.parse_args()
    generateJSONFile(args.datadir)

def generateJSONFile(dirName):
    os.chdir(dirName + "/generated_files")

    entries_dict = process_dir("")

    os.chdir(dirName)
    with open('./psd-dates.json', 'w') as fp:
            json.dump(entries_dict, fp)
    
def process_dir(idroot):
    # must start in root of tree you're processing
    entries = []
    separator = '/'

    for entry in sorted(glob.glob('*')):
        if entry in ['.', '..']:
            continue

        if os.path.isdir(entry):
            items = [idroot, entry]
            nextroot = separator.join( str(x) for x in items)

            # the keys to this dictionary are required by the treeview JS module
            new_entry = {}
            new_entry["id"] = nextroot.lstrip(separator) # trim leading '/'
            dateList = new_entry["id"].split(separator)
            dateStr = ""
            # need to do this more pythonically
            if len(dateList) is 1:
                dateStr = dateList[0]
            if len(dateList) is 2:
                dateStr = dateList[1] + separator + dateList[0]
            if len(dateList) is 3:
                dateStr = dateList[1] + separator + dateList[2] + separator + dateList[0]
            
            new_entry["text"] = dateStr

            os.chdir(entry)
            new_entry['children'] = process_dir(nextroot)
            os.chdir("..")
            
            entries.append(new_entry)
    
    return entries

################################################################################
if __name__ == '__main__':
    main()
