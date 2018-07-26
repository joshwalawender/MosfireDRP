#!/usr/env/python

## Import General Tools
import sys
from os import path
import argparse
import glob
import re
import os
import yaml
import numpy as np

import ccdproc
from astropy.table import Table, Column
from astropy.io import fits

##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description='''
''')
p.add_argument('path', 
               help="Path to files to handle")
args = p.parse_args()
filepath = path.abspath(args.path)
# files = [f for f in args.files if path.splitext(f)[1] != '.original']


##-------------------------------------------------------------------------
## Main Program
##-------------------------------------------------------------------------
def main():
    keywords = ['pwstata7', 'pwstata8', 'pwloca7', 'pwloca8', 'aborted',
                'datafile', 'object', 'truitime', 'maskname', 'mgtname',
                'filter', 'flatspec', 'xoffset', 'yoffset', 'targname',
                'frameid',
                ]
    ifc = ccdproc.ImageFileCollection(filepath, keywords=keywords)
    files = ifc.summary
    photmasks = ['long2pos', 'long2pos_specphot', 'long3pos', 'long3pos_specphot']
    sortnames = Column(['']*len(files), name='sortname', dtype='a90')
    for i,entry in enumerate(files):
        # Remove "(align)" from mask names: treat them the same as normal masks
        if re.search('\s\(align\)', entry['maskname']):
            files[i]['maskname'] = files[i]['maskname'].replace(' (align)', '')

        # Separate targets for photometric masks such as long2pos and LONGSLIT
        maskname = entry['maskname']
        targname = entry['targname'].replace(' ', '_')
        xoffset = float(entry['xoffset'])
        if (maskname in photmasks) or (re.search('LONGSLIT', maskname)):
            sortnames[i] = f"{maskname}_{targname}_{xoffset:.1f}"
        else:
            sortnames[i] = maskname
    files.add_column(sortnames)

    filters = set(files['filter'])
    mask_table = files.group_by('sortname')
    for group in mask_table.groups:
        maskname = group[0]['sortname']
        print(f'Sorting files for: {maskname}')
        for filter in filters:
            mask_files = group[(group['filter'] == filter)
                               & (group['pwstata7'] == 0)
                               & (group['pwstata8'] == 0)
                               & (group['flatspec'] == 0)
                               & (group['aborted'] == False)
                               & (group['mgtname'] != 'mirror')
                               ]
            arc7_files = group[(group['filter'] == filter)
                               & (group['pwstata7'] == 1)
                               & (group['pwstata8'] == 0)
                               & (group['flatspec'] == 0)
                               & (group['aborted'] == False)
                               & (group['mgtname'] != 'mirror')
                               ]
            arc8_files = group[(group['filter'] == filter)
                               & (group['pwstata7'] == 0)
                               & (group['pwstata8'] == 1)
                               & (group['flatspec'] == 0)
                               & (group['aborted'] == False)
                               & (group['mgtname'] != 'mirror')
                               ]
            flat_files = group[(group['filter'] == filter)
                               & (group['pwstata7'] == 0)
                               & (group['pwstata8'] == 0)
                               & (group['flatspec'] == 1)
                               & (group['aborted'] == False)
                               & (group['mgtname'] != 'mirror')
                               ]
            # find lamps off files by parsing OBJECT header keyword
            thermal_flats = mask_files.copy()
            remove_from_thermals = []
            remove_from_mask = []
            for i,entry in enumerate(thermal_flats):
                if re.search('Flat:', entry['object']):
                    remove_from_mask.append(i)
                else:
                    remove_from_thermals.append(i)
            thermal_flats.remove_rows(remove_from_thermals)
            mask_files.remove_rows(remove_from_mask)

            filterpath = path.join(maskname, filter)
            # Write output file for Ne
            if len(arc7_files) > 0:
                print(f'  Found {len(arc7_files)} {filter} Ne files')
                if not path.exists(maskname):
                    os.mkdir(maskname)
                if not path.exists(filterpath):
                    os.mkdir(filterpath)
                with open(path.join(filterpath, 'Ne.txt'), 'w') as Ne_txt:
                    Ne_txt.write(f"{filepath} # Abs. path to files [optional]\n")
                    for entry in arc7_files:
                        Ne_txt.write(f"{entry['file']} # \n")
            # Write output file for Ar
            if len(arc8_files) > 0:
                print(f'  Found {len(arc8_files)} {filter} Ar files')
                if not path.exists(maskname):
                    os.mkdir(maskname)
                if not path.exists(filterpath):
                    os.mkdir(filterpath)
                with open(path.join(filterpath, 'Ar.txt'), 'w') as Ar_txt:
                    Ar_txt.write(f"{filepath} # Abs. path to files [optional]\n")
                    for entry in arc8_files:
                        Ar_txt.write(f"{entry['file']} # \n")
            # Write output file for Flat
            if len(flat_files) > 0:
                print(f'  Found {len(flat_files)} {filter} Flat files')
                if not path.exists(maskname):
                    os.mkdir(maskname)
                if not path.exists(filterpath):
                    os.mkdir(filterpath)
                with open(path.join(filterpath, 'Flat.txt'), 'w') as Flat_txt:
                    Flat_txt.write(f"{filepath} # Abs. path to files [optional]\n")
                    for entry in flat_files:
                        Flat_txt.write(f"{entry['file']} # \n")
            # Write output file for FlatThermal
            if len(thermal_flats) > 0:
                print(f'  Found {len(thermal_flats)} {filter} ThermalFlat files')
                if not path.exists(maskname):
                    os.mkdir(maskname)
                if not path.exists(filterpath):
                    os.mkdir(filterpath)
                with open(path.join(filterpath, 'FlatThermal.txt'), 'w') as FlatThermal_txt:
                    FlatThermal_txt.write(f"{filepath} # Abs. path to files [optional]\n")
                    for entry in thermal_flats:
                        FlatThermal_txt.write(f"{entry['file']} # \n")

            # Sort mask files in to offsets
            if len(mask_files) > 0:
                offset_table = mask_files.group_by('yoffset')
                # Write mask.txt file with mask info
                if not path.exists(maskname):
                    os.mkdir(maskname)
                if not path.exists(filterpath):
                    os.mkdir(filterpath)
                offsets = []
                for group in offset_table.groups:
                    offset = group[0]['yoffset']
                    offsets.append(offset)
                    print(f'  Found {len(group)} {filter} Offset {offset} files')
                    with open(path.join(filterpath, f'Offset_{offset:.1f}.txt'), 'w') as offset_txt:
                        offset_txt.write(f"{filepath} # Abs. path to files [optional]\n")
                        for entry in group:
                            offset_txt.write(f"{entry['file']} # \n")
                # write out YAML file with dictionary of information
                mean_itime = np.mean(mask_files['truitime'])
                stddev_itime = np.std(mask_files['truitime'])
                info = {'maskname': str(maskname),
                        'exptime': float(mean_itime),
                        'exptime_stddev': float(stddev_itime),
                        'filter': str(filter),
                        'offsets': [f"{o:.1f}" for o in offsets],
                        'offset_files': [f'Offset_{o:.1f}.txt' for o in offsets],
                        'mask_files': [str(f) for f in mask_files['file']],
                        'Ne_files': [str(f) for f in arc7_files['file']],
                        'Ar_files': [str(f) for f in arc8_files['file']],
                        'flat_files': [str(f) for f in flat_files['file']],
                        'thermal_flats': [str(f) for f in thermal_flats['file']],
                        }
                with open(path.join(filterpath, 'mask.txt'), 'w') as mask_txt:
                    yaml.dump(info, mask_txt)

            # Hack headers for long2pos
            #   Headers have the wide position with FRAMEID as A which is
            #   identical to the first narrow position.
            #   For each FITS file with zero offset, change FRAMEID to 'C'.
            if maskname[:8] == 'long2pos':
                for entry in mask_files[mask_files['yoffset'] == 0]:
                    filename = entry['file']
                    file = path.join(filepath, filename)
                    hdul = fits.open(file, mode='update')
                    if hdul[0].header['FRAMEID'] == 'A':
                        print(f'Correcting header for {filename}')
                        hdul[0].header.set('FRAMEID', 'C')
                        hdul.flush()


if __name__ == '__main__':
    main()
