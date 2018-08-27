import os
import yaml
import argparse

from MOSFIRE import IO, Wavelength

def nfiles(list_file):
    if os.path.exists(list_file) is False:
        return 0
    else:
        contents = IO.list_file_to_strings(list_file)
        return len(contents)

def main(nf=False):
    with open('mask.txt', 'r') as mask_txt:
        info = yaml.load(mask_txt)

    neon = False
    argon = False
    maskname = info['maskname']
    band = info['filter']
    obsfiles = info['offset_files']

    longexp_threshold = {'Y': 170, 'J': 110, 'H': 110, 'K': 1000}

    longexp = info['exptime'] >= longexp_threshold[band]

    if longexp is True:
        wavfiles = obsfiles
    else:
        wavfiles = []
        # use only Ne if available, otherwise use Ar
        if nfiles('Ne.txt') > 0:
            neon = True
            wavfiles.append('Ne.txt')
        elif nfiles('Ar.txt') > 0:
            argon = True
            wavfiles.append('Ar.txt')
        # if neither Ne nor Ar is available, we must use the skylines
        if len(wavfiles) == 0:
            wavfiles = obsfiles

    redfiles = [f"eps_{file}.fits" for file in obsfiles]
    wavename = Wavelength.filelist_to_wavename(IO.list_file_to_strings(wavfiles), band, maskname, '')
    wavefile = f"lambda_solution_{wavename}"

    template = f'''
maskname = '{maskname}'
band = '{band}'
obsfiles = {obsfiles}
wavfiles = {wavfiles}
redfiles = {redfiles}
Wavelength_file = '{wavefile}'

# User Options:
noninteractiveflag = {str(nf):5s} # Set this flag to autofit wavelength
                           # solution instead of manually verifying the fit.

# Imports
import matplotlib           # Force TkAgg backend for interactivity. This is
matplotlib.use('TkAgg')     # critical to bypass a bug in the MacOSX backend.
from MOSFIRE import Options, Flats, Wavelength, Background, Rectify, Extract

Flats.handle_flats('Flat.txt', maskname, band, Options.flat)

Wavelength.imcombine(wavfiles, maskname, band, Options.wavelength)
Wavelength.fit_lambda_interactively(maskname, band, wavfiles, Options.wavelength,
           neon={neon}, argon={argon}, noninteractive=noninteractiveflag)
Wavelength.fit_lambda(maskname, band, wavfiles, wavfiles, Options.wavelength)
Wavelength.apply_lambda_simple(maskname, band, wavfiles, Options.wavelength,
           smooth={not longexp})

Background.handle_background(obsfiles, Wavelength_file, maskname, band,
           Options.wavelength)

Rectify.handle_rectification(maskname, redfiles, Wavelength_file, band,
        obsfiles, Options.wavelength)

Extract.extract_spectra(maskname, band, width=10,
        interactive=(not noninteractiveflag))

# If you have questions, please submit a ticket on the github issue page:
# https://github.com/Keck-DataReductionPipelines/MosfireDRP/issues
    '''

    with open('Driver.py', 'w') as driver:
        driver.write(template)

if __name__ == '__main__':
    ##-------------------------------------------------------------------------
    ## Parse Command Line Arguments
    ##-------------------------------------------------------------------------
    ## create a parser object for understanding command-line arguments
    p = argparse.ArgumentParser(description='''
    ''')
    p.add_argument('command', help="Command executed (i.e. handle)")
    ## add flags
    p.add_argument("-n", "--ni", dest="noninteractive",
        default=False, action="store_true",
        help="Run in non-interactive mode.")
    args = p.parse_args()

    main(args.noninteractive)
